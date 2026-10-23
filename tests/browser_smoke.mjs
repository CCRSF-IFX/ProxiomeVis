// Run against the demo app and a local Chromium started with --remote-debugging-port=9223.
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';

const debugUrl = process.env.PROXIOME_CDP_URL || 'http://127.0.0.1:9223';
const appUrl = process.env.PROXIOME_SMOKE_URL || 'http://127.0.0.1:3838';
const tabs = await (await fetch(`${debugUrl}/json`)).json();
const socket = new WebSocket(tabs.find(tab => tab.type === 'page').webSocketDebuggerUrl);
await new Promise(resolve => socket.addEventListener('open', resolve, { once: true }));
let serial = 0;
const pending = new Map();
socket.onmessage = ({ data }) => {
  const message = JSON.parse(data);
  pending.get(message.id)?.(message);
  pending.delete(message.id);
};
const call = (method, params = {}) => new Promise((resolve, reject) => {
  const id = ++serial;
  const timeout = setTimeout(() => reject(new Error(`Timed out: ${method}`)), 30000);
  pending.set(id, message => {
    clearTimeout(timeout);
    if (message.error) reject(new Error(JSON.stringify(message.error)));
    else resolve(message.result);
  });
  socket.send(JSON.stringify({ id, method, params }));
});
const evaluate = async expression => {
  const value = await call('Runtime.evaluate', { expression, returnByValue: true, awaitPromise: true });
  if (value.exceptionDetails) throw new Error(JSON.stringify(value.exceptionDetails));
  return value.result.value;
};
const delay = milliseconds => new Promise(resolve => setTimeout(resolve, milliseconds));
async function waitFor(expression) {
  for (let attempt = 0; attempt < 100; attempt++) {
    if (await evaluate(expression)) return;
    await delay(200);
  }
  throw new Error(`Condition not reached: ${expression}`);
}
async function settled() {
  await delay(400);
  await waitFor("!document.documentElement.classList.contains('shiny-busy')");
  const errors = await evaluate("[...document.querySelectorAll('.shiny-output-error:not(.shiny-output-error-validation)')].map(x => x.textContent)");
  assert.deepEqual(errors, []);
}
async function tab(label) {
  await evaluate(`[...document.querySelectorAll('a[data-bs-toggle=tab]')].find(x => x.innerText === ${JSON.stringify(label)}).click()`);
  await settled();
}
async function select(id, value) {
  const selector = `document.getElementById(${JSON.stringify(`report-${id}`)})`;
  await evaluate(`(() => {
    const element = ${selector};
    const radio = [...element.querySelectorAll('input[type=radio]')].find(x => x.value === ${JSON.stringify(value)});
    if (radio) radio.click();
    else if (element.selectize) element.selectize.setValue(${JSON.stringify(value)});
    else { element.value = ${JSON.stringify(value)}; element.dispatchEvent(new Event('change', {bubbles:true})); }
  })()`);
  await settled();
  await waitFor(`Shiny.shinyapp.$inputValues[${JSON.stringify(`report-${id}`)}] === ${JSON.stringify(value)}`);
}
async function download(id, pattern) {
  await waitFor(`document.getElementById('report-${id}')?.getAttribute('href')?.includes('download')`);
  const text = await evaluate(`fetch(document.getElementById('report-${id}').href).then(async response => { if (!response.ok) throw new Error(response.status); return await response.text(); })`);
  assert.match(text, pattern);
}

async function fullscreenZoom(id) {
  const q = JSON.stringify(id);
  await evaluate("document.getElementById("+q+").closest('.card').querySelector('.bslib-full-screen-enter').click()");
  await waitFor("document.getElementById("+q+").closest('.card').dataset.fullScreen === 'true'");
  await settled();
  await evaluate("new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))");
  const before = await evaluate("document.getElementById("+q+")._fullLayout.xaxis.range.slice()");
  const point = await evaluate("(()=>{const r=document.querySelector('#'+"+q+"+' .nsewdrag').getBoundingClientRect();return {x:r.x+r.width/2,y:r.y+r.height/2};})()");
  await call('Input.dispatchMouseEvent', {type:'mouseMoved', ...point});
  await call('Input.dispatchMouseEvent', {type:'mouseWheel', ...point, deltaX:0, deltaY:-120});
  await waitFor("(()=>{const r=document.getElementById("+q+")._fullLayout.xaxis.range;return r[1]-r[0]<"+(before[1]-before[0])+";})()");
  await evaluate("document.querySelector('.bslib-full-screen-exit').click()");
}

try {
  await call('Emulation.setDeviceMetricsOverride', { width: 1600, height: 1100, deviceScaleFactor: 1, mobile: false });
  await call('Page.navigate', { url: appUrl });
  await waitFor("document.getElementById('report-metric_cells')?.textContent.trim() === '4,024'");
  await settled();
  assert.equal(await evaluate("document.getElementById('rds_path') === null"), true);
  assert.equal(await evaluate("!!document.getElementById('report-spatial_report_plot')?.data"), false);
  await download('download_provenance', /source_table_manifest/);
  console.log('PASS: demo load, QC, provenance, locked server path');

  await tab('Embeddings');
  await waitFor("document.getElementById('report-embedding_plot')?.data?.length > 0");
  assert.equal(await evaluate("document.getElementById('report-embedding_reduction').value"), 'umap.cca');
  await select('embedding_split', 'donor');
  await select('embedding_split', 'sample_id');
  await select('embedding_color', 'abundance');
  await download('download_embedding_svg', /<svg/);
  console.log('PASS: single embedding, split controls, marker abundance, SVG download');

  await evaluate("document.getElementById('dataset').click()");
  await waitFor("document.getElementById('shiny-modal')?.innerText.includes('Proximity sidecar')");
  assert.equal(await evaluate("document.getElementById('rds_path') === null"), true);
  await evaluate("document.querySelector('#shiny-modal .modal-footer button').click()");

  await tab('Results');
  await waitFor("document.getElementById('report-effect_detail')?.data?.some(x => x.type === 'box')");
  assert.equal(await evaluate("document.querySelector('#report-results_view .nav-link.active')?.getAttribute('data-value')"), 'Tables');
  assert.equal(await evaluate("document.querySelector('#report-results_view a[data-value=\"Summary\"]') === null"), true);
  assert.equal(await evaluate("document.querySelector('#report-results_view a[data-value=\"Spatial patterns\"]').getAttribute('aria-disabled')"), 'true');
  await download('download_detail_svg', /<svg/);
  await select('donor', 'd1');
  await select('celltype', 'CD4 T');
  for (const readout of ['abundance', 'self_proximity', 'colocalization']) {
    await select('readout', readout);
    await waitFor("document.querySelector('#report-results_table table tbody tr') !== null");
    await evaluate("document.querySelector('#report-results_table table tbody tr').click()");
    await settled();
    await waitFor("document.getElementById('report-effect_detail')?.data?.some(x => x.type === 'box')");
    assert.match(await evaluate("document.getElementById('report-effect_detail').layout.title.text"), /cell-level distribution/);
    await download('download_results_csv', /effect_size/);
  }
  assert.match(await evaluate("document.getElementById('report-result_summary').innerText"), /Donor: d1.*CD4 T/);
  assert.match(await evaluate("document.getElementById('report-result_summary').innerText"), /not donor-population inference/);
  console.log('PASS: shared single-donor scope, all readouts, Dataset modal, stored CSV');

  await tab('Volcano');
  for (const readout of ['abundance', 'self_proximity', 'colocalization']) {
    await select('readout', readout);
    await waitFor("document.getElementById('report-report_volcano_plot')?.data?.some(t => t.mode === 'markers' && t.x.length)");
    assert.equal(await evaluate("document.getElementById('report-report_volcano_plot')._context.scrollZoom"), true);
    await evaluate("(()=>{const p=document.getElementById('report-report_volcano_plot');const i=p.data.findIndex(t=>t.mode==='markers'&&t.x.length);Plotly.Fx.hover(p,[{curveNumber:i,pointNumber:0}]);})()");
    await waitFor("document.querySelector('#report-report_volcano_plot .hoverlayer')?.textContent.includes('Adjusted P value:')");
    const key = await clickFeature('report-report_volcano_plot');
    await waitFor("document.getElementById('report-volcano_distribution')?.data?.some(x => x.type === 'box')");
    await waitFor(JSON.stringify(JSON.parse(key))+".every(marker => document.getElementById('report-volcano_distribution').layout.title.text.includes(marker))");
    assert.equal(await evaluate("(()=>{const v=document.getElementById('report-report_volcano_plot').closest('.card').getBoundingClientRect();const d=document.getElementById('report-volcano_distribution').closest('.card').getBoundingClientRect();return v.right<=d.left+2&&Math.abs(v.top-d.top)<20;})()"), true);
    await download('download_detail_svg', /<svg/);
    await download('download_report_volcano_svg', /<svg/);
    await download('download_report_volcano_csv', /exploratory_component_level/);
  }
  async function clickFeature(id) {
    const key = await evaluate("(()=>{const p=document.getElementById("+JSON.stringify(id)+");const current=Shiny.shinyapp.$inputValues['report-detail_feature'];const i=p.data.findIndex(t=>t.mode==='markers'&&t.x.length&&t.key?.some(k=>k!==current));const d=p.data[i];const n=d.key.findIndex(k=>k!==current);p.emit('plotly_click',{points:[{curveNumber:i,pointNumber:n,x:d.x[n],y:d.y[n],data:d,fullData:p._fullData[i]}]});return d.key[n];})()");
    await waitFor("Shiny.shinyapp.$inputValues['report-detail_feature'] === "+JSON.stringify(key));
    await settled();
    return key;
  }
  const selectedKey = await clickFeature('report-report_volcano_plot');
  await fullscreenZoom('report-report_volcano_plot');
  console.log('PASS: volcano hover/zoom, side-by-side distributions, linked feature clicks for every readout');
  await tab('Tables');
  assert.equal(await evaluate("Shiny.shinyapp.$inputValues['report-detail_feature']"), selectedKey);
  await waitFor("document.querySelector('#report-results_table table') !== null");
  await evaluate("document.querySelector('#report-results_table table tbody tr').click()");
  await settled();
  const fullRows = await evaluate("fetch(document.getElementById('report-download_results_csv').href).then(r=>r.text())");
  await evaluate("$('#report-results_table table').DataTable().search('NO_SUCH_FEATURE_987654').draw()");
  await settled();
  await waitFor("Shiny.shinyapp.$inputValues['report-results_table_rows_all']?.length === 0");
  const filteredRows = await evaluate("fetch(document.getElementById('report-download_filtered_results_csv').href).then(r=>r.text())");
  assert.equal(filteredRows.trim().split('\n').length, 1);
  assert.equal(await evaluate("fetch(document.getElementById('report-download_results_csv').href).then(r=>r.text())"), fullRows);
  await evaluate("$('#report-results_table table').DataTable().search('').draw()");
  await settled();
  console.log('PASS: table click and zero-match CSV leaves all-rows CSV unchanged');

  await tab('Spatial patterns');
  assert.equal(await evaluate("document.querySelector('#report-results_view a[data-value=\"Spatial patterns\"]').getAttribute('aria-disabled')"), 'false');
  await waitFor("document.getElementById('report-spatial_report_plot')?.data?.some(t=>t.mode==='markers'&&t.x.length)");
  assert.equal(await evaluate("document.getElementById('report-spatial_report_plot')._context.scrollZoom"), true);
  await clickFeature('report-spatial_report_plot');
  await fullscreenZoom('report-spatial_report_plot');
  await download('download_spatial_csv', /marker_1/);
  await download('download_spatial_report_svg', /<svg/);
  await select('spatial_kind', 'profiles');
  await select('spatial_panel', 'PHA');
  await download('download_spatial_report_svg', /<svg/);
  await mkdir('output', { recursive: true });
  await evaluate("window.scrollTo(0,0)");
  await writeFile('output/results-explorer.png', Buffer.from((await call('Page.captureScreenshot', {format:'png'})).data,'base64'));
  console.log('PASS: spatial exploration, pair selection, profiles, original SVG and rows');

  await select('spatial_view', 'network');
  await waitFor("document.getElementById('report-network_plot')?.data?.some(t => t.mode === 'markers' && t.key?.length)");
  await evaluate("document.getElementById('report-network_markers').selectize.setValue(['CD19','CD37','CD53','CD81','CD82'])");
  await settled();
  await waitFor("Shiny.shinyapp.$inputValues['report-network_markers']?.length === 5");
  const networkKey = await clickFeature('report-network_plot');
  assert.ok(JSON.parse(networkKey).every(marker => ['CD19','CD37','CD53','CD81','CD82'].includes(marker)));
  await download('download_network_svg', /<svg/);
  await download('download_network_csv', /minimum_detection_fraction/);
  await tab('Volcano');
  assert.equal(await evaluate("Shiny.shinyapp.$inputValues['report-detail_feature']"), networkKey);
  await tab('Spatial patterns');
  await select('readout', 'self_proximity');
  await waitFor("Shiny.shinyapp.$inputValues['report-spatial_view'] === 'matrix'");
  assert.equal(await evaluate("!!document.querySelector('#report-spatial_view input[value=network]')"), false);
  console.log('PASS: network marker list, linked edges, exports, and matrix fallback for self-proximity');

  await select('readout', 'abundance');
  await waitFor("document.querySelector('#report-results_view .nav-link.active')?.getAttribute('data-value') === 'Tables'");
  assert.equal(await evaluate("document.querySelector('#report-results_view a[data-value=\"Spatial patterns\"]').getAttribute('tabindex')"), '-1');
  await evaluate("document.querySelector('#report-results_view a[data-value=\"Spatial patterns\"]').click()");
  await settled();
  assert.equal(await evaluate("document.querySelector('#report-results_view .nav-link.active')?.getAttribute('data-value')"), 'Tables');
  await select('readout', 'colocalization');
  await tab('Spatial patterns');
  await clickFeature('report-spatial_report_plot');
  console.log('PASS: abundance disables Spatial patterns and returns to Tables; spatial readouts enable it');

  await evaluate("document.getElementById('report-explore_3d').click()");
  await waitFor("Shiny.shinyapp.$inputValues['report-report_page'] === 'Cell 3D'");
  await waitFor("Shiny.shinyapp.$inputValues['report-graph-condition'] === ''");
  assert.equal(await evaluate("Shiny.shinyapp.$inputValues['report-graph-donor']"), 'd1');
  assert.equal(await evaluate("Shiny.shinyapp.$inputValues['report-graph-celltype']"), 'CD4 T');
  assert.equal(await evaluate("!!document.querySelector('#report-graph-plot canvas')"), false);
  assert.equal(await evaluate("Shiny.shinyapp.$inputValues['report-graph-markers']"), '');
  console.log('PASS: 3D handoff waits for condition, highlight and Display');

  await call('Emulation.setDeviceMetricsOverride', {width:390,height:844,deviceScaleFactor:1,mobile:true});
  await evaluate("document.querySelector('.navbar-toggle').click()");
  await tab('Results');
  await tab('Volcano');
  await evaluate("window.scrollTo(0,0)");
  assert.ok(await evaluate("document.documentElement.scrollWidth <= window.innerWidth + 2"));
  await writeFile('output/results-mobile.png', Buffer.from((await call('Page.captureScreenshot', {format:'png'})).data,'base64'));
  console.log('PASS: mobile navigation without horizontal overflow');
} finally {
  socket.close();
}
