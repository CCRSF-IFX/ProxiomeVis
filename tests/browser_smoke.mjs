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

try {
  await call('Emulation.setDeviceMetricsOverride', { width: 1600, height: 1100, deviceScaleFactor: 1, mobile: false });
  await call('Page.navigate', { url: appUrl });
  await waitFor("document.getElementById('report-metric_cells')?.textContent.trim() === '4,024'");
  await settled();
  assert.equal(await evaluate("document.getElementById('rds_path') === null"), true);
  await download('download_provenance', /source_table_manifest/);
  console.log('PASS: demo load, QC, provenance, locked server path');

  await tab('Embeddings');
  await waitFor("document.getElementById('report-umap_plot')?.data?.length > 0");
  assert.equal(await evaluate("document.getElementById('report-umap_reduction').value"), 'umap.cca');
  await select('embedding_color', 'abundance');
  await download('download_umap_svg', /<svg/);
  console.log('PASS: default UMAP/PCA, marker abundance, SVG download');

  await tab('Precomputed results');
  for (const readout of ['abundance', 'self_proximity', 'colocalization']) {
    await select('readout', readout);
    await select('donor', 'd1');
    await select('celltype', 'CD4 T');
    await waitFor("document.getElementById('report-effect_heatmap')?.data?.length > 0");
    await download('download_results_csv', /effect_size/);
  }
  await select('donor', '__all__');
  assert.match(await evaluate("document.getElementById('report-fdr_controls').innerText"), /Descriptive result/);
  console.log('PASS: three readouts, donor/cell-type filters, CSV, descriptive FDR state');

  await tab('Report volcano');
  for (const readout of ['abundance', 'self_proximity', 'colocalization']) {
    await select('report_readout', readout);
    await select('report_donor', 'd1');
    await select('report_celltype', 'CD4 T');
    await waitFor("document.getElementById('report-report_volcano_plot')?.data?.some(t => t.mode === 'markers' && t.x.length)");
    assert.equal(await evaluate("document.getElementById('report-report_volcano_plot')._context.scrollZoom"), true);
    assert.equal(await evaluate("[...document.querySelectorAll('#report-report_volcano_plot .annotation-text')].some(x => x.textContent === '-log10(adjusted P)')"), true);
    await evaluate(`(() => {
      const plot = document.getElementById('report-report_volcano_plot');
      const curveNumber = plot.data.findIndex(t => t.mode === 'markers' && t.x.length);
      Plotly.Fx.hover(plot, [{curveNumber, pointNumber: 0}]);
    })()`);
    await waitFor("document.querySelector('#report-report_volcano_plot .hoverlayer')?.textContent.includes('Adjusted P value:')");
    assert.match(await evaluate("document.querySelector('#report-report_volcano_plot .hoverlayer').textContent"), /Marker/);
    if (readout === 'self_proximity') {
      await mkdir('output', { recursive: true });
      const screenshot = await call('Page.captureScreenshot', { format: 'png' });
      await writeFile('output/volcano-interactive.png', Buffer.from(screenshot.data, 'base64'));
    }
    const before = await evaluate("document.getElementById('report-report_volcano_plot').layout.xaxis.range");
    const box = await evaluate(`(() => {
      const r = document.querySelector('#report-report_volcano_plot .nsewdrag').getBoundingClientRect();
      return {x: r.x + r.width * .25, y: r.y + r.height * .25, width: r.width * .5, height: r.height * .5};
    })()`);
    await call('Input.dispatchMouseEvent', { type: 'mousePressed', x: box.x, y: box.y, button: 'left', clickCount: 1 });
    await call('Input.dispatchMouseEvent', { type: 'mouseMoved', x: box.x + box.width, y: box.y + box.height, button: 'left', buttons: 1 });
    await call('Input.dispatchMouseEvent', { type: 'mouseReleased', x: box.x + box.width, y: box.y + box.height, button: 'left', clickCount: 1 });
    await waitFor(`(() => { const r = document.getElementById('report-report_volcano_plot').layout.xaxis.range; return r[1] - r[0] < ${before[1] - before[0]}; })()`);
    await evaluate("document.querySelector('#report-report_volcano_plot [data-title=\"Reset axes\"]').click()");
    await waitFor(`Math.abs(document.getElementById('report-report_volcano_plot').layout.xaxis.range[0] - ${before[0]}) < 1e-6`);
    await call('Input.dispatchMouseEvent', { type: 'mouseWheel', x: box.x, y: box.y, deltaX: 0, deltaY: -100 });
    await waitFor(`(() => { const r = document.getElementById('report-report_volcano_plot').layout.xaxis.range; return r[1] - r[0] < ${before[1] - before[0]}; })()`);
    await download('download_report_volcano_svg', /<svg/);
    await download('download_report_volcano_csv', /exploratory_component_level/);
    if (readout === 'abundance') await download('download_report_abundance_heatmap_svg', /<svg/);
  }
  await select('report_celltype', 'Platelet-rich');
  await waitFor("document.body.innerText.includes('No stored rows match this display scope.')");
  await download('download_report_volcano_svg', /<svg/);
  console.log('PASS: interactive volcano hover, box zoom/reset, frozen abundance heatmap, report exports');

  await tab('Spatial report');
  await select('spatial_celltype', 'CD4 T');
  await select('spatial_donor', 'd2');
  await download('download_spatial_report_svg', /<svg/);
  await select('spatial_kind', 'profiles');
  await select('spatial_panel', 'PHA');
  await download('download_spatial_report_svg', /<svg/);
  await mkdir('output', { recursive: true });
  const screenshot = await call('Page.captureScreenshot', { format: 'png' });
  await writeFile('output/browser-smoke.png', Buffer.from(screenshot.data, 'base64'));
  console.log('PASS: aligned spatial difference/profile views; output/browser-smoke.png');
} finally {
  socket.close();
}
