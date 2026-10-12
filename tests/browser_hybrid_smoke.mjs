// Prepared hybrid demo on :3838 and Chromium with CDP on :9223 (WebGL enabled).
import assert from 'node:assert/strict';
import { mkdir, writeFile } from 'node:fs/promises';

const debugUrl = process.env.PROXIOME_CDP_URL || 'http://127.0.0.1:9223';
const appUrl = process.env.PROXIOME_SMOKE_URL || 'http://127.0.0.1:3838';
const tabs = await (await fetch(`${debugUrl}/json`)).json();
const socket = new WebSocket(tabs.find(t => t.type === 'page').webSocketDebuggerUrl);
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
const delay = ms => new Promise(resolve => setTimeout(resolve, ms));
async function waitFor(expression) {
  for (let n = 0; n < 180; n++) {
    if (await evaluate(expression)) return;
    await delay(250);
  }
  throw new Error(`Condition not reached: ${expression}\n${await evaluate("document.body.innerText.slice(-5000)")}`);
}
async function settled() {
  await delay(400);
  await waitFor("!document.documentElement.classList.contains('shiny-busy')");
  assert.deepEqual(await evaluate("[...document.querySelectorAll('.shiny-output-error:not(.shiny-output-error-validation)')].filter(x => x.offsetParent).map(x => x.textContent)"), []);
}
async function tab(label) {
  await evaluate(`[...document.querySelectorAll('a[data-bs-toggle=tab]')].find(x => x.innerText === ${JSON.stringify(label)}).click()`);
  await settled();
}
async function select(id, value) {
  await evaluate(`document.getElementById(${JSON.stringify(id)}).selectize.setValue(${JSON.stringify(value)})`);
  await settled();
  await waitFor(`JSON.stringify(Shiny.shinyapp.$inputValues[${JSON.stringify(id)}]) === JSON.stringify(${JSON.stringify(value)})`);
}
async function download(id) {
  await waitFor(`document.getElementById(${JSON.stringify(id)})?.getAttribute('href')?.includes('download')`);
  return evaluate(`fetch(document.getElementById(${JSON.stringify(id)}).href).then(async r => { if (!r.ok) throw new Error(r.status); return await r.text(); })`);
}
async function run() {
  await evaluate("document.getElementById('report-analysis-run').click()");
  await waitFor("document.getElementById('report-analysis-run').disabled");
  // A live worker must not prevent navigation in the same session.
  await tab('Overview / QC');
  await tab('Custom analysis');
  await waitFor("document.getElementById('report-analysis-status').textContent.startsWith('Completed.')");
  await settled();
}
async function screenshot(name) {
  await mkdir('output', { recursive: true });
  const shot = await call('Page.captureScreenshot', { format: 'png' });
  await writeFile(`output/${name}.png`, Buffer.from(shot.data, 'base64'));
}
async function checkControlContrast() {
  const failures = await evaluate(`(() => {
    const rgb = value => value.match(/[\\d.]+/g).map(Number);
    const luminance = c => c.slice(0,3).map(v => v/255).map(v => v <= .04045 ? v/12.92 : ((v+.055)/1.055)**2.4)
      .reduce((sum,v,i) => sum + v*[.2126,.7152,.0722][i],0);
    return [...document.querySelectorAll('.btn,.alert')].filter(e => e.offsetParent && !e.disabled).flatMap(e => {
      let bg = [255,255,255];
      for(let parent=e;parent;parent=parent.parentElement) {
        const candidate = rgb(getComputedStyle(parent).backgroundColor);
        if(candidate.length===3 || candidate[3]===1) { bg=candidate; break; }
      }
      const fg=rgb(getComputedStyle(e).color), a=luminance(fg), b=luminance(bg);
      const ratio=(Math.max(a,b)+.05)/(Math.min(a,b)+.05);
      return ratio >= 4.5 ? [] : [{text:e.textContent.trim().slice(0,40),ratio}];
    });
  })()`);
  assert.deepEqual(failures, []);
}

try {
  await call('Page.bringToFront');
  await call('Emulation.setDeviceMetricsOverride', { width: 1600, height: 1100, deviceScaleFactor: 1, mobile: false });
  await call('Page.navigate', { url: appUrl });
  await waitFor("document.getElementById('report-metric_cells')?.textContent.trim() === '4,024'");
  await tab('Custom analysis');
  await select('report-analysis-donor', 'd1');
  await select('report-analysis-celltype', 'CD4 T');
  await select('report-analysis-target', 'PHA');
  await select('report-analysis-reference', 'unstim');
  await select('report-analysis-markers', ['CD45', 'HLA-ABC', 'B2M', 'CD4']);
  await evaluate("document.getElementById('report-analysis-run').click()");
  await waitFor("document.getElementById('report-analysis-status').textContent.includes('Confirm the exploratory')");
  await evaluate("document.getElementById('report-analysis-acknowledge').click()");
  await settled();

  for (const readout of ['abundance', 'self_proximity', 'colocalization']) {
    await select('report-analysis-readout', readout);
    await run();
    const provenance = JSON.parse(await download('report-analysis-provenance_json'));
    assert.equal(provenance.request.readout, readout);
    assert.equal(provenance.request.donor, 'd1');
    assert.equal(provenance.request.celltype, 'CD4 T');
    assert.match(provenance.provenance.bundle_md5, /^[a-f0-9]{32}$/);
    assert.match(await download('report-analysis-results_csv'), /exploratory_component_level/);
    await waitFor("document.getElementById('report-analysis-volcano')?.data?.some(t => t.mode === 'markers' && t.x.length)");
  }
  console.log('PASS: three real background analyses, confirmation, same-session navigation, CSV/provenance');
  await checkControlContrast();

  await evaluate("Plotly.relayout('report-analysis-cells_plot', {dragmode:'select'}); document.getElementById('report-analysis-cells_plot').scrollIntoView({block:'center'})");
  const selectionBox = await evaluate(`(() => {
    const r=document.querySelector('#report-analysis-cells_plot .nsewdrag').getBoundingClientRect();
    return {x:r.x+2,y:r.y+2,width:r.width-4,height:r.height-4};
  })()`);
  await call('Input.dispatchMouseEvent', {type:'mousePressed',x:selectionBox.x,y:selectionBox.y,button:'left',clickCount:1});
  await call('Input.dispatchMouseEvent', {type:'mouseMoved',x:selectionBox.x+selectionBox.width,y:selectionBox.y+selectionBox.height,button:'left',buttons:1});
  await call('Input.dispatchMouseEvent', {type:'mouseReleased',x:selectionBox.x+selectionBox.width,y:selectionBox.y+selectionBox.height,button:'left',clickCount:1});
  await waitFor("Shiny.shinyapp.$inputValues['report-analysis-cells']?.length > 200");
  console.log('PASS: mouse-driven embedding box selection reaches Shiny');

  const selected = await evaluate(`(() => {
    const plot = document.getElementById('report-analysis-cells_plot');
    const points = plot.data.flatMap((trace, curveNumber) =>
      trace.key.slice(0, 100).map((key, pointNumber) => ({
        curveNumber, pointNumber, x: trace.x[pointNumber], y: trace.y[pointNumber]
      })));
    plot.emit('plotly_selected', {points});
    return points.map(point => plot.data[point.curveNumber].key[point.pointNumber]);
  })()`);
  await waitFor(`Shiny.shinyapp.$inputValues['report-analysis-cells']?.length === ${selected.length}`);
  await run();
  const subset = JSON.parse(await download('report-analysis-provenance_json'));
  assert.equal(subset.request.group_a_cells.length, 100);
  assert.equal(subset.request.group_b_cells.length, 100);
  assert.deepEqual([...subset.request.group_a_cells, ...subset.request.group_b_cells].sort(), selected.sort());
  await screenshot('hybrid-analysis');
  console.log('PASS: embedding-selected cells define a new, recorded cohort');

  await evaluate("document.getElementById('report-analysis-run').click()");
  await waitFor("document.getElementById('report-analysis-run').disabled");
  await evaluate("document.getElementById('report-analysis-cancel').click()");
  await waitFor("document.getElementById('report-analysis-status').textContent.includes('cancelled')");
  await run();
  console.log('PASS: cancellation releases the worker and another run succeeds');

  await tab('Cell 3D');
  await waitFor("document.getElementById('report-graph-plot')?.data?.some(t => t.type === 'scatter3d' && t.x.length)");
  await waitFor("Shiny.shinyapp.$inputValues['report-graph-markers']?.length > 0 && !!document.getElementById('report-graph-markers').selectize.options.CD45");
  await select('report-graph-markers', ['CD4', 'CD45']);
  await waitFor("document.querySelector('#report-graph-plot canvas') !== null");
  await waitFor("document.getElementById('report-graph-plot').data.some(t => t.name === 'CD45')");
  const graph = JSON.parse(await download('report-graph-provenance_json'));
  assert.equal(graph.schema_version, 2);
  assert.match(graph.source_pxl_md5, /^[a-f0-9]{32}$/);
  assert.equal(await evaluate("document.getElementById('report-graph-edges') === null"), true);
  const nodeCsv = (await download('report-graph-nodes_csv')).trim().split('\n');
  assert.equal(nodeCsv[0], '"cell_id","node_id","x","y","z","markers"');
  const nodes = nodeCsv.length - 1;
  assert.equal(await evaluate("document.getElementById('report-graph-plot').data.filter(t => t.mode === 'markers').reduce((n,t) => n + t.x.length,0)"), nodes);
  assert.equal(await evaluate("document.getElementById('report-graph-plot').data.every(t => t.mode === 'markers')"), true);
  await checkControlContrast();
  await evaluate("document.getElementById('report-graph-plot').closest('.card').querySelector('.bslib-full-screen-enter').click()");
  await waitFor("document.getElementById('report-graph-plot').closest('.card').dataset.fullScreen === 'true'");
  await settled();
  // Shiny can be idle before Plotly finishes resizing its WebGL interaction area.
  await waitFor(`(() => {
    const plot = document.getElementById('report-graph-plot');
    const canvas = plot.querySelector('canvas');
    if (!canvas || Math.abs(plot._fullLayout.width - plot.clientWidth) > 1 ||
        Math.abs(plot._fullLayout.height - plot.clientHeight) > 1) return false;
    const r = canvas.getBoundingClientRect();
    return r.width > 0 && r.height > 0 &&
      document.elementFromPoint(r.x + r.width * .45, r.y + r.height * .5) === canvas;
  })()`);
  await evaluate("new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)))");
  const cameraBefore = await evaluate("JSON.stringify(document.getElementById('report-graph-plot')._fullLayout.scene.camera)");
  const graphBox = await evaluate(`(() => {
    const r=document.querySelector('#report-graph-plot canvas').getBoundingClientRect();
    return {x:r.x+r.width*.45,y:r.y+r.height*.5};
  })()`);
  await call('Input.dispatchMouseEvent', {type:'mouseMoved',...graphBox});
  await call('Input.dispatchMouseEvent', {type:'mousePressed',...graphBox,button:'left',buttons:1,clickCount:1});
  await delay(100);
  for (let step = 1; step <= 8; step++) {
    await call('Input.dispatchMouseEvent', {type:'mouseMoved',x:graphBox.x+90*step/8,y:graphBox.y+35*step/8,button:'left',buttons:1});
    await delay(75);
  }
  await call('Input.dispatchMouseEvent', {type:'mouseReleased',x:graphBox.x+90,y:graphBox.y+35,button:'left',clickCount:1});
  await waitFor(`JSON.stringify(document.getElementById('report-graph-plot')._fullLayout.scene.camera) !== ${JSON.stringify(cameraBefore)}`);
  const rotatedCamera = await evaluate("JSON.stringify(document.getElementById('report-graph-plot')._fullLayout.scene.camera)");
  await call('Input.dispatchMouseEvent', {type:'mouseWheel',...graphBox,deltaX:0,deltaY:-120});
  await waitFor(`JSON.stringify(document.getElementById('report-graph-plot')._fullLayout.scene.camera) !== ${JSON.stringify(rotatedCamera)}`);
  await screenshot('hybrid-cell-3d');
  await evaluate("document.querySelector('.bslib-full-screen-exit').click()");
  console.log('PASS: prepared 3D graph, all nodes, coloring, exports, full-screen, mouse rotation and zoom');

  await call('Emulation.setDeviceMetricsOverride', { width: 390, height: 844, deviceScaleFactor: 1, mobile: true });
  await tab('Custom analysis');
  assert.ok(await evaluate('document.documentElement.scrollWidth <= window.innerWidth + 2'));
  await screenshot('hybrid-mobile');
  console.log('PASS: mobile navigation/layout without horizontal page overflow');
} finally {
  socket.close();
}
