// Run against the marker demo and Chromium with remote debugging enabled.
import assert from 'node:assert/strict';
import { writeFile } from 'node:fs/promises';

const debugUrl = process.env.PROXIOME_CDP_URL || 'http://127.0.0.1:9249';
const appUrl = process.env.PROXIOME_SMOKE_URL || 'http://127.0.0.1:3849';
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
  const result = await call('Runtime.evaluate', { expression, returnByValue: true, awaitPromise: true });
  if (result.exceptionDetails) throw new Error(JSON.stringify(result.exceptionDetails));
  return result.result.value;
};
const delay = ms => new Promise(resolve => setTimeout(resolve, ms));
async function waitFor(expression, timeoutMs = 60000) {
  const deadline = Date.now() + timeoutMs;
  while (Date.now() < deadline) {
    if (await evaluate(expression)) return;
    await delay(200);
  }
  throw new Error(`Condition not reached: ${expression}`);
}
async function tab(label) {
  await evaluate(`[...document.querySelectorAll('a[data-bs-toggle=tab]')].find(x => x.innerText === ${JSON.stringify(label)}).click()`);
  await delay(300);
}
async function select(id, value) {
  await waitFor(`document.getElementById(${JSON.stringify(`report-${id}`)}) !== null`);
  await evaluate(`(() => {
    const element = document.getElementById(${JSON.stringify(`report-${id}`)});
    const radio = [...element.querySelectorAll('input[type=radio]')].find(x => x.value === ${JSON.stringify(value)});
    if (radio) radio.click();
    else if (element.selectize) element.selectize.setValue(${JSON.stringify(value)});
    else { element.value = ${JSON.stringify(value)}; element.dispatchEvent(new Event('change', {bubbles:true})); }
  })()`);
  await waitFor(`JSON.stringify(Shiny.shinyapp.$inputValues[${JSON.stringify(`report-${id}`)}]) === JSON.stringify(${JSON.stringify(value)})`);
  await delay(300);
}
async function noErrors() {
  await waitFor("!document.documentElement.classList.contains('shiny-busy')");
  assert.deepEqual(await evaluate("[...document.querySelectorAll('.shiny-output-error:not(.shiny-output-error-validation)')].map(x => x.textContent)"), []);
}

async function compactHistory() {
  const gap = await evaluate(`(() => {
    const body = document.getElementById('report-analysis-history_id').closest('.card-body');
    return body.querySelector('p').getBoundingClientRect().top -
      body.closest('.card').querySelector('.card-header').getBoundingClientRect().bottom;
  })()`);
  assert.ok(gap >= 0 && gap < 32, `History content is vertically displaced by ${gap}px`);
}

try {
  await call('Emulation.setDeviceMetricsOverride', { width: 1700, height: 1100, deviceScaleFactor: 1, mobile: false });
  await call('Page.navigate', { url: 'about:blank' });
  await waitFor("location.href === 'about:blank'");
  await call('Page.navigate', { url: appUrl });
  await waitFor("document.getElementById('report-metric_cells')?.textContent.trim().length > 0");
  await tab('Results');
  await select('readout', 'colocalization');
  await select('donor', 'd1');
  await select('celltype', 'B');
  await tab('Spatial patterns');
  await waitFor("document.getElementById('report-spatial_report_plot')?._fullLayout?.xaxis?.ticktext?.length > 0");
  const automaticOrder = await evaluate("document.getElementById('report-spatial_report_plot')._fullLayout.xaxis.ticktext");
  assert.ok(automaticOrder.length >= 30, 'Exercise the full prepared union, not just a small marker subset');
  assert.equal(await evaluate("document.querySelector('.spatial-protein-controls .accordion-button').getAttribute('aria-expanded')"), 'false');
  await select('spatial_layout', 'compare');
  for (const name of ['reference', 'target', 'difference']) {
    await waitFor(`document.getElementById('report-spatial_${name}')?._fullLayout?.xaxis?.ticktext?.length === ${automaticOrder.length}`);
    assert.deepEqual(await evaluate(`document.getElementById('report-spatial_${name}')._fullLayout.xaxis.ticktext`), automaticOrder);
    const pitch = await evaluate(`document.getElementById('report-spatial_${name}')._fullLayout.xaxis._length / ${automaticOrder.length}`);
    assert.ok(pitch >= 14, `Crowded ${name} matrix: ${pitch} pixels per protein`);
  }
  assert.equal(await evaluate("getComputedStyle(document.querySelector('.spatial-comparison')).gridTemplateColumns.split(' ').length"), 1);
  await evaluate("document.querySelector('.spatial-comparison').scrollIntoView({block:'start'})");
  await delay(500);
  const automaticScreenshot = await call('Page.captureScreenshot', { format: 'png' });
  await writeFile('/tmp/proxiomevis-spatial-automatic.png', Buffer.from(automaticScreenshot.data, 'base64'));
  console.log(`PASS: ${automaticOrder.length}-protein automatic union, common order and readable matrix spacing`);
  await call('Emulation.setDeviceMetricsOverride', { width: 700, height: 1100, deviceScaleFactor: 1, mobile: true });
  await delay(700);
  assert.ok(await evaluate("document.documentElement.scrollWidth <= window.innerWidth"), 'Large matrices must scroll inside their cards');
  assert.ok(await evaluate("document.getElementById('report-spatial_reference').closest('.spatial-scroll').scrollWidth > document.getElementById('report-spatial_reference').closest('.spatial-scroll').clientWidth"));
  await call('Emulation.setDeviceMetricsOverride', { width: 1700, height: 1100, deviceScaleFactor: 1, mobile: false });
  await evaluate("document.querySelector('.spatial-protein-controls .accordion-button').click()");
  await select('spatial_reference_donor', 'd2');
  await select('spatial_reference_condition', 'unstim');
  assert.deepEqual(await evaluate("document.getElementById('report-spatial_reference')._fullLayout.xaxis.ticktext"), automaticOrder);
  await evaluate("document.getElementById('report-spatial_apply').click()");
  await waitFor("document.getElementById('report-spatial_selection_status').textContent.includes('Ward.D2 (d2 / unstim)')");
  await noErrors();
  await evaluate("document.getElementById('report-spatial_reset').click()");
  await waitFor("document.getElementById('report-spatial_selection_status').textContent.includes('Prepared reference order (d1 / unstim)')");
  await waitFor(`JSON.stringify(document.getElementById('report-spatial_reference')._fullLayout.xaxis.ticktext) === ${JSON.stringify(JSON.stringify(automaticOrder))}`);
  console.log('PASS: large-matrix mobile scrolling, draft reference change, Ward ordering and prepared-order reset');
  await select('spatial_set', 'manual');
  await select('spatial_markers', ['CD19', 'CD37', 'CD53', 'CD81', 'CD82']);
  await select('spatial_order', 'custom');
  await select('spatial_custom_order', 'CD82, CD81, CD53, CD37, CD19');
  await evaluate("document.getElementById('report-spatial_apply').click()");
  await waitFor("document.getElementById('report-spatial_selection_status').textContent.startsWith('5 proteins applied')");
  await evaluate("document.querySelector('.spatial-protein-controls .accordion-button').click()");
  await select('spatial_layout', 'compare');
  const expected = ['CD82', 'CD81', 'CD53', 'CD37', 'CD19'];
  for (const name of ['reference', 'target', 'difference']) {
    const id = `report-spatial_${name}`;
    await waitFor(`document.getElementById('${id}')?._fullLayout?.xaxis?.ticktext?.length === 5`);
    assert.deepEqual(await evaluate(`document.getElementById('${id}')._fullLayout.xaxis.ticktext`), expected);
    assert.deepEqual(await evaluate(`document.getElementById('${id}')._fullLayout.yaxis.ticktext`), [...expected].reverse());
    assert.equal(await evaluate(`document.getElementById('${id}')._context.scrollZoom`), true);
  }
  await noErrors();
  assert.ok(await evaluate("[...document.querySelectorAll('#report-spatial_difference .legendtext')].some(x=>x.textContent==='4 (capped)')"),
    'The significance-size scale must be visible, not just binary significance labels');
  assert.ok(await evaluate("[...document.querySelectorAll('#report-spatial_reference .legendtext')].some(x=>x.textContent==='100%')"),
    'The condition profile must explain detection-fraction dot sizes');
  const geometry = await evaluate("['reference','target','difference'].map(n=>{const p=document.getElementById('report-spatial_'+n)._fullLayout;return {width:p.xaxis._length,height:p.yaxis._length,top:p.yaxis._offset};})");
  assert.ok(geometry.every(p => Math.abs(p.width - geometry[0].width) < 1 &&
    Math.abs(p.height - geometry[0].height) < 1 && Math.abs(p.top - geometry[0].top) < 1),
    `Matrix plot areas differ: ${JSON.stringify(geometry)}`);
  const rects = await evaluate("['reference','target','difference'].map(n=>{const r=document.getElementById('report-spatial_'+n).closest('.card').getBoundingClientRect();return {x:r.x,y:r.y,w:r.width,h:r.height};})");
  assert.ok(rects[0].x < rects[1].x && rects[1].x < rects[2].x);
  assert.ok(rects.every(r => Math.abs(r.y - rects[0].y) < 2));
  await evaluate("document.querySelector('.spatial-comparison').scrollIntoView({block:'start'})");
  await delay(1200);
  const screenshot = await call('Page.captureScreenshot', { format: 'png' });
  await writeFile('/tmp/proxiomevis-redesign-spatial.png', Buffer.from(screenshot.data, 'base64'));
  await waitFor("document.getElementById('report-download_spatial_compare_csv')?.href?.includes('download')");
  const csv = await evaluate("fetch(document.getElementById('report-download_spatial_compare_csv').href).then(r=>r.text())");
  assert.ok(csv.includes('reference') && csv.includes('target') && csv.includes('difference'));
  assert.equal(csv.trim().split('\n').length, 76);
  console.log('PASS: aligned comparison panels, explicit order, zoom, CSV, desktop layout');

  await select('spatial_layout', 'single');
  await select('spatial_kind', 'profiles');
  await waitFor("document.getElementById('report-spatial_report_plot')._fullLayout.xaxis.ticktext.length === 5");
  assert.deepEqual(await evaluate("document.getElementById('report-spatial_report_plot')._fullLayout.xaxis.ticktext"), expected);
  await call('Emulation.setDeviceMetricsOverride', { width: 700, height: 1100, deviceScaleFactor: 1, mobile: true });
  await select('spatial_layout', 'compare');
  await delay(700);
  assert.equal(await evaluate("getComputedStyle(document.querySelector('.spatial-comparison')).gridTemplateColumns.split(' ').length"), 1);
  await call('Emulation.setDeviceMetricsOverride', { width: 1700, height: 1100, deviceScaleFactor: 1, mobile: false });
  console.log('PASS: single/profile order retained and responsive stacked layout');
  await select('spatial_view', 'network');
  await waitFor("document.getElementById('report-network_plot')?.data?.some(t=>t.key?.length)");
  await noErrors();
  await select('spatial_view', 'matrix');
  console.log('PASS: network view remains available alongside aligned matrices');

  await tab('New comparison');
  await waitFor("document.getElementById('report-analysis-readiness')?.textContent.includes('cells')");
  await compactHistory();
  await select('analysis-donor', 'd1');
  await select('analysis-celltype', 'B');
  await select('analysis-readout', 'abundance');
  await waitFor("document.getElementById('report-analysis-readiness').textContent.includes('Scope is ready')");
  await evaluate("document.getElementById('report-analysis-acknowledge').click(); document.getElementById('report-analysis-run').click()");
  await waitFor("document.getElementById('report-analysis-run').disabled");
  await waitFor("document.getElementById('report-analysis-status')?.textContent.startsWith('Completed')");
  await compactHistory();
  await evaluate("document.getElementById('report-analysis-open_result').click()");
  await waitFor("document.getElementById('report-result_source').value.startsWith('session-')");
  await waitFor("document.getElementById('report-readout').value === 'abundance'");
  await tab('Volcano');
  await waitFor("document.getElementById('report-report_volcano_plot')?.data?.some(t=>t.mode==='markers' && t.x.length)");
  await waitFor("document.getElementById('report-volcano_distribution')?.data?.some(t=>t.type==='box' && t.y.length)");
  await noErrors();
  assert.ok(await evaluate("document.getElementById('report-result_summary').textContent.includes('Session comparison')"));
  await select('result_source', 'prepared');
  await waitFor("document.getElementById('report-result_summary').textContent.includes('Prepared single-donor')");
  console.log('PASS: real background abundance comparison, shared Results viewer, return to prepared results');

  let completedComparisons = 1;
  for (const readout of ['self_proximity', 'colocalization']) {
    await tab('New comparison');
    await select('analysis-readout', readout);
    await waitFor("document.getElementById('report-analysis-readiness').textContent.includes('Scope is ready')");
    const previousId = await evaluate("document.getElementById('report-analysis-history_id').value");
    await evaluate("document.getElementById('report-analysis-run').click()");
    await waitFor("document.getElementById('report-analysis-run').disabled");
    await waitFor(`(document.getElementById('report-analysis-history_id').value !== ${JSON.stringify(previousId)} && document.getElementById('report-analysis-status').textContent.startsWith('Completed')) || document.getElementById('report-analysis-status').textContent.includes('exceeded the time limit')`, 130000);
    const status = await evaluate("document.getElementById('report-analysis-status').textContent");
    if (status.includes('exceeded the time limit')) {
      assert.equal(await evaluate("document.getElementById('report-analysis-history_id').value"), previousId);
      await waitFor("!document.getElementById('report-analysis-run').disabled");
      console.log(`LIMIT: ${readout} reached the existing 120-second worker limit; previous results preserved and Run re-enabled`);
      continue;
    }
    completedComparisons++;
    await evaluate("document.getElementById('report-analysis-open_result').click()");
    await waitFor(`document.getElementById('report-readout').value === '${readout}'`);
    await tab('Spatial patterns');
    for (const name of ['reference', 'target', 'difference']) {
      await waitFor(`document.getElementById('report-spatial_${name}')?._fullLayout?.xaxis?.ticktext?.length === 5`);
      assert.deepEqual(await evaluate(`document.getElementById('report-spatial_${name}')._fullLayout.xaxis.ticktext`), expected);
    }
    await noErrors();
    if (readout === 'colocalization') {
      await select('spatial_view', 'network');
      await waitFor("document.getElementById('report-network_markers').closest('.form-group').offsetParent !== null");
      await noErrors();
      await select('spatial_view', 'matrix');
    }
    await tab('Volcano');
    await waitFor("document.getElementById('report-report_volcano_plot')?.data?.some(t=>t.mode==='markers' && t.x.length)");
    await noErrors();
  }
  await tab('New comparison');
  assert.equal(await evaluate("Object.keys(document.getElementById('report-analysis-history_id').selectize.options).length"), completedComparisons);
  console.log(`PASS: spatial worker outcomes, shared spatial order, ${completedComparisons}-result history`);
} finally {
  socket.close();
}
