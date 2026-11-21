// Start the app with comparison_view_fixture() and allow local dataset loading.
// Set PROXIOME_INDEPENDENT_FIXTURE to a bundle made with comparison_view_fixture(TRUE).
import assert from 'node:assert/strict';

const tabs = await (await fetch(process.env.PROXIOME_CDP_URL || 'http://127.0.0.1:9257/json')).json();
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
    await delay(150);
  }
  throw new Error(`Condition not reached: ${expression}`);
}
async function select(id, value) {
  await waitFor(`document.getElementById(${JSON.stringify(id)}) !== null`);
  await evaluate(`(() => {
    const e = document.getElementById(${JSON.stringify(id)});
    if (e.selectize) e.selectize.setValue(${JSON.stringify(value)});
    else { e.value = ${JSON.stringify(value)}; e.dispatchEvent(new Event('change', {bubbles:true})); }
  })()`);
  await waitFor(`JSON.stringify(Shiny.shinyapp.$inputValues[${JSON.stringify(id)}]) === ${JSON.stringify(JSON.stringify(value))}`);
  await delay(250);
}
async function tab(label) {
  await evaluate(`[...document.querySelectorAll('a[data-bs-toggle=tab]')].find(x=>x.textContent.trim()===${JSON.stringify(label)}).click()`);
  await delay(200);
}
async function run(scope, axis, paired = true) {
  await tab('New comparison');
  await select('report-analysis-analysis_scope', scope);
  await select('report-analysis-comparison_axis', axis);
  await select('report-analysis-readout', 'abundance');
  if (scope === 'across_donors') {
    await select('report-analysis-donors', ['D1', 'D2', 'D3', 'D4']);
    await select('report-analysis-design', paired ? 'paired' : 'independent');
  } else await select('report-analysis-donor', 'D1');
  await select(`report-analysis-${axis === 'condition' ? 'celltype' : 'condition'}`, axis === 'condition' ? 'B' : 'NC');
  await waitFor(`document.getElementById('report-analysis-target').selectize.options[${JSON.stringify(axis === 'condition' ? 'TRT' : 'T')}] !== undefined`);
  await select('report-analysis-reference', axis === 'condition' ? 'NC' : 'B');
  await select('report-analysis-target', axis === 'condition' ? 'TRT' : 'T');
  if (scope === 'within_donor') await evaluate("document.getElementById('report-analysis-acknowledge').checked || document.getElementById('report-analysis-acknowledge').click()");
  await waitFor("document.getElementById('report-analysis-readiness').textContent.includes('Scope is ready')");
  await evaluate("document.getElementById('report-analysis-run').click()");
  await waitFor("document.getElementById('report-analysis-status').textContent.includes('Running in a separate')");
  await waitFor("document.getElementById('report-analysis-status').textContent.includes('Completed.')");
  await evaluate("document.getElementById('report-analysis-open_result').click()");
  await waitFor("Shiny.shinyapp.$inputValues['report-report_page']==='Results'");
  await waitFor(`document.getElementById('report-result_summary').textContent.includes(${JSON.stringify(scope === 'across_donors' ? 'Across donors' : 'Donor: D1')})`);
  await tab('Volcano');
  await waitFor("document.getElementById('report-volcano_distribution')?._fullData?.length > 0 && !document.documentElement.classList.contains('shiny-busy')");
  const errors = await evaluate("[...document.querySelectorAll('.shiny-output-error:not(.shiny-output-error-validation)')].map(x=>x.textContent)");
  assert.deepEqual(errors, []);
  const notice = await evaluate("document.getElementById('report-volcano_notice').textContent");
  assert.ok(notice.includes(scope === 'across_donors' ? 'Across-donor comparison' : 'Exploratory component-level test'));
  if (scope === 'across_donors') {
    const hasLines = await evaluate("document.getElementById('report-volcano_distribution')._fullData.some(x=>x.mode?.includes('lines'))");
    assert.equal(hasLines, paired, 'Only paired donor measurements have connecting lines');
  }
  console.log(`PASS: ${scope}, ${axis}, ${paired ? 'paired' : 'independent'}; real worker, shared viewer and plots`);
}

try {
  await call('Emulation.setDeviceMetricsOverride', { width: 1600, height: 1050, deviceScaleFactor: 1, mobile: false });
  await call('Page.navigate', { url: process.env.PROXIOME_SMOKE_URL || 'http://127.0.0.1:3857' });
  await waitFor("document.getElementById('report-metric_cells')?.textContent.trim().length > 0");
  await run('within_donor', 'condition');
  await run('within_donor', 'celltype');
  await run('across_donors', 'condition');
  await run('across_donors', 'celltype');
  const independentPath = process.env.PROXIOME_INDEPENDENT_FIXTURE;
  assert.ok(independentPath, 'Set PROXIOME_INDEPENDENT_FIXTURE to the independent-donor bundle');
  await evaluate("document.getElementById('dataset').click()");
  await select('rds_path', independentPath);
  await evaluate("document.getElementById('load_rds').click()");
  await waitFor(`document.getElementById('dataset_details').textContent.includes(${JSON.stringify(independentPath)})`);
  await evaluate("document.querySelector('#shiny-modal .modal-footer button').click()");
  await run('across_donors', 'condition', false);
  console.log('PASS: dataset reset, all five comparison scenarios');
} finally {
  socket.close();
}
