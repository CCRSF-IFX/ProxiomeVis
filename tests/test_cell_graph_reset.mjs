import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { runInNewContext } from 'node:vm';

const elements = new Map();
const id = 'report-graph-plot';
const controls = ['_reset', '_png', '_gif'].map(suffix => id + suffix);
const downloads = ['model_glb', 'nodes_csv', 'provenance_json'].map(suffix => id.replace(/plot$/, suffix));
for (const name of [...controls, ...downloads, id + '_gif_status']) {
  const classes = new Set();
  elements.set(name, {
    onclick: () => assert.fail('A stale export callback ran'),
    classList: { add: value => classes.add(value), remove: value => classes.delete(value),
      contains: value => classes.has(value) },
    setAttribute(name, value) { this[name] = value; },
  });
}
let cancelled = 0, disposed = 0, removed = 0;
const el = {
  id, proxiomeCell: 'old-cell', proxiomeCamera: {}, proxiomeInitialized: true,
  proxiomeGifCancel() { cancelled++; },
  proxiomeDispose() { disposed++; },
  querySelectorAll() { return [{ remove() { removed++; } }]; },
};
elements.set(id, el);
const graph = {
  idle: false, controls: { addEventListener() {} },
  pointgroup: {}, linegroup: {}, renderer: {},
};
let reset;
const window = { proxiomeCellGraphGif() { elements.get(id + '_gif').disabled = false; } };
runInNewContext(readFileSync(new URL('../www/cell-graph.js', import.meta.url), 'utf8'), {
  window,
  Shiny: { addCustomMessageHandler(name, callback) {
    assert.equal(name, 'proxiome-cell-graph-reset'); reset = callback;
  } },
  HTMLWidgets: { getInstance: () => ({ widget: graph }) },
  document: { getElementById: name => elements.get(name) },
});
assert.doesNotThrow(() => reset({ id: 'absent' }));
reset({ id });
assert.equal(el.proxiomeCell, null);
assert.equal(el.proxiomeCamera, null);
assert.equal(graph.idle, true);
assert.equal(cancelled, 1);
assert.equal(disposed, 1);
assert.equal(removed, 1);
for (const name of controls) {
  assert.equal(elements.get(name).disabled, true);
  assert.equal(elements.get(name).onclick, null);
}
for (const name of downloads) {
  assert.equal(elements.get(name).inert, true);
  assert.ok(elements.get(name).classList.contains('disabled'));
  assert.equal(elements.get(name)['aria-disabled'], 'true');
}
assert.match(elements.get(id + '_gif_status').textContent, /Display a cell/);
window.proxiomeCellGraph(el, { main: 'new-cell' });
assert.equal(el.proxiomeCell, 'new-cell');
for (const name of controls) assert.equal(elements.get(name).disabled, false);
for (const name of downloads) {
  assert.equal(elements.get(name).inert, false);
  assert.ok(!elements.get(name).classList.contains('disabled'));
  assert.equal(elements.get(name)['aria-disabled'], 'false');
}
console.log('PASS: dataset reset cancels GIF, disposes the old view, blocks exports, and permits a new display');
