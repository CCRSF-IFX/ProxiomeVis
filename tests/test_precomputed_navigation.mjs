import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { runInNewContext } from 'node:vm';

const attributes = new Map();
const classes = new Set();
const tab = {
  classList: { toggle(name, enabled) { enabled ? classes.add(name) : classes.delete(name); } },
  setAttribute(name, value) { attributes.set(name, value); },
  removeAttribute(name) { attributes.delete(name); },
};
const handlers = new Map();
const sizes = new Map();
let resized = 0;
runInNewContext(readFileSync(new URL('../www/precomputed.js', import.meta.url), 'utf8'), {
  $: ready => ready(),
  Event: class { constructor(type) { this.type = type; } },
  window: { dispatchEvent(event) { assert.equal(event.type, 'resize'); resized++; } },
  Shiny: { addCustomMessageHandler(name, callback) {
    handlers.set(name, callback);
  } },
  document: { getElementById: id => id === 'report-results_view' ? {
    querySelector(selector) { assert.equal(selector, 'a[data-value="Spatial patterns"]'); return tab; },
  } : id === 'matrix' ? { style: { setProperty(name, value) { sizes.set(name, value); } } } : null },
});
const handler = handlers.get('precomputed-spatial-availability');
handler({ id: 'report-results_view', disabled: true });
assert.ok(classes.has('disabled'));
assert.equal(attributes.get('aria-disabled'), 'true');
assert.equal(attributes.get('tabindex'), '-1');
assert.match(tab.title, /Self-proximity or Colocalization/);
handler({ id: 'report-results_view', disabled: false });
assert.ok(!classes.has('disabled'));
assert.equal(attributes.get('aria-disabled'), 'false');
assert.ok(!attributes.has('tabindex'));
assert.doesNotThrow(() => handler({ id: 'missing-module', disabled: true }));
handlers.get('precomputed-spatial-size')({ id: 'matrix', proteins: 50 });
assert.equal(sizes.get('--spatial-panel-min'), '1060px');
handlers.get('precomputed-spatial-size')({ id: 'matrix', proteins: 5 });
assert.equal(sizes.get('--spatial-panel-min'), '480px');
assert.equal(resized, 2);
console.log('PASS: spatial tab disabled/enabled accessibility state');
