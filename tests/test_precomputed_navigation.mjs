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
let handler;
runInNewContext(readFileSync(new URL('../www/precomputed.js', import.meta.url), 'utf8'), {
  $: ready => ready(),
  Shiny: { addCustomMessageHandler(name, callback) {
    assert.equal(name, 'precomputed-spatial-availability');
    handler = callback;
  } },
  document: { getElementById: id => id === 'report-results_view' ? {
    querySelector(selector) { assert.equal(selector, 'a[data-value="Spatial patterns"]'); return tab; },
  } : null },
});
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
console.log('PASS: spatial tab disabled/enabled accessibility state');
