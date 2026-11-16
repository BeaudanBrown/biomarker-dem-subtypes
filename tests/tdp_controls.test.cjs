// Synthetic UI-state tests only. No R, report artifacts or participant data.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const test = require('node:test');
const source = fs.readFileSync('R/tdp_controls.js', 'utf8');

class Control {
  constructor(attributes = {}, checked = false) {
    this.attributes = attributes;
    this.checked = checked;
    this.value = 'Overall';
    this.style = {};
    this.handlers = new Map();
  }
  getAttribute(name) { return this.attributes[name]; }
  addEventListener(name, fn) {
    if (!this.handlers.has(name)) this.handlers.set(name, new Set());
    this.handlers.get(name).add(fn);
  }
  removeEventListener(name, fn) { this.handlers.get(name)?.delete(fn); }
  emit(name, target = this) {
    for (const fn of this.handlers.get(name) || []) fn({ target });
  }
}

function fixture(kind, trend = new Control()) {
  const root = new Control();
  const cohorts = ['Controls', 'AD', 'FTD', 'LBD/PD', 'All'].map(name =>
    new Control({ 'data-tdp-cohort': name }, name !== 'All'));
  const all = cohorts[4];
  const panel = new Control();
  const points = new Control({}, true);
  const counts = new Control();
  const tables = ['Overall', 'AD'].map(name => new Control({ 'data-tdp-table-panel': name }));
  root.querySelectorAll = selector => ({
    '[data-tdp-cohort]': kind === 'scatter' ? cohorts : [],
    '[data-tdp-table-panel]': tables
  }[selector] || []);
  root.querySelector = selector => ({
    '[data-tdp-panel-select]': kind === 'scatter' ? null : panel,
    '[data-tdp-points]': kind === 'box' ? points : null,
    '[data-tdp-counts]': kind === 'scatter' ? counts : null
  }[selector]);
  const meta = kind === 'scatter'
    ? cohorts.flatMap(input => ['points', 'trend', 'trend'].map(role =>
      ({ panel: input.getAttribute('data-tdp-cohort'), kind: role })))
    : ['Overall', 'AD'].flatMap(name => (kind === 'box' ? ['box', 'points'] : ['matrix']).map(role =>
      ({ panel: name, kind: role })));
  // Match the reported conversion: box traces omit boxpoints (Plotly defaults
  // to outlier dots), while explicit point traces are scatter/markers.
  const el = { data: meta.map(meta => ({
    meta, type: meta.kind === 'box' ? 'box' : 'scatter',
    ...(meta.kind === 'points' ? { mode: 'markers' } : {})
  })), closest: () => root };
  let visibility;
  const context = vm.createContext({
    document: { getElementById: () => trend },
    Plotly: { restyle: (target, update, indices) => {
      assert.equal(target, el);
      const selected = indices || el.data.map((_, i) => i);
      for (const [key, value] of Object.entries(update)) {
        selected.forEach((index, i) => {
          el.data[index][key] = Array.isArray(value) ? value[i] : value;
        });
      }
      if (update.visible) visibility = Array.from(update.visible);
    } }
  });
  vm.runInContext(source, context);
  const attach = () => context.tdp43Controls(el, {}, {
    kind, counts: { Controls: 10, AD: 20, FTD: 30, 'LBD/PD': 40 }
  });
  attach();
  return { root, cohorts, all, panel, points, counts, tables, trend, attach, el,
    visible: () => visibility };
}

test('box traces cannot draw default outliers when explicit points are hidden', () => {
  const f = fixture('box');
  f.points.checked = false;
  f.root.emit('change', f.points);
  for (const trace of f.el.data.filter(trace => trace.type === 'box')) {
    assert.equal(trace.boxpoints, false, 'Native box outliers must be disabled');
  }
  assert.equal(f.el.data[0].visible, true, 'Keep the box visible');
  assert.equal(f.el.data[1].visible, false, 'Hide the explicit point layer');
  f.panel.value = 'AD';
  f.root.emit('change', f.panel);
  assert.equal(f.el.data[2].boxpoints, false);
  assert.equal(f.el.data[2].visible, true);
  f.points.checked = true;
  f.root.emit('change', f.points);
  assert.equal(f.el.data[3].visible, true, 'Restore only the explicit points');
  assert.equal(f.el.data[2].boxpoints, false, 'Do not duplicate outlier dots');
  assert.equal(f.el.data[3].boxpoints, undefined, 'Do not modify scatter traces');
});

test('All is an independent pooled layer with its own trends and unique participant counts', () => {
  const f = fixture('scatter');
  assert.deepEqual(f.visible(), [true, false, false, true, false, false, true, false, false, true, false, false, false, false, false]);
  assert.match(f.counts.textContent, /Visible total: n=100/);
  f.cohorts[1].checked = false;
  f.root.emit('change', f.cohorts[1]);
  assert.equal(f.all.checked, false);
  assert.match(f.counts.textContent, /Visible total: n=80/);
  f.all.checked = true;
  f.root.emit('change', f.all);
  assert.equal(f.cohorts[1].checked, false, 'All must not select AD');
  assert.deepEqual(f.visible().slice(12), [true, false, false]);
  assert.match(f.counts.textContent, /Visible total: n=100/);
  f.trend.checked = true;
  f.trend.emit('change');
  assert.deepEqual(f.visible().slice(3, 6), [false, false, false]);
  assert.deepEqual(f.visible().slice(12), [true, true, true]);
  f.cohorts.slice(0, 4).forEach(input => { input.checked = false; });
  f.root.emit('change', f.cohorts[0]);
  assert.deepEqual(f.visible().slice(0, 12), Array(12).fill(false));
  assert.deepEqual(f.visible().slice(12), [true, true, true]);
  assert.match(f.counts.textContent, /Visible total: n=100/);
  f.cohorts[1].checked = true;
  f.root.emit('change', f.cohorts[1]);
  assert.match(f.counts.textContent, /Visible total: n=100/, 'Do not double-count AD');
  f.all.checked = false;
  f.root.emit('change', f.all);
  assert.equal(f.cohorts[1].checked, true, 'All must not deselect AD');
  assert.match(f.counts.textContent, /Visible total: n=20/);
  f.cohorts[1].checked = false;
  f.root.emit('change', f.cohorts[1]);
  assert.ok(f.visible().every(value => !value));
  assert.match(f.counts.textContent, /Visible total: n=0/);
});

test('separate scatters share global trends but not cohort selection; rerender cleans handlers', () => {
  const trend = new Control();
  const a = fixture('scatter', trend);
  const b = fixture('scatter', trend);
  a.cohorts[0].checked = false;
  a.root.emit('change', a.cohorts[0]);
  trend.checked = true;
  trend.emit('change');
  assert.equal(a.visible()[0], false);
  assert.equal(b.visible()[0], true);
  assert.equal(b.visible()[1], true);
  a.attach();
  assert.equal(trend.handlers.get('change').size, 2);
  assert.equal(a.root.handlers.get('change').size, 1);
});

test('box selector and point toggle preserve boxes and switch summary tables', () => {
  const f = fixture('box');
  assert.deepEqual(f.visible(), [true, true, false, false]);
  f.points.checked = false;
  f.root.emit('change', f.points);
  assert.deepEqual(f.visible(), [true, false, false, false]);
  f.panel.value = 'AD';
  f.root.emit('change', f.panel);
  assert.deepEqual(f.visible(), [false, false, true, false]);
  assert.equal(f.tables[0].style.display, 'none');
  assert.equal(f.tables[1].style.display, '');
  f.points.checked = true;
  f.root.emit('change', f.points);
  assert.deepEqual(f.visible(), [false, false, true, true]);
});

test('matrix selector shows only the selected matrix and its pairwise n table', () => {
  const f = fixture('matrix');
  assert.deepEqual(f.visible(), [true, false]);
  f.panel.value = 'AD';
  f.root.emit('change', f.panel);
  assert.deepEqual(f.visible(), [false, true]);
  assert.equal(f.tables[0].style.display, 'none');
  assert.equal(f.tables[1].style.display, '');
});
