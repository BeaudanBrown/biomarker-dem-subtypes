function tdp43Controls(el, x, options) {
  "use strict";
  // Display-only controls: every point, fit and summary was computed in R.
  const root = el.closest(".tdp43-widget");
  const cohorts = Array.from(root.querySelectorAll("[data-tdp-cohort]"));
  const panel = root.querySelector("[data-tdp-panel-select]");
  const points = root.querySelector("[data-tdp-points]");
  const trend = document.getElementById("tdp43-trend-toggle");
  const counts = root.querySelector("[data-tdp-counts]");
  const tables = Array.from(root.querySelectorAll("[data-tdp-table-panel]"));

  function update() {
    const selected = new Set(cohorts.filter(input => input.checked)
      .map(input => input.getAttribute("data-tdp-cohort")));
    const selectedPanel = panel ? panel.value : "Overall";
    const visible = el.data.map(trace => {
      const meta = trace.meta;
      if (options.kind === "scatter") {
        return selected.has(meta.panel) &&
          (meta.kind === "points" || Boolean(trend && trend.checked));
      }
      return meta.panel === selectedPanel &&
        (meta.kind !== "points" || !points || points.checked);
    });
    Plotly.restyle(el, { visible: visible });
    tables.forEach(table => {
      table.style.display = table.getAttribute("data-tdp-table-panel") === selectedPanel ? "" : "none";
    });
    if (counts) {
      // All is an independent pooled layer, not a select-all command.
      // Count each participant once even when pooled and subgroup layers overlap.
      const entries = Object.entries(options.counts || {}).filter(([name]) =>
        selected.has("All") || selected.has(name));
      const total = entries.reduce((sum, [, n]) => sum + Number(n), 0);
      counts.textContent = "Visible total: n=" + total +
        (entries.length ? "; " + entries.map(([name, n]) => name + ": n=" + n).join("; ") : "");
    }
  }

  function onChange() {
    update();
  }

  // Avoid duplicate handlers if htmlwidgets renders the same element again.
  if (el.tdp43Dispose) el.tdp43Dispose();
  root.addEventListener("change", onChange);
  if (options.kind === "scatter" && trend) trend.addEventListener("change", update);
  el.tdp43Dispose = function () {
    root.removeEventListener("change", onChange);
    if (trend) trend.removeEventListener("change", update);
  };
  if (options.kind === "box") {
    // ggplotly may omit boxpoints even when ggplot outliers are suppressed.
    // Plotly then draws default outlier dots inside the box trace, independently
    // of the explicit scatter layer. Keep that layer as the only source of dots.
    const boxIndices = el.data.flatMap((trace, index) => trace.type === "box" ? [index] : []);
    if (boxIndices.length) Plotly.restyle(el, { boxpoints: false }, boxIndices);
  }
  update();
}
