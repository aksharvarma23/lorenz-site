/* Lorenz Explorer (simple version): interface. Heavy maths runs in js/worker.js. */
(function () {
  "use strict";
  const L = window.Lorenz;
  const $ = (id) => document.getElementById(id);
  const SIGMA = 10, B = 8 / 3, T_RUN = 60;
  const TR = { pitchfork: 1, homoclinic: 13.9266, attracting: 24.059, hopf: 24.7368 };

  // ---------------------------------------------------------------- worker plumbing
  const worker = new Worker("js/worker.js");
  let reqId = 0;
  const pending = new Map(), latest = {};
  worker.onmessage = (ev) => {
    const { id, result, error } = ev.data;
    const h = pending.get(id);
    if (!h) return;
    pending.delete(id);
    if (error) h.reject(new Error(error)); else h.resolve(result);
  };
  function run(task, p) {
    const id = ++reqId;
    latest[task] = id;
    return new Promise((resolve, reject) => {
      pending.set(id, { resolve, reject });
      worker.postMessage({ id, task, p });
    }).then((res) => ({ res, stale: latest[task] !== id }));
  }

  // ---------------------------------------------------------------- plot styling
  const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  const dark = () => window.matchMedia("(prefers-color-scheme: dark)").matches;
  const axisT = (t) => ({ text: t, font: { color: css("--ink-2"), size: 12 } });
  function layout(extra) {
    const axis = { gridcolor: css("--rule"), zerolinecolor: css("--rule"), linecolor: css("--rule"), tickfont: { color: css("--ink-2"), size: 11 } };
    return Object.assign({
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
      font: { family: "Inter, system-ui, sans-serif", color: css("--ink"), size: 12 },
      margin: { l: 56, r: 14, t: 30, b: 44 }, showlegend: false, hovermode: "closest",
    }, extra, {
      xaxis: Object.assign({}, axis, extra && extra.xaxis),
      yaxis: Object.assign({}, axis, extra && extra.yaxis),
    });
  }
  const CFG = { displaylogo: false, responsive: true, modeBarButtonsToRemove: ["select2d", "lasso2d", "toImage"] };
  const title = (t) => ({ text: t, font: { size: 13, color: css("--ink-2") }, x: 0.02 });

  // ---------------------------------------------------------------- presets and inputs
  const PRESETS = [
    { name: "Rest", r: 0.5, x0: [1, 1, 1], T: 30 },
    { name: "Steady roll", r: 10, x0: [1, 1, 1], T: 30 },
    { name: "Transient chaos", r: 20, x0: [0.3, 18.6, 19.5], T: 60 },
    { name: "Bistable: near C+", r: 24.4, x0: [8.4, 8.4, 23.9], T: 100 },
    { name: "Bistable: far start", r: 24.4, x0: [1, 1, 1], T: 100 },
    { name: "Butterfly", r: 28, x0: [1, 1, 1], T: 60 },
    { name: "Periodic window", r: 100.5, x0: [1, 1, 1], T: 60 },
    { name: "Window, 2 peaks", r: 160, x0: [1, 1, 1], T: 60 },
    { name: "Single cycle", r: 350, x0: [1, 1, 1], T: 40 },
  ];
  let runT = T_RUN;
  const pres = $("presets");
  PRESETS.forEach((p, i) => {
    const bt = document.createElement("button");
    bt.className = "chip"; bt.type = "button";
    bt.innerHTML = `<b>r = ${p.r}</b>${p.name}`;
    bt.addEventListener("click", () => {
      $("r").value = p.r; $("x0").value = p.x0[0]; $("y0").value = p.x0[1]; $("z0").value = p.x0[2];
      runT = p.T; syncSlider(); changed(true);
    });
    pres.appendChild(bt);
  });
  function markPreset(s) {
    [...pres.children].forEach((c, i) => {
      const p = PRESETS[i];
      c.classList.toggle("on", Math.abs(p.r - s.r) < 1e-9 && p.x0.every((v, k) => Math.abs(v - s.x0[k]) < 1e-9));
    });
  }
  // r slider: linear 0..30 on the first 60%, logarithmic 30..400 beyond
  const rFromSlider = (v) => v <= 600 ? +(v / 20).toFixed(2) : +(30 * Math.pow(400 / 30, (v - 600) / 400)).toFixed(1);
  const sliderFromR = (r) => r <= 30 ? Math.max(0, r * 20) : 600 + 400 * Math.log(r / 30) / Math.log(400 / 30);
  function syncSlider() {
    const r = parseFloat($("r").value);
    if (Number.isFinite(r)) $("rRange").value = sliderFromR(Math.min(400, Math.max(0, r)));
  }
  $("rRange").addEventListener("input", () => { $("r").value = rFromSlider(+$("rRange").value); runT = T_RUN; changed(); });
  ["r", "x0", "y0", "z0"].forEach((id) => $(id).addEventListener("input", () => { syncSlider(); runT = T_RUN; changed(); }));

  function readState() {
    return { r: parseFloat($("r").value), x0: [parseFloat($("x0").value), parseFloat($("y0").value), parseFloat($("z0").value)] };
  }

  // ---------------------------------------------------------------- 1. fixed points
  const fmt = (v) => (v >= 0 ? "+" : "−") + Math.abs(v).toFixed(3);
  const fmtC = (e) => Math.abs(e.im) < 1e-10 ? fmt(e.re) : `${fmt(e.re)} ${e.im >= 0 ? "+" : "−"} ${Math.abs(e.im).toFixed(3)}i`;
  function renderFP(s) {
    const fps = L.fixedPoints(s.r, SIGMA, B);
    let h = "<thead><tr><th>Point</th><th>(x, y, z)</th><th>Eigenvalues</th><th>Stability</th></tr></thead><tbody>";
    for (const f of fps) {
      h += `<tr><td><b>${f.label}</b></td><td class="mono num">(${f.p.map(v => v.toFixed(2)).join(", ")})</td>
        <td class="mono num">${f.eig.map(fmtC).join("<br>")}</td>
        <td><span class="pill ${f.stable ? "s" : "u"}">${f.stable ? "stable" : "unstable"}</span><div class="verdict-text">${f.type}</div></td></tr>`;
    }
    $("fpTable").innerHTML = h + "</tbody>";
    return fps;
  }

  // ---------------------------------------------------------------- 2. regime bar
  function locatorSVG(r) {
    const W = 900, H = 150, x0 = 20, x1 = W - 20, tr = TR;
    const hi = Math.max(30, tr.hopf * 1.25, r * 1.08), lo = Math.min(0, r) - 1, linMax = 30;
    const X = (v) => {
      v = Math.max(lo, Math.min(hi, v));
      const split = hi > linMax ? 0.72 : 1;
      if (v <= linMax) return x0 + (x1 - x0) * split * (v - lo) / (Math.min(hi, linMax) - lo);
      return x0 + (x1 - x0) * (split + (1 - split) * Math.log(v / linMax) / Math.log(hi / linMax));
    };
    const bands = [[lo, 1, "--band-1", "rest"], [1, tr.homoclinic, "--band-2", "steady rolls"],
      [tr.homoclinic, tr.attracting, "--band-3", "transient chaos"], [tr.attracting, tr.hopf, "--band-4", "bistable"],
      [tr.hopf, hi, "--band-5", "chaos / windows"]];
    const by = 62, bh = 34;
    let g = "";
    for (const [a, b, col, lab] of bands) {
      const xa = X(a), xb = X(b);
      g += `<rect x="${xa}" y="${by}" width="${Math.max(0, xb - xa)}" height="${bh}" fill="var(${col})"/>`;
      if (xb - xa > 70) g += `<text x="${(xa + xb) / 2}" y="${by + bh / 2 + 4}" text-anchor="middle" font-size="12" fill="var(--ink-2)">${lab}</text>`;
    }
    const marks = [["pitchfork", "pitchfork", "middle", 0], ["homoclinic", "homoclinic", "middle", 0], ["attracting", "attracts", "end", -3], ["hopf", "Hopf", "start", 3]];
    for (const [k, lab, anchor, dx] of marks) {
      const xv = X(tr[k]), val = k === "pitchfork" ? "1" : tr[k].toFixed(2);
      g += `<line x1="${xv}" x2="${xv}" y1="${by - 18}" y2="${by + bh}" stroke="var(--ink-2)" stroke-width="1.2"/>`;
      g += `<text x="${xv + dx}" y="${by - 32}" text-anchor="${anchor}" font-size="11.5" fill="var(--ink-2)">${lab}</text>`;
      g += `<text x="${xv + dx}" y="${by - 21}" text-anchor="${anchor}" font-size="11.5" font-family="var(--mono)" fill="var(--ink)">${val}</text>`;
    }
    const ticks = [0, 5, 10, 15, 20, 25, 30, 50, 100, 200, 400, 1000].filter(t => t >= lo && t <= hi);
    g += `<line x1="${x0}" x2="${x1}" y1="${by + bh}" y2="${by + bh}" stroke="var(--rule)"/>`;
    for (const t of ticks) g += `<text x="${X(t)}" y="${by + bh + 16}" text-anchor="middle" font-size="11" fill="var(--ink-3)">${t}</text>`;
    const xr = X(r), right = xr > W * 0.8;
    g += `<path d="M${xr - 7},${by + bh + 34} L${xr + 7},${by + bh + 34} L${xr},${by + bh + 22} Z" fill="var(--hot)"/>`;
    g += `<line x1="${xr}" x2="${xr}" y1="${by - 4}" y2="${by + bh + 22}" stroke="var(--hot)" stroke-width="2.5"/>`;
    g += `<text x="${xr + (right ? -10 : 10)}" y="${by + bh + 32}" text-anchor="${right ? "end" : "start"}" font-size="12.5" font-weight="600" fill="var(--hot)">your r = ${r}</text>`;
    g += `<text x="${x1}" y="${H - 2}" text-anchor="end" font-size="10.5" fill="var(--ink-3)">linear up to 30, logarithmic beyond</text>`;
    return `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Regime bar with your value of r marked">${g}</svg>`;
  }
  function renderLocator(s) {
    $("locator").innerHTML = locatorSVG(s.r);
    const [name, text] = L.regime(s.r, TR);
    $("regimeText").innerHTML = `<b>${name[0].toUpperCase() + name.slice(1)}.</b> ${text}`;
  }

  // ---------------------------------------------------------------- 3. trajectory
  let lastTraj = null, lastFps = [];
  function renderTrajectory() {
    if (!lastTraj) return;
    const { t, X, verdict } = lastTraj;
    const step = Math.max(1, Math.floor(t.length / 12000));
    const idx = []; for (let i = 0; i < t.length; i += step) idx.push(i);
    const pick = (a) => idx.map(i => a[i]);
    const ax = { gridcolor: css("--rule"), zerolinecolor: css("--rule"), color: css("--ink-2"), backgroundcolor: "rgba(0,0,0,0)", showspikes: false };
    const scale = dark() ? [[0, "#2b4a60"], [1, "#a9dcfa"]] : [[0, "#b8d3e6"], [1, "#0d3550"]];
    Plotly.react("plot3d", [
      { type: "scatter3d", mode: "lines", x: pick(X[0]), y: pick(X[1]), z: pick(X[2]), line: { width: 3, color: pick(t), colorscale: scale }, hoverinfo: "skip" },
      { type: "scatter3d", mode: "markers", hoverinfo: "text", x: lastFps.map(f => f.p[0]), y: lastFps.map(f => f.p[1]), z: lastFps.map(f => f.p[2]),
        text: lastFps.map(f => `${f.label}: ${f.stable ? "stable" : "unstable"}`), marker: { size: 5, symbol: lastFps.map(f => f.stable ? "circle" : "x"), color: css("--ink") } },
      { type: "scatter3d", mode: "markers", x: [X[0][0]], y: [X[1][0]], z: [X[2][0]], marker: { size: 5, color: css("--ok") }, hoverinfo: "skip" },
    ], layout({
      margin: { l: 0, r: 0, t: 26, b: 0 }, title: title("3-D path · green = start · × unstable fixed point"),
      scene: { xaxis: Object.assign({ title: "x" }, ax), yaxis: Object.assign({ title: "y" }, ax), zaxis: Object.assign({ title: "z" }, ax),
        camera: { eye: { x: 1.5, y: -1.5, z: 0.8 } }, aspectmode: "cube" },
    }), CFG);
    Plotly.react("plotXt", [{ type: "scattergl", mode: "lines", x: pick(t), y: pick(X[0]), line: { width: 1.2, color: css("--accent") } }],
      layout({ title: title("x(t)"), xaxis: { title: axisT("time") }, yaxis: { title: axisT("x") } }), CFG);
    const v = $("verdict");
    v.className = "verdict " + verdict.kind;
    v.innerHTML = `<span class="dot"></span><span>Long-run behaviour: <b>${verdict.text}</b></span>`;
  }

  // ---------------------------------------------------------------- CSV download of the trajectory shown
  $("csvBtn").addEventListener("click", () => {
    if (!lastTraj) return;
    const { t, X, r } = lastTraj;
    const rows = ["t,x,y,z"];
    for (let i = 0; i < t.length; i++) rows.push(`${t[i].toFixed(4)},${X[0][i].toFixed(8)},${X[1][i].toFixed(8)},${X[2][i].toFixed(8)}`);
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([rows.join("\n")], { type: "text/csv" }));
    a.download = `lorenz_r${r}_x0_${lastTraj.x0.join("_")}.csv`;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
  });

  // ---------------------------------------------------------------- 4. sensitivity
  let lastSens = null;
  function renderSensitivity() {
    if (!lastSens) return;
    const d = lastSens;
    const t = Array.from(d.t), sep = Array.from(d.sep).map(v => Math.max(v, 1e-300));
    const traces = [{ type: "scattergl", mode: "lines", x: t, y: sep, line: { width: 1.3, color: css("--accent") } }];
    if (Number.isFinite(d.fit.slope)) {
      const xs = [t[d.fit.i0], t[d.fit.i1]];
      traces.push({ type: "scatter", mode: "lines", x: xs, y: xs.map(x => Math.exp(d.fit.icpt + d.fit.slope * x)), line: { dash: "dash", width: 2, color: css("--hot") } });
    }
    traces.push({ type: "scatter", mode: "lines", x: [0, 40], y: [d.L, d.L], line: { width: 1, color: css("--ink-3") } });
    Plotly.react("plotSep", traces, layout({
      title: title("Distance between two runs started 10⁻¹⁰ apart"),
      xaxis: { title: axisT("time") }, yaxis: { type: "log", exponentformat: "power", title: axisT("separation") },
      annotations: [{ x: 40, y: Math.log10(d.L), xref: "x", yref: "y", text: `size of the attractor ≈ ${d.L.toFixed(0)}`, showarrow: false, xanchor: "right", yanchor: "bottom", font: { size: 11, color: css("--ink-2") } }],
    }), CFG);
    const meaning = d.lam > 0.05 ? "positive: chaotic" : d.lam < -0.02 ? "negative: settles down" : "about zero: a repeating cycle";
    $("metrics").innerHTML = `
      <div class="metric"><div class="k">Largest Lyapunov exponent λ₁</div><div class="v">${d.lam.toFixed(3)}</div><div class="s">${meaning}${Math.abs(d.lam - 0.906) < 0.05 && d.lam > 0.5 ? " (published 0.906 at r = 28)" : ""}</div></div>
      <div class="metric"><div class="k">How far ahead it can be predicted</div><div class="v">${d.chaotic ? d.horizonEps.toFixed(1) : "∞"}</div><div class="s">${d.chaotic ? "time units, starting from an error of 0.001" : "nearby runs stay together"}</div></div>`;
    $("horizonText").innerHTML = d.chaotic
      ? `Each <b>10× more accurate</b> starting measurement adds only <b>${d.gain10.toFixed(1)}</b> time units of forecast (ln 10 / λ₁). That is the limit of prediction in a chaotic system.`
      : `λ₁ ≤ 0: small errors do not grow, so the motion is predictable.`;
  }

  // ---------------------------------------------------------------- update loop
  let timer = null;
  function changed(now) {
    markPreset(readState());
    clearTimeout(timer);
    timer = setTimeout(update, now ? 0 : 220);
  }
  async function update() {
    const s = readState();
    if (![s.r, ...s.x0].every(Number.isFinite) || s.r < -50 || s.r > 2000) { $("status").textContent = "Enter a number for r between 0 and 400."; return; }
    $("status").textContent = "Computing…";
    const t0 = performance.now();
    lastFps = renderFP(s);
    renderLocator(s);
    ["pc", "pd"].forEach(id => $(id).classList.add("busy"));
    try {
      const p = { r: s.r, sigma: SIGMA, b: B, x0: s.x0 };
      const [tr, se] = await Promise.all([
        run("trajectory", Object.assign({ T: runT }, p)),
        run("sensitivity", Object.assign({ d0: 1e-10, eps: 1e-3, Tavg: 200 }, p)),
      ]);
      if (!tr.stale) { lastTraj = Object.assign(tr.res, { r: s.r, x0: s.x0 }); renderTrajectory(); $("pc").classList.remove("busy"); }
      if (!se.stale) { lastSens = se.res; renderSensitivity(); $("pd").classList.remove("busy"); }
      $("status").textContent = `Updated in ${Math.round(performance.now() - t0)} ms.`;
    } catch (e) {
      $("status").textContent = "Computation failed: " + e.message;
    }
  }
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => { renderTrajectory(); renderSensitivity(); });

  syncSlider();
  changed(true);
})();
