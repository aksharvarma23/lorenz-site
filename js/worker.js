/* Heavy computations run here so the page stays responsive. */
importScripts("lorenz.js");
const L = self.Lorenz;

function variance(a) {
  let m = 0; for (const v of a) m += v; m /= a.length;
  let s = 0; for (const v of a) s += (v - m) * (v - m);
  return s / a.length;
}

self.onmessage = (ev) => {
  const { id, task, p } = ev.data;
  try {
    if (task === "trajectory") {
      const n = Math.round(p.T * 100) + 1;
      const tj = L.trajectory(p.sigma, p.r, p.b, p.x0, p.T, n);
      const zm = L.zMaxima(p.sigma, p.r, p.b, tj.end, 0, 40);
      const verdict = L.periodicity(zm.zm, p.r);
      self.postMessage({ id, task, result: { t: tj.t, X: tj.X, verdict } });
    } else if (task === "sensitivity") {
      const tw = L.twinSeparation(p.sigma, p.r, p.b, p.x0, p.d0, 40, 4001, 20);
      let Lsize = Math.sqrt(2 * (variance(tw.X[0]) + variance(tw.X[1]) + variance(tw.X[2])));
      Lsize = Math.max(Lsize, 1e-6);
      const ly = L.lyapunovTangent(p.sigma, p.r, p.b, p.x0, p.Tavg, 1, 50);
      const chaotic = ly.lam > 0.05;
      const fit = chaotic ? L.fitSlope(tw.t, tw.sep, Math.max(Lsize, 1)) : { slope: NaN, icpt: NaN, i0: -1, i1: -1 };
      const out = { t: tw.t, sep: tw.sep, fit, L: Lsize, lam: ly.lam, lt: ly.t, running: ly.running, chaotic };
      if (chaotic) {
        out.horizonEps = Math.log(Lsize / p.eps) / ly.lam;
        out.horizonD0 = Math.log(Lsize / p.d0) / ly.lam;
        out.gain10 = Math.log(10) / ly.lam;
        out.gain1000 = Math.log(1000) / ly.lam;
      }
      self.postMessage({ id, task, result: out });
    }
  } catch (err) {
    self.postMessage({ id, task, error: String(err && err.message || err) });
  }
};
