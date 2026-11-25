/* Lorenz system: numerical core (runs in the page and in the Web Worker).
 *
 *   dx/dt = sigma (y - x),   dy/dt = r x - y - x z,   dz/dt = x y - b z
 *
 * Integrator: Dormand-Prince 5(4) with adaptive steps (the method behind
 * MATLAB's ode45 / SciPy's RK45), run at tight tolerances. This is a port of
 * the Python project (lorenz/*.py), which uses SciPy's DOP853.
 */
(function (root) {
  "use strict";

  const SIGMA = 10, B = 8 / 3;
  const LIT = { pitchfork: 1, homoclinic: 13.926, attracting: 24.06, hopf: 24.7368 };

  // ------------------------------------------------------------ vector fields
  function field3(sigma, r, b) {
    return function (y, dy) {
      dy[0] = sigma * (y[1] - y[0]);
      dy[1] = r * y[0] - y[1] - y[0] * y[2];
      dy[2] = y[0] * y[1] - b * y[2];
    };
  }
  // two copies of the system integrated together (twin trajectories)
  function field6(sigma, r, b) {
    return function (y, dy) {
      dy[0] = sigma * (y[1] - y[0]);
      dy[1] = r * y[0] - y[1] - y[0] * y[2];
      dy[2] = y[0] * y[1] - b * y[2];
      dy[3] = sigma * (y[4] - y[3]);
      dy[4] = r * y[3] - y[4] - y[3] * y[5];
      dy[5] = y[3] * y[4] - b * y[5];
    };
  }
  // trajectory + tangent vector v, dv/dt = J(x) v
  function fieldTangent(sigma, r, b) {
    return function (y, dy) {
      const x = y[0], yy = y[1], z = y[2], a = y[3], c = y[4], d = y[5];
      dy[0] = sigma * (yy - x);
      dy[1] = r * x - yy - x * z;
      dy[2] = x * yy - b * z;
      dy[3] = sigma * (c - a);
      dy[4] = (r - z) * a - c - x * d;
      dy[5] = yy * a + x * c - b * d;
    };
  }

  // ------------------------------------------------------------ Dormand-Prince 5(4)
  const A21 = 1 / 5,
    A31 = 3 / 40, A32 = 9 / 40,
    A41 = 44 / 45, A42 = -56 / 15, A43 = 32 / 9,
    A51 = 19372 / 6561, A52 = -25360 / 2187, A53 = 64448 / 6561, A54 = -212 / 729,
    A61 = 9017 / 3168, A62 = -355 / 33, A63 = 46732 / 5247, A64 = 49 / 176, A65 = -5103 / 18656,
    A71 = 35 / 384, A73 = 500 / 1113, A74 = 125 / 192, A75 = -2187 / 6784, A76 = 11 / 84;
  const E1 = 71 / 57600, E3 = -71 / 16695, E4 = 71 / 1920, E5 = -17253 / 339200, E6 = 22 / 525, E7 = -1 / 40;

  /**
   * Integrate y' = f(y) from t0 to t1.
   * opts: rtol, atol, h0, maxStep,
   *       onStep(tOld, yOld, fOld, tNew, yNew, fNew)  called after every accepted step;
   *       return true to stop early.
   * Returns { y, t, steps }.
   */
  function integrate(f, y0, t0, t1, opts) {
    opts = opts || {};
    const n = y0.length;
    const rtol = opts.rtol ?? 1e-10, atol = opts.atol ?? 1e-12;
    const maxStep = opts.maxStep ?? 0.05;
    const onStep = opts.onStep;
    let y = Float64Array.from(y0), t = t0;
    const k1 = new Float64Array(n), k2 = new Float64Array(n), k3 = new Float64Array(n),
      k4 = new Float64Array(n), k5 = new Float64Array(n), k6 = new Float64Array(n),
      k7 = new Float64Array(n), yt = new Float64Array(n), yn = new Float64Array(n);
    f(y, k1);
    let h = opts.h0 ?? Math.min(1e-3, Math.abs(t1 - t0));
    let steps = 0;
    while (t < t1) {
      if (t + h > t1) h = t1 - t;
      if (h < 1e-14) break;
      for (let i = 0; i < n; i++) yt[i] = y[i] + h * A21 * k1[i];
      f(yt, k2);
      for (let i = 0; i < n; i++) yt[i] = y[i] + h * (A31 * k1[i] + A32 * k2[i]);
      f(yt, k3);
      for (let i = 0; i < n; i++) yt[i] = y[i] + h * (A41 * k1[i] + A42 * k2[i] + A43 * k3[i]);
      f(yt, k4);
      for (let i = 0; i < n; i++) yt[i] = y[i] + h * (A51 * k1[i] + A52 * k2[i] + A53 * k3[i] + A54 * k4[i]);
      f(yt, k5);
      for (let i = 0; i < n; i++) yt[i] = y[i] + h * (A61 * k1[i] + A62 * k2[i] + A63 * k3[i] + A64 * k4[i] + A65 * k5[i]);
      f(yt, k6);
      for (let i = 0; i < n; i++) yn[i] = y[i] + h * (A71 * k1[i] + A73 * k3[i] + A74 * k4[i] + A75 * k5[i] + A76 * k6[i]);
      f(yn, k7);
      let err = 0;
      for (let i = 0; i < n; i++) {
        const e = h * (E1 * k1[i] + E3 * k3[i] + E4 * k4[i] + E5 * k5[i] + E6 * k6[i] + E7 * k7[i]);
        const sc = atol + rtol * Math.max(Math.abs(y[i]), Math.abs(yn[i]));
        err += (e / sc) * (e / sc);
      }
      err = Math.sqrt(err / n);
      if (err <= 1 || h < 1e-12) {
        const tOld = t; t += h; steps++;
        let stop = false;
        if (onStep) stop = onStep(tOld, y, k1, t, yn, k7, h);
        y.set(yn); k1.set(k7);   // FSAL
        if (stop) break;
        const fac = err === 0 ? 5 : Math.min(5, Math.max(0.2, 0.9 * Math.pow(err, -0.2)));
        h = Math.min(maxStep, h * fac);
      } else {
        h *= Math.max(0.1, 0.9 * Math.pow(err, -0.25));
      }
    }
    return { y, t, steps };
  }

  // cubic Hermite interpolation of one component on [t0, t0+h]
  function hermite(y0, y1, f0, f1, h, s) {
    const s2 = s * s, s3 = s2 * s;
    return (2 * s3 - 3 * s2 + 1) * y0 + (s3 - 2 * s2 + s) * h * f0 + (-2 * s3 + 3 * s2) * y1 + (s3 - s2) * h * f1;
  }

  /** Uniformly sampled trajectory (n samples over [0, T]) via Hermite interpolation of accepted steps. */
  function trajectory(sigma, r, b, x0, T, n, opts) {
    const f = field3(sigma, r, b);
    const t = new Float64Array(n), X = [new Float64Array(n), new Float64Array(n), new Float64Array(n)];
    const dt = T / (n - 1);
    let k = 0;
    X[0][0] = x0[0]; X[1][0] = x0[1]; X[2][0] = x0[2]; t[0] = 0; k = 1;
    const res = integrate(f, x0, 0, T, Object.assign({ rtol: 1e-10, atol: 1e-12 }, opts, {
      onStep(tO, yO, fO, tN, yN, fN, h) {
        while (k < n && k * dt <= tN + 1e-12) {
          const s = (k * dt - tO) / h;
          for (let c = 0; c < 3; c++) X[c][k] = hermite(yO[c], yN[c], fO[c], fN[c], h, s);
          t[k] = k * dt; k++;
        }
        return false;
      }
    }));
    while (k < n) { for (let c = 0; c < 3; c++) X[c][k] = res.y[c]; t[k] = k * dt; k++; }
    return { t, X, end: Array.from(res.y) };
  }

  /** Successive maxima of z (dz/dt crosses + to -), located by Hermite interpolation. */
  function zMaxima(sigma, r, b, x0, tTrans, tRec, opts) {
    const f = field3(sigma, r, b);
    let y = x0;
    if (tTrans > 0) y = Array.from(integrate(f, y, 0, tTrans, Object.assign({ rtol: 1e-10, atol: 1e-12 }, opts)).y);
    const zm = [];
    const res = integrate(f, y, 0, tRec, Object.assign({ rtol: 1e-10, atol: 1e-12 }, opts, {
      onStep(tO, yO, fO, tN, yN, fN, h) {
        if (fO[2] > 0 && fN[2] <= 0) {
          // z(s) Hermite cubic; its derivative is a quadratic in s; take the root in [0,1]
          const z0 = yO[2], z1 = yN[2], d0 = h * fO[2], d1 = h * fN[2];
          // z'(s) = (6s^2-6s) z0 + (3s^2-4s+1) d0 + (-6s^2+6s) z1 + (3s^2-2s) d1
          const qa = 6 * z0 + 3 * d0 - 6 * z1 + 3 * d1;
          const qb = -6 * z0 - 4 * d0 + 6 * z1 - 2 * d1;
          const qc = d0;
          let s;
          if (Math.abs(qa) < 1e-14) s = -qc / qb;
          else {
            const disc = Math.max(0, qb * qb - 4 * qa * qc);
            const s1 = (-qb + Math.sqrt(disc)) / (2 * qa), s2 = (-qb - Math.sqrt(disc)) / (2 * qa);
            s = (s1 >= -1e-9 && s1 <= 1 + 1e-9) ? s1 : s2;
          }
          s = Math.min(1, Math.max(0, s));
          zm.push(hermite(z0, z1, fO[2], fN[2], h, s));
        }
        return false;
      }
    }));
    return { zm, end: Array.from(res.y) };
  }

  /** Classify recorded z maxima at r: steady, period-k, or chaotic (port of bifurcation.periodicity). */
  function periodicity(zmax, r, tol = 1e-3, maxPeriod = 16) {
    const z = zmax;
    if (z.length < 4) return { kind: "steady", k: 0, text: "steady (trajectory settles on a fixed point)" };
    const half = z.slice(Math.floor(z.length / 2));
    if (half.every(v => Math.abs(v - (r - 1)) < 0.05))
      return { kind: "steady", k: 0, text: "steady (trajectory settles on a fixed point)" };
    // weakly damped spiral onto C± (just below r_H): peaks still creep toward r - 1
    const dist = half.map(v => Math.abs(v - (r - 1)));
    let shrinking = dist[dist.length - 1] < 0.97 * dist[0] && dist[0] < 3;
    for (let i = 1; i < dist.length && shrinking; i++) if (dist[i] > dist[i - 1] * (1 + 1e-6)) shrinking = false;
    if (shrinking) return { kind: "steady", k: 0, text: "steady (spiralling slowly onto C±)" };
    const mean = half.reduce((a, v) => a + v, 0) / half.length;
    const scale = Math.max(1, Math.abs(mean));
    for (let k = 1; k <= maxPeriod; k++) {
      if (half.length > 2 * k) {
        let mx = 0;
        for (let i = k; i < half.length; i++) mx = Math.max(mx, Math.abs(half[i] - half[i - k]));
        if (mx < tol * scale)
          return { kind: "periodic", k, text: `periodic orbit (${k} distinct z-maxima per cycle)` };
      }
    }
    return { kind: "chaotic", k: 0, text: "aperiodic (chaotic)" };
  }

  // ------------------------------------------------------------ linear algebra
  function jacobian(p, sigma, r, b) {
    return [[-sigma, sigma, 0], [r - p[2], -1, -p[0]], [p[1], p[0], -b]];
  }
  // eigenvalues of a real 3x3 matrix: roots of lambda^3 + a lambda^2 + c lambda + d
  function eig3(M) {
    const tr = M[0][0] + M[1][1] + M[2][2];
    const m2 = M[0][0] * M[1][1] - M[0][1] * M[1][0] + M[0][0] * M[2][2] - M[0][2] * M[2][0]
      + M[1][1] * M[2][2] - M[1][2] * M[2][1];
    const det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
      - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
      + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    return cubicRoots(-tr, m2, -det);
  }
  // roots of x^3 + a x^2 + b x + c, returned as [{re, im}], sorted by real part (descending)
  function cubicRoots(a, b, c) {
    const p = b - a * a / 3, q = 2 * a * a * a / 27 - a * b / 3 + c;
    const sh = -a / 3;
    let roots;
    const D = (q * q) / 4 + (p * p * p) / 27;
    if (D > 1e-14 * Math.max(1, Math.abs(q * q))) {
      const sq = Math.sqrt(D);
      const u = Math.cbrt(-q / 2 + sq), v = Math.cbrt(-q / 2 - sq);
      const re1 = u + v + sh;
      const re2 = -(u + v) / 2 + sh, im2 = (Math.sqrt(3) / 2) * (u - v);
      roots = [{ re: re1, im: 0 }, { re: re2, im: Math.abs(im2) }, { re: re2, im: -Math.abs(im2) }];
    } else {
      const m = 2 * Math.sqrt(Math.max(0, -p / 3));
      const arg = m === 0 ? 0 : Math.max(-1, Math.min(1, (3 * q) / (p * m)));
      const th = Math.acos(arg) / 3;
      roots = [0, 1, 2].map(k => ({ re: m * Math.cos(th - (2 * Math.PI * k) / 3) + sh, im: 0 }));
    }
    // one Newton polish per real root
    for (const rt of roots) if (rt.im === 0) {
      for (let it = 0; it < 3; it++) {
        const x = rt.re, fx = ((x + a) * x + b) * x + c, dfx = (3 * x + 2 * a) * x + b;
        if (dfx !== 0) rt.re = x - fx / dfx;
      }
    }
    roots.sort((s, t) => t.re - s.re || t.im - s.im);
    return roots;
  }

  function classify(eigs, tol = 1e-9) {
    const re = eigs.map(e => e.re), im = eigs.map(e => e.im);
    const hasC = im.some(v => Math.abs(v) > tol);
    const nPos = re.filter(v => v > tol).length, nNeg = re.filter(v => v < -tol).length;
    const nZero = 3 - nPos - nNeg;
    if (nZero > 0) return ["non-hyperbolic", "An eigenvalue has zero real part, so this is a bifurcation point and linearisation alone cannot decide stability."];
    const cre = eigs.filter(e => Math.abs(e.im) > tol).map(e => e.re);
    const rotU = hasC && cre.every(v => v > 0), rotS = hasC && cre.every(v => v < 0);
    if (nPos === 0) return hasC ? ["stable spiral", "Nearby trajectories spiral in. Attracting."]
      : ["stable node", "Nearby trajectories decay straight in without rotating. Attracting."];
    if (nNeg === 0) return hasC ? ["unstable spiral", "Everything nearby spirals out. Repelling."]
      : ["unstable node", "Everything nearby is pushed away. Repelling."];
    const kind = `${nPos}-D unstable / ${nNeg}-D stable`;
    if (rotU) return ["saddle-focus", `${kind}: drawn in along one direction, spiralling outward in a plane. Unstable.`];
    if (rotS) return ["saddle-focus", `${kind}: spiralling in within a plane, ejected along one direction. Unstable.`];
    return ["saddle", `${kind}: attracting along some directions, repelling along others. Unstable.`];
  }

  function fixedPoints(r, sigma, b) {
    const pts = [{ label: "origin", p: [0, 0, 0] }];
    const q = b * (r - 1);
    if (q > 0) {
      const s = Math.sqrt(q);
      pts.push({ label: "C+", p: [s, s, r - 1] }, { label: "C−", p: [-s, -s, r - 1] });
    }
    return pts.map(o => {
      const ev = eig3(jacobian(o.p, sigma, r, b));
      const [type, verdict] = classify(ev);
      return Object.assign(o, { eig: ev, type, verdict, stable: ev.every(e => e.re < -1e-9) });
    });
  }

  function hopfR(sigma, b) {
    return sigma - b - 1 <= 0 ? null : (sigma * (sigma + b + 3)) / (sigma - b - 1);
  }
  const divergence = (sigma, b) => -(sigma + 1 + b);

  // ------------------------------------------------------------ shooting (global bifurcations)
  function unstableStart(r, sigma, delta = 1e-9) {
    const lam = 0.5 * (-(sigma + 1) + Math.sqrt((sigma + 1) ** 2 + 4 * sigma * (r - 1)));
    const v = [1, 1 + lam / sigma, 0], nv = Math.hypot(v[0], v[1]);
    return [delta * v[0] / nv, delta * v[1] / nv, 0];
  }
  // number of sign changes of x along the right branch of the origin's unstable manifold
  function crossings(r, sigma, b, T, cap) {
    let count = 0;
    integrate(field3(sigma, r, b), unstableStart(r, sigma), 0, T, {
      rtol: 1e-10, atol: 1e-15,
      onStep(tO, yO, fO, tN, yN) {
        if ((yO[0] > 0) !== (yN[0] > 0) && yN[0] !== 0) count++;
        return count >= cap;
      }
    });
    return count;
  }
  function bisect(pred, lo, hi, tol) {
    const plo = pred(lo), phi = pred(hi);
    if (plo === phi) return null;
    while (hi - lo > tol) {
      const mid = 0.5 * (lo + hi);
      if (pred(mid) === plo) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
  }
  function homoclinicR(sigma, b, tol = 1e-5) {
    const rh = hopfR(sigma, b);
    return bisect(r => crossings(r, sigma, b, 60, 1) > 0, 1.5, rh ?? 60, tol);
  }
  function attractingR(sigma, b, rhom, tol = 1e-4) {
    if (rhom == null) return null;
    const rh = hopfR(sigma, b);
    return bisect(r => crossings(r, sigma, b, 150, 2) > 1, rhom + 0.05, (rh ?? 60) - 1e-6, tol);
  }
  function transitions(sigma, b) {
    const t = { pitchfork: 1, hopf: hopfR(sigma, b) };
    t.homoclinic = homoclinicR(sigma, b);
    t.attracting = attractingR(sigma, b, t.homoclinic);
    return t;
  }

  function regime(r, tr) {
    const rh = tr.hopf, rhom = tr.homoclinic, ra = tr.attracting;
    const f = (v, d) => (v == null ? "?" : v.toFixed(d));
    if (r < 0) return ["unphysical", "Negative heating is unphysical; the origin is the only fixed point."];
    if (r < 1) return ["conduction", "The origin attracts everything: the fluid stays at rest and heat moves by conduction only."];
    if (Math.abs(r - 1) < 1e-9) return ["pitchfork", "Supercritical pitchfork: the origin hands its stability to the new pair C±."];
    if (rhom != null && r < rhom) return ["steady rolls", `1 < r < ${f(rhom, 3)}. C± are stable and the origin is a saddle: every run settles on one of the two steady rolls.`];
    if (ra != null && r < ra) return ["transient chaos", `${f(rhom, 3)} < r < ${f(ra, 2)}. A chaotic set exists but does not attract: runs wander chaotically for a while, then land on C±.`];
    if (rh != null && r < rh) return ["bistable", `${f(ra, 2)} < r < ${f(rh, 4)}. The strange attractor coexists with stable C±: the starting point decides between chaos and a steady roll.`];
    if (rh == null) return ["no Hopf", "σ ≤ b + 1: no Hopf bifurcation, so C± stay stable for all r > 1."];
    return ["chaos", `r > ${f(rh, 4)}. C± have lost stability (subcritical Hopf). No stable fixed point is left: the motion is chaotic, except inside periodic windows (try r = 100.5).`];
  }

  // ------------------------------------------------------------ sensitivity
  function twinSeparation(sigma, r, b, x0, d0, T, n, tTrans) {
    const f3 = field3(sigma, r, b);
    let y = x0;
    if (tTrans > 0) y = Array.from(integrate(f3, y, 0, tTrans, { rtol: 1e-10, atol: 1e-12 }).y);
    const u = 1 / Math.sqrt(3);
    const s0 = [y[0], y[1], y[2], y[0] + d0 * u, y[1] + d0 * u, y[2] + d0 * u];
    const t = new Float64Array(n), sep = new Float64Array(n), X = [new Float64Array(n), new Float64Array(n), new Float64Array(n)];
    const dt = T / (n - 1);
    let k = 1;
    t[0] = 0; sep[0] = d0; for (let c = 0; c < 3; c++) X[c][0] = s0[c];
    const tmp = new Float64Array(6);
    const res = integrate(field6(sigma, r, b), s0, 0, T, {
      rtol: 1e-10, atol: 1e-12,
      onStep(tO, yO, fO, tN, yN, fN, h) {
        while (k < n && k * dt <= tN + 1e-12) {
          const s = (k * dt - tO) / h;
          for (let c = 0; c < 6; c++) tmp[c] = hermite(yO[c], yN[c], fO[c], fN[c], h, s);
          t[k] = k * dt;
          sep[k] = Math.hypot(tmp[3] - tmp[0], tmp[4] - tmp[1], tmp[5] - tmp[2]);
          for (let c = 0; c < 3; c++) X[c][k] = tmp[c];
          k++;
        }
        return false;
      }
    });
    while (k < n) { t[k] = k * dt; sep[k] = sep[k - 1]; for (let c = 0; c < 3; c++) X[c][k] = res.y[c]; k++; }
    return { t, sep, X };
  }

  function fitSlope(t, sep, L, tMin = 1, upperFrac = 1e-3) {
    let stop = sep.length;
    for (let i = 0; i < sep.length; i++) if (sep[i] > upperFrac * L) { stop = i; break; }
    let n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0, i0 = -1, i1 = -1;
    for (let i = 0; i < stop; i++) if (t[i] >= tMin) {
      const ly = Math.log(Math.max(sep[i], 1e-300));
      n++; sx += t[i]; sy += ly; sxx += t[i] * t[i]; sxy += t[i] * ly;
      if (i0 < 0) i0 = i; i1 = i;
    }
    if (n < 10) return { slope: NaN, icpt: NaN, i0, i1 };
    const slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    return { slope, icpt: (sy - slope * sx) / n, i0, i1 };
  }

  /** Largest Lyapunov exponent, tangent-linear Benettin, renormalised every tau. */
  function lyapunovTangent(sigma, r, b, x0, Tavg, tau = 1, tTrans = 20) {
    let y = x0;
    if (tTrans > 0) y = Array.from(integrate(field3(sigma, r, b), y, 0, tTrans, { rtol: 1e-10, atol: 1e-12 }).y);
    const ft = fieldTangent(sigma, r, b);
    let s = [y[0], y[1], y[2], 1 / Math.sqrt(3), 1 / Math.sqrt(3), 1 / Math.sqrt(3)];
    const nren = Math.round(Tavg / tau);
    const running = new Float64Array(nren), tt = new Float64Array(nren);
    let acc = 0;
    for (let k = 0; k < nren; k++) {
      const out = integrate(ft, s, 0, tau, { rtol: 1e-9, atol: 1e-12 }).y;
      const nv = Math.hypot(out[3], out[4], out[5]);
      acc += Math.log(nv);
      s = [out[0], out[1], out[2], out[3] / nv, out[4] / nv, out[5] / nv];
      running[k] = acc / ((k + 1) * tau); tt[k] = (k + 1) * tau;
    }
    return { lam: running[nren - 1], t: tt, running };
  }

  root.Lorenz = {
    SIGMA, B, LIT, integrate, trajectory, zMaxima, periodicity, fixedPoints, eig3, hopfR, divergence,
    transitions, homoclinicR, attractingR, regime, twinSeparation, fitSlope, lyapunovTangent, field3
  };
})(typeof self !== "undefined" ? self : this);
