"""
Lorenz system: regime-by-regime time-series generation and phase portraits.

    dx/dt = sigma (y - x)
    dy/dt = r x - y - x z
    dz/dt = x y - b z          sigma = 10, b = 8/3, r varied

Integrator: scipy solve_ivp, DOP853, rtol 1e-10, atol 1e-12, sampled every dt = 0.01.
Writes one CSV per run (t, x, y, z) to data/, a metadata.json describing every run,
and regime_summary.csv with the computed checks (fixed points, eigenvalues, lambda_1).
"""
import json, os, sys
import numpy as np
from scipy.integrate import solve_ivp

SIGMA, B = 10.0, 8.0 / 3.0
RTOL, ATOL, DT = 1e-10, 1e-12, 0.01
OUT = sys.argv[1] if len(sys.argv) > 1 else "."
os.makedirs(os.path.join(OUT, "data"), exist_ok=True)


def f(t, s, r):
    x, y, z = s
    return [SIGMA * (y - x), r * x - y - x * z, x * y - B * z]


def run(r, x0, T, t_discard=0.0):
    """Integrate, optionally discard a transient, return samples on a uniform grid."""
    y0 = np.asarray(x0, float)
    if t_discard > 0:
        y0 = solve_ivp(f, (0, t_discard), y0, args=(r,), method="DOP853",
                       rtol=RTOL, atol=ATOL).y[:, -1]
    t = np.round(np.arange(0, T + DT / 2, DT), 10)
    sol = solve_ivp(f, (0, T), y0, args=(r,), method="DOP853", rtol=RTOL, atol=ATOL, t_eval=t)
    assert sol.success, sol.message
    return sol.t, sol.y


def fixed_points(r):
    pts = [("origin", np.zeros(3))]
    if r > 1:
        s = np.sqrt(B * (r - 1))
        pts += [("C+", np.array([s, s, r - 1])), ("C-", np.array([-s, -s, r - 1]))]
    return pts


def jac(p, r):
    x, y, z = p
    return np.array([[-SIGMA, SIGMA, 0], [r - z, -1, -x], [y, x, -B]])


def lyap1(r, x0, T=300.0, tau=1.0, t_discard=50.0):
    """Largest Lyapunov exponent: tangent-linear (variational) vector dv/dt = J(x) v,
    renormalised every tau (Benettin). Avoids subtracting two nearly equal states."""
    def fv(t, s):
        x, y, z, a, b_, c = s
        return [SIGMA * (y - x), r * x - y - x * z, x * y - B * z,
                SIGMA * (b_ - a), (r - z) * a - b_ - x * c, y * a + x * b_ - B * c]
    y = np.asarray(x0, float)
    if t_discard > 0:
        y = solve_ivp(f, (0, t_discard), y, args=(r,), method="DOP853", rtol=RTOL, atol=ATOL).y[:, -1]
    v = np.ones(3) / np.sqrt(3); acc = 0.0; n = int(T / tau)
    for _ in range(n):
        s = solve_ivp(fv, (0, tau), np.concatenate([y, v]), method="DOP853", rtol=RTOL, atol=ATOL).y[:, -1]
        y, v = s[:3], s[3:]; d = np.linalg.norm(v); acc += np.log(d); v = v / d
    return acc / (n * tau)


def observed(r, x0, t_discard=300.0, T=200.0):
    """Classify the long-run behaviour from exact z-maxima (solver events)."""
    def zdot(t, s, r): return s[0] * s[1] - B * s[2]
    zdot.direction = -1
    y = solve_ivp(f, (0, t_discard), np.asarray(x0, float), args=(r,), method="DOP853",
                  rtol=RTOL, atol=ATOL).y[:, -1]
    sol = solve_ivp(f, (0, T), y, args=(r,), method="DOP853", rtol=RTOL, atol=ATOL, events=zdot)
    zm = sol.y_events[0][:, 2] if len(sol.t_events[0]) else np.array([])
    if zm.size < 4 or np.all(np.abs(zm[-10:] - (r - 1)) < 0.5):
        return "steady (settles on a fixed point)"
    zm = zm[len(zm) // 2:]
    for k in range(1, 17):
        if zm.size > 2 * k and np.max(np.abs(zm[k:] - zm[:-k])) < 1e-3 * max(1, abs(zm.mean())):
            return f"periodic ({k} z-peak{'s' if k > 1 else ''} per cycle)"
    return "chaotic (aperiodic)"


# ------------------------------------------------------------------ the runs
# (id, r, initial condition, record length, transient discarded, expected behaviour)
RUNS = [
    ("r00.5_conduction",     0.5, (1, 1, 1),            30,  0,  "decays to the origin: fluid at rest"),
    ("r10_steady_roll",     10.0, (1, 1, 1),            30,  0,  "damped spiral onto C+ or C-: steady convection roll"),
    ("r20_transient_chaos", 20.0, (0.3, 18.6, 19.5),   100,  0,  "chaotic-looking wandering, then settles on C+/-"),
    ("r24.4_start_near_C",  24.4, None,                100,  0,  "bistable window: starting near C+ -> settles on C+"),
    ("r24.4_start_far",     24.4, (1, 1, 1),           100,  0,  "bistable window: starting far away -> sustained chaos"),
    ("r28_chaos",           28.0, (1, 1, 1),           200, 20,  "strange attractor (butterfly): aperiodic, never settles"),
    ("r100.5_periodic",    100.5, (1, 1, 1),           50, 200,  "periodic window: a closed loop"),
    ("r160_periodic",      160.0, (1, 1, 1),           50, 200,  "periodic window: a closed loop"),
    ("r350_periodic",      350.0, (1, 1, 1),           20, 500,  "stable symmetric limit cycle (r > 313)"),
]

meta = {"system": "Lorenz 1963", "sigma": SIGMA, "b": B,
        "solver": {"method": "DOP853 (scipy.integrate.solve_ivp)", "rtol": RTOL, "atol": ATOL},
        "sampling_dt": DT, "columns": ["t", "x", "y", "z"], "runs": []}
summary = []

for rid, r, x0, T, tdisc, expected in RUNS:
    if x0 is None:                                   # just off C+ for the bistable test
        s = np.sqrt(B * (r - 1)); x0 = (s + 0.5, s + 0.5, r - 1 + 0.5)
    t, X = run(r, x0, T, tdisc)
    path = os.path.join(OUT, "data", f"lorenz_{rid}.csv")
    np.savetxt(path, np.column_stack([t, X.T]), delimiter=",", header="t,x,y,z", comments="", fmt="%.10g")

    fps = []
    for lab, p in fixed_points(r):
        ev = np.linalg.eigvals(jac(p, r)); ev = ev[np.argsort(-ev.real)]
        fps.append({"label": lab, "point": [round(float(v), 4) for v in p],
                    "eigenvalues": [f"{e.real:+.4f}{e.imag:+.4f}j" for e in ev],
                    "stable": bool(np.all(ev.real < 0))})
    lam = lyap1(r, x0, T=300.0, t_discard=tdisc if tdisc else 50.0)
    obs = observed(r, x0)
    meta["runs"].append({"id": rid, "file": f"data/lorenz_{rid}.csv", "r": r,
                         "initial_condition": [float(v) for v in x0], "record_length": T,
                         "transient_discarded": tdisc, "samples": int(len(t)),
                         "expected": expected, "observed_long_run": obs,
                         "fixed_points": fps, "lambda1_estimate": round(lam, 4)})
    summary.append((rid, r, len(t), lam, expected, obs))
    print(f"{rid:22s} r={r:6.1f}  n={len(t):6d}  lambda1={lam:+.4f}  {obs}", flush=True)

# twin trajectories at r = 28, 1e-10 apart: the sensitivity data set
r = 28.0
y0 = solve_ivp(f, (0, 20), [1, 1, 1], args=(r,), method="DOP853", rtol=RTOL, atol=ATOL).y[:, -1]
d = 1e-10 * np.ones(3) / np.sqrt(3)
def f2(t, s): return f(t, s[:3], r) + f(t, s[3:], r)
t = np.round(np.arange(0, 40 + DT / 2, DT), 10)
sol = solve_ivp(f2, (0, 40), np.concatenate([y0, y0 + d]), method="DOP853", rtol=RTOL, atol=ATOL, t_eval=t)
sep = np.linalg.norm(sol.y[3:] - sol.y[:3], axis=0)
np.savetxt(os.path.join(OUT, "data", "lorenz_r28_twin_trajectories.csv"),
           np.column_stack([sol.t, sol.y[:3].T, sol.y[3:].T, sep]), delimiter=",",
           header="t,x1,y1,z1,x2,y2,z2,separation", comments="", fmt="%.12g")
meta["runs"].append({"id": "r28_twin_trajectories", "file": "data/lorenz_r28_twin_trajectories.csv",
                     "r": 28.0, "initial_separation": 1e-10, "record_length": 40,
                     "expected": "separation grows like exp(0.906 t) until it saturates at the attractor size"})

with open(os.path.join(OUT, "metadata.json"), "w") as fh:
    json.dump(meta, fh, indent=1)
with open(os.path.join(OUT, "regime_summary.csv"), "w") as fh:
    fh.write("run,r,samples,lambda1_estimate,expected,observed_long_run\n")
    for rid, r_, n, lam, ex, ob in summary:
        fh.write(f"{rid},{r_},{n},{lam:.4f},\"{ex}\",\"{ob}\"\n")
print("done")
