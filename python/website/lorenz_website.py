"""
Lorenz Explorer, Python version of the website's calculations.

The website (index.html + js/*.js) was ported from the Python project. This script
does the same four calculations in Python, so every number shown on the site can be
reproduced and checked:

    Panel 1  fixed points, Jacobian eigenvalues, stable / unstable
    Panel 2  where r sits among the bifurcations (1, 13.93, 24.06, 24.74)
    Panel 3  the trajectory, its long-run label (steady / periodic / chaotic), CSV export
    Panel 4  twin runs 1e-10 apart, largest Lyapunov exponent, prediction horizon

    dx/dt = sigma (y - x),   dy/dt = r x - y - x z,   dz/dt = x y - b z,   sigma = 10, b = 8/3

Solver: scipy solve_ivp, DOP853, rtol 1e-10, atol 1e-12 (the website uses a
Dormand-Prince 5(4) solver at the same tolerances).

Usage
    python lorenz_website.py --preset butterfly          one of the nine presets
    python lorenz_website.py --r 24.4 --x0 1 1 1 --T 100  any r and starting point
    python lorenz_website.py --preset butterfly --csv run.csv --plot run.png
    python lorenz_website.py --all                       all nine presets in one table
    python lorenz_website.py --shoot                     recompute 13.93 and 24.06 by shooting
"""
import argparse
import numpy as np
from scipy.integrate import solve_ivp

SIGMA, B = 10.0, 8.0 / 3.0
RTOL, ATOL = 1e-10, 1e-12

# Transition values used by panel 2. Pitchfork and Hopf are exact formulas; the two
# global ones were found by shooting in the Python project (run with --shoot to redo it).
TRANSITIONS = dict(pitchfork=1.0, homoclinic=13.9266, attracting=24.059,
                   hopf=SIGMA * (SIGMA + B + 3) / (SIGMA - B - 1))

# The nine presets on the website = the nine runs of the 11 November time-series data.
PRESETS = {
    "rest":             dict(name="Rest",                r=0.5,   x0=(1, 1, 1),          T=30),
    "steady":           dict(name="Steady roll",         r=10.0,  x0=(1, 1, 1),          T=30),
    "transient":        dict(name="Transient chaos",     r=20.0,  x0=(0.3, 18.6, 19.5),  T=60),
    "bistable-near":    dict(name="Bistable: near C+",   r=24.4,  x0=(8.4, 8.4, 23.9),   T=100),
    "bistable-far":     dict(name="Bistable: far start", r=24.4,  x0=(1, 1, 1),          T=100),
    "butterfly":        dict(name="Butterfly",           r=28.0,  x0=(1, 1, 1),          T=60),
    "window-3":         dict(name="Periodic window",     r=100.5, x0=(1, 1, 1),          T=60),
    "window-2":         dict(name="Window, 2 peaks",     r=160.0, x0=(1, 1, 1),          T=60),
    "single-cycle":     dict(name="Single cycle",        r=350.0, x0=(1, 1, 1),          T=40),
}


# ---------------------------------------------------------------- equations
def field(r, sigma=SIGMA, b=B):
    def f(t, s):
        x, y, z = s
        return [sigma * (y - x), r * x - y - x * z, x * y - b * z]
    return f


def field_twin(r, sigma=SIGMA, b=B):
    def f(t, s):
        x, y, z, X, Y, Z = s
        return [sigma * (y - x), r * x - y - x * z, x * y - b * z,
                sigma * (Y - X), r * X - Y - X * Z, X * Y - b * Z]
    return f


def field_tangent(r, sigma=SIGMA, b=B):
    """Trajectory plus a tangent vector v with dv/dt = J(x) v."""
    def f(t, s):
        x, y, z, a, c, d = s
        return [sigma * (y - x), r * x - y - x * z, x * y - b * z,
                sigma * (c - a), (r - z) * a - c - x * d, y * a + x * c - b * d]
    return f


def solve(f, y0, T, **kw):
    return solve_ivp(f, (0, T), y0, method="DOP853", rtol=RTOL, atol=ATOL, **kw)


# ---------------------------------------------------------------- panel 1
def jacobian(p, r, sigma=SIGMA, b=B):
    x, y, z = p
    return np.array([[-sigma, sigma, 0.0], [r - z, -1.0, -x], [y, x, -b]])


def classify(eig, tol=1e-9):
    re, im = eig.real, eig.imag
    cplx = np.any(np.abs(im) > tol)
    npos, nneg = int(np.sum(re > tol)), int(np.sum(re < -tol))
    if npos + nneg < 3:
        return "non-hyperbolic"
    if npos == 0:
        return "stable spiral" if cplx else "stable node"
    if nneg == 0:
        return "unstable spiral" if cplx else "unstable node"
    return "saddle-focus" if cplx else "saddle"


def fixed_points(r, sigma=SIGMA, b=B):
    pts = [("origin", np.zeros(3))]
    q = b * (r - 1)
    if q > 0:
        s = np.sqrt(q)
        pts += [("C+", np.array([s, s, r - 1])), ("C-", np.array([-s, -s, r - 1]))]
    out = []
    for label, p in pts:
        ev = np.linalg.eigvals(jacobian(p, r, sigma, b))
        ev = ev[np.argsort(-ev.real)]
        out.append(dict(label=label, point=p, eig=ev, type=classify(ev), stable=bool(np.all(ev.real < -1e-9))))
    return out


# ---------------------------------------------------------------- panel 2
def regime(r, tr=TRANSITIONS):
    if r < 1:
        return "rest: the origin attracts everything; heat moves by conduction only"
    if r < tr["homoclinic"]:
        return "steady rolls: C+- stable, every run settles on one roll"
    if r < tr["attracting"]:
        return "transient chaos: a chaotic set exists but does not attract; runs settle on C+-"
    if r < tr["hopf"]:
        return "bistable: chaos coexists with stable C+-; the starting point decides"
    return "chaos: C+- unstable (subcritical Hopf); chaotic except inside periodic windows"


def _unstable_start(r, sigma=SIGMA, delta=1e-9):
    lam = 0.5 * (-(sigma + 1) + np.sqrt((sigma + 1) ** 2 + 4 * sigma * (r - 1)))
    v = np.array([1.0, 1.0 + lam / sigma, 0.0])
    return delta * v / np.linalg.norm(v)


def _sign_changes(r, T, sigma=SIGMA, b=B):
    ev = lambda t, u: u[0]
    sol = solve_ivp(field(r, sigma, b), (0, T), _unstable_start(r, sigma), method="DOP853",
                    rtol=RTOL, atol=1e-15, events=ev)
    return len(sol.t_events[0])


def _bisect(pred, lo, hi, tol):
    plo = pred(lo)
    if plo == pred(hi):
        return None
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if pred(mid) == plo:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def shoot_transitions(sigma=SIGMA, b=B):
    """Homoclinic explosion and onset of attracting chaos, by shooting along the
    origin's unstable direction (same method as the Python project)."""
    rh = sigma * (sigma + b + 3) / (sigma - b - 1)
    hom = _bisect(lambda r: _sign_changes(r, 60.0, sigma, b) > 0, 1.5, rh, 1e-5)
    att = _bisect(lambda r: _sign_changes(r, 150.0, sigma, b) > 1, hom + 0.05, rh - 1e-6, 1e-4)
    return dict(pitchfork=1.0, homoclinic=hom, attracting=att, hopf=rh)


# ---------------------------------------------------------------- panel 3
def trajectory(r, x0, T, dt=0.01):
    t = np.linspace(0, T, int(round(T / dt)) + 1)
    sol = solve(field(r), np.asarray(x0, float), T, t_eval=t)
    return sol.t, sol.y


def z_maxima(r, x0, T):
    ev = lambda t, u: u[0] * u[1] - B * u[2]       # dz/dt; + to - is a maximum of z
    ev.direction = -1
    sol = solve(field(r), np.asarray(x0, float), T, events=ev)
    return sol.y_events[0][:, 2] if len(sol.t_events[0]) else np.empty(0)


def periodicity(zm, r, tol=1e-3, max_period=16):
    """The long-run label shown under panel 3."""
    if zm.size < 4:
        return "steady (trajectory settles on a fixed point)"
    half = zm[len(zm) // 2:]
    if np.all(np.abs(half - (r - 1)) < 0.05) or np.all(np.abs(half) < 0.05):   # on C+- or at the origin
        return "steady (trajectory settles on a fixed point)"
    dist = np.abs(half - (r - 1))
    if dist[-1] < 0.97 * dist[0] and dist[0] < 3 and np.all(np.diff(dist) <= dist[:-1] * 1e-6):
        return "steady (spiralling slowly onto C+-)"
    scale = max(1.0, abs(half.mean()))
    for k in range(1, max_period + 1):
        if half.size > 2 * k and np.max(np.abs(half[k:] - half[:-k])) < tol * scale:
            return f"periodic orbit ({k} distinct z-maxima per cycle)"
    return "aperiodic (chaotic)"


# ---------------------------------------------------------------- panel 4
def twin_separation(r, x0, d0=1e-10, T=40.0, t_transient=20.0):
    y = solve(field(r), np.asarray(x0, float), t_transient).y[:, -1]
    u = np.ones(3) / np.sqrt(3)
    t = np.linspace(0, T, 4001)
    sol = solve(field_twin(r), np.concatenate([y, y + d0 * u]), T, t_eval=t)
    sep = np.linalg.norm(sol.y[3:] - sol.y[:3], axis=0)
    L = float(np.sqrt(2 * np.sum(np.var(sol.y[:3], axis=1))))      # attractor size
    return sol.t, sep, max(L, 1e-6)


def fit_slope(t, sep, L, t_min=1.0, upper_frac=1e-3):
    big = np.nonzero(sep > upper_frac * L)[0]
    stop = big[0] if big.size else len(t)
    m = (t >= t_min) & (np.arange(len(t)) < stop)
    return np.polyfit(t[m], np.log(sep[m]), 1)[0] if m.sum() >= 10 else np.nan


def lyapunov(r, x0, T_avg=200.0, tau=1.0, t_transient=50.0):
    """Largest exponent, tangent-linear Benettin method, renormalised every tau."""
    y = solve(field(r), np.asarray(x0, float), t_transient).y[:, -1]
    s = np.concatenate([y, np.ones(3) / np.sqrt(3)])
    acc = 0.0
    n = int(round(T_avg / tau))
    for _ in range(n):
        s = solve(field_tangent(r), s, tau).y[:, -1]
        nv = np.linalg.norm(s[3:])
        acc += np.log(nv)
        s[3:] /= nv
    return acc / (n * tau)


# ---------------------------------------------------------------- run one case
def run(r, x0, T, eps=1e-3, csv=None, plot=None, quiet=False):
    fps = fixed_points(r)
    t, X = trajectory(r, x0, T)
    label = periodicity(z_maxima(r, X[:, -1], 40.0), r)
    ts, sep, L = twin_separation(r, x0)
    lam = lyapunov(r, x0)
    chaotic = lam > 0.05
    slope = fit_slope(ts, sep, max(L, 1.0)) if chaotic else np.nan
    horizon = np.log(L / eps) / lam if chaotic else np.inf
    if not quiet:
        print(f"\nr = {r:g}, start = {tuple(x0)}, run length = {T:g}")
        print("Panel 1  fixed points")
        for f in fps:
            ev = ", ".join(f"{e.real:+.3f}" + (f" {'+' if e.imag >= 0 else '-'} {abs(e.imag):.3f}i" if abs(e.imag) > 1e-10 else "")
                           for e in f["eig"])
            print(f"   {f['label']:7s} {np.round(f['point'], 2)}  eig: {ev}   {'stable' if f['stable'] else 'unstable'} ({f['type']})")
        print("Panel 2 ", regime(r))
        print("Panel 3 ", label)
        print(f"Panel 4  lambda1 = {lam:.3f}" + (f", slope fit = {slope:.3f}" if chaotic else "")
              + (f", horizon from eps = {eps:g}: {horizon:.1f} time units, +{np.log(10) / lam:.2f} per 10x better data"
                 if chaotic else ", horizon: infinite (not chaotic)"))
    if csv:
        np.savetxt(csv, np.column_stack([t, X.T]), delimiter=",", header="t,x,y,z", comments="", fmt="%.8f")
        if not quiet:
            print(f"saved {csv}  ({len(t)} rows)")
    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(14, 4.2))
        a = fig.add_subplot(1, 3, 1, projection="3d"); a.plot(*X, lw=0.4); a.set_title("3-D path")
        a = fig.add_subplot(1, 3, 2); a.plot(t, X[0], lw=0.8); a.set_xlabel("time"); a.set_ylabel("x"); a.set_title("x(t)")
        a = fig.add_subplot(1, 3, 3); a.semilogy(ts, sep, lw=0.8); a.axhline(L, color="grey", lw=0.8)
        a.set_xlabel("time"); a.set_ylabel("separation"); a.set_title(f"twin runs, lambda1 = {lam:.3f}")
        fig.tight_layout(); fig.savefig(plot, dpi=140)
        if not quiet:
            print(f"saved {plot}")
    return dict(label=label, lam=lam, horizon=horizon)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--preset", choices=PRESETS)
    ap.add_argument("--r", type=float)
    ap.add_argument("--x0", type=float, nargs=3, default=(1.0, 1.0, 1.0))
    ap.add_argument("--T", type=float, default=60.0, help="run length for panel 3 (default 60)")
    ap.add_argument("--csv", help="save the panel-3 time series (t,x,y,z) to this file")
    ap.add_argument("--plot", help="save a PNG of the 3-D path, x(t) and the twin separation")
    ap.add_argument("--all", action="store_true", help="table of all nine presets")
    ap.add_argument("--shoot", action="store_true", help="recompute the transition values by shooting")
    a = ap.parse_args()

    if a.shoot:
        tr = shoot_transitions()
        print("by shooting:", {k: round(v, 4) for k, v in tr.items()}, " published: 1, 13.926, 24.06, 24.7368")
    if a.all:
        print(f"{'preset':22s} {'r':>6s}  {'panel 3 label':48s} {'lambda1':>8s} {'horizon':>8s}")
        for p in PRESETS.values():
            res = run(p["r"], p["x0"], p["T"], quiet=True)
            h = f"{res['horizon']:.1f}" if np.isfinite(res["horizon"]) else "inf"
            print(f"{p['name']:22s} {p['r']:6g}  {res['label']:48s} {res['lam']:8.3f} {h:>8s}")
    if a.preset or a.r is not None:
        if a.preset:
            p = PRESETS[a.preset]
            r, x0, T = p["r"], p["x0"], p["T"]
        else:
            r, x0, T = a.r, a.x0, a.T
        run(r, x0, T, csv=a.csv, plot=a.plot)
    if not (a.shoot or a.all or a.preset or a.r is not None):
        ap.print_help()


if __name__ == "__main__":
    main()
