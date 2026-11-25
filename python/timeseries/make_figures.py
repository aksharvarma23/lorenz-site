import json, os, sys
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

OUT = sys.argv[1]; FIG = os.path.join(OUT, "figures"); os.makedirs(FIG, exist_ok=True)
meta = json.load(open(os.path.join(OUT, "metadata.json")))
runs = {r["id"]: r for r in meta["runs"]}
B = meta["b"]
DARK="#10243B"; ACC="#1E5F86"; RED="#C0392B"; GREY="#93A3B2"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"savefig.dpi":200,
                     "axes.spines.top":False,"axes.spines.right":False})

def load(rid):
    d = np.loadtxt(os.path.join(OUT, runs[rid]["file"]), delimiter=",", skiprows=1)
    return d[:, 0], d[:, 1], d[:, 2], d[:, 3]

def fp_marks(ax, r):
    ax.plot(0, 0, marker="x" if r > 1 else "o", color="k", ms=6, mew=1.5, ls="none", zorder=5)
    if r > 1:
        s = np.sqrt(B * (r - 1))
        stable = r < 24.7368
        for sg in (1, -1):
            ax.plot(sg * s, r - 1, marker="o" if stable else "x", color="k", ms=6, mew=1.5,
                    mfc="k" if stable else "none", ls="none", zorder=5)

def timeline(ax, x, z, cmap="Blues", lw=0.7):
    pts = np.column_stack([x, z]).reshape(-1, 1, 2)
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
    lc = LineCollection(seg, cmap=cmap, norm=plt.Normalize(-0.35, 1.0), lw=lw)
    lc.set_array(np.linspace(0, 1, len(seg))); ax.add_collection(lc); ax.autoscale()

# ---------------------------------------------------------------- Fig 1: phase portraits
panels = [("r00.5_conduction", "r = 0.5  ·  rest"), ("r10_steady_roll", "r = 10  ·  steady roll"),
          ("r20_transient_chaos", "r = 20  ·  transient chaos"), (None, "r = 24.4  ·  two starts, two fates"),
          ("r28_chaos", "r = 28  ·  strange attractor"), ("r100.5_periodic", "r = 100.5  ·  periodic, 3 peaks"),
          ("r160_periodic", "r = 160  ·  periodic, 2 peaks"), ("r350_periodic", "r = 350  ·  single limit cycle"),
          ("hero", "r = 28 in 3-D")]
fig = plt.figure(figsize=(11, 10.5))
for k, (rid, title) in enumerate(panels):
    if rid == "hero":
        ax = fig.add_subplot(3, 3, k + 1, projection="3d")
        t, x, y, z = load("r28_chaos"); n = 8000
        ax.plot(x[:n], y[:n], z[:n], lw=0.35, color=ACC)
        ax.set_xlabel("x", labelpad=-8); ax.set_ylabel("y", labelpad=-8); ax.set_zlabel("z", labelpad=-8)
        ax.tick_params(labelsize=6, pad=-3); ax.set_title(title, fontsize=10, color=DARK)
        continue
    ax = fig.add_subplot(3, 3, k + 1)
    if rid is None:
        t, x, y, z = load("r24.4_start_far");     timeline(ax, x, z, "Reds", 0.45)
        t, x, y, z = load("r24.4_start_near_C");  ax.plot(x, z, color=ACC, lw=2.6, zorder=7)
        r = 24.4
        sC = np.sqrt(B * (r - 1))
        ax.annotate("start near C+:\nspirals in, stays", xy=(sC, r - 1), xytext=(9.0, 3.0), fontsize=8,
                    color=ACC, ha="center", zorder=8,
                    arrowprops=dict(arrowstyle="->", color=ACC, lw=1.2),
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
        ax.text(0.02, 0.98, "start far away → chaos forever", color=RED, transform=ax.transAxes,
                va="top", fontsize=8, zorder=8,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
    else:
        t, x, y, z = load(rid); r = runs[rid]["r"]
        timeline(ax, x, z, "Reds" if r >= 24.06 and "periodic" not in rid else "Blues",
                 0.45 if rid == "r28_chaos" else 0.8)
        ax.plot(x[0], z[0], "o", color="#2E7D52", ms=4, zorder=6)
    fp_marks(ax, r)
    ax.set_title(title, fontsize=10, color=DARK); ax.set_xlabel("x"); ax.set_ylabel("z")
    ax.margins(0.08)
fig.suptitle("Phase portraits (x–z projection) across the r-axis  ·  σ = 10, b = 8/3\n"
             "colour runs light → dark with time; green dot = start; ● stable / × unstable fixed point",
             fontsize=11, color=DARK, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.965))
fig.savefig(os.path.join(FIG, "fig1_phase_portraits.png"), facecolor="white"); plt.close(fig)

# ---------------------------------------------------------------- Fig 2: time series
rows = [("r00.5_conduction", "r = 0.5"), ("r10_steady_roll", "r = 10"), ("r20_transient_chaos", "r = 20"),
        ("r24.4_start_near_C", "r = 24.4 (near C+)"), ("r24.4_start_far", "r = 24.4 (far start)"),
        ("r28_chaos", "r = 28"), ("r100.5_periodic", "r = 100.5"), ("r160_periodic", "r = 160"),
        ("r350_periodic", "r = 350")]
fig, axs = plt.subplots(len(rows), 1, figsize=(10, 12.5))
for ax, (rid, lab) in zip(axs, rows):
    t, x, y, z = load(rid)
    tmax = min(t[-1], 60.0); m = t <= tmax
    col = RED if ("far" in rid or rid == "r28_chaos") else ACC
    ax.plot(t[m], x[m], lw=0.8, color=col)
    ax.axhline(0, color=GREY, lw=0.5)
    ax.set_ylabel("x", rotation=0, labelpad=8)
    ax.text(1.005, 0.5, lab + "\n" + runs[rid]["observed_long_run"].split(" (")[0],
            transform=ax.transAxes, va="center", fontsize=8.5, color=DARK)
    ax.set_xlim(0, tmax)
axs[-1].set_xlabel("time (dimensionless)")
axs[2].annotate("roll reverses 15 times,\nthen settles at t ≈ 33", xy=(33, 0), xytext=(42, 12),
                fontsize=8, color=DARK, arrowprops=dict(arrowstyle="->", color=DARK, lw=0.8))
fig.suptitle("Generated time series, x(t), for each regime (first 60 time units shown)", fontsize=11, color=DARK)
fig.tight_layout(rect=(0, 0, 0.88, 0.98))
fig.savefig(os.path.join(FIG, "fig2_time_series.png"), facecolor="white"); plt.close(fig)

# ---------------------------------------------------------------- Fig 3: twin trajectories
d = np.loadtxt(os.path.join(OUT, "data", "lorenz_r28_twin_trajectories.csv"), delimiter=",", skiprows=1)
t, x1, x2, sep = d[:, 0], d[:, 1], d[:, 4], d[:, 7]
fig, (a, b) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
a.plot(t, x1, color=ACC, lw=1.1, label="run 1")
a.plot(t, x2, color=RED, lw=0.9, ls="--", label="run 2 (started 10⁻¹⁰ away)")
a.set_ylabel("x"); a.legend(fontsize=8, loc="upper left", ncol=2, frameon=False)
a.set_title("r = 28: two runs that start 10⁻¹⁰ apart agree, then part company", fontsize=10.5, color=DARK)
b.semilogy(t, sep, color=DARK, lw=1)
m = (t > 1) & (sep < 1e-2)
p = np.polyfit(t[m], np.log(sep[m]), 1)
b.semilogy(t[m], np.exp(np.polyval(p, t[m])), "--", color=RED, lw=1.4, label=f"fitted slope {p[0]:.3f}")
b.semilogy(t, 1e-10 * np.exp(0.906 * t), ":", color=GREY, lw=1.2, label="published 0.906")
b.set_ylim(1e-11, 1e3); b.set_xlabel("time"); b.set_ylabel("separation")
b.legend(fontsize=8, frameon=False, loc="lower right")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig3_twin_trajectories.png"), facecolor="white"); plt.close(fig)
json.dump({"twin_fit_slope": float(p[0])}, open(os.path.join(OUT, "figures", "twin_fit.json"), "w"))
print("figures done; twin slope", round(p[0], 4))
