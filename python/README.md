# Python scripts

The Python code behind this project's time-series data and the website.

```
pip install -r python/requirements.txt
```

## `timeseries/`: the 11 November time-series progress update

These are the exact scripts used to make the progress-update PDF and its data.

| File | What it does |
|---|---|
| `lorenz_timeseries.py` | Solves the Lorenz equations for the nine runs (r = 0.5, 10, 20, 24.4 near C+, 24.4 far, 28, 100.5, 160, 350) plus a twin pair at r = 28. Writes one CSV per run (`t,x,y,z`, every 0.01), `metadata.json` (start, fixed points, eigenvalues, lambda_1, observed behaviour) and `regime_summary.csv` |
| `make_figures.py` | Draws the three figures from that data: phase portraits, x(t) for each run, and the twin-trajectory plot |

```
cd python/timeseries
python lorenz_timeseries.py out      # data into out/data, out/metadata.json, out/regime_summary.csv
python make_figures.py out           # figures into out/figures
```

The output reproduces the PDF's data files byte for byte.

## `website/`: the website's calculations in Python

`lorenz_website.py` does the same four calculations as the site, so any number on the page can be checked:

| Site panel | Function |
|---|---|
| 1. Fixed points and their stability | `fixed_points()` |
| 2. Where your r sits among the bifurcations | `regime()`, `shoot_transitions()` |
| 3. Trajectory, label, CSV | `trajectory()`, `z_maxima()`, `periodicity()` |
| 4. Is it chaotic? | `twin_separation()`, `fit_slope()`, `lyapunov()` |

```
cd python/website
python lorenz_website.py --all                          # all nine presets, as on the site
python lorenz_website.py --preset butterfly --csv run.csv --plot run.png
python lorenz_website.py --r 24.4 --x0 1 1 1 --T 100     # any r and starting point
python lorenz_website.py --shoot                        # recompute 13.9266 and 24.059 by shooting
```

`--all` gives:

| Preset | r | Label | lambda_1 |
|---|---|---|---|
| Rest | 0.5 | steady | -0.476 |
| Steady roll | 10 | steady | -0.595 |
| Transient chaos | 20 | steady | -0.151 |
| Bistable: near C+ | 24.4 | steady (spiralling slowly onto C+-) | -0.010 |
| Bistable: far start | 24.4 | aperiodic (chaotic) | 0.784 |
| Butterfly | 28 | aperiodic (chaotic) | 0.908 |
| Periodic window | 100.5 | periodic, 3 peaks | -0.009 |
| Window, 2 peaks | 160 | periodic, 2 peaks | 0.001 |
| Single cycle | 350 | periodic, 1 peak | 0.001 |

The website runs the same methods in JavaScript (Dormand-Prince 5(4) instead of DOP853, same tolerances). The labels match exactly. The 200-unit lambda_1 averages agree to about 0.01, which is normal scatter for a finite average. The long-run published value at r = 28 is 0.906.

This folder is not deployed: `.vercelignore` keeps it off the website.
