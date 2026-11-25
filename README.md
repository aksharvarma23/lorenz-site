# Lorenz System Explorer (website)

A simple simulator for the BITS F316 project "Bifurcations and the Limits of Prediction in the Lorenz System"
(Akshar Ananya Varma, 2023B5A80799H).

    dx/dt = sigma (y - x),   dy/dt = r x - y - x z,   dz/dt = x y - b z,   sigma = 10, b = 8/3

Pick a heating level r (or one of nine presets) and a starting point. Four panels update, one per course learning outcome:

| Panel | Shows | Course outcome |
|---|---|---|
| 1. Fixed points and their stability | origin, C+, C-, Jacobian eigenvalues, stable / unstable | A: classifying fixed points |
| 2. Where your r sits among the bifurcations | bands with transitions at 1, 13.93, 24.06, 24.74 | E: bifurcation analysis |
| 3. Trajectory | 3-D path, x(t), steady / periodic / chaotic label, CSV download | C: solving ODEs on a computer |
| 4. Is it chaotic? | twin runs 1e-10 apart, lambda_1, prediction horizon | F: conditions for chaos |

The nine presets are the nine runs of the 11 November time-series data (r = 0.5, 10, 20, 24.4 near C+, 24.4 far, 28, 100.5, 160, 350).

It is a static site: plain HTML, CSS and JavaScript, no build step. All the maths runs in the visitor's browser.

## Deploying

This folder is a Git repository (GitHub `aksharvarma23/lorenz-site`) connected to Vercel: every push to `main`
redeploys the site. In GitHub Desktop: commit, then Push origin.

## Run locally

```bash
python -m http.server 8000     # inside this folder, then open http://localhost:8000
```

Opening index.html directly from the file system will not work: browsers block Web Workers on `file://` pages.

## What is inside

| Path | Contents |
|---|---|
| `index.html` | the page |
| `css/style.css` | styles |
| `js/lorenz.js` | numerical core: Dormand-Prince 5(4) solver (rtol 1e-10, atol 1e-12), fixed points and eigenvalues, z-maxima classifier, twin runs, tangent-linear Lyapunov exponent |
| `js/worker.js` | runs the heavy computations off the main thread |
| `js/app.js` | interface, presets, plots, CSV download |
| `vendor/plotly.min.js` | Plotly.js 2.35.2 (bundled so the site does not depend on a CDN) |
| `assets/` | the favicon |
| `python/` | the Python scripts (not deployed): see below |

## Python scripts

- `python/timeseries/`: the exact scripts used for the 11 November time-series progress update (data and figures).
- `python/website/lorenz_website.py`: the website's four calculations in Python, so every number on the page can be
  checked (`python lorenz_website.py --all`).

See `python/README.md`. `.vercelignore` keeps this folder off the website.
