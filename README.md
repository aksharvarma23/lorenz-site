# Lorenz System Explorer (website)

Interactive simulator for the BITS F316 project "Bifurcations and the Limits of Prediction in the Lorenz System"
(Akshar Ananya Varma, 2023B5A80799H).

It is a static site: plain HTML, CSS and JavaScript, with no build step. All the maths runs in the visitor's browser
(`js/lorenz.js`, in a Web Worker `js/worker.js`), so it needs no server.

## Deploy on Vercel

**Option 1: Vercel CLI (quickest)**
```bash
npm i -g vercel
cd lorenz-site
vercel            # first time: log in, accept the defaults
vercel --prod     # publishes to https://<project-name>.vercel.app
```

**Option 2: GitHub + Vercel dashboard**
1. Create a GitHub repository and upload the contents of this folder (index.html must be at the top level).
2. On vercel.com: Add New -> Project -> Import the repository.
3. Framework Preset: **Other**. Build Command: leave empty. Output Directory: leave empty (the root).
4. Deploy. Every push to the repository redeploys automatically.

## Run locally
```bash
cd lorenz-site
python3 -m http.server 8000     # then open http://localhost:8000
```
Opening index.html directly from the file system will not work, because browsers block Web Workers and
`fetch` on `file://` pages.

## What is inside
| Path | Contents |
|---|---|
| `index.html` | the page |
| `css/style.css` | styles (light and dark themes) |
| `js/lorenz.js` | numerical core: Dormand-Prince 5(4) integrator (rtol 1e-10, atol 1e-12), fixed points and eigenvalues, shooting for the homoclinic (13.9266) and attracting (24.059) values, tangent-linear Lyapunov exponent, z-maxima classifier |
| `js/worker.js` | runs the heavy computations off the main thread |
| `js/app.js` | interface and plots |
| `vendor/plotly.min.js` | Plotly.js 2.35.2 (bundled so the site does not depend on a CDN) |
| `assets/` | the favicon |

## Checks (browser results vs the Python project)
| Quantity | Browser | Python |
|---|---|---|
| Homoclinic explosion (shooting) | 13.9266 | 13.9266 |
| Chaos becomes attracting (shooting) | 24.059 | 24.059 |
| Same for sigma = 16, b = 4 | 19.377 / 32.360 | 19.38 / 32.36 |
| lambda_1 at r = 28 (200 time units) | 0.915 | 0.915 |
| lambda_1 at r = 20 (settles on C+) | -0.154 | -0.154 (eigenvalue -0.155) |
| r = 100.5 / 160 / 350 | period 3 / 2 / 1 | period 3 / 2 / 1 |

A shareable link keeps the settings in the address, for example `/#r=24.4&sigma=10&b=2.666667&x0=1,1,1&T=100`.
