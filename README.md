# ChemSphere

An interactive chemistry education explorer built with Next.js, React, TypeScript and CSS. It includes a compound library, browser-side calculations, reaction balancing demo, periodic table, molecular geometry studio and learning modules.

## Run locally

```bash
npm install
npm run dev
```

Open http://localhost:3000.

## Optional Python scientific service

The frontend remains a Next.js app. The companion FastAPI service adds validated
RDKit structure parsing and 2D coordinates, SymPy equation balancing, NumPy/SciPy
kinetics, and Pint unit conversion.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r api/requirements.txt
python -m uvicorn api.main:app --reload --port 8000
```

The service endpoints are `/api/structure`, `/api/balance`, `/api/units`, and
`/api/kinetics/first-order`. RDKit is used only for validated molecular inputs;
educational UI fallbacks remain available when the service is offline.

## Notes

This is an educational application using local sample data. It is not a substitute for laboratory safety procedures or professional chemical analysis.
