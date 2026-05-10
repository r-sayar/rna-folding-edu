# Educational tools — Hölzer RNA lectures × research code

Small, self-contained demos that bridge the FU Berlin
"Algorithmische Bioinformatik / RNA analysis" lectures (Hölzer, WS 2024)
with research-grade RNA tools — the Rouskin Lab
(<https://github.com/rouskinlab>) for prediction / probing-ensembles, and
the Backofen Lab (<https://github.com/BackofenLab>) for visualization.

The goal is *pedagogical*: each demo implements the central algorithm
from scratch, then connects it to a production-grade tool that does the
same thing on real data.

## What's here

| File | Lecture concept | Bridge to research code |
|---|---|---|
| `nussinov_vs_efold.py` | RNA secondary-structure prediction, Nussinov DP (Lecture I, "RNA folding") | [`eFold`](https://github.com/rouskinlab/efold) — Evoformer-style deep learning predictor |
| `dreem_em_demo.py` | Mixture interpretation of probing data (extension of the lecture's "single structure" view) | [`DREEM`](https://github.com/rouskinlab/dreem) / [`SEISMIC-RNA`](https://github.com/rouskinlab/seismic-rna) — EM-based ensemble deconvolution |
| `visualizer/index.html` | Same Nussinov DP, *animated* | Uses [`fornac`](https://github.com/ViennaRNA/fornac) (Backofen Lab's [`vaRRI`](https://github.com/BackofenLab/vaRRI) ships it); pedagogy modelled on [`RNA-Playground`](https://github.com/BackofenLab/RNA-Playground) |

### Browser visualizer

`visualizer/index.html` is a single-page tool that walks through the
Nussinov DP step-by-step. The DP matrix fills diagonally as a heatmap,
then the traceback commits one base pair at a time — and after each
commit the structure on the right is re-rendered with fornac, so you
literally watch the strand fold from a straight line into a stem-loop.

![visualizer screenshot](visualizer/preview.png)

Open the file directly in a browser, or serve the directory:

```bash
cd visualizer
python3 -m http.server 8765
# then open http://localhost:8765/index.html
```

URL parameters: `?seq=GGGAAACCC&minloop=3&speed=20&autoplay=1`.

## Why this connection

The lectures end at "one sequence → one optimal structure" via Nussinov
(max base pairs) and Zuker (min free energy). The Rouskin Lab works on
the same primary task — RNA 2D structure — but pushes it in two modern
directions that the lectures don't cover:

1. **Replace the hand-engineered objective with a learned one.** That's
   eFold: an AlphaFold-Evoformer-inspired network trained on >300k
   structures. `nussinov_vs_efold.py` runs both on the same RNAs (a
   tiny hairpin, a hammerhead ribozyme from Lecture II, and a tRNA) and
   shows the dot-bracket strings side by side.

2. **Drop the "one structure" assumption entirely.** That's DREEM:
   an EM algorithm that treats each chemical-probing read as coming
   from one of K latent conformations and recovers all of them.
   `dreem_em_demo.py` builds a 2-conformation toy ensemble, simulates
   noisy DMS-MaPseq reads from a 50/50 mixture, throws away the
   population average, and lets EM pull both structures back out.

## Running

Both scripts are stdlib + numpy. eFold and matplotlib are optional —
the scripts degrade gracefully and tell you what to install.

```bash
# minimum
pip install numpy

# add eFold output to nussinov_vs_efold.py (pulls torch — heavy)
pip install efold

# add the convergence + recovered-profile plot to dreem_em_demo.py
pip install matplotlib

python nussinov_vs_efold.py
python dreem_em_demo.py
```

`dreem_em_demo.py` writes `dreem_em_demo.png` next to itself when
matplotlib is available.

## Reading order

1. Skim Lecture I sections "RNA secondary structure" and "RNA folding"
   (Nussinov / Zuker / mutual information).
2. Open `nussinov_vs_efold.py`. The `nussinov_fill` /
   `nussinov_traceback` pair is the lecture's algorithm in ~50 lines.
3. Run it. The output is two dot-bracket strings per RNA — read them
   against the structures the lecture draws.
4. Skim Lecture II sections that show real ribozyme structures from
   RNAfold / LocARNA. That's where the "one structure per sequence"
   assumption starts to break down.
5. Open `dreem_em_demo.py`. The `em()` function is DREEM's E-step /
   M-step in ~30 lines.
6. Run it. The plot shows that the population-average DMS profile is a
   blurred average of two distinct conformations, and that EM recovers
   both.

## What's *not* here (deliberately)

- A re-trained tiny eFold. The eFold model architecture is in
  `rouskinlab/efold` under `efold/models/`; training on a small subset
  needs a GPU and a few hours and is a separate exercise.
- Pseudoknots, multi-loop energetics, suboptimal structures. Nussinov
  doesn't model them; that's why Zuker / ViennaRNA exist. The lecture
  flags this.
- Real probing data. DREEM-on-real-data lives at
  <https://rnadreem.org> / `seismic-rna`; this demo uses synthetic
  reads so the ground truth is known and the EM behaviour is visible.

## Files

```
educational_tools/
├── README.md
├── nussinov_vs_efold.py      # ~210 LOC, stdlib + optional efold
├── dreem_em_demo.py          # ~210 LOC, numpy + optional matplotlib
└── visualizer/
    ├── index.html            # the page (own JS, ~400 lines inline)
    ├── preview.png           # screenshot of a folded 12-nt hairpin
    ├── NOTICE.md             # third-party attribution (fornac, RNA-Playground)
    └── vendor/fornac/        # Apache-2.0 — fornac.js + d3.js + LICENSE
```

## Credits

- [`fornac`](https://github.com/ViennaRNA/fornac) (Peter Kerpedjiev,
  Apache-2.0) — the RNA-structure renderer the visualizer uses; vendored
  from the Backofen Lab's [`vaRRI`](https://github.com/BackofenLab/vaRRI)
  distribution.
- [`RNA-Playground`](https://github.com/BackofenLab/RNA-Playground)
  (Backofen Lab) — the "DP matrix + structure side-by-side" teaching
  layout was the inspiration; the visualizer's JavaScript is original.
