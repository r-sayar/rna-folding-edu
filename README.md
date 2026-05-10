# Educational tools — Hölzer RNA lectures × Rouskin Lab

Two small, self-contained Python scripts that bridge the FU Berlin
"Algorithmische Bioinformatik / RNA analysis" lectures (Hölzer, WS 2024)
with the Rouskin Lab's research code at <https://github.com/rouskinlab>.

The goal is *pedagogical*: each script implements the central algorithm
from scratch in a few hundred lines, then connects it to the
production-grade Rouskin Lab tool that does the same thing on real data.

## What's here

| File | Lecture concept | Rouskin Lab tool it mirrors |
|---|---|---|
| `nussinov_vs_efold.py` | RNA secondary-structure prediction, Nussinov DP (Lecture I, "RNA folding") | [`eFold`](https://github.com/rouskinlab/efold) — Evoformer-style deep learning predictor |
| `dreem_em_demo.py` | Mixture interpretation of probing data (extension of the lecture's "single structure" view) | [`DREEM`](https://github.com/rouskinlab/dreem) / [`SEISMIC-RNA`](https://github.com/rouskinlab/seismic-rna) — EM-based ensemble deconvolution |

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
└── dreem_em_demo.py          # ~210 LOC, numpy + optional matplotlib
```
