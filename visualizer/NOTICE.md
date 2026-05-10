# NOTICE — third-party assets

This directory bundles third-party code so the visualizer works offline:

## fornac (RNA secondary-structure renderer)

Files: `vendor/fornac/fornac.js`, `vendor/fornac/fornac.css`, `vendor/fornac/d3.js`,
`vendor/fornac/LICENSE`, `vendor/fornac/README.md`.

- Copyright © 2015 Peter Kerpedjiev
- License: Apache License 2.0 (see `vendor/fornac/LICENSE`)
- Upstream: <https://github.com/ViennaRNA/fornac>
- Snapshot taken from the Backofen Lab's vaRRI distribution:
  <https://github.com/BackofenLab/vaRRI/tree/main/fornac>

These files are unmodified copies of the upstream releases.

## 3Dmol.js (WebGL molecular viewer for the experimental 3D panel)

Files: `vendor/3dmol/3Dmol-min.js`, `vendor/3dmol/LICENSE`.

- Copyright © 2014 University of Pittsburgh and contributors
- License: BSD-3-Clause (see `vendor/3dmol/LICENSE`)
- 3Dmol.js bundles code from GLmol, Three.js, and jQuery — all
  permissive licenses, attributions inside the same LICENSE file.
- Upstream: <https://github.com/3dmol/3Dmol.js>
- Snapshot: v2.5.4 from <https://unpkg.com/3dmol@2.5.4/build/3Dmol-min.js>

Used to fetch and render structures from the RCSB PDB
(<https://www.rcsb.org>) directly in the browser.

## Pedagogical inspiration

The "DP matrix + structure side by side" teaching layout follows the
Backofen Lab's [RNA-Playground](https://github.com/BackofenLab/RNA-Playground)
(`Nussinov.html`). The JavaScript in `index.html` is original — written
from scratch — and does not derive from RNA-Playground's source.
