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

## Pedagogical inspiration

The "DP matrix + structure side by side" teaching layout follows the
Backofen Lab's [RNA-Playground](https://github.com/BackofenLab/RNA-Playground)
(`Nussinov.html`). The JavaScript in `index.html` is original — written
from scratch — and does not derive from RNA-Playground's source.
