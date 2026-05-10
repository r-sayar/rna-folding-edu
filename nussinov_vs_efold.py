"""
nussinov_vs_efold.py
====================

Bridges the FU Berlin RNA-analysis lectures (Hölzer, WS 2024) with the
Rouskin Lab's eFold predictor (https://github.com/rouskinlab/efold).

What it does
------------
1. Implements Nussinov's classic O(N^3) base-pair-maximisation algorithm
   from scratch (lecture I, "Comparative RNA analysis" / "RNA folding").
2. (Optional) Runs the Rouskin Lab's pretrained eFold model on the same
   sequences via `pip install efold`.
3. Prints both predicted dot-bracket structures side by side and computes
   simple base-pair F1 between them.

The point is *pedagogical*: the tractable lecture algorithm and a modern
deep-learning predictor are doing the same task (RNA sequence → 2D
structure), so seeing both outputs on a few small RNAs makes the
comparison concrete.

Run
---
    python nussinov_vs_efold.py                    # Nussinov only
    pip install efold && python nussinov_vs_efold.py   # adds eFold output

Dependencies
------------
- Python 3.9+. Stdlib only for the Nussinov part.
- Optional: `efold` (pulls torch).
"""

from __future__ import annotations

import textwrap
from dataclasses import dataclass

# ---------------------------------------------------------------------------
# 1. Nussinov algorithm  (lecture I: "RNA folding via maximising base pairs")
# ---------------------------------------------------------------------------

# Canonical base pairs from the lecture: Watson-Crick A-U, G-C, plus G-U
# wobble.  Anything else is treated as non-pairing.
CANONICAL_PAIRS = {("A", "U"), ("U", "A"),
                   ("G", "C"), ("C", "G"),
                   ("G", "U"), ("U", "G")}

# Minimum hairpin loop length: at least 3 unpaired bases between i and j.
MIN_LOOP = 3


def can_pair(a: str, b: str) -> bool:
    return (a, b) in CANONICAL_PAIRS


def nussinov_fill(seq: str) -> list[list[int]]:
    """Fill the N×N DP matrix M where M[i][j] = max # base pairs in seq[i..j]."""
    n = len(seq)
    M = [[0] * n for _ in range(n)]

    # diagonal & sub-diagonal stay 0 (can't pair within a window <= MIN_LOOP)
    for length in range(MIN_LOOP + 1, n + 1):       # window length
        for i in range(0, n - length + 1):
            j = i + length - 1
            # case 1: j unpaired
            best = M[i][j - 1]
            # case 2: j pairs with some k in [i, j-MIN_LOOP-1]
            for k in range(i, j - MIN_LOOP):
                if can_pair(seq[k], seq[j]):
                    left  = M[i][k - 1] if k > i else 0
                    right = M[k + 1][j - 1]
                    best  = max(best, left + 1 + right)
            M[i][j] = best
    return M


def nussinov_traceback(seq: str, M: list[list[int]]) -> list[tuple[int, int]]:
    """Recover one optimal pair list from the filled matrix."""
    pairs: list[tuple[int, int]] = []

    def trace(i: int, j: int) -> None:
        if j - i <= MIN_LOOP:
            return
        # Did j stay unpaired?
        if M[i][j] == M[i][j - 1]:
            trace(i, j - 1)
            return
        # Otherwise some k pairs with j.
        for k in range(i, j - MIN_LOOP):
            if not can_pair(seq[k], seq[j]):
                continue
            left  = M[i][k - 1] if k > i else 0
            right = M[k + 1][j - 1]
            if M[i][j] == left + 1 + right:
                pairs.append((k, j))
                if k > i:
                    trace(i, k - 1)
                trace(k + 1, j - 1)
                return

    trace(0, len(seq) - 1)
    return sorted(pairs)


def pairs_to_dotbracket(n: int, pairs: list[tuple[int, int]]) -> str:
    s = ["."] * n
    for i, j in pairs:
        s[i], s[j] = "(", ")"
    return "".join(s)


def nussinov(seq: str) -> tuple[str, list[tuple[int, int]]]:
    seq = seq.upper().replace("T", "U")
    M = nussinov_fill(seq)
    pairs = nussinov_traceback(seq, M)
    return pairs_to_dotbracket(len(seq), pairs), pairs


# ---------------------------------------------------------------------------
# 2. eFold wrapper  (Rouskin Lab — drop-in modern alternative)
# ---------------------------------------------------------------------------

def efold_predict(seq: str) -> str | None:
    """Return eFold's dot-bracket prediction, or None if eFold isn't installed."""
    try:
        from efold import inference          # type: ignore
    except Exception as exc:                 # ImportError or torch missing
        print(f"  [eFold not available: {exc}]")
        print(f"  [install with: pip install efold]")
        return None
    return inference(seq.upper().replace("T", "U"), fmt="dotbracket")


# ---------------------------------------------------------------------------
# 3. Comparison helpers
# ---------------------------------------------------------------------------

def dotbracket_to_pairs(db: str) -> set[tuple[int, int]]:
    stack: list[int] = []
    pairs: set[tuple[int, int]] = set()
    for i, c in enumerate(db):
        if c == "(":
            stack.append(i)
        elif c == ")":
            if stack:
                pairs.add((stack.pop(), i))
    return pairs


def pair_f1(db1: str, db2: str) -> float:
    p1, p2 = dotbracket_to_pairs(db1), dotbracket_to_pairs(db2)
    if not p1 and not p2:
        return 1.0
    tp = len(p1 & p2)
    if tp == 0:
        return 0.0
    precision = tp / len(p1)
    recall    = tp / len(p2)
    return 2 * precision * recall / (precision + recall)


# ---------------------------------------------------------------------------
# 4. Demo sequences
# ---------------------------------------------------------------------------

@dataclass
class Demo:
    name: str
    seq: str
    note: str


DEMOS = [
    Demo(
        name="tiny hairpin",
        seq="GGGAAACCC",
        note="textbook 9-mer, single stem-loop — sanity check",
    ),
    Demo(
        name="GC hairpin (lecture-style)",
        seq="CGCGAAUUCGCG",
        note="symmetric stem — should give a clean 4-bp hairpin",
    ),
    Demo(
        name="hammerhead ribozyme HQ880434.1 (lecture II, slide 4-5)",
        seq=("GCUUAUUUGGAAACGUCCAUGACGCUUAUCAUGUGUGAAGCCUGG"
             "CAAGGAGCCAAACGGUGUUAACACGUG"),
        note="74 nt env. hammerhead ribozyme, RNAfold output shown in lecture",
    ),
    Demo(
        name="yeast tRNA-Phe (a classic)",
        seq=("GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAG"
             "GUCCUGUGUUCGAUCCACAGAAUUCGCACCA"),
        note="76 nt tRNA, the cloverleaf the lecture promises on slide 6",
    ),
]


# ---------------------------------------------------------------------------
# 5. Pretty printer
# ---------------------------------------------------------------------------

def show(seq: str, db: str, label: str) -> None:
    print(f"  {label:>10}: {db}")


def main() -> None:
    print(textwrap.dedent("""\
        ============================================================
        Nussinov  vs.  eFold  (Rouskin Lab)
        Same task — sequence → 2D structure — two very different
        approaches: combinatorial DP vs deep learning.
        ============================================================
    """))

    for d in DEMOS:
        seq = d.seq.upper().replace("T", "U")
        print(f"--- {d.name}  ({len(seq)} nt) ---")
        print(f"  note: {d.note}")
        print(f"       seq: {seq}")

        nuss_db, nuss_pairs = nussinov(seq)
        show(seq, nuss_db, "Nussinov")
        print(f"  Nussinov base pairs: {len(nuss_pairs)}")

        efold_db = efold_predict(seq)
        if efold_db is not None:
            show(seq, efold_db, "eFold")
            f1 = pair_f1(nuss_db, efold_db)
            print(f"  base-pair F1 (Nussinov vs eFold) = {f1:.2f}")
        print()

    print(textwrap.dedent("""\
        Take-away
        ---------
        Nussinov optimises a single objective (max canonical pairs) and
        knows nothing about loop energetics, pseudoknots, or evolutionary
        signal.  eFold has *learned* a function on >300k structures so it
        can pick up patterns Nussinov misses (and sometimes hallucinate
        ones the data don't justify).  Both are doing the lecture's
        sequence → structure mapping; comparing them on the same input is
        the cheapest way to feel the difference.
    """))


if __name__ == "__main__":
    main()
