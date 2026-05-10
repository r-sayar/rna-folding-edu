"""
dreem_em_demo.py
================

A from-scratch, fully synthetic demo of DREEM
(Detection of RNA folding Ensembles using Expectation-Maximisation,
Tomezsko et al. 2020, Rouskin Lab — https://github.com/rouskinlab/dreem).

Why this exists
---------------
The FU Berlin lectures (Hölzer, WS 2024) cover RNA secondary-structure
*prediction* — Nussinov, Zuker, comparative analysis.  They assume each
sequence has *one* structure.  In reality the same RNA often folds into
several distinct conformations and chemical-probing reads come from a
mixture.  DREEM is the algorithm that pulls the mixture apart.

This script:
  1. Simulates two RNA conformations and the per-base mutation rates
     they would produce in a DMS-MaPseq experiment (paired bases mutate
     rarely; unpaired bases mutate often).
  2. Generates N noisy reads from a 50/50 mixture of those two
     conformations.
  3. Throws the population-average profile away and runs 2-cluster EM
     from random init to recover both conformations and their
     proportions.
  4. Plots log-likelihood convergence and recovered vs true profiles.

If matplotlib isn't installed, the script still prints the recovered
profiles; install matplotlib to also get the figure.

Run
---
    python dreem_em_demo.py
    pip install matplotlib   # optional, for the plot

Concepts mapped
---------------
  - "structure ensemble"            -> two true mu vectors
  - "DMS-MaPseq read"               -> Bernoulli sample over bases
  - "population-average profile"    -> mean of all reads (loses ensemble info)
  - "EM E-step"                     -> posterior cluster weight per read
  - "EM M-step"                     -> mixing prop & mu per cluster
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass

import numpy as np

# ---------------------------------------------------------------------------
# 1. Toy ground truth: two structures over the same 30-nt RNA
# ---------------------------------------------------------------------------
#
# Each structure is encoded as a binary mask of length L:
#     1 -> paired (DMS protected, low mutation rate)
#     0 -> unpaired (DMS reactive, high mutation rate)
#
# The two structures share some bases but differ in the middle stem,
# which is exactly the situation DREEM was built for.

L = 30

#                                111111111122222222223
#                      0123456789012345678901234567890
STRUCT_A = np.array([0,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,0])
STRUCT_B = np.array([0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0])

assert STRUCT_A.shape == (L,) and STRUCT_B.shape == (L,)

# DMS mutation rate per base, conditional on conformation:
MU_PAIRED   = 0.01      # ~1% — base is protected
MU_UNPAIRED = 0.10      # ~10% — base is exposed and reactive

def struct_to_mu(struct: np.ndarray) -> np.ndarray:
    return np.where(struct == 1, MU_PAIRED, MU_UNPAIRED)

MU_A_TRUE = struct_to_mu(STRUCT_A)
MU_B_TRUE = struct_to_mu(STRUCT_B)


# ---------------------------------------------------------------------------
# 2. Simulate a sequencing experiment: 50/50 mixture, N reads
# ---------------------------------------------------------------------------

@dataclass
class Sim:
    reads: np.ndarray   # shape (N, L), 0/1
    truth: np.ndarray   # shape (N,), 0 or 1 — which conformation each read came from


def simulate(n_reads: int = 2000, frac_a: float = 0.5, seed: int = 0) -> Sim:
    rng = np.random.default_rng(seed)
    truth = (rng.random(n_reads) >= frac_a).astype(int)   # 0 = A, 1 = B
    mu = np.where(truth[:, None] == 0, MU_A_TRUE, MU_B_TRUE)  # (N, L)
    reads = (rng.random((n_reads, L)) < mu).astype(np.uint8)
    return Sim(reads=reads, truth=truth)


# ---------------------------------------------------------------------------
# 3. EM for a mixture of K independent-Bernoulli "profiles"
# ---------------------------------------------------------------------------
#
# Per cluster k we have a vector mu_k in [0,1]^L.
# Likelihood of read r under cluster k:
#     P(r | k) = prod_i mu_k[i]^r[i] * (1-mu_k[i])^(1-r[i])
#
# E-step: posterior weight w[r,k] = pi_k * P(r|k) / sum_j pi_j * P(r|j)
# M-step: pi_k        = mean_r w[r,k]
#         mu_k[i]     = sum_r w[r,k] * r[i] / sum_r w[r,k]
#
# This is exactly what DREEM does for the structure-ensemble problem,
# minus the bells & whistles (per-read coverage masks, EM restarts,
# BIC for choosing K, base-quality filters, etc.).

def em_log_likelihood(reads: np.ndarray, mu: np.ndarray, pi: np.ndarray) -> float:
    """Total log-likelihood under the current params (used for convergence)."""
    # log P(r | k) computed in log space to avoid underflow
    eps = 1e-9
    log_mu      = np.log(mu      + eps)         # (K, L)
    log_one_mu  = np.log(1 - mu  + eps)
    # (N, K) = reads (N,L) @ log_mu.T  + (1-reads) @ log_one_mu.T
    log_pr = reads @ log_mu.T + (1 - reads) @ log_one_mu.T
    log_pr += np.log(pi + eps)                  # add log mixing prop
    # logsumexp over K, then sum over N
    m = log_pr.max(axis=1, keepdims=True)
    return float(np.sum(m.squeeze(1) + np.log(np.exp(log_pr - m).sum(axis=1))))


def em(reads: np.ndarray, k: int = 2, max_iter: int = 100,
       tol: float = 1e-4, seed: int = 1) -> dict:
    rng = np.random.default_rng(seed)
    n, length = reads.shape
    # init: random mu in [0.02, 0.15], uniform pi
    mu = rng.uniform(0.02, 0.15, size=(k, length))
    pi = np.full(k, 1.0 / k)

    history = []
    prev_ll = -math.inf
    for it in range(max_iter):
        # E-step  ---------------------------------------------------------
        eps = 1e-9
        log_mu     = np.log(mu     + eps)
        log_one_mu = np.log(1 - mu + eps)
        log_pr = reads @ log_mu.T + (1 - reads) @ log_one_mu.T   # (N, K)
        log_pr += np.log(pi + eps)
        m = log_pr.max(axis=1, keepdims=True)
        log_pr -= m
        w = np.exp(log_pr)
        w /= w.sum(axis=1, keepdims=True)                        # (N, K)

        # M-step  ---------------------------------------------------------
        nk = w.sum(axis=0)              # (K,)
        pi = nk / n
        mu = (w.T @ reads) / nk[:, None]
        mu = np.clip(mu, 1e-4, 1 - 1e-4)

        # convergence check
        ll = em_log_likelihood(reads, mu, pi)
        history.append(ll)
        if abs(ll - prev_ll) < tol:
            break
        prev_ll = ll

    return dict(mu=mu, pi=pi, w=w, ll_history=history, iterations=len(history))


# ---------------------------------------------------------------------------
# 4. Match recovered clusters back to the true ones (label permutation)
# ---------------------------------------------------------------------------

def best_permutation(mu_recovered: np.ndarray) -> tuple[int, int]:
    """Pick the assignment of (recovered cluster -> true A/B) that minimises L2 distance."""
    a, b = mu_recovered
    d_aA = np.linalg.norm(a - MU_A_TRUE) + np.linalg.norm(b - MU_B_TRUE)
    d_aB = np.linalg.norm(a - MU_B_TRUE) + np.linalg.norm(b - MU_A_TRUE)
    return (0, 1) if d_aA <= d_aB else (1, 0)


# ---------------------------------------------------------------------------
# 5. Reporting
# ---------------------------------------------------------------------------

def print_profile(name: str, mu: np.ndarray) -> None:
    bar = "".join("#" if v > 0.05 else "." for v in mu)
    nums = " ".join(f"{v:.2f}" for v in mu)
    print(f"  {name:>14}: |{bar}|")
    print(f"  {'rates':>14}: {nums}")


def maybe_plot(sim: Sim, result: dict, perm: tuple[int, int]) -> str | None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"  [matplotlib not available: {exc}; skipping plot]")
        return None

    pop_avg = sim.reads.mean(axis=0)
    rec_a = result["mu"][perm[0]]
    rec_b = result["mu"][perm[1]]

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), constrained_layout=True)
    x = np.arange(L)

    axes[0].plot(result["ll_history"], "-o")
    axes[0].set_title(f"EM log-likelihood ({result['iterations']} iters)")
    axes[0].set_xlabel("iteration")
    axes[0].set_ylabel("log-likelihood")

    width = 0.18
    axes[1].bar(x - 1.5 * width, MU_A_TRUE, width, label="true A",       color="#1f77b4")
    axes[1].bar(x - 0.5 * width, rec_a,     width, label="recovered A",  color="#aec7e8")
    axes[1].bar(x + 0.5 * width, MU_B_TRUE, width, label="true B",       color="#d62728")
    axes[1].bar(x + 1.5 * width, rec_b,     width, label="recovered B",  color="#ff9896")
    axes[1].plot(x, pop_avg, "k--", lw=1, label="population avg (what you'd get without DREEM)")
    axes[1].set_title("Per-base DMS mutation rate: truth vs EM-recovered")
    axes[1].set_xlabel("base index")
    axes[1].set_ylabel("mutation rate")
    axes[1].legend(fontsize=8, ncols=3)

    out = "dreem_em_demo.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# 6. Driver
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 60)
    print(" DREEM mini: recovering 2 RNA conformations from a")
    print(" mixture of DMS-MaPseq reads, using Expectation-Maximisation.")
    print("=" * 60)
    print()

    sim = simulate(n_reads=2000, frac_a=0.5, seed=42)
    print(f"Simulated {sim.reads.shape[0]} reads of length {L}.")
    print(f"True mixing: A = {1 - sim.truth.mean():.2f}, "
          f"B = {sim.truth.mean():.2f}\n")

    print_profile("true A",  MU_A_TRUE)
    print_profile("true B",  MU_B_TRUE)
    print_profile("pop. avg.", sim.reads.mean(axis=0))
    print("\n(notice: the population average looks like a single, blurred profile.")
    print(" Without EM you cannot see that it's two structures.)\n")

    result = em(sim.reads, k=2, max_iter=200, seed=1)
    perm = best_permutation(result["mu"])

    print(f"EM converged in {result['iterations']} iterations.")
    print(f"Recovered mixing: A = {result['pi'][perm[0]]:.2f}, "
          f"B = {result['pi'][perm[1]]:.2f}\n")
    print_profile("recovered A", result["mu"][perm[0]])
    print_profile("recovered B", result["mu"][perm[1]])

    err_a = np.linalg.norm(result["mu"][perm[0]] - MU_A_TRUE)
    err_b = np.linalg.norm(result["mu"][perm[1]] - MU_B_TRUE)
    print(f"\nL2 error vs truth: A = {err_a:.3f}, B = {err_b:.3f}")

    out = maybe_plot(sim, result, perm)
    if out:
        print(f"\nPlot written to: {out}")

    print("\nTake-away")
    print("---------")
    print("Population-average DMS rates hide structural heterogeneity.")
    print("EM (DREEM) treats each read as coming from one of K latent")
    print("conformations and alternates between guessing assignments and")
    print("re-fitting per-conformation rates.  With enough reads it pulls")
    print("the original structures back out — which is how the Rouskin lab")
    print("discovered alternative folds in HIV, SARS-CoV-2, pri-miRNAs, etc.")


if __name__ == "__main__":
    main()
