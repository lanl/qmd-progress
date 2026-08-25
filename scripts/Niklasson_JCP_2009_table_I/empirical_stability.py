#!/usr/bin/env python
"""
Empirical (Option B) stability validator for the XLBO kernelTimesRes recurrence
under QUASI-RANDOM adaptive timestep histories.

Why not periodic-orbit eigenvalues: computing a stable alpha cap per 5-bit pattern
by treating that pattern as an infinite periodic cycle is pathological. Several
"one split every 5 steps forever" orbits (idx 15,23,27,29,30 - cyclic rotations of
the same FFFFh orbit) show a ~0.9%/cycle variable-Verlet PARAMETRIC RESONANCE even
at alpha=0, so an alpha cap search returns 0 -- yet the committed alpha runs stably
in production. The resolution: real adaptive trajectories are quasi-random and never
lock into those exact periodic cadences. So we validate the way the physics actually
happens: propagate the exact recurrence over random dt histories at a realistic split
rate and measure the ACTUAL growth rate of the homogeneous solution.

Recurrence (kernelTimesRes path, scalar linear-response test, Niklasson Eq 23):
  coeffs/alpha looked up from the FIVE PRIOR steps (rolling window, pre-push order):
    n_{t+1} = P_n n_t - P_n1 n_{t-1}
              + ka_scale*kappa*(gamma-1)*n_t
              + alpha_idx*(C0 n_t + C1 n_{t-1} + ... + C5 n_{t-5})
  with P_n=1+r, P_n1=r, r=dt_t/dt_{t-1}, ka_scale=0.5*dt_t*(dt_t+dt_{t-1}).

Growth metric: renormalized product of per-step Jacobians (finite-sample largest
Lyapunov exponent) over a long random history. lambda > 0 => energy grows => unstable
for that (alpha table, gamma, split rate). We report growth per step and per unit
physical time, maximized over gamma in the physical band.

Run:  python scripts/Niklasson_JCP_2009_table_I/empirical_stability.py
"""

import numpy as np
import stability_faithful as sf
from stability_faithful import history_index, C0, C1, C2, C3, C4, C5

# committed (halved-even) alpha, for baseline comparison
ALPHA_COMMITTED = np.array([
    0.036000, 0.186207, 0.069231, 0.129496, 0.081818, 0.126162, 0.076378, 0.121484,
    0.004415, 0.012298, 0.012000, 0.015600, 0.002267, 0.005701, 0.005744, 0.019081,
    0.003450, 0.018000, 0.005520, 0.021261, 0.003273, 0.014477, 0.004672, 0.017476,
    0.002298, 0.011817, 0.003857, 0.016216, 0.002531, 0.012420, 0.004320, 0.018000])

dK = np.array([
    0.75, 0.28571429, 0.39285714, 0.17063492, 0.33333333, 0.06785714, 0.015625, 0.078125,
    2.0, 1.14285714, 2.25, 1.38095238, 5.08333333, 2.71428571, 4.703125, 2.83203125,
    7.0, 3.0, 4.89285714, 2.53968254, 8.25, 3.72857143, 5.78125, 3.08984375,
    11.75, 4.57142857, 7.0, 3.33333333, 10.66666667, 4.325, 6.25, 3.0])

DT_CUR = np.array([0.5 if (p & 1) == 0 else 1.0 for p in range(32)])
KAPPA = sf.KAPPA
REF_RATE = 0.054   # full-step reference dissipation per unit physical time


def make_history(n_steps, p_split, seed):
    """Quasi-random dt-ratio stream. With probability p_split a 'full step' slot is
    emitted as two half steps instead; otherwise a full step. p_split is per-slot,
    so the fraction of half-STEPS ~ 2*p_split/(1+p_split)."""
    rng = np.random.RandomState(seed)
    seq = []
    for _ in range(n_steps):
        if rng.random() < p_split:
            seq += [0.5, 0.5]
        else:
            seq += [1.0]
    return np.array(seq)


def step_jacobian(dt_n, dt_prev, gamma, idx, alpha_table):
    """6x6 companion Jacobian of one step of the recurrence (same as the transfer
    matrix in stability_faithful, parameterized by an arbitrary alpha table)."""
    if dt_prev > 1e-12:
        r = dt_n / dt_prev
        P_n = 1.0 + r; P_n1 = r; ka = 0.5 * dt_n * (dt_n + dt_prev)
    else:
        P_n = 2.0; P_n1 = 1.0; ka = dt_n * dt_n
    a = alpha_table[idx]
    A = np.zeros((6, 6))
    A[0, 0] = P_n + ka * KAPPA * (gamma - 1.0) + a * C0[idx]
    A[0, 1] = -P_n1 + a * C1[idx]
    A[0, 2] = a * C2[idx]; A[0, 3] = a * C3[idx]; A[0, 4] = a * C4[idx]; A[0, 5] = a * C5[idx]
    for i in range(5):
        A[i + 1, i] = 1.0
    return A


def lyapunov(dt_seq, gamma, alpha_table, warmup=50):
    """Finite-sample largest Lyapunov exponent (growth per step) of the homogeneous
    recurrence along dt_seq. Uses periodic renormalization of a single vector to
    avoid overflow; returns average log-growth per step after a warmup."""
    v = np.ones(6) / np.sqrt(6.0)
    total = 0.0
    count = 0
    win = [1.0] * 5   # rolling window of prior dt ratios, most-recent first
    N = len(dt_seq)
    for t in range(N):
        dt_n = dt_seq[t]
        dt_prev = win[0]
        idx = history_index(win)                 # coeffs from prior 5 steps
        A = step_jacobian(dt_n, dt_prev, gamma, idx, alpha_table)
        v = A @ v
        nrm = np.linalg.norm(v)
        if nrm > 0:
            v = v / nrm
            if t >= warmup:
                total += np.log(nrm)
                count += 1
        win = [dt_n] + win[:4]
    return total / max(count, 1)


def max_growth(alpha_table, p_split, n_steps=6000, seed=0,
               gammas=np.linspace(0.7, 1.0, 16), per_time=True):
    """Worst-case growth over the physical gamma band, for one random history.
    Returns (growth_per_step, growth_per_phys_time)."""
    dt_seq = make_history(n_steps, p_split, seed)
    T = dt_seq.sum()
    worst_step = -1e9
    for g in gammas:
        lam = lyapunov(dt_seq, g, alpha_table)
        if lam > worst_step:
            worst_step = lam
    # per physical time: multiply by (#steps / total physical time) = avg steps/fs
    per_t = worst_step * (len(dt_seq) / T) if per_time else worst_step
    return worst_step, per_t


def rule_alpha(ceiling):
    """General history-independent rule: alpha = min(rate-match, ceiling)."""
    ideal = REF_RATE * DT_CUR / dK
    return np.minimum(ideal, ceiling)


def main():
    print("=" * 76)
    print("EMPIRICAL STABILITY (Option B) - quasi-random histories, kappa=%.2f" % KAPPA)
    print("growth per step (Lyapunov); >0 means energy grows. avg over seeds.")
    print("=" * 76)

    SEEDS = range(6)
    # Production split rate: ~4500/20000 half-STEPS. fraction_half = 2p/(1+p) => p~0.29.
    rates = {"rare 1/20 (p=.05)": 0.05, "1/10 (p=.10)": 0.10,
             "production ~1/4.4 (p=.29)": 0.29, "heavy 1/2 (p=.50)": 0.50}

    # ---- 1. committed alpha baseline: is it stable across split rates? ----
    print("\n[1] Committed alpha table (baseline):")
    print("%-28s %14s %16s" % ("split rate", "growth/step", "growth/phys-time"))
    for name, p in rates.items():
        gs, gt = np.mean([max_growth(ALPHA_COMMITTED, p, seed=s) for s in SEEDS], axis=0)
        tag = "  STABLE" if gs <= 1e-6 else "  <-- GROWS"
        print("%-28s %14.2e %16.2e%s" % (name, gs, gt, tag))

    # ---- 2. sweep the global ceiling under the general rule at production rate ----
    print("\n[2] General rule alpha=min(0.054*dt/dK, ceiling), production split rate:")
    print("%-12s %14s %16s %s" % ("ceiling", "growth/step", "growth/phys-time", "status"))
    p_prod = 0.29
    best_ceiling = None
    for ceil in [0.02, 0.03, 0.05, 0.075, 0.10, 0.12, 0.15, 0.18, 0.20]:
        at = rule_alpha(ceil)
        gs, gt = np.mean([max_growth(at, p_prod, seed=s) for s in SEEDS], axis=0)
        status = "stable" if gs <= 1e-6 else "GROWS"
        if gs <= 1e-6:
            best_ceiling = ceil
        print("%-12.3f %14.2e %16.2e  %s" % (ceil, gs, gt, status))
    print("\nLargest stable ceiling at production rate: %s" %
          (("%.3f" % best_ceiling) if best_ceiling else "none in scan"))

    # ---- 3. does that ceiling stay stable across ALL split rates? ----
    if best_ceiling:
        print("\n[3] General rule with ceiling=%.3f across split rates:" % best_ceiling)
        print("%-28s %14s %16s" % ("split rate", "growth/step", "growth/phys-time"))
        at = rule_alpha(best_ceiling)
        for name, p in rates.items():
            gs, gt = np.mean([max_growth(at, p, seed=s) for s in SEEDS], axis=0)
            tag = "  STABLE" if gs <= 1e-6 else "  <-- GROWS"
            print("%-28s %14.2e %16.2e%s" % (name, gs, gt, tag))


if __name__ == "__main__":
    main()
