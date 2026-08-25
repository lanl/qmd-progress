#!/usr/bin/env python
"""
Numerical time-reversibility test for the XLBO backbone under (a) constant kappa
and (b) per-pattern kappa.

Motivation: energy conservation in XLBOMD rests on the near time-reversibility of
the Verlet backbone `n_{t+1} = P_n n_t - P_n1 n_{t-1} + ka_scale*kappa*(q-n_t)`.
The intentional alpha-dissipation aside (set alpha=0 here so we test the backbone
alone), a genuinely reversible integrator, run forward N steps and then backward
over the reversed timestep sequence, returns to its starting state to ~machine
precision. Any *return error* is irreversibility injected by the scheme.

We compare three setups, all with the SAME boundary handling so the common
start-of-trajectory artifact cancels in the B-vs-A and C-vs-B comparisons:

  A. uniform full steps, constant kappa           -> harness baseline
  B. variable steps (with splits), constant kappa -> cost of variable stepping
  C. variable steps, kappa depends on dt-pattern  -> the proposed per-pattern kappa

The scalar linear-response model q[n]-n = (gamma-1)*n is used (Niklasson Eq. 23);
reversibility is a structural property, so a linear surrogate is sufficient and
lets the return error be read cleanly.

Run:  python scripts/Niklasson_JCP_2009_table_I/reversibility_test.py
"""

import numpy as np

KAPPA_BASE = 1.82
GAMMA = 0.9          # a converged-SCF-like response (near 1); result is insensitive


def verlet_coeffs(dt_n, dt_prev):
    if dt_prev > 1e-12:
        r = dt_n / dt_prev
        return 1.0 + r, r, 0.5 * dt_n * (dt_n + dt_prev)
    return 2.0, 1.0, dt_n * dt_n


def history_index(window5):
    """bit k set if the k-th prior step is a full step (dt ~ 1.0)."""
    idx = 0
    for k in range(5):
        if abs(window5[k] - 1.0) < 0.1:
            idx |= (1 << k)
    return idx


def kappa_pattern(idx):
    """Example per-pattern kappa: soften the shadow oscillator for half-ending
    patterns (even index = most recent step was a half step). This is exactly the
    'reduce kappa across a split' idea whose reversibility cost we want to expose."""
    most_recent_is_full = (idx & 1) == 1
    return KAPPA_BASE if most_recent_is_full else 1.40


def propagate(dt_seq, kappa_mode):
    """Forward-propagate n over dt_seq, then run the same rule backward over the
    reversed sequence, and return the L2 return error ||n_backward_end - n_start||.

    kappa_mode: 'const' or 'pattern'. State history window is tracked so 'pattern'
    can look kappa up per the rolling 5-step pattern (as the code would)."""
    M = len(dt_seq)
    # --- forward pass ---
    x = np.zeros(M + 1)
    x[0] = 1.0
    # small perturbation so there is genuine dynamics to reverse
    x[0] = 1.0
    x[1] = 1.0 + 0.05        # a "kick": n_1 differs from n_0
    # rolling window of the last 5 dt ratios (most-recent first); unknown -> full
    win = [1.0, 1.0, 1.0, 1.0, 1.0]

    def kap(window):
        return KAPPA_BASE if kappa_mode == 'const' else kappa_pattern(history_index(window))

    for k in range(2, M + 1):
        dt_n = dt_seq[k - 1]
        dt_prev = dt_seq[k - 2]
        Pn, Pn1, ka = verlet_coeffs(dt_n, dt_prev)
        kk = kap(win)                       # window = the five steps prior to this one
        x[k] = Pn * x[k - 1] - Pn1 * x[k - 2] + ka * kk * (GAMMA - 1.0) * x[k - 1]
        win = [dt_n] + win[:4]              # push current step onto the window

    # --- backward pass: same update rule, reversed timestep ordering ---
    E = dt_seq[::-1]
    y = np.zeros(M + 1)
    y[0] = x[M]
    y[1] = x[M - 1]
    # reversed rolling window: the steps prior (in backward time) to the current one
    winb = [dt_seq[0]] * 5
    for j in range(2, M + 1):
        dt_n = E[j - 1]
        dt_prev = E[j - 2]
        Pn, Pn1, ka = verlet_coeffs(dt_n, dt_prev)
        kk = kap(winb)
        y[j] = Pn * y[j - 1] - Pn1 * y[j - 2] + ka * kk * (GAMMA - 1.0) * y[j - 1]
        winb = [dt_n] + winb[:4]

    return abs(y[M] - x[0]), x


def make_split_sequence(n_blocks=40, split_every=5):
    """A stream of full steps with a full step replaced by two halves periodically."""
    seq = []
    for b in range(n_blocks):
        if (b + 1) % split_every == 0:
            seq += [0.5, 0.5]     # split
        else:
            seq += [1.0]
    return seq


def main():
    print("=" * 74)
    print("XLBO backbone time-reversibility  (alpha=0, gamma=%.2f)" % GAMMA)
    print("return error = ||n after forward+backward pass  -  n_start||")
    print("=" * 74)

    N = 60
    uniform = [1.0] * N
    variable = make_split_sequence(n_blocks=48, split_every=5)   # ~9 splits

    errA, _ = propagate(uniform, 'const')
    errB, _ = propagate(variable, 'const')
    errC, _ = propagate(variable, 'pattern')

    print("\n%-42s %14s" % ("setup", "return error"))
    print("-" * 58)
    print("%-42s %14.3e" % ("A. uniform steps,   constant kappa", errA))
    print("%-42s %14.3e" % ("B. variable steps,  constant kappa", errB))
    print("%-42s %14.3e" % ("C. variable steps,  per-pattern kappa", errC))

    print("\nInterpretation:")
    print("  A  -> harness baseline (only the single start-of-trajectory")
    print("        boundary is asymmetric; ~machine precision if reversible).")
    print("  B/A ratio = %.3e : irreversibility from variable stepping alone." %
          (errB / errA if errA > 0 else float('inf')))
    print("  C/B ratio = %.3e : EXTRA irreversibility from making kappa" %
          (errC / errB if errB > 0 else float('inf')))
    print("        depend on the dt-pattern (the reversed trajectory sees a")
    print("        different pattern at each split, so a different kappa is")
    print("        applied -> the reversible backbone no longer inverts).")

    # sensitivity: how the gap grows with the number of splits
    print("\n" + "-" * 58)
    print("Return error vs number of splits (constant vs per-pattern kappa):")
    print("%8s %16s %16s" % ("#splits", "B const-kappa", "C pattern-kappa"))
    for split_every in [12, 8, 6, 4, 3]:
        seq = make_split_sequence(n_blocks=48, split_every=split_every)
        nsp = sum(1 for i in range(0, len(seq) - 1) if seq[i] == 0.5 and seq[i + 1] == 0.5)
        eB, _ = propagate(seq, 'const')
        eC, _ = propagate(seq, 'pattern')
        print("%8d %16.3e %16.3e" % (nsp, eB, eC))


if __name__ == "__main__":
    main()
