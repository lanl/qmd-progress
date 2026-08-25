#!/usr/bin/env python
"""
Faithful XLBO stability model for the kernelTimesRes integrator path.

This script reproduces, as exactly as possible, the recurrence actually executed
by `prg_xlbo_nint_kernelTimesRes` in `src/prg_xlbo_mod.F90` under the settings in
`examples/gpmdk/run/water/input.in` (XLBOLevel1=T, KernelType=ByParts). Earlier
stability scripts did NOT match the code:

  * variable_timestep_stability_corrected.py: fixed Verlet (2P_n - P_{n-1}) and
    plain kappa -- code uses variable Verlet (P_n=1+r, P_n1=r) and ka_scale*kappa.
  * find_max_alpha_per_pattern.py (source of alpha_max_table.txt): applied a
    1/(1+cc*kappa) denominator and a single coefficient set per cycle. The
    kernelTimesRes path has NO cc, NO denominator (restoring term is explicit in
    n_0 because KK0Res is computed from the current n before integration), and it
    ROLLS the 5-bit history window every step.

What the code does each step (src/prg_xlbo_mod.F90:475-483):

    call xlbo_compute_coeffs(...)   ! uses dt_history(1:5) = the FIVE PRIOR steps
    n = P_n*n_0 - P_n1*n_1 - ka_scale*kappa*KK0Res
          + alpha_use*(C0u*n_0 + C1u*n_1 + ... + C5u*n_5)
    shift history: n_5<-n_4 ... n_1<-n_0 n_0<-n
    push_dt_history(dt)             ! only AFTER integration

Scalar linear-response test (Niklasson Eq. 23): q[n_0]-n_0 = (gamma-1)*n_0, and
KK0Res -> -(q-n_0) in the scalar limit, so the restoring term contributes
ka_scale*kappa*(gamma-1) to the n_0 coefficient. gamma is swept over [-1,1]; a
pattern is stable iff max|eig(cycle matrix)| <= 1 across that sweep.

Run (use `python`, not python3):
    cd /Users/mewall/packages/gpmd/qmd-progress
    python scripts/Niklasson_JCP_2009_table_I/stability_faithful.py
"""

import numpy as np

# --------------------------------------------------------------------------
# Tables hard-coded from src/prg_xlbo_mod.F90 (branch split_step_alpha_ratefix).
# Indexed by 5-bit pattern: bit 0 = most recent dt, bit 4 = oldest dt.
# --------------------------------------------------------------------------
KAPPA = 1.82

C0 = np.array([
    -6.00000000, -0.85714286, -1.92857143, -0.32539683,
    -1.33333333,  0.08571429,  0.87500000,  0.40625000,
    22.00000000,  4.92857143, 16.00000000,  4.23809524,
    44.33333333,  9.57142857, 28.37500000,  7.47656250,
   -62.00000000,-10.50000000,-29.42857143, -6.58730159,
   -63.00000000,-11.34285714,-30.75000000, -7.11718750,
   -98.00000000,-14.71428571,-39.00000000, -7.83333333,
   -75.66666667,-11.80000000,-30.00000000, -6.00000000])

C1 = np.array([
    14.00000000,  3.00000000,  3.00000000,  0.52380952,
     1.33333333, -1.94285714, -2.12500000, -1.64062500,
   -64.00000000,-31.00000000,-31.00000000,-14.28571429,
  -110.33333333,-46.28571429,-50.62500000,-21.72265625,
   160.00000000, 53.00000000, 53.00000000, 19.23809524,
   147.00000000, 47.77142857, 52.25000000, 18.63671875,
   245.00000000, 68.00000000, 68.00000000, 21.00000000,
   170.66666667, 44.80000000, 49.00000000, 14.00000000])

C2 = np.array([
    -8.00000000, -0.57142857,  0.71428571,  2.05555556,
     1.66666667,  3.70000000,  3.06250000,  3.06250000,
    67.00000000, 47.28571429, 34.00000000, 26.66666667,
    79.33333333, 48.00000000, 32.81250000, 23.51562500,
  -133.00000000,-63.00000000,-40.28571429,-23.77777778,
   -94.00000000,-42.80000000,-27.12500000,-15.42187500,
  -192.00000000,-73.14285714,-44.00000000,-19.83333333,
  -102.66666667,-35.70000000,-21.00000000, -8.00000000])

C3 = np.array([
    -3.00000000, -4.57142857, -4.78571429, -5.25396825,
    -4.66666667, -4.84285714, -4.81250000, -4.82812500,
   -28.00000000,-24.21428571,-22.00000000,-19.61904762,
   -16.33333333,-14.28571429,-13.56250000,-12.26953125,
    32.00000000, 17.50000000, 13.71428571,  8.12698413,
     7.00000000,  3.37142857,  2.62500000,  0.90234375,
    42.00000000, 16.85714286, 12.00000000,  3.66666667,
     4.66666667, -0.30000000, -1.00000000, -3.00000000])

C4 = np.full(32, 4.0)
C5 = np.full(32, -1.0)

# Committed alpha (halved-even version) on branch split_step_alpha_ratefix.
ALPHA = np.array([
    0.036000, 0.186207, 0.069231, 0.129496,
    0.081818, 0.126162, 0.076378, 0.121484,
    0.004415, 0.012298, 0.012000, 0.015600,
    0.002267, 0.005701, 0.005744, 0.019081,
    0.003450, 0.018000, 0.005520, 0.021261,
    0.003273, 0.014477, 0.004672, 0.017476,
    0.002298, 0.011817, 0.003857, 0.016216,
    0.002531, 0.012420, 0.004320, 0.018000])

# Pre-halving alpha (before the rate fix) for comparison.
ALPHA_PREFIX = np.array([
    0.072000, 0.186207, 0.138462, 0.129496,
    0.163636, 0.126162, 0.152756, 0.121485,
    0.008830, 0.012305, 0.024000, 0.015603,
    0.004534, 0.005701, 0.011489, 0.019081,
    0.007714, 0.018000, 0.011043, 0.021260,
    0.006545, 0.014477, 0.009343, 0.017476,
    0.004596, 0.011816, 0.007714, 0.016216,
    0.005061, 0.012471, 0.008640, 0.018000])

GAMMA = np.linspace(-1.0, 1.0, 201)


def history_index(window5):
    """Port of get_K5_history_index. window5[k] = dt_history(k+1), i.e.
    window5[0] is the most recent prior step. Bit k set if dt ~ 1.0."""
    idx = 0
    for k in range(5):
        if abs(window5[k] - 1.0) < 0.1:
            idx |= (1 << k)
    return idx


def step_matrix(dt_n, dt_prev, gamma, idx, alpha_table):
    """6x6 transfer matrix for state [n_0, n_1, n_2, n_3, n_4, n_5].

    Verlet coefficients come from (dt_n, dt_prev); dissipation coefficients and
    alpha come from `idx` (the five-prior-step window), exactly as the code
    splits them in xlbo_compute_coeffs / prg_xlbo_nint_kernelTimesRes."""
    if dt_prev > 1e-12:
        r = dt_n / dt_prev
        P_n = 1.0 + r
        P_n1 = r
        ka_scale = 0.5 * dt_n * (dt_n + dt_prev)
    else:
        P_n = 2.0
        P_n1 = 1.0
        ka_scale = dt_n * dt_n

    a = alpha_table[idx]
    A = np.zeros((6, 6))
    A[0, 0] = P_n + ka_scale * KAPPA * (gamma - 1.0) + a * C0[idx]
    A[0, 1] = -P_n1 + a * C1[idx]
    A[0, 2] = a * C2[idx]
    A[0, 3] = a * C3[idx]
    A[0, 4] = a * C4[idx]
    A[0, 5] = a * C5[idx]
    for i in range(5):
        A[i + 1, i] = 1.0
    return A


def cycle_max_eig(dt_cycle, gamma, alpha_table):
    """Spectral radius of the product of per-step matrices over one period.

    dt_cycle is the periodic sequence of dt ratios. At cycle position j the
    current step is dt_cycle[j]; the five-prior window used for the index is
    dt_cycle[j-1], dt_cycle[j-2], ... (mod L), matching the code (coeffs looked
    up before the dt is pushed)."""
    L = len(dt_cycle)
    M = np.eye(6)
    for j in range(L):
        dt_n = dt_cycle[j]
        dt_prev = dt_cycle[(j - 1) % L]
        window5 = [dt_cycle[(j - 1 - k) % L] for k in range(5)]
        idx = history_index(window5)
        A = step_matrix(dt_n, dt_prev, gamma, idx, alpha_table)
        M = A @ M
    return np.max(np.abs(np.linalg.eigvals(M)))


def max_over_gamma(dt_cycle, alpha_table):
    return max(cycle_max_eig(dt_cycle, g, alpha_table) for g in GAMMA)


def pattern_to_cycle(pattern):
    """Length-5 periodic dt-ratio sequence implied by a 5-bit pattern.
    bit 0 (most recent) -> position 0. bit=1 -> full (1.0), bit=0 -> half (0.5)."""
    return [1.0 if (pattern >> k) & 1 else 0.5 for k in range(5)]


def find_alpha_margin(dt_cycle, alpha_table, hi=200.0):
    """Largest scalar m such that scaling ALL alpha entries by m keeps the cycle
    stable (max|eig| <= 1 over the gamma sweep). Bisection on m."""
    def rho(m):
        return max_over_gamma(dt_cycle, alpha_table * m) - 1.0
    if rho(1.0) > 0.0:
        # Already unstable at m=1: find m<1 boundary.
        lo, hiv = 0.0, 1.0
    else:
        lo, hiv = 1.0, hi
        if rho(hiv) <= 0.0:
            return hiv  # stable even at the cap; margin >= hi
    for _ in range(60):
        mid = 0.5 * (lo + hiv)
        if rho(mid) > 0.0:
            hiv = mid
        else:
            lo = mid
    return lo


# --------------------------------------------------------------------------
# 1. Sanity check: uniform full steps with the ORIGINAL fixed coefficients
#    (c = [-6,14,-8,-3,4,-1], alpha=0.018, Verlet 2/1) must be stable across
#    gamma in [-1,1]. Cross-checks the transfer-matrix wiring against the
#    classic test_original_coefficients result.
# --------------------------------------------------------------------------
def sanity_original():
    c = np.array([-6.0, 14.0, -8.0, -3.0, 4.0, -1.0])
    a = 0.018
    worst = 0.0
    for g in GAMMA:
        A = np.zeros((6, 6))
        A[0, 0] = 2.0 + KAPPA * (g - 1.0) + a * c[0]
        A[0, 1] = -1.0 + a * c[1]
        A[0, 2] = a * c[2]
        A[0, 3] = a * c[3]
        A[0, 4] = a * c[4]
        A[0, 5] = a * c[5]
        for i in range(5):
            A[i + 1, i] = 1.0
        worst = max(worst, np.max(np.abs(np.linalg.eigvals(A))))
    return worst


def main():
    print("=" * 78)
    print("FAITHFUL XLBO STABILITY MODEL  (kernelTimesRes path, kappa=%.2f)" % KAPPA)
    print("=" * 78)

    worst = sanity_original()
    status = "STABLE" if worst <= 1.0 else "UNSTABLE"
    print("\n[Sanity] Uniform full steps, original fixed coeffs c=[-6,14,-8,-3,4,-1],")
    print("         alpha=0.018:  max|eig| over gamma in [-1,1] = %.6f  -> %s" %
          (worst, status))
    print("         (must be STABLE to trust the transfer-matrix wiring)")

    # ----------------------------------------------------------------------
    # 2. All 32 patterns evaluated as genuine length-5 periodic cycles under
    #    the committed (halved-even) alpha, with the code's rolling index.
    # ----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("Committed alpha, each 5-bit pattern as a length-5 periodic cycle")
    print("(rolling window lookup, exactly as the code does):")
    print("-" * 78)
    print("%4s  %-9s  %10s  %8s  %8s   %s" %
          ("idx", "pattern", "alpha[idx]", "max|eig|", "margin", "status"))
    n_unstable = 0
    for p in range(32):
        cyc = pattern_to_cycle(p)
        rho = max_over_gamma(cyc, ALPHA)
        margin = find_alpha_margin(cyc, ALPHA)
        st = "OK" if rho <= 1.0 else "**UNSTABLE**"
        if rho > 1.0:
            n_unstable += 1
        bits = "".join(str((p >> k) & 1) for k in range(5))  # b0..b4
        print("%4d  %-9s  %10.6f  %8.5f  %8.2fx   %s" %
              (p, bits, ALPHA[p], rho, margin, st))

    # ----------------------------------------------------------------------
    # 3. Representative REAL adaptive cycles the mdloop actually produces:
    #    a stream of full steps with an occasional split (one full replaced by
    #    two half steps). These are the cases that matter in production.
    # ----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("Representative real adaptive cycles (committed alpha):")
    print("-" * 78)
    real_cycles = {
        "all-full [1]*1":          [1.0],
        "all-half [0.5]*1":        [0.5],
        "one split in 6 full":     [1.0]*6 + [0.5, 0.5],
        "one split in 10 full":    [1.0]*10 + [0.5, 0.5],
        "two splits in 10 full":   [1.0]*5 + [0.5, 0.5] + [1.0]*5 + [0.5, 0.5],
        "isolated split (5 full)": [1.0]*5 + [0.5, 0.5],
    }
    print("%-26s  %8s  %8s   %s" % ("cycle", "max|eig|", "margin", "status"))
    for name, cyc in real_cycles.items():
        rho = max_over_gamma(cyc, ALPHA)
        margin = find_alpha_margin(cyc, ALPHA)
        st = "OK" if rho <= 1.0 else "**UNSTABLE**"
        if rho > 1.0:
            n_unstable += 1
        print("%-26s  %8.5f  %8.2fx   %s" % (name, rho, margin, st))

    # ----------------------------------------------------------------------
    # 4. Compare committed (halved) vs pre-halving alpha on the split cycles,
    #    to see whether the rate fix changed the stability margin.
    # ----------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("Committed (halved-even) vs pre-halving alpha, split-containing cycles:")
    print("-" * 78)
    print("%-26s  %10s  %10s" % ("cycle", "rho(halved)", "rho(prefix)"))
    for name, cyc in real_cycles.items():
        rho_h = max_over_gamma(cyc, ALPHA)
        rho_p = max_over_gamma(cyc, ALPHA_PREFIX)
        print("%-26s  %10.5f  %10.5f" % (name, rho_h, rho_p))

    print("\n" + "=" * 78)
    if n_unstable == 0:
        print("RESULT: all committed-alpha cycles STABLE (max|eig| <= 1) over gamma[-1,1].")
    else:
        print("RESULT: %d cycle(s) UNSTABLE under committed alpha -- see ** above." %
              n_unstable)
    print("=" * 78)


if __name__ == "__main__":
    main()
