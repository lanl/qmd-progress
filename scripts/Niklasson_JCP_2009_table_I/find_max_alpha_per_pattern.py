#!/usr/bin/env python
"""
Find maximum stable alpha for each pattern using pattern-specific coefficients.
"""

import numpy as np
from scipy.optimize import brentq

def construct_transfer_matrix(gamma, alpha, kappa_use, c, cc=1.0):
    """Construct transfer matrix for one timestep."""
    A = np.zeros((6, 6))
    denom = 1.0 + cc * kappa_use
    A[0, 0] = (2.0 + cc * kappa_use * (gamma - 1.0) + alpha * c[0]) / denom
    A[0, 1] = (-1.0 + alpha * c[1]) / denom
    A[0, 2:6] = alpha * c[2:6] / denom
    for i in range(5):
        A[i+1, i] = 1.0
    return A

def get_kappa_sequence(kappa_base, dt_pattern):
    """Compute kappa_use for each step using Verlet-consistent scaling."""
    kappa_seq = []
    for i in range(len(dt_pattern)):
        dt_n = dt_pattern[i]
        dt_prev = dt_pattern[i-1] if i > 0 else dt_pattern[-1]
        kappa_use = kappa_base * 0.5 * dt_n * (dt_n + dt_prev)
        kappa_seq.append(kappa_use)
    return kappa_seq

def max_rho(alpha, kappa_base, c, dt_pattern, cc=1.0):
    """Find maximum spectral radius over γ ∈ [-1, 1]."""
    gamma_vals = np.linspace(-1.0, 1.0, 201)
    kappa_seq = get_kappa_sequence(kappa_base, dt_pattern)

    max_eig = 0
    for gamma in gamma_vals:
        M = np.eye(6)
        for i, dt_n in enumerate(dt_pattern):
            M = construct_transfer_matrix(gamma, alpha, kappa_seq[i], c, cc) @ M
        max_eig = max(max_eig, np.max(np.abs(np.linalg.eigvals(M))))
    return max_eig

def find_max_alpha(kappa_base, c, dt_pattern, cc=1.0, alpha_low=1e-6, alpha_high=1.0):
    """Find maximum stable alpha for given pattern."""
    def objective(alpha):
        return max_rho(alpha, kappa_base, c, dt_pattern, cc) - 1.0

    obj_low = objective(alpha_low)
    obj_high = objective(alpha_high)

    if obj_low > 0:
        # Even minimum alpha is unstable
        return None, max_rho(alpha_low, kappa_base, c, dt_pattern, cc)
    if obj_high < 0:
        # Maximum alpha is still stable
        return alpha_high, max_rho(alpha_high, kappa_base, c, dt_pattern, cc)

    try:
        alpha_opt = brentq(objective, alpha_low, alpha_high, xtol=1e-8)
        return alpha_opt, 1.0
    except ValueError:
        return None, None

def pattern_to_dt(pattern_idx):
    """Convert pattern index (0-31) to timestep sequence."""
    dt_seq = []
    for i in range(5):
        bit = (pattern_idx >> i) & 1
        dt_seq.append(1.0 if bit else 0.5)
    return dt_seq

# Read coefficient table
coefficient_table = {}
with open('scripts/Niklasson_JCP_2009_table_I/K5_variable_timestep_coefficients_table.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('==') or line.startswith('Idx') or line.startswith('--'):
            continue
        if line.startswith('ENCODING') or line.startswith('Bit') or line.startswith('Pattern'):
            break

        parts = line.split()
        if len(parts) >= 10:
            try:
                idx = int(parts[0])
                c0 = float(parts[3])
                c1 = float(parts[4])
                c2 = float(parts[5])
                c3 = float(parts[6])
                c4 = float(parts[7])
                c5 = float(parts[8])
                d_K = float(parts[9])

                coefficient_table[idx] = {
                    'c': np.array([c0, c1, c2, c3, c4, c5]),
                    'd_K': d_K
                }
            except (ValueError, IndexError):
                continue

# Parameters
kappa_base = 1.82
cc = 0.99

print("=" * 80)
print("Maximum Stable α for Each Pattern with Pattern-Specific Coefficients")
print("=" * 80)
print()
print(f"κ_base = {kappa_base}")
print(f"cc = {cc}")
print()
print("Pattern | dt_history            | d_K      | α_max     | α_scaled  | Ratio")
print("--------|------------------------|----------|-----------|-----------|-------")
print("                                                            0.018*3/d_K   α_max/α_scaled")

alpha_values = []
d_K_values = []
paper_alpha = 0.018
d_K_31 = coefficient_table[31]['d_K']

for pattern_idx in range(32):
    dt_pattern = pattern_to_dt(pattern_idx)
    coeff_data = coefficient_table[pattern_idx]
    c = coeff_data['c']
    d_K = coeff_data['d_K']

    alpha_opt, max_eig = find_max_alpha(kappa_base, c, dt_pattern, cc)
    alpha_scaled = paper_alpha * d_K_31 / d_K

    dt_str = str(dt_pattern).replace("'", "")

    if alpha_opt is not None:
        ratio = alpha_opt / alpha_scaled
        alpha_values.append(alpha_opt)
        d_K_values.append(d_K)
        print(f"  {pattern_idx:2d}    | {dt_str:22s} | {d_K:8.4f} | {alpha_opt:9.6f} | {alpha_scaled:9.6f} | {ratio:6.2f}×")
    else:
        print(f"  {pattern_idx:2d}    | {dt_str:22s} | {d_K:8.4f} | {'UNSTABLE':>9s} | {alpha_scaled:9.6f} | {'N/A':>6s}")

print()
print("=" * 80)
print("Statistics")
print("=" * 80)

if alpha_values:
    print(f"Minimum α_max: {min(alpha_values):.6f} (pattern {alpha_values.index(min(alpha_values))})")
    print(f"Maximum α_max: {max(alpha_values):.6f} (pattern {alpha_values.index(max(alpha_values))})")
    print(f"Average α_max: {np.mean(alpha_values):.6f}")
    print(f"Std dev:       {np.std(alpha_values):.6f}")
    print()
    print(f"Paper uses α = {paper_alpha}")
    print(f"Minimum safety margin: {min(alpha_values)/paper_alpha:.1f}× (pattern {alpha_values.index(min(alpha_values))})")
    print()

    # Check correlation with d_K
    if len(d_K_values) > 1:
        correlation = np.corrcoef(d_K_values, alpha_values)[0, 1]
        print(f"Correlation between d_K and α_max: {correlation:.3f}")
        if abs(correlation) > 0.5:
            print(f"  {'Strong' if abs(correlation) > 0.7 else 'Moderate'} {'positive' if correlation > 0 else 'negative'} correlation")
        else:
            print(f"  Weak correlation - simple d_K scaling won't work")
else:
    print("No patterns are stable!")

print()
print("=" * 80)
print("Alpha scaling analysis")
print("=" * 80)

# Try to find a scaling that works
if alpha_values:
    # What constant would give us alpha = constant / d_K?
    # We want the constant such that min(constant / d_K) matches min(alpha_values)
    # That means: constant = min(alpha_values) * d_K[argmin(alpha_values)]
    min_idx = alpha_values.index(min(alpha_values))
    constant = alpha_values[min_idx] * d_K_values[min_idx]

    print(f"If we use α = {constant:.4f} / d_K:")

    # Check how many would be stable
    stable_count = 0
    for i, (alpha_max, d_K) in enumerate(zip(alpha_values, d_K_values)):
        alpha_scaled = constant / d_K
        if alpha_scaled <= alpha_max:
            stable_count += 1

    print(f"  {stable_count}/{len(alpha_values)} patterns would be stable")

    if stable_count == len(alpha_values):
        print(f"  ✓ This scaling works for all patterns!")
    else:
        print(f"  ✗ This scaling does not work for all patterns")
