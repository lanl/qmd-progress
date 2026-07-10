#!/usr/bin/env python
"""Compute alpha values using formula: alpha = 0.018 * 3.0 / d_K, capped at alpha_max/2"""

import numpy as np

# d_K values from XLBO_K5_dK table
d_K = [
    0.75000000, 0.28947368, 0.38888889, 0.41666667,
    0.32926829, 0.42765957, 0.35294118, 0.44444444,
    6.11111111, 4.39090909, 2.25000000, 3.46153846,
    11.90476190, 6.13636364, 4.70000000, 2.82692308,
    7.82608696, 3.00000000, 4.89130435, 2.53846154,
    8.24137931, 3.72972973, 5.77777778, 3.08571429,
    11.74193548, 4.56756757, 7.00000000, 3.32954545,
    10.66666667, 4.34782609, 6.25000000, 3.00000000
]

# Alpha_max values from stability analysis
alpha_max = [
    0.144000, 0.372414, 0.276923, 0.258992,
    0.327273, 0.252323, 0.305512, 0.242969,
    0.017660, 0.024610, 0.048000, 0.031206,
    0.009067, 0.011402, 0.022979, 0.038162,
    0.015428, 0.036000, 0.022087, 0.042521,
    0.013091, 0.028954, 0.018685, 0.034951,
    0.009192, 0.023633, 0.015428, 0.032432,
    0.010122, 0.024941, 0.017280, 0.036000
]

base_alpha = 0.018
d_K_base = 3.0  # Pattern 31

print("Pattern | d_K      | formula_alpha | alpha_max | alpha_max/2 | Use (min)  | Capped?")
print("--------|----------|---------------|-----------|-------------|------------|--------")

alpha_values = []
for i in range(32):
    formula = base_alpha * d_K_base / d_K[i]
    half_max = alpha_max[i] / 2.0
    use_alpha = min(formula, half_max)
    capped = "YES" if formula > half_max else "NO"
    alpha_values.append(use_alpha)
    print(f"{i:7d} | {d_K[i]:8.5f} | {formula:13.6f} | {alpha_max[i]:9.6f} | {half_max:11.6f} | {use_alpha:10.6f} | {capped:7s}")

print("\n" + "="*80)
print("Fortran array for prg_xlbo_mod.F90:")
print("="*80)
print()
print("  real(dp), parameter :: XLBO_K5_alpha(0:31) = [ &")
for i in range(0, 32, 4):
    values = ", ".join([f"{alpha_values[j]:14.6f}_dp" for j in range(i, min(i+4, 32))])
    if i < 28:
        print(f"    {values}, &")
    else:
        print(f"    {values}]")
