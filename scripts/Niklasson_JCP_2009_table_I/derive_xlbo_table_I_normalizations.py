#!/usr/bin/env python3
"""
Derive XLBO dissipation coefficients for K=3 through K=10
Using the normalizations from Table I in Niklasson 2009

Normalization pattern:
  c_K = (-1)^K
  c_{K-1} = 2*(K-3)*(-1)^(K-1)

This gives:
K=3: c_2=2*(3-3)*(-1)^2=0,    c_3=(-1)^3=-1
K=4: c_3=2*(4-3)*(-1)^3=-2,   c_4=(-1)^4=1
K=5: c_4=2*(5-3)*(-1)^4=4,    c_5=(-1)^5=-1
K=6: c_5=2*(6-3)*(-1)^5=-6,   c_6=(-1)^6=1
K=7: c_6=2*(7-3)*(-1)^6=8,    c_7=(-1)^7=-1
K=8: c_7=2*(8-3)*(-1)^7=-10,  c_8=(-1)^8=1
K=9: c_8=2*(9-3)*(-1)^8=12,   c_9=(-1)^9=-1
K=10: c_9=2*(10-3)*(-1)^9=-14, c_10=(-1)^10=1
"""

import numpy as np

# Table I normalizations (last two coefficients fixed)
# Pattern: c_K = (-1)^K, c_{K-1} = 2*(K-3)*(-1)^(K-1)
TABLE_I_NORMALIZATIONS = {
    3: (0, -1),      # c_2=2*(3-3)*(-1)^2=0, c_3=(-1)^3=-1
    4: (-2, 1),      # c_3=2*(4-3)*(-1)^3=-2, c_4=(-1)^4=1
    5: (4, -1),      # c_4=2*(5-3)*(-1)^4=4, c_5=(-1)^5=-1
    6: (-6, 1),      # c_5=2*(6-3)*(-1)^5=-6, c_6=(-1)^6=1
    7: (8, -1),      # c_6=2*(7-3)*(-1)^6=8, c_7=(-1)^7=-1
    8: (-10, 1),     # c_7=2*(8-3)*(-1)^7=-10, c_8=(-1)^8=1
    9: (12, -1),     # c_8=2*(9-3)*(-1)^8=12, c_9=(-1)^9=-1
    10: (-14, 1),    # c_9=2*(10-3)*(-1)^9=-14, c_10=(-1)^10=1
}

# Table I reference values for verification
TABLE_I_REFERENCE = {
    3: np.array([-2, 3, 0, -1]),
    4: np.array([-3, 6, -2, -2, 1]),
    5: np.array([-6, 14, -8, -3, 4, -1]),
    6: np.array([-14, 36, -27, -2, 12, -6, 1]),
    7: np.array([-36, 99, -88, 11, 32, -25, 8, -1]),
    8: np.array([-99, 286, -286, 78, 78, -90, 42, -10, 1]),
    9: np.array([-286, 858, -936, 364, 168, -300, 184, -63, 12, -1]),
}

def derive_K(K):
    """
    Derive coefficients for given K using Table I normalization

    Returns: coeffs (array of K+1 values), d_K
    """
    # Get normalization: c_K = (-1)^K, c_{K-1} = 2*(K-3)*(-1)^(K-1)
    if K not in TABLE_I_NORMALIZATIONS:
        print(f"Warning: No Table I normalization for K={K}, using pattern")
        c_Kminus1_fixed = float(2 * (K - 3) * ((-1) ** (K - 1)))
        c_K_fixed = float((-1) ** K)
    else:
        c_Kminus1_fixed, c_K_fixed = TABLE_I_NORMALIZATIONS[K]
        c_Kminus1_fixed = float(c_Kminus1_fixed)
        c_K_fixed = float(c_K_fixed)

    # Number of unknowns: c_0 through c_{K-2}
    n_unknowns = K - 1

    # Build constraint equations
    A = []
    b = []

    # Moment 0: Σc_k = 0
    row = np.ones(n_unknowns)
    A.append(row)
    b.append(-(c_Kminus1_fixed + c_K_fixed))

    # Moment 1: Σk·c_k = 0
    row = np.array([k for k in range(n_unknowns)])
    A.append(row)
    b.append(-((K-1) * c_Kminus1_fixed + K * c_K_fixed))

    # Odd moments: 3, 5, 7, ..., until we have enough equations
    m = 3
    while len(A) < n_unknowns:
        row = np.array([k**m for k in range(n_unknowns)])
        A.append(row)
        b.append(-((K-1)**m * c_Kminus1_fixed + K**m * c_K_fixed))
        m += 2  # Next odd moment

    A = np.array(A)
    b = np.array(b)

    # Solve
    solution = np.linalg.solve(A, b)

    # Assemble full coefficients
    coeffs = np.zeros(K+1)
    coeffs[:n_unknowns] = solution
    coeffs[K-1] = c_Kminus1_fixed
    coeffs[K] = c_K_fixed

    # Calculate d_K from moment 2
    moment2 = sum(k**2 * coeffs[k] for k in range(K+1))
    d_K = moment2 / 2

    return coeffs, d_K

def verify_moments(K, coeffs, d_K):
    """Verify that coefficients satisfy moment constraints"""
    print(f"\n  Moment verification:")

    # Moment 0
    m0 = sum(coeffs)
    status = "✓" if abs(m0) < 1e-10 else "✗"
    print(f"    {status} Moment 0: {m0:.6e} (should be 0)")

    # Moment 1
    m1 = sum(k * coeffs[k] for k in range(K+1))
    status = "✓" if abs(m1) < 1e-10 else "✗"
    print(f"    {status} Moment 1: {m1:.6e} (should be 0)")

    # Moment 2
    m2 = sum(k**2 * coeffs[k] for k in range(K+1))
    expected_m2 = 2 * d_K
    status = "✓" if abs(m2 - expected_m2) < 1e-10 else "✗"
    print(f"    {status} Moment 2: {m2:.6e} (should be {expected_m2:.6e})")

    # Odd moments up to 2K-3
    for m in range(3, 2*K-2, 2):
        moment = sum(k**m * coeffs[k] for k in range(K+1))
        status = "✓" if abs(moment) < 1e-6 else "✗"
        print(f"    {status} Moment {m}: {moment:.6e} (should be 0)")

print("="*80)
print("XLBO Dissipation Coefficients: K=3 through K=10")
print("Using Table I normalizations from Niklasson 2009")
print("="*80)

all_coeffs = {}

for K in [3, 4, 5, 6, 7, 8, 9, 10]:
    print(f"\n{'='*80}")
    print(f"K = {K} ({K+1} points)")
    print(f"{'='*80}")

    coeffs, d_K = derive_K(K)
    all_coeffs[K] = (coeffs, abs(d_K))  # Use absolute value of d_K

    if K in TABLE_I_NORMALIZATIONS:
        c_Kminus1, c_K = TABLE_I_NORMALIZATIONS[K]
        print(f"\nTable I normalization: c_{K-1}={c_Kminus1}, c_{K}={c_K}")

    print(f"\nCoefficients:")
    for k, c in enumerate(coeffs):
        marker = " *" if k in [K-1, K] else ""
        print(f"  c_{k} = {c:8.1f}{marker}")
    print(f"\n  d_K = {abs(d_K):.1f}")

    verify_moments(K, coeffs, d_K)

    # Check against Table I if available
    if K in TABLE_I_REFERENCE:
        ref = TABLE_I_REFERENCE[K]
        if np.allclose(coeffs, ref, atol=0.01):
            print(f"\n  ✓✓✓ MATCHES Table I! ✓✓✓")
        else:
            print(f"\n  ✗ DIFFERS from Table I")
            print(f"\n  Table I: c = {ref}")
            print(f"  Ours:    c = {coeffs.astype(int)}")
            diff = coeffs - ref
            max_diff_idx = np.argmax(np.abs(diff))
            print(f"  Max diff: c_{max_diff_idx} = {diff[max_diff_idx]:.6f}")

# Print Fortran implementation
print("\n" + "="*80)
print("FORTRAN IMPLEMENTATION")
print("="*80)

for K in [3, 4, 5, 6, 7, 8, 9, 10]:
    coeffs, d_K = all_coeffs[K]

    print(f"\n  ! K={K} ({K+1} points)")
    print(f"  real(dp), parameter :: C_K{K}(0:{K}) = [", end="")
    for k in range(K+1):
        if k > 0:
            print(", ", end="")
        print(f"{coeffs[k]:.1f}_dp", end="")
    print("]")
    print(f"  real(dp), parameter :: d_K{K} = {d_K:.1f}_dp")

# Summary table
print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"\n{'K':>3s} {'Points':>7s} {'d_K':>10s} {'c_0':>10s} {'c_K':>10s} {'Match Table I':>15s}")
print("-" * 60)
for K in [3, 4, 5, 6, 7, 8, 9, 10]:
    coeffs, d_K = all_coeffs[K]
    if K in TABLE_I_REFERENCE:
        match = "✓" if np.allclose(coeffs, TABLE_I_REFERENCE[K], atol=0.01) else "✗"
    else:
        match = "N/A"
    print(f"{K:3d} {K+1:7d} {d_K:10.1f} {coeffs[0]:10.1f} {coeffs[K]:10.1f} {match:>15s}")
