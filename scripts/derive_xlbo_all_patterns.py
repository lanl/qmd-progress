#!/usr/bin/env python3
"""
Derive XLBO dissipation coefficients for all 32 timestep patterns (K=5).

Based on Niklasson et al., J. Chem. Phys. 130, 214109 (2009), Eq. (18-19).

For K=5, we remove odd-order terms up to δt⁵, giving 5 constraints:
  1. δt⁰: Σc_i = 0
  2. δt¹: Σα_i·c_i = 0
  3. δt²: Σα_i²·c_i = 2d₅
  4. δt³: Σα_i³·c_i = 0
  5. δt⁵: Σα_i⁵·c_i = 0

where α_i are cumulative time distances from P_n.

Pattern encoding: 5-bit binary where bit i = 0 for δt/2, 1 for δt
"""

import numpy as np

def pattern_to_alpha(pattern_bits):
    """
    Convert a 5-bit pattern to cumulative time distances.

    Args:
        pattern_bits: list of 5 bits [b0, b1, b2, b3, b4]
                     where 0 = δt/2, 1 = δt

    Returns:
        alpha: array [α₁, α₂, α₃, α₄, α₅] in units of δt
    """
    alpha = np.zeros(5)
    cumulative = 0.0

    for i in range(5):
        if pattern_bits[i] == 0:
            cumulative += 0.5
        else:
            cumulative += 1.0
        alpha[i] = cumulative

    return alpha


def solve_xlbo_coefficients(alpha):
    """
    Solve for XLBO coefficients given cumulative timesteps.

    Constraints (with c₅ = -1):
      1. c₀ + c₁ + c₂ + c₃ + c₄ = 1
      2. α₁c₁ + α₂c₂ + α₃c₃ + α₄c₄ = α₅
      3. α₁²c₁ + α₂²c₂ + α₃²c₃ + α₄²c₄ = 2d₅ + α₅²
      4. α₁³c₁ + α₂³c₂ + α₃³c₃ + α₄³c₄ = α₅³
      5. α₁⁵c₁ + α₂⁵c₂ + α₃⁵c₃ + α₄⁵c₄ = α₅⁵

    This is 5 equations in 5 unknowns (c₀, c₁, c₂, c₃, c₄).

    Returns:
        c: array [c₀, c₁, c₂, c₃, c₄, c₅]
        d5: the d₅ parameter
    """
    # Build constraint matrix A and right-hand side b
    # We'll solve the system Ax = b where x = [c₀, c₁, c₂, c₃, c₄]

    A = np.array([
        [1,       1,       1,       1,       1      ],  # Constraint 1
        [0,       alpha[0], alpha[1], alpha[2], alpha[3]],  # Constraint 2
        [0,       alpha[0]**2, alpha[1]**2, alpha[2]**2, alpha[3]**2],  # Constraint 3 (temp)
        [0,       alpha[0]**3, alpha[1]**3, alpha[2]**3, alpha[3]**3],  # Constraint 4
        [0,       alpha[0]**5, alpha[1]**5, alpha[2]**5, alpha[3]**5]   # Constraint 5
    ])

    b = np.array([
        1,           # From c₅ = -1
        alpha[4],    # From constraint 2
        0,           # Placeholder for constraint 3 (will solve for d₅)
        alpha[4]**3, # From constraint 4
        alpha[4]**5  # From constraint 5
    ])

    # We have 5 equations but constraint 3 has d₅ as unknown
    # We'll solve constraints 1, 2, 4, 5 first (4 equations, 5 unknowns)
    # Then use constraint 3 to get d₅

    # Actually, we have 5 constraints total. Let's use all of them.
    # We need to set b[2] = 2*d5 + alpha[4]**2, but d5 is unknown.
    #
    # Alternative: We have 5 constraints, 6 unknowns (c₀-c₅, d₅).
    # With c₅ = -1 fixed, we have 5 unknowns (c₀-c₄, d₅).
    # But constraint 3 is the only one with d₅.
    # So we can solve constraints 1,2,4,5 for 4 of the c's in terms of one c,
    # then use constraint 3 to determine both the remaining c and d₅.

    # Better approach: eliminate d₅ from constraint 3 by using it as the definition.
    # Solve constraints 1,2,4,5 as a 4x5 system (underdetermined),
    # pick a normalization, then compute d₅ from constraint 3.

    # Even better: Use constraints 1,2,4,5 (without constraint 3) to get c's,
    # treating this as overdetermined. Then compute d₅ from constraint 3.

    # Actually the cleanest: we have 4 independent constraints (1,2,4,5) and 5 unknowns.
    # We need one more constraint. Let's fix c₄ = 4 (from Table I optimization).

    # Let me reconsider. We actually have 5 constraints and 5 unknowns (c₀-c₄).
    # Constraint 3 tells us what d₅ must be to satisfy the system.
    # So we solve constraints 1,2,4,5 as 4 equations in 5 unknowns (underdetermined),
    # then we have 1 degree of freedom. We pick c₄ based on optimization (minimize |c₀|?),
    # then compute d₅ from constraint 3.

    # Hmm, let me think more carefully. Looking at our earlier derivation:
    # - We have 5 constraints
    # - We have 6 unknowns: c₀, c₁, c₂, c₃, c₄, c₅
    # - We fix c₅ = -1, leaving 5 unknowns
    # - This gives us exactly 5 equations in 5 unknowns - fully determined!
    # - d₅ is NOT an unknown; it's determined by constraint 3 after we solve for c's

    # Wait no. Let's read the paper again. In Eq. 19:
    # δt² P̈_n = (1/d_K) Σc_k P_{n-k} + O(δt⁴) + O(δt^{2K-3})
    #
    # So d_K appears on the LHS. The constraint is that the RHS should equal δt² P̈_n.
    # So constraint 3 is: Σα_i² c_i = d_K × (something from the LHS scaling)

    # Actually, from Eq. 19: δt² P̈_n = (1/d_K) Σc_k P_{n-k}
    # From Eq. 20: F^diss = μ P̈^diss = (μ / d_K δt²) Σc_k P_{n-k}
    #
    # The constraint is that we want: Σc_k P_{n-k} = d_K × δt² × P̈_n × (scaling)
    #
    # From our earlier work, constraint 3 was:
    # (c₁/2 + 2c₂ + 9c₃/2 + 8c₄ + 25c₅/2) / (1/2) = d₅
    # Simplifying: c₁ + 4c₂ + 9c₃ + 16c₄ + 25c₅ = d₅
    # Wait that doesn't match α² either.

    # Let me re-derive from the Taylor expansion. For uniform δt:
    # P_{n-k} = P_n - k δt P'_n + (k² δt²/2) P''_n - ...
    #
    # For non-uniform with α_k:
    # P_{n-k} = P_n - α_k δt P'_n + (α_k² δt²/2) P''_n - ...
    #
    # Summing: Σc_k P_{n-k} = (Σc_k) P_n - (Σα_k c_k) δt P'_n + (Σα_k² c_k)(δt²/2) P''_n - ...
    #
    # We want:
    # - Σc_k = 0 (remove constant term)
    # - Σα_k c_k = 0 (remove δt term)
    # - Σα_k² c_k = 2 d₅ (coefficient of δt² P''_n should be d₅ δt²)

    # So constraint 3 is: Σα_i² c_i = 2 d₅
    # This means d₅ is whatever value makes this constraint work.
    # It's not a free parameter; it's determined by the solution!

    # OK so the approach is:
    # 1. Solve constraints 1, 2, 4, 5 (4 equations, 5 unknowns) → 1 DOF
    # 2. Pick the remaining DOF (e.g., minimize |λ_max| or some other criterion)
    # 3. Compute d₅ from constraint 3

    # For simplicity, let's fix c₄ = 4 as in the uniform case (from optimization in paper)

    # System: 4 equations in 4 unknowns (c₀, c₁, c₂, c₃) with c₄ = 4, c₅ = -1
    A_reduced = np.array([
        [1,       1,       1,       1      ],  # Constraint 1
        [0,       alpha[0], alpha[1], alpha[2]],  # Constraint 2
        [0,       alpha[0]**3, alpha[1]**3, alpha[2]**3],  # Constraint 4
        [0,       alpha[0]**5, alpha[1]**5, alpha[2]**5]   # Constraint 5
    ])

    b_reduced = np.array([
        1 - 4,                            # c₀ + c₁ + c₂ + c₃ = 1 - c₄ = -3
        alpha[4] - 4*alpha[3],            # Σα c = α₅ - c₄α₄
        alpha[4]**3 - 4*alpha[3]**3,      # Σα³c = α₅³ - c₄α₄³
        alpha[4]**5 - 4*alpha[3]**5       # Σα⁵c = α₅⁵ - c₄α₄⁵
    ])

    # Solve the system
    try:
        c_partial = np.linalg.solve(A_reduced, b_reduced)
    except np.linalg.LinAlgError:
        print(f"Singular matrix for alpha = {alpha}")
        return None, None

    c = np.append(c_partial, [4.0, -1.0])

    # Compute d₅ from constraint 3
    sum_alpha2_c = sum(alpha[i]**2 * c[i+1] for i in range(5))
    d5 = sum_alpha2_c / 2.0

    return c, d5


def main():
    print("=" * 80)
    print("XLBO Dissipation Coefficients for All 32 Timestep Patterns (K=5)")
    print("=" * 80)
    print()

    # Store results
    results = []

    for pattern_idx in range(32):
        # Decode pattern index to bits
        pattern_bits = [(pattern_idx >> i) & 1 for i in range(5)]

        # Convert to cumulative alpha values
        alpha = pattern_to_alpha(pattern_bits)

        # Solve for coefficients
        c, d5 = solve_xlbo_coefficients(alpha)

        if c is None:
            print(f"Pattern {pattern_idx:2d}: FAILED (singular matrix)")
            continue

        # Build pattern string
        pattern_str = ','.join(['1' if b == 1 else '½' for b in pattern_bits])

        results.append({
            'idx': pattern_idx,
            'bits': pattern_bits,
            'pattern': pattern_str,
            'alpha': alpha,
            'c': c,
            'd5': d5
        })

        print(f"Pattern {pattern_idx:2d}: [{pattern_str}]")
        print(f"  α = [{', '.join(f'{a:.1f}' for a in alpha)}]")
        print(f"  c = [{', '.join(f'{ci:7.2f}' for ci in c)}]")
        print(f"  d₅ = {d5:7.3f}")
        print()

    print("=" * 80)
    print("Fortran array initialization:")
    print("=" * 80)
    print()

    for r in results:
        comment = f"! Pattern {r['idx']:2d}: [{r['pattern']}]"
        coeff_str = f"xlbo_coeff_set({r['c'][0]:.1f}_dp, {r['c'][1]:.1f}_dp, {r['c'][2]:.1f}_dp, {r['c'][3]:.1f}_dp, {r['c'][4]:.1f}_dp, {r['c'][5]:.1f}_dp, {r['d5']:.3f}_dp)"

        if r['idx'] < 31:
            print(f"    {coeff_str}, & {comment}")
        else:
            print(f"    {coeff_str} &  {comment}")

    print()
    print("=" * 80)
    print("Verification: Check pattern 30 and 31 match our manual derivation")
    print("=" * 80)

    # Pattern 30: [½,1,1,1,1]
    r30 = results[30]
    print(f"Pattern 30 (½,1,1,1,1):")
    print(f"  Expected: c = [-30, 49, -21, -1, 4, -1], d₅ = -6.25")
    print(f"  Computed: c = {r30['c']}, d₅ = {r30['d5']:.3f}")

    # Pattern 31: [1,1,1,1,1]
    r31 = results[31]
    print(f"Pattern 31 (1,1,1,1,1):")
    print(f"  Expected: c = [-6, 14, -8, -3, 4, -1], d₅ = -3.0")
    print(f"  Computed: c = {r31['c']}, d₅ = {r31['d5']:.3f}")


if __name__ == '__main__':
    main()
