#!/usr/bin/env python
"""
Derive K=5 XLBO dissipation coefficients for variable timesteps

For K=5, we have 6 points (including current): n, n-1, n-2, n-3, n-4, n-5
The last 5 timesteps can each be either dt (full step = 1.0) or dt/2 (half step = 0.5)
This gives 2^5 = 32 possible histories.

For each history, we calculate τ_k (time elapsed from step n to step n-k in units of Δt)
and derive coefficients c_0 through c_5 using the same normalization as uniform K=5:
  c_4 = 4
  c_5 = -1

Modified equations:
  P_{n-k} = Σ(m=0 to M) [(-τ_k·Δt)^m / m!] · P_n^(m) + O(Δt^(M+1))

  Δt² · P̈_n = (1/d_K) · Σ(k=0 to K) c_k · P_{n-k}

Moment constraints (using τ_k instead of k):
  - Moment 0: Σ c_k = 0
  - Moment 1: Σ τ_k·c_k = 0
  - Moment 3: Σ τ_k³·c_k = 0
  - Moment 5: Σ τ_k⁵·c_k = 0
"""

import numpy as np
import itertools

# K=5 normalization (same as uniform case)
K = 5
c_Kminus1_fixed = 4.0   # c_4 = 2*(K-3)*(-1)^(K-1) = 4
c_K_fixed = -1.0         # c_5 = (-1)^K = -1

def calculate_tau_values(timestep_history):
    """
    Calculate τ_k values from timestep history.

    timestep_history: list of 5 timesteps (most recent first)
                     Each element is either 0.5 (half step) or 1.0 (full step)

    Returns: array of 6 τ values [τ_0, τ_1, τ_2, τ_3, τ_4, τ_5]
             where τ_k is cumulative time back from step n
    """
    tau = np.zeros(6)
    tau[0] = 0.0  # Current step

    for k in range(1, 6):
        tau[k] = tau[k-1] + timestep_history[k-1]

    return tau

def derive_coefficients_for_history(tau_values):
    """
    Derive c_0 through c_3 using moment constraints with given τ_k values.
    c_4 and c_5 are fixed by normalization.

    Returns: (coeffs, d_K)
    """
    # Number of unknowns: c_0 through c_3
    n_unknowns = 4

    # Build constraint equations using τ_k values
    A = []
    b = []

    # Moment 0: Σc_k = 0
    row = np.ones(n_unknowns)
    A.append(row)
    b.append(-(c_Kminus1_fixed + c_K_fixed))

    # Moment 1: Σ τ_k·c_k = 0
    row = np.array([tau_values[k] for k in range(n_unknowns)])
    A.append(row)
    b.append(-(tau_values[4] * c_Kminus1_fixed + tau_values[5] * c_K_fixed))

    # Moment 3: Σ τ_k³·c_k = 0
    row = np.array([tau_values[k]**3 for k in range(n_unknowns)])
    A.append(row)
    b.append(-(tau_values[4]**3 * c_Kminus1_fixed + tau_values[5]**3 * c_K_fixed))

    # Moment 5: Σ τ_k⁵·c_k = 0
    row = np.array([tau_values[k]**5 for k in range(n_unknowns)])
    A.append(row)
    b.append(-(tau_values[4]**5 * c_Kminus1_fixed + tau_values[5]**5 * c_K_fixed))

    A = np.array(A)
    b = np.array(b)

    # Solve for c_0, c_1, c_2, c_3
    try:
        solution = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        print(f"Warning: Singular matrix for tau = {tau_values}")
        return None, None

    # Assemble full coefficients
    coeffs = np.zeros(6)
    coeffs[:4] = solution
    coeffs[4] = c_Kminus1_fixed
    coeffs[5] = c_K_fixed

    # Calculate d_K from moment 2
    moment2 = sum(tau_values[k]**2 * coeffs[k] for k in range(6))
    d_K = abs(moment2 / 2.0)  # Use absolute value as in the uniform K=5 case

    return coeffs, d_K

def verify_moments(tau_values, coeffs, d_K):
    """Verify that coefficients satisfy moment constraints"""
    results = {}

    # Moment 0
    m0 = sum(coeffs)
    results['m0'] = m0

    # Moment 1
    m1 = sum(tau_values[k] * coeffs[k] for k in range(6))
    results['m1'] = m1

    # Moment 2
    m2 = sum(tau_values[k]**2 * coeffs[k] for k in range(6))
    results['m2'] = m2
    results['d_K'] = d_K

    # Moment 3
    m3 = sum(tau_values[k]**3 * coeffs[k] for k in range(6))
    results['m3'] = m3

    # Moment 5
    m5 = sum(tau_values[k]**5 * coeffs[k] for k in range(6))
    results['m5'] = m5

    return results

# Generate all 32 possible timestep histories
# Each of the 5 most recent timesteps can be 0.5 or 1.0
timestep_options = [0.5, 1.0]
all_histories = list(itertools.product(timestep_options, repeat=5))

print("="*80)
print("K=5 XLBO Coefficients for Variable Timesteps")
print("="*80)
print(f"\nNormalization: c_4 = {c_Kminus1_fixed}, c_5 = {c_K_fixed}")
print(f"Total histories: {len(all_histories)}")
print()

# Store results
results_table = []

for i, history in enumerate(all_histories):
    tau = calculate_tau_values(history)
    coeffs, d_K = derive_coefficients_for_history(tau)

    if coeffs is None:
        continue

    # Verify moments
    verification = verify_moments(tau, coeffs, d_K)

    # Store for summary
    results_table.append({
        'index': i,
        'history': history,
        'tau': tau,
        'coeffs': coeffs,
        'd_K': d_K,
        'verification': verification
    })

# Print detailed results for a few examples
print("\nExample 1: All full steps [1.0, 1.0, 1.0, 1.0, 1.0]")
print("-" * 80)
ex1 = [r for r in results_table if r['history'] == (1.0, 1.0, 1.0, 1.0, 1.0)][0]
print(f"τ values: {ex1['tau']}")
print(f"Coefficients: {ex1['coeffs']}")
print(f"d_K = {ex1['d_K']:.6f}")
print(f"Moment verification: m0={ex1['verification']['m0']:.2e}, m1={ex1['verification']['m1']:.2e}, m3={ex1['verification']['m3']:.2e}")

print("\nExample 2: All half steps [0.5, 0.5, 0.5, 0.5, 0.5]")
print("-" * 80)
ex2 = [r for r in results_table if r['history'] == (0.5, 0.5, 0.5, 0.5, 0.5)][0]
print(f"τ values: {ex2['tau']}")
print(f"Coefficients: {ex2['coeffs']}")
print(f"d_K = {ex2['d_K']:.6f}")
print(f"Moment verification: m0={ex2['verification']['m0']:.2e}, m1={ex2['verification']['m1']:.2e}, m3={ex2['verification']['m3']:.2e}")

print("\nExample 3: Mixed [0.5, 1.0, 1.0, 1.0, 0.5]")
print("-" * 80)
ex3 = [r for r in results_table if r['history'] == (0.5, 1.0, 1.0, 1.0, 0.5)][0]
print(f"τ values: {ex3['tau']}")
print(f"Coefficients: {ex3['coeffs']}")
print(f"d_K = {ex3['d_K']:.6f}")
print(f"Moment verification: m0={ex3['verification']['m0']:.2e}, m1={ex3['verification']['m1']:.2e}, m3={ex3['verification']['m3']:.2e}")

# Summary table
print("\n" + "="*80)
print("SUMMARY: All 32 Histories")
print("="*80)
print(f"{'ID':>3s} {'History':>30s} {'τ_5':>6s} {'c_0':>10s} {'c_1':>10s} {'c_2':>10s} {'c_3':>10s} {'d_K':>10s}")
print("-" * 80)

for r in results_table:
    history_str = str([h for h in r['history']])
    print(f"{r['index']:3d} {history_str:>30s} {r['tau'][5]:6.1f} {r['coeffs'][0]:10.4f} {r['coeffs'][1]:10.4f} {r['coeffs'][2]:10.4f} {r['coeffs'][3]:10.4f} {r['d_K']:10.6f}")

# Check for any patterns in d_K values
d_K_values = [r['d_K'] for r in results_table]
print(f"\nd_K statistics:")
print(f"  Min: {min(d_K_values):.6f}")
print(f"  Max: {max(d_K_values):.6f}")
print(f"  Mean: {np.mean(d_K_values):.6f}")
print(f"  Std: {np.std(d_K_values):.6f}")

# Save results to file
output_file = "k5_variable_timestep_coefficients.txt"
with open(output_file, 'w') as f:
    f.write("K=5 XLBO Coefficients for Variable Timesteps\n")
    f.write("=" * 80 + "\n")
    f.write(f"Normalization: c_4 = {c_Kminus1_fixed}, c_5 = {c_K_fixed}\n\n")

    for r in results_table:
        f.write(f"History {r['index']}: {[h for h in r['history']]}\n")
        f.write(f"  τ: {r['tau']}\n")
        f.write(f"  c: {r['coeffs']}\n")
        f.write(f"  d_K: {r['d_K']:.6f}\n\n")

print(f"\nResults saved to {output_file}")
