#!/usr/bin/env python
"""
Generate Fortran lookup table for K=5 variable timestep coefficients

The lookup table maps a 5-bit integer (representing the timestep history)
to the corresponding coefficients and d_K value.
"""

import numpy as np
import itertools

# Run the derivation
from derive_k5_variable_timesteps import (
    calculate_tau_values, derive_coefficients_for_history,
    c_Kminus1_fixed, c_K_fixed
)

# Generate all 32 possible timestep histories
timestep_options = [0.5, 1.0]
all_histories = list(itertools.product(timestep_options, repeat=5))

# Store results
results = []

for i, history in enumerate(all_histories):
    tau = calculate_tau_values(history)
    coeffs, d_K = derive_coefficients_for_history(tau)

    # Convert history to bit pattern
    # bit_pattern: 1 = full step (1.0), 0 = half step (0.5)
    # bit 0 = most recent step, bit 4 = oldest step
    bit_pattern = 0
    for j, dt in enumerate(history):
        if dt == 1.0:
            bit_pattern |= (1 << j)

    results.append({
        'index': i,
        'bit_pattern': bit_pattern,
        'history': history,
        'tau': tau,
        'coeffs': coeffs,
        'd_K': d_K
    })

# Sort by bit pattern for easier lookup
results.sort(key=lambda x: x['bit_pattern'])

# Generate Fortran code
print("!"*80)
print("! K=5 XLBO Variable Timestep Coefficient Lookup Table")
print("!"*80)
print("!")
print("! Lookup based on 5-bit timestep history:")
print("! - Bit 0 (LSB): most recent timestep (n-1 to n)")
print("! - Bit 4 (MSB): oldest timestep (n-5 to n-4)")
print("! - Bit value: 0 = half step (dt/2), 1 = full step (dt)")
print("!")
print("! Example: history [0.5, 1.0, 1.0, 1.0, 0.5] -> bits 01110 -> index 14")
print("!")

print("\n! Coefficient arrays indexed by bit pattern (0 to 31)")
print("real(dp), parameter :: XLBO_K5_C0(0:31) = [&")
for i, r in enumerate(results):
    comma = "," if i < 31 else " "
    print(f"  {r['coeffs'][0]:12.8f}_dp{comma}  ! {r['bit_pattern']:2d} {str(list(r['history']))}")
print("  ]")

for coeff_idx in range(1, 6):
    print(f"\nreal(dp), parameter :: XLBO_K5_C{coeff_idx}(0:31) = [&")
    for i, r in enumerate(results):
        comma = "," if i < 31 else " "
        print(f"  {r['coeffs'][coeff_idx]:12.8f}_dp{comma}  ! {r['bit_pattern']:2d}")
    print("  ]")

print("\nreal(dp), parameter :: XLBO_K5_dK(0:31) = [&")
for i, r in enumerate(results):
    comma = "," if i < 31 else " "
    print(f"  {r['d_K']:12.8f}_dp{comma}  ! {r['bit_pattern']:2d} {str(list(r['history']))}")
print("  ]")

print("\n!" + "="*78)
print("! Usage example:")
print("!" + "="*78)
print("!")
print("! integer :: history_bits, k")
print("! real(dp) :: C_K5(0:5), d_K")
print("!")
print("! ! Build bit pattern from last 5 timesteps")
print("! history_bits = 0")
print("! do k = 1, 5")
print("!   if (abs(dt_history(k) - 1.0_dp) < 0.1_dp) then  ! Full step")
print("!     history_bits = ibset(history_bits, k-1)")
print("!   endif")
print("! end do")
print("!")
print("! ! Lookup coefficients")
print("! C_K5(0) = XLBO_K5_C0(history_bits)")
print("! C_K5(1) = XLBO_K5_C1(history_bits)")
print("! C_K5(2) = XLBO_K5_C2(history_bits)")
print("! C_K5(3) = XLBO_K5_C3(history_bits)")
print("! C_K5(4) = XLBO_K5_C4(history_bits)")
print("! C_K5(5) = XLBO_K5_C5(history_bits)")
print("! d_K = XLBO_K5_dK(history_bits)")
print("!")

# Summary statistics
print("\n!" + "="*78)
print("! Statistics")
print("!" + "="*78)
d_K_values = [r['d_K'] for r in results]
print(f"! d_K range: [{min(d_K_values):.6f}, {max(d_K_values):.6f}]")
print(f"! d_K mean:  {np.mean(d_K_values):.6f}")
print(f"! d_K std:   {np.std(d_K_values):.6f}")
print("!")
print(f"! Standard K=5 (all full steps): index={results[31]['bit_pattern']}, d_K={results[31]['d_K']:.6f}")
print(f"! All half steps: index={results[0]['bit_pattern']}, d_K={results[0]['d_K']:.6f}")
print("!")
