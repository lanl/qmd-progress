# K=5 XLBO Coefficients for Variable Timesteps

## Overview

Generalization of K=5 XLBO dissipation to handle variable timesteps of dt/2 or dt.

## Mathematical Foundation

### Modified Equations

**Original (uniform Δt):**
```
P_{n-k} = Σ(m=0 to M) [(-k·Δt)^m / m!] · P_n^(m) + O(Δt^(M+1))
```

**Generalized (variable Δt):**
```
P_{n-k} = Σ(m=0 to M) [(-τ_k·Δt)^m / m!] · P_n^(m) + O(Δt^(M+1))
```

where **τ_k** is the cumulative time elapsed from step n to step n-k, in units of Δt.

### Example

For timestep history `[0.5, 1.0, 1.0, 1.0, 0.5]` (most recent first):
- τ_0 = 0.0
- τ_1 = 0.5 (half step back)
- τ_2 = 0.5 + 1.0 = 1.5
- τ_3 = 1.5 + 1.0 = 2.5
- τ_4 = 2.5 + 1.0 = 3.5
- τ_5 = 3.5 + 0.5 = 4.0

### Normalization

Same as uniform K=5:
- c_4 = 4
- c_5 = -1

### Moment Constraints

With τ_k replacing k:
- Moment 0: Σ c_k = 0
- Moment 1: Σ τ_k·c_k = 0
- Moment 3: Σ τ_k³·c_k = 0
- Moment 5: Σ τ_k⁵·c_k = 0

From these 4 equations, we solve for c_0, c_1, c_2, c_3.

The dissipation parameter d_K is calculated from moment 2:
- d_K = |Σ τ_k²·c_k / 2|

## Results

### All 32 Possible Histories

With 5 historical timesteps, each either dt/2 or dt, there are 2^5 = 32 possibilities.

**Key Examples:**

| History | τ_5 | c_0 | c_1 | c_2 | c_3 | c_4 | c_5 | d_K |
|---------|-----|-----|-----|-----|-----|-----|-----|-----|
| [1, 1, 1, 1, 1] | 5.0 | -6.0 | 14.0 | -8.0 | -3.0 | 4.0 | -1.0 | 3.00 |
| [0.5, 0.5, 0.5, 0.5, 0.5] | 2.5 | -6.0 | 14.0 | -8.0 | -3.0 | 4.0 | -1.0 | 0.75 |
| [0.5, 1, 1, 1, 0.5] | 4.0 | 28.4 | -50.6 | 32.8 | -13.6 | 4.0 | -1.0 | 4.70 |

**Observations:**
1. All full steps and all half steps use the **same coefficients** (uniform case)
2. But d_K differs by factor of 4 (scales with τ²)
3. Mixed histories have **different coefficients** specific to their pattern
4. d_K ranges from 0.016 to 11.75

### Statistics

- **d_K range**: [0.016, 11.75]
- **d_K mean**: 3.54
- **d_K std**: 3.00

## Implementation

### Lookup Table Approach

The coefficients are stored in Fortran arrays indexed by a 5-bit integer representing the timestep history:

```fortran
! Bit pattern encoding:
! Bit 0 (LSB): most recent timestep (n-1 to n)
! Bit 4 (MSB): oldest timestep (n-5 to n-4)
! Bit value: 0 = half step (dt/2), 1 = full step (dt)

! Example: [0.5, 1.0, 1.0, 1.0, 0.5] -> bits 01110 -> index 14

real(dp), parameter :: XLBO_K5_C0(0:31) = [...]
real(dp), parameter :: XLBO_K5_C1(0:31) = [...]
real(dp), parameter :: XLBO_K5_C2(0:31) = [...]
real(dp), parameter :: XLBO_K5_C3(0:31) = [...]
real(dp), parameter :: XLBO_K5_C4(0:31) = [...]
real(dp), parameter :: XLBO_K5_C5(0:31) = [...]
real(dp), parameter :: XLBO_K5_dK(0:31) = [...]
```

### Usage

```fortran
integer :: history_bits, k
real(dp) :: C_K5(0:5), d_K

! Build bit pattern from last 5 timesteps
history_bits = 0
do k = 1, 5
  if (abs(dt_history(k) - 1.0_dp) < 0.1_dp) then  ! Full step
    history_bits = ibset(history_bits, k-1)
  endif
end do

! Lookup coefficients
C_K5(0) = XLBO_K5_C0(history_bits)
C_K5(1) = XLBO_K5_C1(history_bits)
C_K5(2) = XLBO_K5_C2(history_bits)
C_K5(3) = XLBO_K5_C3(history_bits)
C_K5(4) = XLBO_K5_C4(history_bits)
C_K5(5) = XLBO_K5_C5(history_bits)
d_K = XLBO_K5_dK(history_bits)
```

## Files

- `../derive_k5_variable_timesteps.py` - Derives all 32 coefficient sets
- `../generate_k5_fortran_lookup.py` - Generates Fortran lookup tables
- `../k5_variable_timestep_coefficients.txt` - Detailed results for all histories
- `../k5_variable_timestep_fortran.f90` - Fortran parameter declarations

## Validation

All coefficient sets satisfy the moment constraints to machine precision:
- |m0| < 1e-10
- |m1| < 1e-10
- |m3| < 1e-6
- |m5| < 1e-6

## References

Based on:
> Niklasson, A. M. N., et al. (2009). Extended Lagrangian Born–Oppenheimer molecular dynamics with dissipation. *The Journal of Chemical Physics*, 130(21), 214109.

Equations (18) and (19), generalized to variable τ_k.
