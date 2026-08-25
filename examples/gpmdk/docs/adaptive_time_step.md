# Adaptive Time Step (Split-Step XLBO)

This document describes the adaptive time-step ("split-step") feature of the gpmdk QMD
driver: what it does, how to enable it, the extended-Lagrangian (XLBO) dissipation rule that
makes it work across a mixed sequence of full and half steps, its measured energy-conservation
behavior on small and large systems, and guidance on when it is appropriate (NVT vs NVE).

## What the adaptive time step does

The MD loop (`examples/gpmdk/src/gpmdcov_mdloop.F90`) runs the user time step `dt` normally,
but **splits a full step into two half-steps** (`dt/2`) whenever the projected maximum atomic
displacement over the next step would exceed a threshold:

```fortran
this_maxdisp = maxval(user_timestep * sy%velocity)      ! projected max displacement
if (gpmdt%adaptive_timestep .and. &
    (first_substep_taken .or. (this_maxdisp > maxdist)) .and. &
    mdstep > gpmdt%minimization_steps) then
  ! take two dt/2 half-steps instead of one dt step
```

with `maxdist = 0.02` (Angstrom) as the split threshold. The purpose is to **bound the
per-step atomic displacement**: if an atom would move too far in a single step, the step is
halved so the fast degrees of freedom stay well integrated. For large *reactive* systems this
also prevents atoms from jumping far enough in one step to trigger unphysical reactions.

During the initial XLBO history build the first few user steps are always split (the K=5
integrator needs at least 6 half-steps of history), after which normal adaptive splitting
applies. Splitting is suppressed during the minimization/annealing phase
(`mdstep <= minimization_steps`).

Each time a new split begins, the driver prints:

```
Rank <r>  Splitting print_mdstep <n>
```

Trajectory output and `Mdstep, Energy, ...` lines remain on the user-step schedule (keyed to
`print_mdstep`), independent of split decisions.

## Enabling it

Set the input option in the `GPMD{ ... }` block:

```
GPMD{
  AdaptiveTimeStep= T
  ...
}
```

Parsed in `gpmdcov_parser.F90` into `gpmdt%adaptive_timestep`. Default is off (fixed `dt`).

## How the electronic integrator handles the mixed step sequence

The XLBO integrator (`src/prg_xlbo_mod.F90`) is a K=5 six-point multistep predictor whose
coefficients were derived assuming **uniform** step spacing. A split sequence places the
history points `n_0..n_5` at non-uniform intervals, so the code selects the coefficients from
**per-pattern lookup tables** (`XLBO_K5_C0..C5(0:31)`, `XLBO_K5_dK(0:31)`) indexed by a 5-bit
integer encoding the recent full/half pattern (`get_K5_history_index`, most-recent-first). All
32 patterns satisfy the moment constraints (m0 = m1 = m3 = m5 = 0), so the *conservative*
backbone stays correct under splitting.

Two structural facts make this work and are worth stating explicitly, because both were
non-obvious and both were arrived at only after rejected alternatives (see "Rejected
approaches" below):

- **kappa is constant (`kappa = 1.82`)**, not pattern-dependent. The physical `omega²` is a
  constant; the `dt²` dependence already lives in the Verlet scaling factor
  `ka_scale = 0.5 * dt_n * (dt_n + dt_prev)`. That factor is **symmetric** under
  `dt_n <-> dt_prev`, so each split interface returns the identical coefficient when the
  trajectory is reversed — the Verlet backbone remains time-reversible to machine precision
  (~2e-16) even under variable stepping. Making kappa pattern-dependent *destroys* this
  (reversibility error ~8e-4); do not do it.
- **Within the electronic integrator, the split-related contribution is in the dissipation
  term (alpha), not the conservative backbone.** The backbone is reversible for any fixed
  schedule (previous point) and the coefficient tables preserve the moment constraints, so
  the only XLBO-internal split asymmetry was the *dissipation rate* (addressed by the alpha
  rate-match below). This scopes what kappa/coefficients can and cannot affect — it does
  **not** mean the alpha rule eliminates the observed drift. The dominant residual drift
  comes from a level *above* the integrator: the state-dependent choice of *when* to split
  (see "Fixed-schedule control experiment" below), which no XLBO coefficient controls.

## The general (history-independent) alpha rule

The dissipation coefficient `alpha` must match the full-step dissipation *rate per unit
physical time*. The subtlety: two `dt/2` half-steps span the same physical interval as one
`dt` full step, so a naive constant `alpha·d_K` per *step* makes half-ending patterns dissipate
at ~2× the reference rate per femtosecond — a split/unsplit asymmetry that drives secular
energy drift. The rule below removes that asymmetry by scaling alpha with `dt/d_K`:

```
alpha_idx = min(0.054 * dt_current / d_K(idx), 0.15)
```

where `0.054 = 0.018 * 3.0` is the full-step reference dissipation rate per unit physical
time. This is a single **history-independent** rule applied to all 32 patterns, so it supports
**any** split history rather than a table tuned to one run's split statistics. Patterns
`idx 1, 3, 5, 6, 7` hit the 0.15 stability ceiling; the rest are rate-matched.

The ceiling (0.15) was chosen empirically, not from periodic-orbit eigenvalues (which are
pathological here — several "1 split every N steps forever" orbits show parametric resonance
even at alpha = 0, yet run stably in practice). Validation propagates the exact recurrence over
*quasi-random* histories (`scripts/Niklasson_JCP_2009_table_I/empirical_stability.py`). For a
converged SCF the operating point is `gamma = 1 - |q-n|/|n| >= ~0.999`, where the scheme is
stable (peak amplitude ~1.3 over 20000 steps); 0.15 keeps a safety margin as gamma degrades,
while 0.20 was rejected as more fragile once gamma slips below 0.999.

The honest finding: **alpha alone cannot fully close the friction deficit.** Over a
split-and-recovery transient the aggregate dissipation falls below the pure-full-step baseline,
and the deficit is owned almost entirely by `idx 7` (the post-split tail — three full steps
after a split), which would need alpha ≈ 0.69 to rate-match but is stability-locked near 0.15.
Raising the ceiling barely changes aggregate friction. The general rule is therefore a
robustness/code-quality improvement — one clean rule, no pathological per-pattern caps — and,
as shown below, carries no drift penalty relative to the earlier tuned table. The real levers
for the residual deficit are **structural** (gentler splitting, or a longer-history K=10
kernel), not alpha or kappa.

## Energy-conservation methodology (equilibrate-then-restart)

To give every test a common initial condition:

1. **Equilibrate:** `ts=0.25`, 2000 steps, `GPMD{ DumpEach=100, RestartFromDump=F }` →
   equilibrated to ~330 K. Writes `restart.dmp`.
2. **Pin it:** `cp restart.dmp restart_equil.dmp; chmod 444`. Before each test:
   `cp -f restart_equil.dmp restart.dmp; chmod 644`, run with `RestartFromDump=T`.
3. **Caveat:** the restart is not numerically bit-consistent — the XLBO history reload does
   not truly continue the integrator (mdstep resets to 0), so each run re-settles over a short
   transient. **Discard the first ~50 MD records** before fitting.
4. **Drift metric:** OLS slope of `Energy Total [eV]` vs time, ×1000 → eV/ps (a fitted trend,
   not an endpoint difference).
5. **Convergence test:** split the post-transient window in half and fit each half. If the
   halves disagree in sign or by >3×, the full-window slope is oscillation-phase noise, not
   drift.

Analysis scripts: `examples/gpmdk/run/water/drift_analysis.py`, `toggle_count.py`.

## Results — 300-atom water

System: `coords_300.dat`, `Method=DiagEfFull`, `kBT=0.025`.

Short runs (1000 steps from pinned equil, skip 50):

| run                  | t_end fs | ptp eV | rms eV | drift eV/ps | 1st / 2nd half | converged |
|----------------------|---------:|-------:|-------:|------------:|:--------------:|:---------:|
| 0.5 FIXED (no adapt)  | 499.5   | 0.035  | 0.0074 | +0.0005     | +.003 / −.005  | YES       |
| 0.45 adaptive         | 449.8   | 0.065  | 0.014  | +0.079      | +.038 / +.191  | NO        |
| 0.5  adaptive         | 499.8   | 0.070  | 0.018  | −0.092      | −.274 / +.046  | NO        |
| 0.6  adaptive         | 599.7   | 0.013  | 0.0025 | −0.001      | +.003 / −.004  | YES       |

At 1000 steps the intermittent-split cases (0.45, 0.5) fail the half-window test — the "drift"
is the same order as the energy oscillation it rides on.

Long run (10000 steps, ts=0.5 adaptive, ~5000 fs) — resolves the convergence question:

| window   | n     | span (fs)     | drift eV/ps | ptp eV | rms eV |
|----------|------:|---------------|------------:|-------:|-------:|
| **FULL** | 19669 | 12.5 – 4999.8 | **+0.006**  | 0.100  | 0.021  |
| 1st half |  9834 | 12.5 – 2501   | +0.026      | 0.100  | 0.025  |
| 2nd half |  9835 | 2501 – 5000   | +0.015      | 0.052  | 0.014  |

**Converged: YES** (halves agree in sign, within 3×; stable at skip = 50/200/500). Total energy
wander over the 5 ps run is ~0.03 eV — buried well inside the 0.1 eV physical oscillation, and
two orders of magnitude below the short-run headline estimate. No secular leak.

## Result — large system (Mac1)

System: `Mac1_gmx_kernel` (total energy ≈ −148 422 eV), `dt = 0.4 fs`, adaptive split active,
graph partitioning on. Log format: `Mdstep, Energy, Egap, Resnorm, Temp` (time = step × dt).

| window   | n     | steps         | drift eV/ps | ptp eV | rms eV |
|----------|------:|---------------|------------:|-------:|-------:|
| **FULL** | 28130 | 51 – 28180    | **+0.596**  | 6.83   | 1.94   |
| 1st half | 14065 | 51 – 14115    | +0.576      | —      | —      |
| 2nd half | 14065 | 14115 – 28180 | +0.523      | —      | —      |

**Converged: YES.** Split fraction 21.5% (6047/28180 post-warmup steps), duration ≈ 11.3 ps.

The absolute drift (+0.60 eV/ps) is ~100× the water value, but the drift is **extensive**
(scales with system size), not a new large-system pathology — same mechanism, more atoms. For
this ~30k-atom system that is ≈ **2×10⁻⁵ eV/atom/ps**, small and converged.

Unlike the water case, this drift is **not** buried inside the physical energy oscillation:
over the 11.3 ps run it accumulates ~6.7 eV, comparable to the ptp (6.83 eV) itself. It is a
genuine secular trend — small per atom, but real and cumulative at scale, so it matters for
long NVE runs (see NVE guidance below).

## Fixed-schedule control experiment (drift mechanism)

To separate "how *often* the step size toggles" from "what *fraction* of steps are split," the
state-dependent trigger was replaced with a deterministic cadence hitting the same 30.6% split
fraction, on the 300-atom water restart:

| scheme (300-atom water) | split frac | F↔S toggles | drift eV/ps |
|-------------------------|-----------:|------------:|------------:|
| velocity/displacement trigger | ~30.6% | baseline | ~ +0.006 (10k, converged) |
| **fixed cadence, same fraction** | 30.6% | ~8× more | **−0.018** (quietest) |

**The drift is not caused by toggle frequency or split fraction.** The fixed-cadence run has
~8× *more* full↔half switches at the *same* split fraction, yet is the quietest with no secular
drift. What distinguishes it is that the split decision is **state-independent** — it never
looks at the trajectory, so it is trivially time-reversible (forward and backward passes make
identical choices). The displacement-triggered scheme decides from the *entering* endpoint
only, so at threshold crossings the forward and backward passes disagree — a phase-locked,
one-signed energy injection correlated with the fast mode. That phase correlation, not the
switching rate, is the drift source. This is Hairer's "trap of variable step size":
**state-dependent step-size selection breaks time-reversibility.**

## When to use it — NVT vs NVE

- **NVT:** appropriate. The thermostat absorbs the small phase-locked injection; the
  displacement bound (its actual purpose) is delivered; there is no reversibility requirement
  to violate.
- **NVE, moderate duration:** usable with care. The drift is a small, converged trend that
  does not blow up (water: +0.006 eV/ps over 5 ps — buried inside the 0.1 eV oscillation;
  large system: +0.60 eV/ps ≈ 2×10⁻⁵ eV/atom/ps over 11 ps). Note the large-system drift is
  *not* hidden inside the oscillation — it accumulates to ~6.7 eV, comparable to the ptp — so
  it is a genuine secular trend. Judge acceptability against the per-atom rate and the
  intended run length, not against the oscillation amplitude.
- **NVE, long duration / strict energy conservation:** not recommended with the current
  state-dependent trigger. Either use a fixed-cadence schedule (drift-free, but cannot
  *guarantee* the displacement bound) or await the reversible backtracking scheme below.

## Planned improvement — reversible backtracking

Exact reversibility with a displacement-bounding (state-dependent) trigger requires making the
split decision a **symmetric function of both interval endpoints**, realized by backtracking:

```
fast(x) ≡ [ dt*max|v(x)| > maxdist ]
split iff fast(x_n) OR fast(x_{n+1})
  - fast(x_n) true  -> split immediately (reversible fast path; v_n already in hand)
  - fast(x_n) false -> take a tentative full step, test fast(x_{n+1}^full);
                       if true, discard and redo as two half-steps
```

Extra cost is confined to threshold-straddling steps (one tentative full-step force
evaluation), so it does not broadly add split steps — preserving the point of the feature
(longer average steps → longer reachable durations). Keeping the displacement (`dt*|v|`)
criterion is preferred over a force/acceleration criterion: the residual O(h³)
"probe-endpoint ≠ kept-endpoint" mismatch is identical for either, but displacement gets the
reversible fast path for free (`v_n` is available at loop top; force is not) and it bounds
displacement directly rather than as a curvature proxy.

## Rejected approaches (do not revisit)

The path to the current design passed through several approaches that were tried and found not
to work. They are recorded here so they are not re-attempted:

- **Charge interpolation (Lagrange or cubic spline) onto a uniform dt/2 grid.** Early
  split-step branches interpolated the history `n_0..n_5` to a uniform grid and applied the
  fixed uniform coefficients. It is mathematically equivalent to adaptive coefficients and gave
  0.0 eV agreement on short tests, but the energy-drift problem persists with it, and it was
  removed in favor of the per-pattern coefficient tables. Do not reintroduce interpolation.
- **Per-pattern (dt-pattern-dependent) kappa.** Softening kappa for split/half-ending patterns
  to suppress a formal stability mode breaks time-reversibility of the Verlet backbone
  (reversibility error ~8e-4 vs ~2e-16 for constant kappa). Keep kappa constant.
- **Removing kappa scaling while extending history to 11 points but keeping only C0–C5.** This
  was ~19× *worse* on drift. An 11-point (dt/2) history needs the extended C0–C10 (K=10)
  coefficient set to preserve the 5·dt_base temporal window; using C0–C5 with 11 points
  underutilizes the history.
- **Periodic-orbit eigenvalue caps for the alpha ceiling.** Ill-posed here (the neutral
  charge-conservation mode sits at exactly 1.0; several split cadences show parametric
  resonance at unphysical gamma even at alpha = 0). Use empirical quasi-random-history
  propagation instead.

A K=10 (11-point history) kernel is a plausible *future* structural improvement for matching
arbitrary split histories — its coefficients (`C0_K10..C10_K10`, `kappa_K10 = 1.88`) are
present in `prg_xlbo_mod.F90` but the K=10 predictor is not the active adaptive scheme, and its
derivation for the mixed-history case was left unfinished.

## Reproduction

- 300-atom water: pinned `restart_equil.dmp`, `input_restart_ts05_10k.in`,
  `drift_analysis.py`, `toggle_count.py` under `examples/gpmdk/run/water/`.
- Large system: `examples/gpmdk/run/Mac1_gmx_kernel/energy_split.txt`, analyzed with
  `drift_analysis.py` adapted for `Mdstep, Energy` records (time = step × 0.4 fs).
