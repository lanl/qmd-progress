# gpmdk optimizations since the Venado hackathon (manuscript notes)

Anchor: commit `e89b45c` "Add gpmdk for the Venado hackathon" (2024-07-02), extended back to
the **early-2024 vectorization campaign** (from `fee2dc5`, 2024-04-03), which is the natural
"before" state for the manuscript. GPU *enablement* groundwork is older still (crusher offload
build `c0b9a8b`, 2022; GPU-run edits `09d11dc`, 2023-07) but is infrastructure, not algorithmic
optimization.
Span: 2024-04 → 2026-08-25, ~330 commits, ~140 performance-related.
This is a reconstruction from commit history (subjects, the 51 commits with bodies, and diffs
of the major commits), organized by theme. Commit hashes are given so each claim is traceable;
quantitative figures are cited only where a commit body recorded them.

> Caveat for the manuscript: most HPC-optimization commits recorded *what* changed but not
> a measured speedup. The numbers below that carry a `%`/`×` come from commit bodies; the
> rest are qualitative and should be backed by fresh benchmarks before publication.

---

## Cross-cutting optimization principles

Three recurring principles organize most of the HPC work below. Each is a design rule that
shows up repeatedly across otherwise-unrelated subsystems (Ewald, neighbor list, kernel, DM
build), so they are worth stating once as themes and then tracing through the per-subsystem
catalog (§0–§6) where the concrete commits live.

### P1. Minimize allocation cost — reuse and resize existing buffers instead of re-allocating

The hot path runs every MD step, so per-step `allocate`/`deallocate` churn (especially of BML
matrices and neighbor-list arrays) is pure overhead. The recurring fix is to allocate once and
then *reuse* or *resize in place*:

- **Allocate-once / guarded reuse of BML work matrices** (`9033070` "Reuse some matrices",
  `8e9d037` "protections against double alloc"): the kernel's `ptham/ptrho/zq/zqt/ptaux/…`
  matrices are created behind an `if(.not.bml_allocated(...))` guard and kept across steps
  rather than rebuilt each call.
- **Resize on the existing buffer, avoid a transpose allocation** (`a20cb00`): replaced
  `bml_transpose(zq_bml,zqt_bml)` (allocates a result) with `bml_copy` + the new
  `bml_transpose_inplace(zqt_bml)`; `zqt` is resized in place rather than reallocated.
- **`move_alloc` to hand off a buffer with zero copy** (`7e8e21b` "Use move_alloc to eliminate
  CPU copy"): neighbor-list arrays (`nnType/nnStruct/nrnnStruct/nrnnlist`) are built in local
  temporaries and *moved* into the `nll%` structure instead of allocated-then-copied.
- **Smaller / right-sized work arrays** (`599e723` "Allocate smaller array for work in kernel",
  `b3d5c84` "Reduce allocations in hcsf method", `9a7f5bc`/`3fdc774` resize-and-init).
- **Persistent pairwise/Ewald parameters** (`0b654b6`→`285a541` "Working persistent ewald
  params"): the per-pair Slater-Koster/Ewald parameter arrays are computed once and kept
  resident rather than rebuilt each step (the CPU-side analogue of P2's persistent GPU data).

### P2. Prioritize GPU offload of whole kernels — keep data resident, minimize host↔device traffic

Once individual kernels are on the GPU, the bottleneck shifts from compute to PCIe/NVLink data
movement. The strategy is to offload *enough of the per-step chain* that large arrays (BML
matrices, neighbor lists, force/potential buffers) can stay device-resident across the whole
SCF/MD step, so nothing is copied back to the host between kernels:

- **Persistent GPU neighbor list** (`60f0ad3` "Preserve nlist allocation. Use existing GPU
  pointers in Ewald real", `9f8a208`): the nlist is allocated/freed on the device at
  creation/destruction and its device pointers are reused in Ewald real — no per-step copyin.
- **Fewer transfers by resizing device buffers in place** (`baec46c` "Fewer data transfers in
  offload Ewald real"): instead of copyin/delete every step, the code tracks `maxnn` vs
  `maxnn_old` and only issues `!$acc exit data delete(... maxnn_old ...)` + `enter data
  create(... maxnn ...)` when the neighbor count actually grows, otherwise a cheap `!$acc
  update device`.
- **BML matrices passed by device pointer** (`9989a03`): `!$acc parallel loop
  deviceptr(aux_bml_ptr)` runs the charge kernel directly on the resident sparse matrix,
  avoiding host↔device copies of the large matrices entirely.
- **Offload the whole per-step chain** so residency pays off — Ewald real+recip, `get_hsmat`
  and H–S derivative, DM build / `prg_get_charges` / eval-vector, Pulay+SK forces, graph
  update (catalog in §1–§3). Each was ported specifically so the intermediate arrays never
  round-trip to the host.

### P3. Tune the OpenACC pragmas — pick the right loop level, enable collapse, drop redundant clauses

Getting a kernel *onto* the GPU is only the first step; the pragma structure then decides
occupancy and memory-access efficiency. Recurring tuning moves:

- **`collapse(2)` on the SK/Hamiltonian nested loops** (`c692490` "Opts for nvidia build":
  *"Change get_hsmat nested loop to enable collapse(2) … the collapse() clause improves
  performance"*) — flattens the atom×orbital nest into one larger parallel iteration space.
- **`worker` → `vector`** (`c052c49` "Change worker to vector in ACC", `e890f53` "Remove worker
  clause from openacc pragmas"): dropped the middle `worker` level in the Coulomb/Ewald loops
  in favor of two-level gang/vector, which mapped better to the hardware.
- **Explicit gang/vector sizing** (`a05ee53`): `!$acc parallel loop independent gang
  vector_length(64) num_gangs(1024)` on Ewald real.
- **Move `present`/`private` to match residency** (`285a541`): once params are persistent, the
  inner loop switches from `!$acc private(...)` (per-iteration temporaries) to
  `!$acc present(...)` (resident arrays), removing redundant allocation on the device.
- Gating offloaded pragmas behind `USE_OFFLOAD` when they regressed CPU builds (`8972061`,
  `dbce9d9`, `3cd01d1`) — a portability discipline rather than a speedup.

---

## 0. Pre-hackathon: CPU vectorization (early 2024, the "before" baseline)

Before any GPU offload, the hot kernels were restructured for vectorization — the groundwork
that later offload built on.

- **Vectorized H and S construction** (`85c0d00`, 2024-05-13, +515 lines in `ham_latte_mod.F90`):
  added `get_hsmat_vect` / `get_SKBlock_vect` alongside the scalar versions. The Slater–Koster
  block construction was reworked to fill a pre-allocated batched array
  `block(maxnorbi,maxnorbi,nats)` and dense `ham_vect/over_vect`, replacing per-pair scalar
  assembly — a memory-layout change that exposes SIMD parallelism over atoms/orbitals.
  Preceded by a `get_SKBlock` clean-up in preparation (`d610afc`).
- **Vectorized Hamiltonian derivative** (`fafc12e`, `13ce3da`, 2024-05-16/17;
  `hsderivative_latte_mod.F90` rewritten ±200 lines) and a follow-up reducing memory access in
  the vectorized path (`c5744a0`).
- **Vectorized Ewald/Coulomb real** (`0adbf75`→`f135398`, 2024-04-10–12): forces and potentials
  re-derived correctly in vector form, matrices replaced with vectors, and the slower scalar
  path removed (`d4cd4bd`). Earliest "GPU performance optimizations" commit (`fee2dc5`,
  2024-04-03) already trimmed `nonorthocoulombforces` (−69 lines, net −29).
- Removed memory contention in `prg_canon_response_p1_dpdmu` (`a3ac1d8`).

These `_vect` routines are the ones later offloaded (§1) — the CPU vectorization and the GPU
port share the same restructured, batched data layout.

## 1. GPU offload (OpenACC / OpenMP-offload, NVHPC) — the dominant thread

The single largest body of work: moving the per-MD-step hot path onto the GPU. Rolled out
routine-by-routine over ~2025, each verified numerically before the next.

- **Build/infra:** OpenMP-offload build for NVHPC (`95d0c14`, 2025-08-04); NVTX
  instrumentation for profiling (`905c009`, `89fd181`, 2024-07); NVHPC <25.0 back-compat
  (`dcdefd5`, `7833080`).
- **Ewald / Coulomb (real + reciprocal):** initial offload of Ewald real (`a05ee53`,
  `070e5c2`, 2025-08-11), shown faster than CPU (`0718b22`); nonortho-Coulomb offload and
  speedups (`305119e`, `8b4a120`, `cc3e19c`); reciprocal-Ewald memory footprint limited
  (`1f2de34`); fewer host↔device transfers (`baec46c`, 2025-10-07).
- **Density-matrix / electronic structure:** OpenACC kernel for `prg_get_charges`
  (`9989a03`, 2025-08-20; improved `6793462`); offload of `prg_get_evalsDvalsEvects`
  (`b25f46d`), Pulay forces (`4a94d6a`), SK-force (`a6c45c5`), `prg_buildzdiag`
  (`00408fe`); trace calcs (`d97e893`); atomic-density calc (`3828eb3`); elec-energy calc
  (`a3fba17`, `28f717f`).
- **Hamiltonian build / SK matrix:** offload of `get_hsmat` and H–S derivative
  (`0dc5329`, `e7a9d0d`, 2026-01-13); `get_nonorthocoul` optimizations (`bd146a2`).
- **Response / kernel:** OMP-offload `gpmdcov_response` (`fd55477`, 2024-07-29); MAGMA-pointer
  offload response kernel (`c96c9c8`, `9ce5a9d`, 2024-09).

**Concrete offload technique (from the diffs), for the methods section:**
- **OpenACC as the offload model** (`!$acc`), not raw CUDA. Kernels use explicit
  gang/vector tiling — e.g. Ewald real: `!$acc parallel loop independent gang
  vector_length(64) num_gangs(1024)` (`a05ee53`).
- **BML matrices passed by device pointer:** `!$acc parallel loop deviceptr(aux_bml_ptr)`
  in `prg_get_charges` (`9989a03`), avoiding host↔device copies of the large sparse matrices;
  inner charge accumulation uses `!$acc loop reduction(+:this_charge)`.
- **Explicit data regions** decouple transfer from compute: `!$acc enter data copyin(...)` /
  `!$acc exit data copyout(...)` around the force/potential arrays (`a05ee53`), so persistent
  arrays stay resident across the SCF loop.
- A documented serialization constraint worth citing as a limitation: the Ewald-real force/
  potential accumulation loop could not be fully parallelized without enlarging `raboff/droff`
  and adding reduction clauses (commit comment in `a05ee53`) — an honest note on remaining
  headroom.

Design pattern worth highlighting in the manuscript: a persistent-GPU-data strategy
(allocate once, reuse device pointers across steps) rather than per-call transfer — see
"persistent GPU nlist" and neighbor-list memory management below. This is principle **P2**
(data residency) realized in concrete pragmas; the pragma-level tuning is **P3**.

## 2. Neighbor list — new data structure + GPU residency

A ground-up neighbor-list rewrite, then made GPU-resident.

- New/faster neighbor list (`d7e1462`, `a1bfb1a`, `dea5981`, 2026-02-19–23); "only use new
  nlist for offload build" (`2d9f728`).
- **GPU-resident:** memory managed at creation/destruction on device (`41170fa`, 2026-02-27);
  allocation preserved across steps + existing GPU pointers reused in Ewald real (`60f0ad3`,
  `9f8a208`, 2026-03-02) — "persistent GPU nlist."
- Sparse neighbor list: compute-condition tuning and fixes (`11fd8d2`, `5108c3c` 2025-05;
  `1241294`, `5f0064c`, `256ae35`); SEDACS neighbor list adopted (`096cf95`, 2026-05-21).
- Max-density tuning to cut memory (`7a21ef4`, 2026-02-25) and later made a user knob:
  `GPMD{ MaxDensity= }` (`d456637`, 2026-08-17).

## 3. Graph partitioning / subgraph update — MPI communication + GPU

The O(N) graph machinery: less communication, faster updates, then offloaded.

- **MPI communication reduction:** subgraph graph-update method with less MPI comm
  (`46a1396`, 2024-08-08; groundwork `372c8ca`, `b2f6264`); more-efficient graph-update
  communication (`d828f83`, 2025-09-09).
- **Faster subgraph update:** `f8474d0`, `68719a8`, `93a6c9d` (2025-09-08–09); multiple
  parts per rank (`2fd9a86`, 2025-09-10).
- **New partitioners:** Box partitioning (`9391d45`, 2025-09-11; works with small subgraphs
  on Perlmutter `ec7cd9a`); SEDACS partition (`c365036`, 2024-07-16; fix `aec78d0`).
- **GPU:** offload of graph partitioning in progress (`3353115`, 2025-09-05); GPU-accelerated
  graph update (`ce78bdb`, 2026-02-26) and a numerically-correct offloaded version
  (`4d2b406`, 2026-03-10); several revert/guard commits show the offload was gated behind
  `USE_OFFLOAD` when it regressed CPU builds (`dbce9d9`, `3cd01d1`).

## 4. Memory footprint & single precision

Enabling larger systems per GPU. (Overlaps principle **P1**: several of these are
allocation-reduction/reuse commits viewed through the "fit more atoms per device" lens.)

- Single-precision build option (`8ac0fc0`, 2026-01-05) extended to Ewald real + graph
  partitioning (`5b74788`, 2026-01-09); lower-memory kernel method (`6b03dfc`).
- Decrease memory footprint / access (`d265e20`, `9a4e3bf` 2026-02-23; `d44b400`, `0b49975`
  rankN-update); reduce allocations in hcsf (`b3d5c84`); smaller kernel work array
  (`599e723`); GPU-deallocation forcing + memory reports for leak-hunting (`7e6c61f`,
  `7513382`, `006a1b7`, `981fe9b`, 2025-10).

## 5. Core numerics speedups (CPU-side)

- `dH` optimization (`533e37d`, 2025-01-31); vectorized Hamiltonian-derivative work predates
  the anchor but continued.
- Optimizations in deortho and DM build (`e7ad895`, `bd9a638`, 2026-03-05).
- Optimize µ (chemical potential) bisection search (`2e8e70d`, 2026-03-04).
- Threading additions (`db77553`, 2026-03-10; thread charge calc `89fd181`).

## 6. MPI restart & correctness enablers

Not speed per se, but enable the large/long runs the optimizations target.

- Working MPI-friendly restart (`1c0e538`, `791c8f0`, 2026-03-04).
- Rescale-restart-velocities option (`c096641`, 2026-03-30).
- MPI-rank correctness fixes that gate the adaptive scheme: sync XLBO charges across ranks
  (`67a98ce`), update `print_mdstep` on all ranks (`32be3f8`), off-by-one in replicate cell
  offsets (`4bbac0c`).

---

## 7. Adaptive time step + XLBO (the feature this manuscript centers on)

Separate from the HPC optimizations; this is the scientific contribution. Full detail lives
in `adaptive_time_step.md`. Arc:

- Adaptive "split-step" introduced and made compatible with minimization (`be1897d`,
  2026-05-05) and with the kernel/Langevin path (`0dc37ba`).
- K=5 variable-timestep XLBO coefficients derived (per-pattern tables), with several
  correctness fixes: missing C4/C5 lookups (`2cbae24`), Verlet coefficients from current dt
  not `dt_history(1)` (`7b2152b`), kappa-scaling fixes (`7f62eec`, `e1daea9`, `588af7a`).
- Dissipation (alpha) rule: pattern-specific → general history-independent rule matching
  full-step dissipation *rate per unit physical time* (`a279420`, `296f24a`, `ea7df42`);
  conservative alpha for large-system stability (`38439a3`, body: "~30% slower dissipation
  for much better stability", stable on 2088+ atoms).
- K=10 (11-point) kernel explored and shelved: coefficients kept for future, extended-history
  option removed (`259e75f`, 2026-07-11).

### 7a. LATE-BREAKING (Aug 2026, uncommitted analysis — verify before manuscript)

Two findings from the current session that revise the drift story and are **not yet in the
committed docs**:

1. **`abs()` split-trigger fix (`8a9c3b4`, on `split_step_clean`).** The adaptive split trigger
   used `maxval(user_timestep*sy%velocity)` — max over *signed* components — which missed fast
   atoms moving in −x/−y/−z, under-splitting them into repulsive-wall overshoots and a one-signed
   energy leak. Fixed to `user_timestep*maxval(abs(sy%velocity))`. On 300-atom water (ts=0.38,
   plateau-offset metric): fixed-step −0.001 eV, **legacy adaptive dropped to +0.055 eV** (was in
   the drift regime before). Strong evidence the previously-documented adaptive drift was
   *substantially an abs-bug artifact*, not the reversibility mechanism. Large-system (Mac1)
   confirmation run pending — doc baseline was +0.596 eV/ps.
2. **Reversible-backtracking scheme is worse, not better.** A symmetric-split (both-endpoint,
   backtracking) variant was built (`ReversibleSplit=T`, branch `reversible_split`) to remove the
   phase-locked drift. Measured plateau offset **+0.558 eV — ~10× legacy.** Rollback audit
   confirmed the checkpoint is complete (not a coding bug); the culprit is the intrinsic O(h³)
   "probe-endpoint ≠ kept-endpoint" leak, now measured to be large. Conclusion so far: the abs
   fix makes plain legacy adaptive conserve well, and the reversible scheme's motivation is
   largely gone. **Do not feature the reversible scheme as a win in the manuscript** unless the
   large-system data overturns this.

---

## Suggested manuscript framing

- **Methods-section spine — three optimization principles** (see "Cross-cutting optimization
  principles" above): **(P1)** minimize allocation cost by reusing/resizing existing buffers
  rather than re-allocating each step; **(P2)** prioritize whole-kernel GPU offload so large
  arrays stay device-resident and host↔device traffic is minimized; **(P3)** tune the OpenACC
  pragmas (loop level, `collapse`, gang/vector sizing, drop redundant clauses). Organizing the
  optimization narrative around these three is cleaner than a flat per-routine list — each
  subsystem (§1–§6) is then an *instance* of the principles.
- **Headline HPC contribution:** end-to-end GPU offload of an O(N) graph-based DFTB QMD driver
  (Ewald, DM build, H/S derivatives, forces, graph update), with persistent-GPU data structures
  and reduced MPI communication in the subgraph update — enabling [X]-atom systems at [Y] ns/day
  on [Venado/Perlmutter]. *(Fill X/Y from fresh benchmarks; history lacks them.)*
- **Scientific contribution:** adaptive split-step timestep with a history-independent XLBO
  dissipation rule that bounds per-step atomic displacement while preserving the conservative
  Verlet backbone to machine precision under variable stepping.
- **Honest energy-conservation section:** fixed-step baseline, the abs-trigger fix, and the
  measured limits of state-dependent splitting (and why backtracking does not help).
- **Benchmarks to (re)generate:** per-routine GPU speedup vs CPU; strong/weak scaling; ns/day vs
  system size; memory vs single/double precision. The commit history proves the *changes* but
  rarely the *magnitudes*.
