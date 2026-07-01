# Branch Comparison Script

## Overview

The `compare_branches.sh` script automatically compares MD simulation behavior between two branches, building each branch, running simulations with specified parameters, and generating a detailed comparison report.

## Location

This script is located in `examples/gpmdk/` and works with any run directory under `examples/gpmdk/`.

## Usage

```bash
./compare_branches.sh <branch1> <branch2> <mdsteps> <timestep> [run_dir]
```

**Arguments:**
- `branch1`: First branch name (e.g., split_step)
- `branch2`: Second branch name (e.g., xlbo_adapt)
- `mdsteps`: Number of MD steps to run
- `timestep`: Timestep in femtoseconds
- `run_dir`: Directory containing input.in, relative to examples/gpmdk/ (default: run/water)

**Examples:**
```bash
# Use default run/water directory
./compare_branches.sh split_step xlbo_adapt 50 0.6

# Specify a different run directory
./compare_branches.sh split_step xlbo_adapt 100 0.35 run/ammonia

# Run from examples/gpmdk/ directory
cd examples/gpmdk
./compare_branches.sh split_step xlbo_adapt 50 0.4
```

## What the Script Does

1. **Checks out and builds each branch**
   - Automatically switches between branches
   - Builds each branch using `make -j4 install`
   - Returns to original branch when complete

2. **Modifies input parameters**
   - Backs up original `input.in`
   - Sets requested `TimeStep` and `MDSteps`
   - Restores original after completion

3. **Runs simulations**
   - Executes GPMD with `OMP_NUM_THREADS=4`
   - Captures full output for analysis

4. **Analyzes and compares results**
   - Extracts energy data from both runs
   - Counts split-step occurrences
   - Calculates energy conservation metrics
   - Computes differences between branches

5. **Generates reports**
   - Prints comparison table to console
   - Saves detailed results to file
   - Preserves output files for inspection

## Output Files

The script generates files in the specified run directory:
- `out_<branch1>_comparison` - Full simulation output for branch 1
- `out_<branch2>_comparison` - Full simulation output for branch 2
- `comparison_<branch1>_vs_<branch2>_ts<timestep>_md<mdsteps>.txt` - Detailed comparison report
- `input.in.backup` - Backup of original input.in (created if needed)

## Comparison Metrics

The script reports:
- **Total MD Steps**: Number of output steps produced
- **Split-steps triggered**: How many times timestep was split
- **Split-step range**: Which steps had splits
- **Mean Energy**: Average total energy
- **Std Dev**: Energy fluctuation magnitude
- **Total Drift**: Maximum energy change
- **Drift (% of E)**: Relative energy drift
- **Linear Drift**: Systematic energy drift rate
- **Maximum/Mean energy difference**: How much branches differ

## Interpretation

**Energy Conservation Quality:**
- < 0.01% drift = Excellent
- 0.01-0.1% drift = Good
- \> 0.1% drift = Poor (investigate)

**Branch Differences:**
- < 1e-10 eV: Identical (within machine precision)
- < 1e-6 eV: Essentially identical (numerical noise)
- < 0.001 eV: Small differences
- \> 0.001 eV: Measurable differences

## Example Output

```
===========================================================================
   COMPARISON at TimeStep=0.6 fs, MDSteps=50
===========================================================================

Metric                                      split_step         xlbo_adapt
---------------------------------------------------------------------------
Total MD Steps                                      50                 50
Split-steps triggered                               96                 96
Split-step range                                  3-98               3-98

Mean Energy (eV)                          -1360.058794       -1360.060119
Std Dev (eV)                                  0.009284           0.008485
Total Drift (eV)                              0.045340           0.043950
Drift (% of E)                                  0.0033             0.0032
Linear Drift (eV/step)                    1.604643e-04       1.332754e-04

Maximum energy difference                 3.390000e-03 eV
Mean absolute difference                  1.535800e-03 eV
```

## Requirements

- **Python**: Must have numpy installed
- **Git**: Repository must be a git repo
- **Build system**: CMake/Make setup must work
- **Both branches**: Must exist and be buildable

## Notes

- Script automatically backs up and restores `input.in`
- Returns to original branch on completion
- Handles build failures gracefully
- Safe to run multiple times
- Uses conda/system Python (not python3)

## Troubleshooting

**"Build failed"**: Check that both branches compile successfully manually first

**"No energy data found"**: Simulation may have crashed - check output files directly

**"operands could not be broadcast"**: Branches produced different numbers of output steps - ensure both have MDSteps fix applied
