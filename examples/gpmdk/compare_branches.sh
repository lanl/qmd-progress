#!/bin/bash
#
# Compare MD simulation behavior between two branches
#
# Usage: ./compare_branches.sh <branch1> <branch2> <mdsteps> <timestep> [run_dir]
#
# Example: ./compare_branches.sh split_step xlbo_adapt 50 0.4 run/water
#

set -e

# Check arguments
if [ "$#" -lt 4 ] || [ "$#" -gt 5 ]; then
    echo "Usage: $0 <branch1> <branch2> <mdsteps> <timestep> [run_dir]"
    echo "Example: $0 split_step xlbo_adapt 50 0.4 run/water"
    echo ""
    echo "Arguments:"
    echo "  branch1, branch2: Git branch names to compare"
    echo "  mdsteps:          Number of MD steps to run"
    echo "  timestep:         Timestep in femtoseconds"
    echo "  run_dir:          Directory containing input.in (default: run/water)"
    exit 1
fi

BRANCH1="$1"
BRANCH2="$2"
MDSTEPS="$3"
TIMESTEP="$4"
RUN_SUBDIR="${5:-run/water}"

# Script is in examples/gpmdk/, so repo root is ../..
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build"
RUN_DIR="${SCRIPT_DIR}/${RUN_SUBDIR}"
INPUT_FILE="${RUN_DIR}/input.in"

# Check that run directory exists
if [ ! -d "${RUN_DIR}" ]; then
    echo "ERROR: Run directory not found: ${RUN_DIR}"
    exit 1
fi

if [ ! -f "${INPUT_FILE}" ]; then
    echo "ERROR: input.in not found in: ${RUN_DIR}"
    exit 1
fi

# Save original branch
ORIGINAL_BRANCH=$(git rev-parse --abbrev-ref HEAD)

echo "=========================================="
echo "   Branch Comparison Tool"
echo "=========================================="
echo "Branch 1:   ${BRANCH1}"
echo "Branch 2:   ${BRANCH2}"
echo "MD Steps:   ${MDSTEPS}"
echo "TimeStep:   ${TIMESTEP} fs"
echo "Run Dir:    ${RUN_SUBDIR}"
echo ""

# Function to run simulation on a branch
run_branch() {
    local BRANCH="$1"
    local OUTPUT_FILE="$2"

    echo "=========================================="
    echo "Running ${BRANCH}"
    echo "=========================================="

    # Checkout branch
    echo "Checking out ${BRANCH}..."
    cd "${REPO_ROOT}"
    git checkout "${BRANCH}" 2>&1 | grep -v "^M\s" || true

    # Build
    echo "Building ${BRANCH}..."
    cd "${BUILD_DIR}"
    make -j4 install > /dev/null 2>&1 || {
        echo "ERROR: Build failed on ${BRANCH}"
        cd "${REPO_ROOT}"
        git checkout "${ORIGINAL_BRANCH}" 2>&1 | grep -v "^M\s" || true
        exit 1
    }

    # Modify input.in
    echo "Setting TimeStep=${TIMESTEP} and MDSteps=${MDSTEPS} in input.in..."
    cd "${RUN_DIR}"

    # Backup original input.in if not already backed up
    if [ ! -f input.in.backup ]; then
        cp input.in input.in.backup
    fi

    # Restore from backup and modify
    cp input.in.backup input.in
    sed -i.tmp "s/TimeStep=.*/TimeStep= ${TIMESTEP}/" input.in
    sed -i.tmp "s/MDSteps=.*/MDSteps= ${MDSTEPS}/" input.in
    rm -f input.in.tmp

    # Run simulation
    echo "Running simulation on ${BRANCH}..."
    OMP_NUM_THREADS=4 "${BUILD_DIR}/gpmdk" input.in > "${OUTPUT_FILE}" 2>&1 || {
        echo "ERROR: Simulation failed on ${BRANCH}"
        cd "${REPO_ROOT}"
        git checkout "${ORIGINAL_BRANCH}" 2>&1 | grep -v "^M\s" || true
        exit 1
    }

    echo "${BRANCH} complete!"
    echo ""
}

# Run both branches
cd "${REPO_ROOT}"
OUTPUT1="${RUN_DIR}/out_${BRANCH1}_comparison"
OUTPUT2="${RUN_DIR}/out_${BRANCH2}_comparison"

run_branch "${BRANCH1}" "${OUTPUT1}"
run_branch "${BRANCH2}" "${OUTPUT2}"

# Return to original branch
echo "Returning to ${ORIGINAL_BRANCH}..."
cd "${REPO_ROOT}"
git checkout "${ORIGINAL_BRANCH}" 2>&1 | grep -v "^M\s" || true

# Restore original input.in
if [ -f "${INPUT_FILE}.backup" ]; then
    cp "${INPUT_FILE}.backup" "${INPUT_FILE}"
fi

# Analyze results
echo "=========================================="
echo "   Analyzing Results"
echo "=========================================="
echo ""

cd "${RUN_DIR}"

# Export variables for Python
export BRANCH1="${BRANCH1}"
export BRANCH2="${BRANCH2}"
export TIMESTEP="${TIMESTEP}"
export MDSTEPS="${MDSTEPS}"

python << PYTHON_ANALYSIS
import numpy as np
import sys
import os

def analyze_with_splits(filename, label):
    """Analyze energy data and split-steps from output file"""
    energies = []
    split_steps = []

    if not os.path.exists(filename):
        print(f"ERROR: Output file not found: {filename}")
        sys.exit(1)

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("Mdstep, Energy"):
                parts = line.split()
                step = int(parts[5])
                energy = float(parts[6])
                energies.append((step, energy))
            elif "Splitting mdstep" in line:
                split_step = int(line.split()[-1])
                split_steps.append(split_step)

    if len(energies) == 0:
        print(f"ERROR: No energy data found in {filename}")
        sys.exit(1)

    energy_arr = np.array([e[1] for e in energies])

    e_mean = np.mean(energy_arr)
    e_std = np.std(energy_arr)
    e_min = np.min(energy_arr)
    e_max = np.max(energy_arr)
    e_drift = e_max - e_min
    e_drift_pct = (e_drift / abs(e_mean)) * 100

    timesteps = np.arange(len(energy_arr))
    coeffs = np.polyfit(timesteps, energy_arr, 1)
    slope = coeffs[0]

    return {
        'label': label,
        'steps': len(energy_arr),
        'mean': e_mean,
        'std': e_std,
        'drift': e_drift,
        'drift_pct': e_drift_pct,
        'slope': slope,
        'energies': energy_arr,
        'split_count': len(split_steps),
        'split_range': (min(split_steps), max(split_steps)) if split_steps else (None, None)
    }

# Get branch names from environment
branch1 = os.environ.get('BRANCH1', 'branch1')
branch2 = os.environ.get('BRANCH2', 'branch2')
timestep = os.environ.get('TIMESTEP', '?')
mdsteps = os.environ.get('MDSTEPS', '?')

# Analyze both outputs
data1 = analyze_with_splits(f'out_{branch1}_comparison', branch1)
data2 = analyze_with_splits(f'out_{branch2}_comparison', branch2)

# Print comparison table
print("="*75)
print(f"   COMPARISON at TimeStep={timestep} fs, MDSteps={mdsteps}")
print("="*75)
print("")
print(f"{'Metric':<35} {branch1:>18} {branch2:>18}")
print("-"*75)
print(f"{'Total MD Steps':<35} {data1['steps']:>18} {data2['steps']:>18}")
print(f"{'Split-steps triggered':<35} {data1['split_count']:>18} {data2['split_count']:>18}")
if data1['split_count'] > 0:
    range1 = f"{data1['split_range'][0]}-{data1['split_range'][1]}"
    range2 = f"{data2['split_range'][0]}-{data2['split_range'][1]}"
    print(f"{'Split-step range':<35} {range1:>18} {range2:>18}")
print("")
print(f"{'Mean Energy (eV)':<35} {data1['mean']:>18.6f} {data2['mean']:>18.6f}")
print(f"{'Std Dev (eV)':<35} {data1['std']:>18.6f} {data2['std']:>18.6f}")
print(f"{'Total Drift (eV)':<35} {data1['drift']:>18.6f} {data2['drift']:>18.6f}")
print(f"{'Drift (% of E)':<35} {data1['drift_pct']:>18.4f} {data2['drift_pct']:>18.4f}")
print(f"{'Linear Drift (eV/step)':<35} {data1['slope']:>18.6e} {data2['slope']:>18.6e}")
print("")

# Calculate differences
max_diff = np.max(np.abs(data1['energies'] - data2['energies']))
mean_diff = np.mean(np.abs(data1['energies'] - data2['energies']))
print(f"{'Maximum energy difference':<35} {max_diff:>18.6e} eV")
print(f"{'Mean absolute difference':<35} {mean_diff:>18.6e} eV")
print("")

# Conclusion
print("="*75)
print("                            CONCLUSION")
print("="*75)
print("")

if data1['split_count'] > 0:
    print(f"Split-steps triggered: {data1['split_count']} times", end="")
    if data1['split_range'][0]:
        print(f" (steps {data1['split_range'][0]}-{data1['split_range'][1]})")
    else:
        print()
    print("")

if max_diff < 1e-10:
    print("✓ Both methods produce IDENTICAL results")
elif max_diff < 1e-6:
    print("✓ Both methods produce essentially identical results")
    print(f"  (max difference {max_diff:.2e} eV - likely numerical noise)")
elif max_diff < 0.001:
    print(f"~ Methods show small differences (max {max_diff:.5f} eV)")
    print(f"  Mean difference: {mean_diff:.2e} eV")
else:
    print(f"⚠️  Methods show measurable differences:")
    print(f"  Max difference:  {max_diff:.5f} eV")
    print(f"  Mean difference: {mean_diff:.5f} eV")

print("")
print(f"Both branches maintain {'excellent' if max(data1['drift_pct'], data2['drift_pct']) < 0.01 else 'good'} energy conservation")
print(f"({branch1}: {data1['drift_pct']:.4f}%, {branch2}: {data2['drift_pct']:.4f}%)")
print("")

# Save results to file
output_file = f"comparison_{branch1}_vs_{branch2}_ts{timestep}_md{mdsteps}.txt"
with open(output_file, 'w') as f:
    f.write("="*75 + "\n")
    f.write(f"   COMPARISON: {branch1} vs {branch2}\n")
    f.write("="*75 + "\n")
    f.write(f"TimeStep: {timestep} fs\n")
    f.write(f"MDSteps:  {mdsteps}\n")
    f.write(f"Date:     {os.popen('date').read().strip()}\n")
    f.write("\n")
    f.write(f"{'Metric':<35} {branch1:>18} {branch2:>18}\n")
    f.write("-"*75 + "\n")
    f.write(f"{'Total MD Steps':<35} {data1['steps']:>18} {data2['steps']:>18}\n")
    f.write(f"{'Split-steps triggered':<35} {data1['split_count']:>18} {data2['split_count']:>18}\n")
    if data1['split_count'] > 0:
        range1 = f"{data1['split_range'][0]}-{data1['split_range'][1]}"
        range2 = f"{data2['split_range'][0]}-{data2['split_range'][1]}"
        f.write(f"{'Split-step range':<35} {range1:>18} {range2:>18}\n")
    f.write("\n")
    f.write(f"{'Mean Energy (eV)':<35} {data1['mean']:>18.6f} {data2['mean']:>18.6f}\n")
    f.write(f"{'Std Dev (eV)':<35} {data1['std']:>18.6f} {data2['std']:>18.6f}\n")
    f.write(f"{'Total Drift (eV)':<35} {data1['drift']:>18.6f} {data2['drift']:>18.6f}\n")
    f.write(f"{'Drift (% of E)':<35} {data1['drift_pct']:>18.4f} {data2['drift_pct']:>18.4f}\n")
    f.write(f"{'Linear Drift (eV/step)':<35} {data1['slope']:>18.6e} {data2['slope']:>18.6e}\n")
    f.write("\n")
    f.write(f"{'Maximum energy difference':<35} {max_diff:>18.6e} eV\n")
    f.write(f"{'Mean absolute difference':<35} {mean_diff:>18.6e} eV\n")

print(f"Results saved to: {output_file}")
print("")
PYTHON_ANALYSIS

echo "=========================================="
echo "   Comparison Complete"
echo "=========================================="
echo ""
echo "Output files:"
echo "  ${OUTPUT1}"
echo "  ${OUTPUT2}"
echo "  comparison_${BRANCH1}_vs_${BRANCH2}_ts${TIMESTEP}_md${MDSTEPS}.txt"
echo ""
    """Analyze energy data and split-steps from output file"""
    energies = []
    split_steps = []

    if not os.path.exists(filename):
        print(f"ERROR: Output file not found: {filename}")
        sys.exit(1)

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("Mdstep, Energy"):
                parts = line.split()
                step = int(parts[5])
                energy = float(parts[6])
                energies.append((step, energy))
            elif "Splitting mdstep" in line:
                split_step = int(line.split()[-1])
                split_steps.append(split_step)

    if len(energies) == 0:
        print(f"ERROR: No energy data found in {filename}")
        sys.exit(1)

    energy_arr = np.array([e[1] for e in energies])

    e_mean = np.mean(energy_arr)
    e_std = np.std(energy_arr)
    e_min = np.min(energy_arr)
    e_max = np.max(energy_arr)
    e_drift = e_max - e_min
    e_drift_pct = (e_drift / abs(e_mean)) * 100

    timesteps = np.arange(len(energy_arr))
    coeffs = np.polyfit(timesteps, energy_arr, 1)
    slope = coeffs[0]

    return {
        'label': label,
        'steps': len(energy_arr),
        'mean': e_mean,
        'std': e_std,
        'drift': e_drift,
        'drift_pct': e_drift_pct,
        'slope': slope,
        'energies': energy_arr,
        'split_count': len(split_steps),
        'split_range': (min(split_steps), max(split_steps)) if split_steps else (None, None)
    }

# Get branch names from environment
branch1 = os.environ.get('BRANCH1', 'branch1')
branch2 = os.environ.get('BRANCH2', 'branch2')
timestep = os.environ.get('TIMESTEP', '?')
mdsteps = os.environ.get('MDSTEPS', '?')

# Analyze both outputs
data1 = analyze_with_splits(f'out_{branch1}_comparison', branch1)
data2 = analyze_with_splits(f'out_{branch2}_comparison', branch2)

# Print comparison table
print("="*75)
print(f"   COMPARISON at TimeStep={timestep} fs, MDSteps={mdsteps}")
print("="*75)
print("")
print(f"{'Metric':<35} {branch1:>18} {branch2:>18}")
print("-"*75)
print(f"{'Total MD Steps':<35} {data1['steps']:>18} {data2['steps']:>18}")
print(f"{'Split-steps triggered':<35} {data1['split_count']:>18} {data2['split_count']:>18}")
if data1['split_count'] > 0:
    range1 = f"{data1['split_range'][0]}-{data1['split_range'][1]}"
    range2 = f"{data2['split_range'][0]}-{data2['split_range'][1]}"
    print(f"{'Split-step range':<35} {range1:>18} {range2:>18}")
print("")
print(f"{'Mean Energy (eV)':<35} {data1['mean']:>18.6f} {data2['mean']:>18.6f}")
print(f"{'Std Dev (eV)':<35} {data1['std']:>18.6f} {data2['std']:>18.6f}")
print(f"{'Total Drift (eV)':<35} {data1['drift']:>18.6f} {data2['drift']:>18.6f}")
print(f"{'Drift (% of E)':<35} {data1['drift_pct']:>18.4f} {data2['drift_pct']:>18.4f}")
print(f"{'Linear Drift (eV/step)':<35} {data1['slope']:>18.6e} {data2['slope']:>18.6e}")
print("")

# Calculate differences
max_diff = np.max(np.abs(data1['energies'] - data2['energies']))
mean_diff = np.mean(np.abs(data1['energies'] - data2['energies']))
print(f"{'Maximum energy difference':<35} {max_diff:>18.6e} eV")
print(f"{'Mean absolute difference':<35} {mean_diff:>18.6e} eV")
print("")

# Conclusion
print("="*75)
print("                            CONCLUSION")
print("="*75)
print("")

if data1['split_count'] > 0:
    print(f"Split-steps triggered: {data1['split_count']} times", end="")
    if data1['split_range'][0]:
        print(f" (steps {data1['split_range'][0]}-{data1['split_range'][1]})")
    else:
        print()
    print("")

if max_diff < 1e-10:
    print("✓ Both methods produce IDENTICAL results")
elif max_diff < 1e-6:
    print("✓ Both methods produce essentially identical results")
    print(f"  (max difference {max_diff:.2e} eV - likely numerical noise)")
elif max_diff < 0.001:
    print(f"~ Methods show small differences (max {max_diff:.5f} eV)")
    print(f"  Mean difference: {mean_diff:.2e} eV")
else:
    print(f"⚠️  Methods show measurable differences:")
    print(f"  Max difference:  {max_diff:.5f} eV")
    print(f"  Mean difference: {mean_diff:.5f} eV")

print("")
print(f"Both branches maintain {'excellent' if max(data1['drift_pct'], data2['drift_pct']) < 0.01 else 'good'} energy conservation")
print(f"({branch1}: {data1['drift_pct']:.4f}%, {branch2}: {data2['drift_pct']:.4f}%)")
print("")
PYTHON_ANALYSIS

echo "=========================================="
echo "   Comparison Complete"
echo "=========================================="
echo ""
echo "Output files:"
echo "  ${OUTPUT1}"
echo "  ${OUTPUT2}"
echo "  comparison_${BRANCH1}_vs_${BRANCH2}_ts${TIMESTEP}_md${MDSTEPS}.txt"
echo ""
