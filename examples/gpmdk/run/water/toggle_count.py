import sys, re

log = sys.argv[1]
# Reconstruct per-user-step split state from "Splitting print_mdstep N" markers
# and the "Mdstep, Energy" lines (one per completed user step).
split_steps = set()
mdsteps = []
for line in open(log):
    m = re.search(r"Splitting print_mdstep\s+(\d+)", line)
    if m:
        n = int(m.group(1))
        if n >= 1:
            split_steps.add(n)
        continue
    m = re.search(r"Mdstep, Energy.*?\s+(\d+)\s+([-0-9.]+)", line)
    if m:
        mdsteps.append(int(m.group(1)))

steps = sorted(set(mdsteps))
n = len(steps)
nsplit = sum(1 for s in steps if s in split_steps)
# count full<->split toggles across consecutive user steps
toggles = 0
for i in range(1, len(steps)):
    if steps[i] == steps[i-1] + 1:
        if (steps[i] in split_steps) != (steps[i-1] in split_steps):
            toggles += 1
print(f"{log}")
print(f"  user steps          : {n}")
print(f"  split steps         : {nsplit}  ({100*nsplit/max(n,1):.1f}%)")
print(f"  F<->S toggles       : {toggles}  ({toggles/max(n,1):.4f} per step)")
