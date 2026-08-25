import sys, re

log = sys.argv[1] if len(sys.argv) > 1 else "out_restart_ts05_10k.log"
skip = int(sys.argv[2]) if len(sys.argv) > 2 else 50

t, e = [], []
pend_t = None
with open(log) as f:
    for line in f:
        m = re.search(r"Time \[fs\] =\s+([-0-9.Ee+]+)", line)
        if m:
            pend_t = float(m.group(1))
            continue
        m = re.search(r"Energy Total \[eV\] =\s+([-0-9.Ee+]+)", line)
        if m and pend_t is not None:
            t.append(pend_t)
            e.append(float(m.group(1)))
            pend_t = None

# de-duplicate by time (keep last energy at each unique time), keep monotone order
pairs = list(zip(t, e))
print(f"parsed {len(pairs)} (Time,Energy) records; time span {pairs[0][0]:.1f}..{pairs[-1][0]:.1f} fs")

pairs = pairs[skip:]
t = [p[0] for p in pairs]
e = [p[1] for p in pairs]
n = len(pairs)
print(f"after skipping first {skip}: {n} records, {t[0]:.1f}..{t[-1]:.1f} fs")

# split-fraction proxy: each user step emits one record if full, two if split (one per
# half), so consecutive-time gaps cluster at ~ts (full) and ~ts/2 (half). Classify each
# gap against the midpoint of the observed gap range and report the half fraction. This
# flags whether a run is in the discriminating INTERMITTENT regime (~20-40% split) vs the
# near-all-split (>~80%) or near-fixed (<~5%) limits, where the two schemes coincide.
if n > 2:
    gaps = sorted(t[i] - t[i-1] for i in range(1, n))
    gmin = gaps[len(gaps)//20]        # 5th percentile (robust min)
    gmax = gaps[-(len(gaps)//20)-1]   # 95th percentile (robust max)
    if gmax - gmin > 1e-6:            # two distinct gap sizes present -> mixed schedule
        thr = 0.5 * (gmin + gmax)
        n_half = sum(1 for g in (t[i]-t[i-1] for i in range(1, n)) if g < thr)
        # two half-records per split step, one full-record per full step:
        # split_steps = n_half/2 ; full_steps = (n-1) - n_half ; frac = split/(split+full)
        split_steps = n_half / 2.0
        full_steps = (n - 1) - n_half
        frac = split_steps / max(split_steps + full_steps, 1)
        regime = ("INTERMITTENT (discriminating)" if 0.15 <= frac <= 0.55
                  else "near-all-split (schemes coincide)" if frac > 0.55
                  else "near-fixed (schemes coincide)")
        print(f"split-fraction proxy: ~{frac*100:.0f}% of user steps split  "
              f"[half-gap~{gmin:.3f} full-gap~{gmax:.3f} fs]  -> {regime}")
    else:
        print(f"split-fraction proxy: uniform gap ~{gmin:.3f} fs (fixed schedule, no splits)")

def ols_slope(tt, ee):
    m = len(tt)
    mt = sum(tt)/m; me = sum(ee)/m
    num = sum((tt[i]-mt)*(ee[i]-me) for i in range(m))
    den = sum((tt[i]-mt)**2 for i in range(m))
    return num/den  # eV/fs

def stats(tt, ee, label):
    s = ols_slope(tt, ee) * 1000.0  # eV/ps
    ptp = max(ee) - min(ee)
    mean = sum(ee)/len(ee)
    rms = (sum((x-mean)**2 for x in ee)/len(ee))**0.5
    print(f"{label:22s} n={len(ee):5d}  {tt[0]:7.1f}-{tt[-1]:7.1f}fs  drift={s:+.4f} eV/ps  ptp={ptp:.4f}  rms={rms:.5f}")
    return s

full = stats(t, e, "FULL window")
# partition halves by PHYSICAL TIME, not record index: with variable/split steps the
# record cadence follows mdstep (one record per half), so record-index halves are skewed
# toward the split-dense stretch. Time-midpoint halves are cadence-independent.
t_mid = 0.5 * (t[0] + t[-1])
h = next((i for i in range(n) if t[i] >= t_mid), n//2)
h1 = stats(t[:h], e[:h], "1st half")
h2 = stats(t[h:], e[h:], "2nd half")
conv = "YES" if (h1*h2 > 0 and 0.33 < abs(h1/h2) < 3.0) else "NO"
print(f"convergence (halves agree in sign & within 3x): {conv}   (1st={h1:+.4f}, 2nd={h2:+.4f})")
