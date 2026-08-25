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
h = n//2
h1 = stats(t[:h], e[:h], "1st half")
h2 = stats(t[h:], e[h:], "2nd half")
conv = "YES" if (h1*h2 > 0 and 0.33 < abs(h1/h2) < 3.0) else "NO"
print(f"convergence (halves agree in sign & within 3x): {conv}   (1st={h1:+.4f}, 2nd={h2:+.4f})")
