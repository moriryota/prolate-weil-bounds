"""Paper IV, Section 8: deficit of the zeta zeros relative to the density (1+eps) a/pi, and deficit after thinning.
Diagnostic only: zeros from mpmath as floats, finite range; not part of any proof.
Usage: python numerics-IV/zero_deficit.py mu1 mu2 ..."""
import sys, math
import mpmath as mp
mp.mp.dps = 20
mus = [float(x) for x in sys.argv[1:]] or [5, 10, 20, 50]
Tcap = 1.6 * 2 * math.pi * math.e * max(mus)
g = []; n = 1
while True:
    t = float(mp.im(mp.zetazero(n)))
    if t > Tcap: break
    g.append(t); n += 1
print(f"# {len(g)} positive zeros up to {Tcap:.1f} (float; not interval-certified)")
def deficit(points, d, L):
    """max over intervals within [-L, L] of d*|I| - #points strictly inside; endpoints at points or at +-L."""
    P = [p for p in points if -L < p < L]; m = len(P); best = (0.0, None, None)
    ends = [-L] + P + [L]
    for i in range(len(ends)):
        for j in range(i + 1, len(ends)):
            val = d * (ends[j] - ends[i]) - (j - i - 1)
            if val > best[0]: best = (val, ends[i], ends[j])
    return best
def thin(points, a, q):
    out = []; last = -1e300
    for p in sorted(points):
        x = a * p / math.pi
        if x - last >= q: out.append(p); last = x
    return out
pts = sorted([-x for x in g] + g)
for mu in mus:
    a = 0.5 * math.log(mu)
    print(f"mu={mu}")
    for name, eps in (("0", 0.0), ("0.05", 0.05), ("1/log mu", 1 / math.log(mu))):
        d = (1 + eps) * a / math.pi; Te = 2 * math.pi * mu ** (1 + eps); L = min(1.5 * Te, Tcap)
        D, lo, hi = deficit(pts, d, L)
        print(f"  exp3 eps={name:8s} D={D:9.3f}  2mu^(1+eps)={2*mu**(1+eps):9.3f}  ratio={D/(2*mu**(1+eps)):.4f}  hi/T_eps={hi/Te:.3f}  (L={L:.0f})", flush=True)
    eps = 1 / math.log(mu); d = (1 + eps) * a / math.pi; Te = 2 * math.pi * mu ** (1 + eps); L = min(1.5 * Te, Tcap)
    base = [p for p in pts if -L < p < L]
    for qn, q in (("0.05", 0.05), ("0.01", 0.01), ("1/mu", 1 / mu), ("1/mu^2", 1 / mu ** 2)):
        th = thin(base, a, q); D, lo, hi = deficit(th, d, L)
        print(f"  exp4 q={qn:7s} removed={len(base)-len(th):4d}/{len(base)}  D={D:9.3f}  D/mu={D/mu:.4f}", flush=True)
