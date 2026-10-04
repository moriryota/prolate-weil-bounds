"""Paper IV, Section 7: largest ratio P_h(J)/(h l L^2 + l + L^2) of Hypothesis (LP) over scanned windows, for the first NZ zeros.
Diagnostic only (floats, finite height); not part of any proof.
Usage: python numerics-IV/local_pairs.py NZ"""
import sys, math, bisect
import mpmath as mp
mp.mp.dps = 20
NZ = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
g = [float(mp.im(mp.zetazero(n))) for n in range(1, NZ + 1)]
print(f"# {NZ} zeros up to t = {g[-1]:.1f} (floats; diagnostic)")
gaps = sorted((g[k + 1] - g[k], g[k]) for k in range(NZ - 1))
print("# closest pairs (gap, t):", [(round(a, 5), round(b, 2)) for a, b in gaps[:5]])
for h in (0.4, 0.1, 0.03, 0.01):
    best = (0, None)
    for l in (0.5, 2, 10, 50, 250, 1000):
        step = max(l / 4, 0.25)
        t0 = 0.0
        while t0 + l <= g[-1]:
            i = bisect.bisect_left(g, t0); j = bisect.bisect_right(g, t0 + l)
            P = 0
            for k in range(i, j):
                m = k + 1
                while m < j and g[m] - g[k] < h: P += 1; m += 1
            L = math.log(2 + t0 + l)
            r = P / (h * l * L * L + l + L * L)
            if r > best[0]: best = (r, (round(t0, 2), l, P))
            t0 += step
    print(f"h={h:5.2f}: max P_h(J)/(h l L^2 + l + L^2) = {best[0]:.4f}  at (t0, l, P) = {best[1]}", flush=True)
