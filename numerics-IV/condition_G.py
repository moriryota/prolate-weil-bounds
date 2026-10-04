"""Paper IV, Section 8: deficit of the zeros after thinning (packing lower bound max(N - P, N^2/(N + 2P))), eps = 1/log mu,
for several separations q in the normalized variable x = a t/pi. Diagnostic only; not part of any proof.
Usage: python numerics-IV/condition_G.py mu1 mu2 ..."""
import sys, math
import mpmath as mp
mp.mp.dps = 20
mus = [float(x) for x in sys.argv[1:]] or [5, 10, 20]
Tcap = 1.6 * 2 * math.pi * math.e * max(mus)
g = []; n = 1
while True:
    t = float(mp.im(mp.zetazero(n)))
    if t > Tcap: break
    g.append(t); n += 1
print(f"# {len(g)} positive zeros up to {Tcap:.1f} (floats, finite range; diagnostic only)")
for mu in mus:
    a = 0.5 * math.log(mu); eps = 1 / math.log(mu)
    Te = 2 * math.pi * mu ** (1 + eps); L = min(1.5 * Te, Tcap)
    X = sorted([a * s * t / math.pi for t in g if t < L for s in (-1, 1)]); XL = a * L / math.pi
    m = len(X); ends = [-XL] + X + [XL]
    for qn, q in (("0", 0.0), ("1/log mu", 1 / math.log(mu)), ("(log mu)^-1/3", math.log(mu) ** (-1 / 3)), ("0.9", 0.9)):
        close = [1 if (X[k + 1] - X[k]) < q else 0 for k in range(m - 1)]   # adjacent close pairs (P counts pairs < q; adjacent dominate for q < typical gap)
        # cumulative count of all close pairs (not only adjacent) via two-pointer per start
        bestG = (0, None, None); bestN = 0
        for i in range(len(ends)):
            P = 0; jptr = None
            for j in range(i + 1, len(ends)):
                N = j - i - 1
                if N >= 1:
                    k = i + N - 1            # index in X of newest point inside: X[i .. i+N-1] (ends offset 1)
                    # pairs (k', k) with X[k]-X[k'] < q, k' >= i
                    kk = k - 1
                    while kk >= i and X[k] - X[kk] < q: P += 1; kk -= 1
                length = ends[j] - ends[i]
                p = max(N - P, (N * N / (N + 2 * P)) if N else 0)
                vG = (1 + eps) * length - p; vN = (1 + eps) * length - N
                if vG > bestG[0]: bestG = (vG, ends[i], ends[j])
                if vN > bestN: bestN = vN
        print(f"mu={mu:5.1f} q={qn:14s} D_G={bestG[0]:9.3f}  D_unthinned={bestN:9.3f}  D_G/mu={bestG[0]/mu:.3f}  worst I=({bestG[1]:.1f},{bestG[2]:.1f})", flush=True)
