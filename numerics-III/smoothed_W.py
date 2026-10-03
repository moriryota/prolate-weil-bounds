"""Paper III, Section 8 (second table): Galerkin values of W/||k||^2 for the CCM vector (s = 0) and the smoothed vectors
(s = 1, 2), with h_n~ = psi_n(x/lam) (1 - exp(-c^2(1-x^2/lam^2)/bl))^s / sqrt(lam). Not part of any proof.
Usage: python numerics-III/smoothed_W.py mu bl [N] [nq]"""
import sys
sys.path.insert(0, 'common')
import mpmath as mp
from flint import arb, ctx
from weil_spectral import WeilSpectral, gauss_legendre
from prolate_kvector import pswf_even, legendre_norm_vals
mu = int(sys.argv[1]); bl = mp.mpf(sys.argv[2])
N = int(sys.argv[3]) if len(sys.argv) > 3 else 140
nq = int(sys.argv[4]) if len(sys.argv) > 4 else 60
prec = 480
mp.mp.dps = int(prec * 0.30103) + 10; ctx.prec = prec
lam = mp.sqrt(mu); a = mp.log(lam); c = 2 * mp.pi * mu
K = int(float(c) / 2) + 60
ks, ps = pswf_even(c, K)
kmax = ks[-1]
lam_n = []
P0 = legendre_norm_vals(mp.mpf(0), kmax)
for m in (0, 2):
    co = ps[m][1]; psi0 = sum(co[i] * P0[k] for i, k in enumerate(ks)); ipsi = co[0] * mp.sqrt(2)
    lam_n.append(c / (2 * mp.pi) * (ipsi / psi0)**2)
xs, ws = gauss_legendre(nq)
xs = [mp.mpf(x.mid().str(mp.mp.dps, radius=False)) for x in xs]
ws = [mp.mpf(w.mid().str(mp.mp.dps, radius=False)) for w in ws]
W = WeilSpectral(arb(mu).sqrt().log(), N, prec=prec)
M = W.matrix()
Mm = [[mp.mpf(M[i, j].mid().str(mp.mp.dps, radius=False)) for j in range(N)] for i in range(N)]
wl = 40 * bl / c**2                     # layer zone in u = x/lam: [1 - wl, 1]
def pieces(lo, hi, cuts):
    pts = sorted(set([lo, hi] + [p for p in cuts if lo < p < hi]))
    return list(zip(pts[:-1], pts[1:]))
def integ(f, segs):
    tot = mp.mpf(0)
    for lo, hi in segs:
        half, mid = (hi - lo) / 2, (hi + lo) / 2
        tot += sum(w * half * f(mid + half * x) for x, w in zip(xs, ws))
    return tot
for s in (0, 1, 2):
    def layer(u): return (1 - mp.e**(-c**2 * (1 - u**2) / bl))**s if s else mp.mpf(1)
    def psis(u):
        P = legendre_norm_vals(u, kmax)
        return [sum(ps[m][1][i] * P[k] for i, k in enumerate(ks)) for m in (0, 2)]
    # integrals and norms of hn~ (in u-variable, then rescale)
    usegs = pieces(mp.mpf(-1), mp.mpf(1), [-1 + wl, 1 - wl, -1 + wl / 40, 1 - wl / 40])
    def vals(u):
        p = psis(u); L = layer(u); return [p[0] * L, p[1] * L]
    I = [integ(lambda u, m=m: vals(u)[m], usegs) * mp.sqrt(lam) for m in (0, 1)]
    G = [[integ(lambda u, i=i, j=j: vals(u)[i] * vals(u)[j], usegs) for j in (0, 1)] for i in (0, 1)]
    al, be = I[1], -I[0]
    hnorm2 = al**2 * G[0][0] + 2 * al * be * G[0][1] + be**2 * G[1][1]
    def h(x):
        if abs(x) >= lam: return mp.mpf(0)
        v = vals(x / lam); return (al * v[0] + be * v[1]) / mp.sqrt(lam)
    def g(y):
        u = mp.e**y
        return mp.e**(y / 2) * sum(h(n * u) for n in range(1, int(mp.floor(lam / u)) + 1))
    cuts = [a - mp.log(n) for n in range(1, int(mu) + 1)] + [a + mp.log(1 - wl) - mp.log(n) for n in range(1, int(mu) + 1)] \
         + [a + mp.log(1 - wl / 40) - mp.log(n) for n in range(1, int(mu) + 1)]
    ysegs = pieces(-a, a, cuts)
    coef = [mp.mpf(0)] * N; norm2 = mp.mpf(0)
    for lo, hi in ysegs:
        half, mid = (hi - lo) / 2, (hi + lo) / 2
        for x, w in zip(xs, ws):
            y = mid + half * x; gy = g(y); P = legendre_norm_vals(y / a, N - 1)
            for j in range(N): coef[j] += w * half * gy * P[j] / mp.sqrt(a)
            norm2 += w * half * gy**2
    proj2 = sum(x**2 for x in coef)
    Qv = sum(coef[i] * sum(Mm[i][j] * coef[j] for j in range(N)) for i in range(N))
    rho_out_ref = (al**2 * (1 - lam_n[0]) + be**2 * (1 - lam_n[1]))   # out-of-band mass of the UNSMOOTHED combination with these coefficients
    print(f"mu={mu} bl={bl} s={s}: W/||k||^2 = {mp.nstr(Qv / proj2, 8)}  ||k||^2/||h||^2 = {mp.nstr(norm2 / hnorm2, 6)}  "
          f"W/(rho_out_ref) = {mp.nstr(Qv / rho_out_ref, 6)}  proj_loss = {mp.nstr(1 - proj2 / norm2, 3)}", flush=True)
