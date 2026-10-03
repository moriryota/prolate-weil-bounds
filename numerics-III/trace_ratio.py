"""Paper III, Remark 5.2: J^2/B for the smoothed vector and its band-edge Poisson term J1 = lam f^(c). Not part of any proof.
Usage: python numerics-III/trace_ratio.py mu [beta] [dps]"""
import sys
sys.path.insert(0, 'common')
import mpmath as mp
from prolate_kvector import pswf_even, legendre_norm_vals
mu = mp.mpf(sys.argv[1]); beta = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf(4)
mp.mp.dps = int(sys.argv[3]) if len(sys.argv) > 3 else 80
c = 2 * mp.pi * mu; lam = mp.sqrt(mu)
K = int(float(c) / 2) + 40
ks, ps = pswf_even(c, K, which=(0, 2))
L = 2 * K + 80
xg, wg = mp.gauss_quadrature(2 * L, 'legendre')
Pg = [legendre_norm_vals(x, L) for x in xg]
def lay(u): return (1 - mp.e**(-c**2 * (1 - u**2) / beta))**2
def fn(u, Pv=None):
    Pv = Pv or legendre_norm_vals(u, ks[-1])
    l = lay(u)
    return [sum(ps[m][1][i] * Pv[k] for i, k in enumerate(ks)) * l for m in (0, 2)]
vals = [fn(x, Pv) for x, Pv in zip(xg, Pg)]
I = [sum(w * v[m] for w, v in zip(wg, vals)) for m in (0, 1)]
al, be = I[1], -I[0]
fv = [al * v[0] + be * v[1] for v in vals]
co = [sum(w * v * Pv[k] for w, v, Pv in zip(wg, fv, Pg)) for k in range(L + 1)]
nrm = sum(x**2 for x in co)
X, Wq = mp.gauss_quadrature(24, 'legendre')
kk = [k for k in range(len(co)) if k % 2 == 0]
def fh(w):
    return sum(co[k] * mp.sqrt((2 * k + 1) / mp.mpf(2)) * 2 * (-1)**(k // 2) * mp.sqrt(mp.pi / (2 * w)) * mp.besselj(k + mp.mpf(1) / 2, w) for k in kk)
tot = mp.mpf(0); nseg = int(mp.ceil(c / (mp.pi / 4)))
for j in range(nseg):
    lo, hi = c * j / nseg, c * (j + 1) / nseg
    tot += sum(wq * fh((hi + lo) / 2 + (hi - lo) / 2 * x) ** 2 for x, wq in zip(X, Wq)) * (hi - lo) / 2
B = nrm - 2 * tot / (2 * mp.pi)
J = sum(al * v[0] + be * v[1] for v in (fn(n / mu) for n in range(1, int(mp.ceil(mu)))) ) / lam
f0 = al * fn(mp.mpf(0))[0] + be * fn(mp.mpf(0))[1]
J1 = lam * fh(c)   # m = 1 Poisson term at X = lam (exactly the band edge)
print(f"   J1 = lam f^(c): J1^2/B = {mp.nstr(J1**2 / B, 8)}  J1/J = {mp.nstr(J1 / J, 6)}  J^2/(cB) = {mp.nstr(J**2 / (c * B), 6)}", flush=True)
print(f"mu={mp.nstr(mu, 6)} beta={beta}: B/||f||^2 = {mp.nstr(B / nrm, 8)}  J^2/||f||^2 = {mp.nstr(J**2 / nrm, 8)}  J^2/B = {mp.nstr(J**2 / B, 8)}  "
      f"h(0)^2/(lam B) = {mp.nstr(f0**2 / lam / B, 6)}  int f check = {mp.nstr(sum(w * v for w, v in zip(wg, fv)), 3)}", flush=True)
