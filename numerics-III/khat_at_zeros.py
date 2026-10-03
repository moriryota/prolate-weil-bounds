"""Paper III, Section 8: |k^(gamma)|^2 at the first zeros versus the boundary-term prediction J^2/(1/4+gamma^2)
(CCM vector). Not part of any proof.  Usage: python numerics-III/khat_at_zeros.py mu NZ"""
import sys
sys.path.insert(0, 'common')
import mpmath as mp
from prolate_kvector import pswf_even, legendre_norm_vals
mu = int(sys.argv[1]); NZ = int(sys.argv[2]) if len(sys.argv) > 2 else 30
mp.mp.dps = 50
lam = mp.sqrt(mu); a = mp.log(lam); c = 2 * mp.pi * mu
K = int(float(c) / 2) + 60
ks, ps = pswf_even(c, K)
kmax = ks[-1]
P0 = legendre_norm_vals(mp.mpf(0), kmax)
lam_n = []
for m in (0, 2):
    co = ps[m][1]; psi0 = sum(co[i] * P0[k] for i, k in enumerate(ks)); ipsi = co[0] * mp.sqrt(2)
    lam_n.append(c / (2 * mp.pi) * (ipsi / psi0)**2)
I0 = mp.sqrt(lam) * ps[0][1][0] * mp.sqrt(2); I4 = mp.sqrt(lam) * ps[2][1][0] * mp.sqrt(2)
al, be = I4, -I0
B = al**2 * (1 - lam_n[0]) + be**2 * (1 - lam_n[1])
def h(x):
    if abs(x) >= lam: return mp.mpf(0)
    Pz = legendre_norm_vals(x / lam, kmax)
    return (al * sum(ps[0][1][i] * Pz[k] for i, k in enumerate(ks)) + be * sum(ps[2][1][i] * Pz[k] for i, k in enumerate(ks))) / mp.sqrt(lam)
def g(y):
    u = mp.e**y
    return mp.e**(y / 2) * sum(h(n * u) for n in range(1, int(mp.floor(lam / u)) + 1))
J = mp.e**(-a / 2) * sum(h(n / lam) for n in range(1, int(mu)))   # left limit at y = -a (n = mu term: h(lam) excluded, jump)
J = mp.e**(-a / 2) * (sum(h(n / lam) for n in range(1, int(mu))) + h(lam * (1 - mp.mpf(10)**-30)))  # include n = mu from inside
br = sorted(set([-a, a] + [a - mp.log(n) for n in range(2, int(mu) + 1)]))
X, W = mp.gauss_quadrature(60, 'legendre')
pts = []
for lo, hi in zip(br[:-1], br[1:]):
    for j in range(8):
        l2, h2 = lo + (hi - lo) * j / 8, lo + (hi - lo) * (j + 1) / 8
        for x, w in zip(X, W):
            y = (l2 + h2) / 2 + (h2 - l2) / 2 * x; pts.append((y, w * (h2 - l2) / 2, g(y)))
print(f"mu={mu}: B={mp.nstr(B, 6)}  J^2/B={mp.nstr(J**2 / B, 6)}  ||k||^2/||h||^2={mp.nstr(sum(w * v**2 for _, w, v in pts) / (al**2 + be**2), 6)}", flush=True)
Sk = Sb = mp.mpf(0)
for j in range(1, NZ + 1):
    gam = mp.im(mp.zetazero(j))
    kh = sum(w * v * mp.expj(gam * y) for y, w, v in pts)
    bh2 = J**2 / (mp.mpf(1) / 4 + gam**2)
    Sk += 2 * abs(kh)**2; Sb += 2 * bh2       # +- gamma
    if j <= 8 or j % 10 == 0:
        print(f"  j={j:3d} gamma={mp.nstr(gam, 8)}  |k^|^2/B={mp.nstr(abs(kh)**2 / B, 6)}  J^2/(1/4+g^2)/B={mp.nstr(bh2 / B, 6)}  ratio={mp.nstr(abs(kh)**2 / bh2, 6)}  "
              f"partial sums /B: k {mp.nstr(Sk / B, 6)}  b {mp.nstr(Sb / B, 6)}", flush=True)
