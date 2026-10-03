"""Paper III, Section 8: direct quadrature of g^_{n,L}(lam t) compared with the bound of Lemma 4.1.
Not part of any proof.  Usage: python numerics-III/cancellation_direct.py mu n dps [n_bulk_segments]"""
import sys
sys.path.insert(0, 'common')
import mpmath as mp
from prolate_kvector import pswf_even, legendre_norm_vals
mu = int(sys.argv[1]); nn = int(sys.argv[2]); mp.mp.dps = int(sys.argv[3])
c = 2 * mp.pi * mu; lam = mp.sqrt(mu); beta = 4
K = int(float(c) / 2) + 60
ks, ps = pswf_even(c, K, which=(nn // 2,))
co = ps[nn // 2][1]; kmax = ks[-1]
def psi_dpsi(u):
    # P-bar_k and derivative via recurrence: P'_k from (1-u^2)P'_k = k(P_{k-1} - u P_k)
    P = [mp.mpf(1), u]
    for n in range(1, kmax): P.append(((2 * n + 1) * u * P[n] - n * P[n - 1]) / (n + 1))
    dP = [mp.mpf(0)] + [k * (P[k - 1] - u * P[k]) / (1 - u**2) if abs(1 - u**2) > 0 else mp.mpf(k * (k + 1)) / 2 for k in range(1, kmax + 1)]
    s = sum(co[i] * mp.sqrt((2 * k + 1) / mp.mpf(2)) * P[k] for i, k in enumerate(ks))
    ds = sum(co[i] * mp.sqrt((2 * k + 1) / mp.mpf(2)) * dP[k] for i, k in enumerate(ks))
    return s, ds
def integrand_val(u):
    s, ds = psi_dpsi(u)
    E = mp.e**(-c**2 * (1 - u**2) / beta); L = (1 - E)**2
    dL = 2 * (1 - E) * (-E) * (c**2 * 2 * u / beta)       # d/du of (1-E)^2, E' = E * (c^2 2u / beta)
    return u * (ds * L + s * dL)
A = abs(psi_dpsi(mp.mpf(1) - mp.mpf(10)**(-mp.mp.dps + 5))[0])
# nodes: bulk [0, 1-w] and layer [1-w, 1], w = 60 beta / c^2, subdivided; tabulate integrand once
w = 60 * mp.mpf(beta) / c**2
X, W = mp.gauss_quadrature(40, 'legendre')
def segs(lo, hi, n): return [(lo + (hi - lo) * j / n, lo + (hi - lo) * (j + 1) / n) for j in range(n)]
NB = int(sys.argv[4]) if len(sys.argv) > 4 else 400
S = segs(mp.mpf(0), 1 - w, NB) + segs(1 - w, mp.mpf(1), 60 if NB == 400 else NB)
nodes = []
for lo, hi in S:
    for x, wq in zip(X, W):
        u = (lo + hi) / 2 + (hi - lo) / 2 * x; nodes.append((u, wq * (hi - lo) / 2, integrand_val(u)))
Cf = 20 + 321 * mp.e**8
print(f"mu={mu} n={nn} dps={mp.mp.dps} c={mp.nstr(c, 8)} A=|psi(1)|={mp.nstr(A, 8)}  nodes={len(nodes)}", flush=True)
for tm in ((1, 2, 5, 10, 20, 40) if NB == 400 else (20, 40, 80)):
    t = tm * c
    val = 2 * mp.sqrt(lam) * sum(wq * f * mp.cos(c * t * u) for u, wq, f in nodes)
    lead = 2 * mp.sqrt(lam) * A * c**2 / (2 * c * t)      # |2 sqrt(lam) psi'(1)/(ct)| ~ sqrt(lam) A c (1-s)/t, crude size
    bnd = Cf * mp.sqrt(lam) * A * c**2 / t**2
    print(f"  t={tm:3d}c: |g^|={mp.nstr(abs(val), 8)}  bound={mp.nstr(bnd, 8)}  ratio={mp.nstr(abs(val) / bnd, 6)}  (size of leading term ~{mp.nstr(lead, 4)})", flush=True)
