"""Numerical check of Slepian's identity in our normalisation:
  d/dc log(1 - lambda_n(c)) = -(lambda/(1-lambda)) * 2 psi_n(1)^2 / c   (equivalently d/dc log lambda = 2 psi(1)^2/c),
psi_n normalised in L^2(-1,1); lambda_n from the centre eigen-relation (as in kvector_endpoint.py)."""
import sys
sys.path.insert(0, 'common')
import mpmath as mp
from prolate_kvector import pswf_even, legendre_norm_vals
mp.mp.dps = 200
def lam_psi1(c, m):
    K = int(float(c) / 2) + 70
    ks, ps = pswf_even(c, K)
    coef = ps[m][1]
    P0 = legendre_norm_vals(mp.mpf(0), ks[-1]); P1 = legendre_norm_vals(mp.mpf(1), ks[-1])
    psi0 = sum(coef[i] * P0[k] for i, k in enumerate(ks)); psi1 = sum(coef[i] * P1[k] for i, k in enumerate(ks))
    chi = coef[0] * mp.sqrt(2) / psi0       # = int psi / psi(0) = mu_n (c=... ) -> chi_(n) = sqrt(lambda)*... see below
    # finite Fourier eigenvalue at 0: int psi = mu_n psi(0) with F_c psi = mu psi; lambda = (c/2pi) mu^2
    lam = c / (2 * mp.pi) * chi**2
    return lam, psi1
for mu in [5, 7]:
    c = 2 * mp.pi * mu
    for m, n in [(0, 0), (2, 4)]:
        h = mp.mpf('1e-8')
        lp, _ = lam_psi1(c + h, m); lm, _ = lam_psi1(c - h, m); l0, p1 = lam_psi1(c, m)
        lhs = (mp.log(1 - lp) - mp.log(1 - lm)) / (2 * h)
        rhs = -(l0 / (1 - l0)) * 2 * p1**2 / c
        print(f"mu={mu} n={n}: 1-lambda={mp.nstr(1 - l0, 6)}  d/dc log(1-lambda) = {mp.nstr(lhs, 10)}  vs  -(lam/(1-lam)) 2 psi(1)^2/c = {mp.nstr(rhs, 10)}  ratio {mp.nstr(lhs / rhs, 12)}", flush=True)
