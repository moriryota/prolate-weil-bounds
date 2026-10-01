"""Theorem C: numerical sanity checks of the intermediate inequalities (not part of the proof).
  Psi^2 <= 8 c rho_out ||h||^2           (Lemma L5b)
  h(0)^2 <= 2.0408 lam rho_out ||h||^2   (Lemma H0)
  lambda_n >= 395/399                     (Lemma S1e)"""
import sys
sys.path.insert(0, 'common')
import mpmath as mp
from prolate_kvector import pswf_even, legendre_norm_vals
for mu in [5, 7, 11]:
    mp.mp.dps = 160
    lam = mp.sqrt(mu); c = 2 * mp.pi * mu
    K = int(float(c) / 2) + 60
    ks, ps = pswf_even(c, K)
    P0 = legendre_norm_vals(mp.mpf(0), ks[-1]); P1 = legendre_norm_vals(mp.mpf(1), ks[-1])
    psi = lambda m, P: sum(ps[m][1][i] * P[k] for i, k in enumerate(ks))
    I0 = mp.sqrt(lam) * ps[0][1][0] * mp.sqrt(2); I4 = mp.sqrt(lam) * ps[2][1][0] * mp.sqrt(2)
    chi0 = I0 / (psi(0, P0) / mp.sqrt(lam)); chi2 = I4 / (psi(2, P0) / mp.sqrt(lam))
    al, be = I4, -I0
    h2 = al**2 + be**2
    rho = (al**2 * (1 - chi0**2) + be**2 * (1 - chi2**2)) / h2
    Psi = abs(al) * chi0 * abs(psi(0, P1)) + abs(be) * chi2 * abs(psi(2, P1))
    h0 = I0 * I4 * (chi2 - chi0) / (chi0 * chi2)
    h0d = (al * psi(0, P0) + be * psi(2, P0)) / mp.sqrt(lam)
    print(f"mu={mu:2d}: rho_out={mp.nstr(rho, 5)}  lambda_0={mp.nstr(chi0**2, 12)}  lambda_4={mp.nstr(chi2**2, 12)}  (>= 395/399={mp.nstr(mp.mpf(395)/399, 6)})")
    print(f"        Psi^2/(8c rho||h||^2) = {mp.nstr(Psi**2 / (8 * c * rho * h2), 5)}   "
          f"h(0)^2/(2.0408 lam rho||h||^2) = {mp.nstr(h0**2 / (mp.mpf('2.0408') * lam * rho * h2), 5)}   "
          f"[h(0) formula vs direct: {mp.nstr(h0, 8)} / {mp.nstr(h0d, 8)}]", flush=True)
