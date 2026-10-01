"""Endpoint value h_λ(λ^-) of the CCM vector versus its out-of-band mass:
heuristic  ĥ_λ(ξ) ≈ h_λ(λ) sin(2πλξ)/(πξ) for |ξ| > λ  gives  ‖ĥ_λ 1_{|ξ|>λ}‖² ≈ h_λ(λ)²/(π² λ),
i.e. the ratio  h_λ(λ)² / (π² λ ρ_out ‖h_λ‖²)  should be close to 1."""
import sys
sys.path.insert(0, 'common')
import mpmath as mp
from prolate_kvector import pswf_even, legendre_norm_vals
for mu in [5, 7, 11, 13, 15]:
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
    hl = (al * psi(0, P1) + be * psi(2, P1)) / mp.sqrt(lam)
    print(f"mu={mu:2d}: h(lam-) = {mp.nstr(hl, 6)};  h(lam)^2/(pi^2 lam rho ||h||^2) = {mp.nstr(hl**2 / (mp.pi**2 * lam * rho * h2), 6)};  "
          f"psi4(1)^2/(c(1-chi2^2)) = {mp.nstr(psi(2, P1)**2 / (c * (1 - chi2**2)), 6)}", flush=True)
