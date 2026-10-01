"""Out-of-band mass of the CCM vector h_λ = α h_{0,λ} + β h_{4,λ} (∫h_λ = 0), from exact prolate identities:
  F(h_{2m,λ})|_{[-λ,λ]} = χ_m h_{2m,λ}  =>  χ_m = ∫h_{2m,λ} / h_{2m,λ}(0),
  ρ_out = ‖F(h_λ)·1_{|ξ|>λ}‖² / ‖h_λ‖² = (α²(1-χ_0²) + β²(1-χ_2²)) / (α² + β²)   (h_0, h_4 orthonormal).
Also ‖k_λ‖² (k_λ = E(h_λ) on the window) by piecewise Gauss quadrature, and h_λ(0).
Compared with W(k_λ)/‖k_λ‖² from the Galerkin runs of prolate_kvector.py (largest N computed)."""
import sys
sys.path.insert(0, 'common')
import mpmath as mp
from prolate_kvector import pswf_even, legendre_norm_vals
W_GALERKIN = {5: '1.39668e-17', 7: '9.32768e-28', 11: '1.28507e-48', 13: '3.49221e-59', 15: '8.71732e-70', 17: '1.86581e-80'}
def run(mu, dps=160):
    mp.mp.dps = dps
    lam = mp.sqrt(mu); a = mp.log(lam); c = 2 * mp.pi * mu
    K = int(float(c) / 2) + 60
    ks, ps = pswf_even(c, K)
    def val0(m):  # ψ_m(0)
        return sum(ps[m][1][i] * mp.sqrt((2 * k + 1) / mp.mpf(2)) * mp.legendre(k, 0) for i, k in enumerate(ks))
    # h_{n,λ}(x) = ψ_n(x/λ)/√λ:  ∫h = √λ coef_0 √2,  h(0) = ψ_n(0)/√λ
    I0 = mp.sqrt(lam) * ps[0][1][0] * mp.sqrt(2); I4 = mp.sqrt(lam) * ps[2][1][0] * mp.sqrt(2)
    chi0 = I0 / (val0(0) / mp.sqrt(lam)); chi2 = I4 / (val0(2) / mp.sqrt(lam))
    alpha, beta = I4, -I0
    rho = (alpha**2 * (1 - chi0**2) + beta**2 * (1 - chi2**2)) / (alpha**2 + beta**2)
    h0 = (alpha * val0(0) + beta * val0(2)) / mp.sqrt(lam)
    kmax = ks[-1]
    def h(x):
        if abs(x) >= lam: return mp.mpf(0)
        Pz = legendre_norm_vals(x / lam, kmax)
        return (alpha * sum(ps[0][1][i] * Pz[k] for i, k in enumerate(ks)) + beta * sum(ps[2][1][i] * Pz[k] for i, k in enumerate(ks))) / mp.sqrt(lam)
    def g(y):
        u = mp.e**y
        return mp.e**(y / 2) * sum(h(n * u) for n in range(1, int(mp.floor(lam / u)) + 1))
    breaks = sorted(set([-a, a] + [a - mp.log(n) for n in range(2, int(mu) + 1) if a - mp.log(n) > -a]))
    k2 = mp.mpf(0)
    for lo, hi in zip(breaks[:-1], breaks[1:]):
        k2 += mp.quad(lambda y: g(y)**2, [lo, hi])
    return dict(chi0=chi0, chi2=chi2, alpha=alpha, beta=beta, rho=rho, h2=alpha**2 + beta**2, k2=k2, h0=h0)
if __name__ == '__main__':
    fuchs = lambda mu: (mp.mpf(2)**14 / 3) * mp.sqrt(2) * mp.pi**5 * mp.e**(-4 * mp.pi * mu) * mp.mpf(mu)**4.5
    print("mu | 1-chi0 | 1-chi2 (exact) | Fuchs/exact | beta^2/(a^2+b^2) | rho_out | W/||k||^2 | (W/||k||^2)/rho_out | W/||h||^2 / rho_out | ||k||^2/||h||^2 | h(0)")
    for mu in [5, 7, 11, 13, 15, 17]:
        r = run(mu)
        W = mp.mpf(W_GALERKIN[mu]); Wn = W  # W(k)/||k||^2
        print(f"{mu:2d} | {mp.nstr(1-r['chi0'],4)} | {mp.nstr(1-r['chi2'],6)} | {mp.nstr(fuchs(mu)/(1-r['chi2']),5)} | {mp.nstr(r['beta']**2/r['h2'],4)} | "
              f"{mp.nstr(r['rho'],5)} | {mp.nstr(Wn,5)} | {mp.nstr(Wn/r['rho'],5)} | {mp.nstr(Wn*r['k2']/r['h2']/r['rho'],5)} | {mp.nstr(r['k2']/r['h2'],5)} | {mp.nstr(r['h0'],4)}", flush=True)
