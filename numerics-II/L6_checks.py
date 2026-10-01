"""Checks for the Hermite limit (paper II, Section 7): (1) share of ||e_inf||^2 on |y| <= Y0;  (2) E_n(mu) = chi_n/mu vs 2pi(2n+1) (n = 0, 4);
(3) L2 distance ||h_{n,lam} - g_n|| and sup on |x| <= 2 (mu = 5, 7, 11, 15)."""
import sys
sys.path.insert(0, 'common')
import mpmath as mp
from prolate_kvector import pswf_even, legendre_norm_vals
mp.mp.dps = 30
def Ht(n, u): return mp.hermite(n, u) * mp.exp(-u**2 / 2) / mp.sqrt(2**n * mp.factorial(n) * mp.sqrt(mp.pi))
def g(n, x): return (2 * mp.pi) ** mp.mpf(0.25) * Ht(n, mp.sqrt(2 * mp.pi) * x)
I0, I4 = g(0, 0), g(4, 0)
hinf = lambda x: I4 * g(0, x) - I0 * g(4, x)
def einf(y):
    x = mp.exp(y); s = 0; n = 1
    while n * x < 9: s += hinf(n * x); n += 1
    return mp.exp(y / 2) * s
tot = 2 * mp.quad(lambda y: einf(y)**2, mp.linspace(0, 3, 13))
for Y0 in [0.5, 0.8, 1.0, 1.2]:
    print(f"share of ||e_inf||^2 on |y|<={Y0}: {mp.nstr(2 * mp.quad(lambda y: einf(y)**2, mp.linspace(0, Y0, 9)) / tot, 8)}", flush=True)
for mu in [5, 7, 11, 15]:
    mp.mp.dps = 160
    lam = mp.sqrt(mu); c = 2 * mp.pi * mu
    K = int(float(c) / 2) + 60
    ks, ps = pswf_even(c, K)
    for m, n in [(0, 0), (2, 4)]:
        coef = ps[m][1]
        def psi(s, coef=coef):
            P = legendre_norm_vals(s, ks[-1]); return sum(coef[i] * P[k] for i, k in enumerate(ks))
        sgn = 1 if psi(mp.mpf(0)) * g(n, 0) > 0 else -1
        h = lambda x: sgn * psi(x / lam) / mp.sqrt(lam) if abs(x) < lam else mp.mpf(0)
        mp.mp.dps = 30
        d2 = 2 * mp.quad(lambda x: (h(x) - g(n, x))**2, mp.linspace(0, lam, 12)) + 2 * mp.quad(lambda x: g(n, x)**2, [lam, mp.inf])
        sup = max(abs(h(x) - g(n, x)) for x in mp.linspace(0, 2, 81))
        chi = ps[m][0]
        print(f"mu={mu} n={n}: E_n=chi/mu={mp.nstr(mp.mpf(chi) / mu, 8)} (limit {mp.nstr(2 * mp.pi * (2 * n + 1), 8)})  ||h-g||_L2={mp.nstr(mp.sqrt(d2), 4)}  sup_|x|<=2|h-g|={mp.nstr(sup, 4)}", flush=True)
        mp.mp.dps = 160
