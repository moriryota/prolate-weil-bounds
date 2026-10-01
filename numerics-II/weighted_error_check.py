"""Numerics for the weighted Davis-Kahan bound: ||x (h_{n,lam} - g_n)|| and ||(1+|x|)(h - g)|| for mu = 5, 7, 11, 15 (n = 0, 4)."""
import sys
sys.path.insert(0, 'common')
import mpmath as mp
from prolate_kvector import pswf_even, legendre_norm_vals
def Ht(n, u): return mp.hermite(n, u) * mp.exp(-u**2 / 2) / mp.sqrt(2**n * mp.factorial(n) * mp.sqrt(mp.pi))
def g(n, x): return (2 * mp.pi) ** mp.mpf(0.25) * Ht(n, mp.sqrt(2 * mp.pi) * x)
for mu in [5, 7, 11, 15]:
    mp.mp.dps = 160
    lam = mp.sqrt(mu); c = 2 * mp.pi * mu
    ks, ps = pswf_even(c, int(float(c) / 2) + 60)
    for m, n in [(0, 0), (2, 4)]:
        coef = ps[m][1]
        def psi(s, coef=coef):
            P = legendre_norm_vals(s, ks[-1]); return sum(coef[i] * P[k] for i, k in enumerate(ks))
        sgn = 1 if psi(mp.mpf(0)) * g(n, 0) > 0 else -1
        h = lambda x: sgn * psi(x / lam) / mp.sqrt(lam) if abs(x) < lam else mp.mpf(0)
        mp.mp.dps = 30
        pts = mp.linspace(0, lam, 14)
        v2 = 2 * mp.quad(lambda x: (h(x) - g(n, x))**2, pts) + 2 * mp.quad(lambda x: g(n, x)**2, [lam, mp.inf])
        xv2 = 2 * mp.quad(lambda x: x**2 * (h(x) - g(n, x))**2, pts) + 2 * mp.quad(lambda x: x**2 * g(n, x)**2, [lam, mp.inf])
        w2 = 2 * mp.quad(lambda x: (1 + x)**2 * (h(x) - g(n, x))**2, pts) + 2 * mp.quad(lambda x: (1 + x)**2 * g(n, x)**2, [lam, mp.inf])
        print(f"mu={mu} n={n}: ||v||={mp.nstr(mp.sqrt(v2), 4)} (x mu = {mp.nstr(mu * mp.sqrt(v2), 4)})  ||xv||={mp.nstr(mp.sqrt(xv2), 4)} (x mu = {mp.nstr(mu * mp.sqrt(xv2), 4)})  ||(1+|x|)v||={mp.nstr(mp.sqrt(w2), 4)}", flush=True)
        mp.mp.dps = 160
