"""Residual of the Hermite function g_n for the rescaled prolate operator
A_mu = -(d/dx)(1 - x^2/mu)(d/dx) + 4 pi^2 x^2 :  (A_mu - E_n^inf) g_n = (1/mu) (x^2 g_n')'  on (-lam, lam).
rho_n := ||(x^2 g_n')'||_{L^2(R)};  Davis-Kahan predicts ||h_n - g_n|| <~ rho_n/(mu * gap), gap = 8 pi (even sector)."""
import mpmath as mp
mp.mp.dps = 30
def Ht(n, u): return mp.hermite(n, u) * mp.exp(-u**2 / 2) / mp.sqrt(2**n * mp.factorial(n) * mp.sqrt(mp.pi))
def g(n, x): return (2 * mp.pi) ** mp.mpf(0.25) * Ht(n, mp.sqrt(2 * mp.pi) * x)
for n in [0, 2, 4, 6]:
    f = lambda x: mp.diff(lambda y: y**2 * mp.diff(lambda z: g(n, z), y), x)
    rho = mp.sqrt(2 * mp.quad(lambda x: f(x)**2, [0, 1, 2, 4, 8]))
    print(f"n={n}: rho_n = ||(x^2 g')'|| = {mp.nstr(rho, 6)};  rho/(8 pi) = {mp.nstr(rho / (8 * mp.pi), 5)}", flush=True)
print("observed (L6_checks.py): ||h_0-g_0|| ~ 0.025/mu, ||h_4-g_4|| ~ 0.22/mu ; eigenvalue error ~ 0.76/mu (n=0), 11/mu (n=4)")
