"""The mu -> infinity limit of the CCM vector.
Rescaled PSWF  h_{n,lam}(x) = psi_n(x/lam)/sqrt(lam)  ->  g_n(x) = (2pi)^{1/4} Htilde_n(sqrt(2pi) x)  (Hermite function, Fourier eigenvalue i^n = 1 for n = 0, 4).
h_inf = I4 g0 - I0 g4,  I_n = int g_n = g_n(0);  e_inf(y) = e^{y/2} sum_{n>=1} h_inf(n e^y).
Checks: h_inf(0) = 0, int h_inf = 0, e_inf(y) = e_inf(-y) (Poisson + Fourier invariance), kappa_inf = ||e_inf||^2/||h_inf||^2."""
import mpmath as mp
mp.mp.dps = 30
def Ht(n, u): return mp.hermite(n, u) * mp.exp(-u**2 / 2) / mp.sqrt(2**n * mp.factorial(n) * mp.sqrt(mp.pi))
def g(n, x): return (2 * mp.pi) ** mp.mpf(0.25) * Ht(n, mp.sqrt(2 * mp.pi) * x)
I0, I4 = g(0, 0), g(4, 0)
h = lambda x: I4 * g(0, x) - I0 * g(4, x)
print("I0, I4 =", mp.nstr(I0, 12), mp.nstr(I4, 12), "; h(0) =", mp.nstr(h(0), 5), "; int h =", mp.nstr(mp.quad(h, [-mp.inf, 0, mp.inf]), 5))
print("||h||^2 =", mp.nstr(I0**2 + I4**2, 12), "(check:", mp.nstr(mp.quad(lambda x: h(x)**2, [-mp.inf, 0, mp.inf]), 12), ")")
def e(y):
    x = mp.exp(y); s = mp.mpf(0); n = 1
    while True:
        t = h(n * x); s += t
        if n * x > 8 and abs(t) < mp.mpf(10)**-28: break
        n += 1
    return mp.exp(y / 2) * s
for y in [0.3, 1, 2, 3]:
    print(f"e({y}) = {mp.nstr(e(mp.mpf(y)), 12)},  e(-{y}) = {mp.nstr(e(mp.mpf(-y)), 12)}")
E = 2 * mp.quad(lambda y: e(y)**2, mp.linspace(0, 6, 25))
tail = mp.quad(lambda y: e(y)**2, [6, 7])
print("||e_inf||^2 (2*int_0^6) =", mp.nstr(E, 12), "; int_6^7 e^2 =", mp.nstr(tail, 3))
print("kappa_inf =", mp.nstr(E / (I0**2 + I4**2), 12))
