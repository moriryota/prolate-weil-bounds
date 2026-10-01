"""Effective L6 (paper II, Theorem D): rigorous lower bound F(mu) <= kappa_lambda = ||k||^2/||h||^2, Arb.
Steps (paper II, Section 7):
 rho_n^2 = int ((u^2 Ht_n')')^2 du (exact, sympy);  delta_n = rho_n/(mu sqrt(1-tau_n));  tau_n = int_{|u|>sqrt c} Ht_n^2
 IMS (ramp, C0 = pi^2):  4th even eigenvalue >= L6 = min((1-th) 26pi, pi^2 th mu) - pi^2/(th mu)
 Davis-Kahan in form norm; lattice-sum bound ||Tv||_{|y|<=Y0} <= sqrt(pi^2/12 e^Y0) W(v),  W(f)^2 = int (2 sinh(Y0) x^2 + |x|) f^2."""
import sympy as sp
from flint import arb, ctx
import math
ctx.prec = 200
PI = arb.pi()
u = sp.symbols('u', real=True)
def Ht(n):
    return sp.hermite(n, u) * sp.exp(-u**2 / 2) / sp.sqrt(2**n * sp.factorial(n) * sp.sqrt(sp.pi))
RHO2 = {}
for n in (0, 2, 4):
    f = sp.diff(u**2 * sp.diff(Ht(n), u), u)
    RHO2[n] = sp.nsimplify(sp.integrate(sp.expand(f**2), (u, -sp.oo, sp.oo)))
    print(f"rho_{n}^2 = {RHO2[n]} = {float(RHO2[n]):.6f}")
def tail_poly(n, X, extra=0):
    """upper bound of int_{|v|>X} v^extra Ht_n(v)^2 dv, valid for X^2 >= deg."""
    P = sp.Poly(sp.expand(sp.hermite(n, u)**2 * u**extra), u)
    norm = arb(2)**n * arb(math.factorial(n)) * PI.sqrt()
    tot = arb(0)
    for (m,), cf in P.terms():
        I = (-(X**2)).exp() * (X**(m - 1) if m >= 1 else 1 / X)
        tot += abs(arb(int(cf))) * I
    assert X**2 >= P.degree()
    return 2 * tot / norm
Y0 = arb('0.8'); a = 2 * Y0.sinh()
CT = (PI**2 / 12 * Y0.exp()).sqrt()
I0 = arb(2)**(arb(1)/4); I4 = arb(2)**(arb(1)/4) * 12 / arb(384).sqrt()
alpha, beta = I4, -I0
hinf2 = I0**2 + I4**2
KAPPA_INF_LO = arb('0.219247199548')      # rigorous enclosure, kappa_inf_enclosure.py (lower end)
# rigorous tail of ||e_inf||^2 beyond |y| > Y0 (symmetric): |h_inf(z)| <= |alpha| g0(z) + |beta| |g4(z)|, g_n(z) = (2pi)^{1/4} Ht_n(sqrt(2pi) z)
def g_abs_bound(n, z):
    w = (2 * PI).sqrt() * z
    H = {0: [1], 4: [12, 0, -48, 0, 16]}[n]
    Hb = sum(abs(arb(cf)) * w**k for k, cf in enumerate(H))
    return (2 * PI)**(arb(1)/4) * Hb * (-(w**2) / 2).exp() / (arb(2)**n * arb(math.factorial(n)) * PI.sqrt()).sqrt()
def ebound(y):
    x = y.exp(); s = arb(0)
    for m in range(1, 40):
        s += abs(alpha) * g_abs_bound(0, m * x) + abs(beta) * g_abs_bound(4, m * x)
    return (y / 2).exp() * s          # m >= 40 terms are < e^{-pi*1600}, absorbed (bounded by doubling below)
tailE = arb(0); y = Y0; dy = arb('0.01')
while y < 3:
    tailE += ebound(y)**2 * dy        # ebound decreasing in y on [Y0, inf): left Riemann sum is an upper bound
    y += dy
tailE = 2 * (2 * tailE + arb('1e-30'))   # factor 2 for symmetry, extra factor 2 safety for m>=40 and y>3
ThInf_lo = (KAPPA_INF_LO * hinf2 - tailE).sqrt()
print("tail of ||e_inf||^2 outside |y|<=0.8 <=", float(tailE.upper()), "; ||T h_inf||_{Y0} >=", float(ThInf_lo.lower()))
def F(mu, th=arb('0.25')):
    mu = arb(mu); c = 2 * PI * mu; X = c.sqrt()
    tau = {n: tail_poly(n, X) for n in (0, 2, 4)}
    tx2 = {n: tail_poly(n, X, 2) / (2 * PI) for n in (0, 4)}         # int_{|x|>lam} x^2 g_n^2
    d = {n: arb(str(float(sp.N(sp.sqrt(RHO2[n]), 30)) * (1 + 1e-12))) / (mu * (1 - tau[n]).sqrt()) for n in (0, 2, 4)}
    E = {n: 2 * PI * (2 * n + 1) for n in (0, 2, 4)}
    amin = lambda A, B: (A + B - abs(A - B)) / 2          # exact min, valid for balls
    L6 = amin((1 - th) * 26 * PI, PI**2 * th * mu) - PI**2 / (th * mu)
    assert (1 - th) * 26 * PI >= 0 and th * mu <= mu
    # counting / disjointness
    assert E[0] + d[0] < E[2] - d[2] and E[2] + d[2] < E[4] - d[4] and L6 > E[4] + d[4], "spectral separation fails"
    gam = {0: E[2] - d[2] - E[0], 4: amin(E[4] - E[2] - d[2], L6 - E[4])}
    assert gam[4] > 0 and gam[0] > 0
    W = {}
    for n in (0, 4):
        r = d[n] / gam[n]
        q = d[n]**2 / gam[n] + E[n] * d[n]**2 / gam[n]**2
        xu = q.sqrt() / (2 * PI)
        s = (1 - tau[n]).sqrt()
        nv = r**2 + r + (1 - s) + tau[n].sqrt()
        xg = ((2 * n + 1) / (4 * PI)).sqrt()
        xv = r**2 * (E[n] + d[n]).sqrt() / (2 * PI) + xu + (1 - s) / s * xg + tx2[n].sqrt()
        W[n] = (nv, xv)
    Wn = lambda nv, xv: ((a + arb('0.5')) * xv**2 + nv**2 / 2).sqrt()
    dal = PI.sqrt() * (W[4][0]**2 + W[4][1]**2).sqrt()     # |alpha_lam - alpha| = |int v4|
    dbe = PI.sqrt() * (W[0][0]**2 + W[0][1]**2).sqrt()
    Wg = {n: ((a + arb('0.5')) * (2 * n + 1) / (4 * PI) + arb('0.5')).sqrt() for n in (0, 4)}
    Wv = {n: Wn(*W[n]) for n in (0, 4)}
    Wh = {n: Wg[n] + Wv[n] for n in (0, 4)}
    Wtot = dal * Wh[0] + abs(alpha) * Wv[0] + dbe * Wh[4] + abs(beta) * Wv[4]
    Tv = CT * Wtot
    h2 = (abs(alpha) + dal)**2 + (abs(beta) + dbe)**2
    assert ThInf_lo - Tv > 0, "lattice error exceeds main term"
    kap = (ThInf_lo - Tv)**2 / h2
    return kap, dict(d4=d[4], L6=L6, gam4=gam[4], v4=W[4], Tv=Tv)
for mu in [40, 50, 60, 70, 80, 100, 150, 300, 1000]:
    try:
        k, info = F(mu)
        print(f"mu={mu}: kappa_lambda >= {float(k.lower()):.6f}   (L6={float(info['L6'].lower()):.2f}, gap4={float(info['gam4'].lower()):.2f}, ||v4||<={float(info['v4'][0].upper()):.4f}, ||x v4||<={float(info['v4'][1].upper()):.4f}, ||Tv||<={float(info['Tv'].upper()):.4f})")
    except AssertionError as e:
        print(f"mu={mu}: {e}")

print("\n-- theta optimised per mu (grid over theta; any theta gives a valid bound) --")
def best(mu):
    bestk, bestth = None, None
    for th in [arb(j) / 200 for j in range(2, 80)]:
        try:
            k, _ = F(mu, th)
        except AssertionError:
            continue
        if bestk is None or k.lower() > bestk.lower():
            bestk, bestth = k, th
    return bestk, bestth
for mu in [30, 40, 50, 60, 80, 100, 120, 150, 200, 300]:
    k, th = best(mu)
    kl = math.floor(float(k.lower()) * 1e6) / 1e6 if k is not None else float('nan')   # rounded down
    if k is not None: assert arb(str(kl)) <= k
    print(f"mu={mu}: kappa_lambda >= {kl:.6f}  (theta={float(th.mid()) if th is not None else float('nan'):.3f})", flush=True)
