"""Theorem C: rigorous evaluation of every numerical constant (python-flint Arb),
and numerical sanity checks of the intermediate quantities at mu = 5, 7, 11.
Upper bounds are printed from ball.upper() rounded up, lower bounds from ball.lower() rounded down;
each printed constant is asserted against the ball in the safe direction."""
import sys
sys.path.insert(0, 'python')
from flint import arb, ctx
ctx.prec = 200
PI = arb.pi()

def up_sci(x, nd=3):
    """scientific string with nd decimals in the mantissa, >= every point of the ball x"""
    from decimal import Decimal, ROUND_CEILING, getcontext
    getcontext().prec = 60
    d = Decimal(x.upper().str(45, radius=False).replace(' ', ''))
    q = Decimal(1).scaleb(d.adjusted() - nd)
    r = (d / q).to_integral_value(ROUND_CEILING) * q
    while not arb(str(r)) >= x:
        r += q
    return f"{r:.{nd}e}"

def up(x, nd=4):
    """decimal string >= every point of the ball x"""
    import math
    v = float(x.upper().mid()) if hasattr(x.upper(), 'mid') else float(x.upper())
    s = math.ceil(v * 10**nd) / 10**nd
    if not arb(str(s)) >= x:
        s = round(s + 10**-nd, nd)
    assert arb(str(s)) >= x, (s, x)
    return s

c0 = 10 * PI                      # c >= 10 pi; all constants below are decreasing in c
s_max = arb('0.3')

# (S1a'-1) eta: (pq)^{-1/4}(ct^2)^{1/2} - 1 = ((1-u)(1-su))^{-1/4} - 1 <= 4 f(1/4) u, u = 1/t^2 <= 1/4 (convex, f(0)=0)
f14 = ((1 - arb(1) / 4) * (1 - s_max / 4)) ** (arb(-1) / 4) - 1
print("eta coefficient 4 f(1/4) <=", up(4 * f14))

# (S1a'-2) LG tail: int_t^inf |delta| dz <= 0.504*sqrt(4/3)/(2 c t^2)
d_int = arb('0.5038') * (arb(4) / 3).sqrt() / 2
print("int |delta| dz <= D/(c t^2), D <=", up(d_int))
eps1 = d_int / c0 / 4 * 4        # times 1/t^2 coefficient: D/c
eps2 = (d_int / 2 / c0 / 4).exp() - 1   # at t = 2; |eps2| <= (e^{D/(2 c t^2)} - 1) <= coefficient/t^2
eps2_coef = eps2 * 4
print("eps1 coefficient (per 1/t^2) <=", up(eps1, 5), "; eps2 coefficient <=", up(eps2_coef, 5))
rho_coef = 4 * f14 + eps1 + eps2_coef + (f14) * eps2_coef
print("rho remainder coefficient (per 1/t^2) <=", up(rho_coef))
assert rho_coef <= arb('0.4')

# (S1a'-3) amplitude: E0 = w^2 + w_z^2 at t0 = 1 + 1/(2c^2)  (endpoint ODE bound, Lemma L5b proof, |psi| <= A0)
t0 = 1 + 1 / (2 * c0**2)
c2p = 1 + 1 / (4 * c0**2)          # c^2 p(t0)
r_lo, r_hi = t0**2 - s_max, t0**2
w2 = t0 * c2p.sqrt()               # sqrt(p q) <= c sqrt(p r) <= t0 sqrt(c^2 p)
pq14_lo = (c2p * r_lo) ** (arb(1) / 4)
term1 = t0**2 / 2 / pq14_lo
term2 = t0 / 2 * (1 / pq14_lo + c2p ** (arb(3) / 4) * r_lo ** (arb(-5) / 4) / c0**2)
E0 = w2 + (term1 + term2) ** 2
print("E0/A0^2 <=", up(E0), "; |w_z|/A0 <=", up(term1 + term2))
Ainf2 = arb('0.35').exp() * E0      # Gronwall with J <= 7/20 (Liouville-Green error J <= 7/20, Lemma L5b proof)
print("A_inf^2/A0^2 <=", up(Ainf2), "; |A_inf|/A0 <=", up(Ainf2.sqrt(), 3))
AMP = arb('1.78'); assert Ainf2.sqrt() <= AMP
AMP2 = AMP**2
near = AMP2 / arb('0.7').sqrt()     # |psi~(u)|^2 <= A0^2 * near/(c sqrt(u^2-1)), 1<u<2
print("near-field coefficient 1.78^2/sqrt(0.7) <=", up(near, 3))
assert near <= arb('3.79')

# (S1b) constants
lint = PI**2 / 12 + PI**3 / (3 * c0)
print("per-unit-interval int l(cu)^2 du <=", up(lint, 3)); assert lint <= arb('1.152')
acosh2 = (2 + arb(3).sqrt()).log()
bint = 2 * arb('3.79') * acosh2
print("int_1^2 b^2 coefficient (per Psi^2/c) <=", up(bint, 3)); assert bint <= arb('9.99')
A_ = 3 * (AMP2 * 2 * (2 * arb('1.152')) + arb('9.99'))
B_ = 3 * AMP2 * 2
print("M = (A + B K1^2) Psi^2/c + 0.75 h(0)^2/lam with A <=", up(A_, 2), ", B <=", up(B_, 2))
assert A_ <= arb('73.8') and B_ <= arb('19.02')
h0c = 2 * (arb(399) / 395) ** 2 * arb('0.75')
print("0.75 h(0)^2/lam <= (coef) rho||h||^2, coef <=", up(h0c, 3)); assert h0c <= arb('1.531')
MA = 8 * arb('73.8') + arb('1.531'); MB = 8 * arb('19.02')
print("M/(rho||h||^2) <= MA + MB K1^2, MA <=", up(MA, 1), ", MB <=", up(MB, 2))
assert MA <= arb('592') and MB <= arb('152.2')

# L: sum_m [min(2, kappa/m) + 0.4/m^2]/m, kappa = 2c/3; exact value vs closed form 2 log(c/3) + 5.49
def L_exact(c):
    kap = 2 * c / 3
    m0 = int((kap / 2).upper().floor().mid().__float__())
    s = arb(0)
    for m in range(1, m0 + 1): s += arb(2) / m
    N = 200000
    for m in range(m0 + 1, N + 1): s += kap / m**2
    s += kap / N                       # tail sum_{m>N} 1/m^2 <= 1/N
    return s + arb('0.4') * PI**2 / 6
print("L(c): exact sum vs closed form 2 log(c/3) + 5.49")
for mu in [5, 7, 11]:
    c = 2 * PI * mu
    print(f"  mu={mu}: exact <= {up(L_exact(c), 3)},  closed = {up(2 * (c / 3).log() + arb('5.49'), 3)}")
    assert L_exact(c) <= 2 * (c / 3).log() + arb('5.49')

# Theorem C factors
def Rtr(lT):  # Trudgian 2014 Cor. 1, T0 = e
    return arb('0.112') * lT + arb('0.278') * lT.log() + arb('2.510') + arb('0.2') / arb(1).exp()

def P_of(mu):
    mu = arb(mu); lam = mu.sqrt(); c = 2 * PI * mu
    lT = 20 * mu                                  # T = e^{20 mu}
    lT2 = lT + 1                                  # log(T+2) <= 20 mu + 1 (used for n* and r_T)
    nstar = (lT2 - (2 * PI).log()) / (4 * PI) + 2 * Rtr(lT2)
    rT = 1 / (arb('11.15') * lT2)
    pref = 2 * nstar / rT**2 * (arb('0.5') + lam / 2 * (1 / (2 * rT)).log())
    K1 = 2 * (c / 3).log() + arb('8.08')
    assert arb('5.49') + PI / 2 + 1 <= arb('8.08')
    P = pref * (arb('592') + arb('152.2') * K1**2)
    H = 2 * lam * (9 * mu**3)**2 * (lT + 1) / PI * (-lT).exp()
    return P, H, nstar, rT, pref, K1

obs = {5: 1.6232, 7: 1.8048, 11: 2.0471, 13: 2.1169, 15: 2.2512}
print("\n mu |   P(mu) <=   | P/observed | H(mu)e^{16mu} <= | n* <=  | K1 <=")
for mu in [5, 7, 11, 13, 15, 20, 50, 100, 1000]:
    P, H, ns, rT, pref, K1 = P_of(mu)
    o = obs.get(mu)
    He = H * (16 * arb(mu)).exp()
    assert He <= 1
    print(f" {mu:4d} | {up_sci(P):>10} | {up_sci(P / o) if o else '     -':>10} | {up_sci(He):>12}     | {up_sci(ns, 4):>10} | {up_sci(K1, 4)}")
import math
P1, *_ = P_of(1000); P2, *_ = P_of(2000)
print("local exponent d log P / d log mu at mu ~ 1400:", round(math.log(float(P2.mid()) / float(P1.mid())) / math.log(2), 3))

print("\nP(mu)/(mu^{7/2} (log mu)^3) on a grid (not a proof of monotonicity):")
for mu in [5, 6, 8, 10, 15, 20, 30, 50, 100, 300, 1000, 10**4, 10**6]:
    P, *_ = P_of(mu)
    print(f"  mu={mu}: {up_sci(P / (arb(mu)**arb('3.5') * arb(mu).log()**3))}")
