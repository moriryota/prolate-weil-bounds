"""Rigorous bound for the Liouville–Green error term of Lemma L5a(ii).
D := δ·c²t⁴ written in u = 1/t² ∈ (0, 1/4], s = χ/c² ∈ [0, 1/2], ε = 1/c² ∈ (0, 1/(10π)²]  (c ≥ 10π, i.e. μ ≥ 5).
We verify with Arb on a grid of boxes that |D| ≤ K, and that the denominator stays away from 0."""
import sympy as sp
from flint import arb, ctx
ctx.prec = 120
t, c, chi = sp.symbols('t c chi', positive=True)
u, s, e = sp.symbols('u s e', positive=True)
num = -2*c**4*t**4 + 3*c**4*t**2 + 6*c**2*chi*t**4 - 10*c**2*chi*t**2 + 2*c**2*chi - chi**2*t**2 + 2*chi**2
den = 4*(c**6*t**8 - c**6*t**6 - 3*c**4*chi*t**6 + 3*c**4*chi*t**4 + 3*c**2*chi**2*t**4 - 3*c**2*chi**2*t**2 - chi**3*t**2 + chi**3)
D = sp.simplify((num / den * c**2 * t**4).subs({chi: s * c**2}).subs({t: 1 / sp.sqrt(u), c: 1 / sp.sqrt(e)}))
D = sp.factor(sp.simplify(D))
print("D(u,s,e) =", D)
N, Dd = sp.fraction(sp.together(D))
N = sp.expand(N); Dd = sp.expand(Dd)
fN = sp.lambdify((u, s, e), N, modules=[{'sqrt': lambda x: x.sqrt()}]); fD = sp.lambdify((u, s, e), Dd, modules=[{'sqrt': lambda x: x.sqrt()}])
EPS = 1 / (100 * arb.pi()**2)
nu, ns, ne = 40, 20, 4
worst = arb(0); mind = None
for i in range(nu):
    U = arb((arb(i) + arb(i + 1)) / (8 * nu)).union(arb(i) / (4 * nu)).union(arb(i + 1) / (4 * nu))
    for j in range(ns):
        S = (arb(j) / (2 * ns)).union(arb(j + 1) / (2 * ns))
        for k in range(ne):
            E = (EPS * k / ne).union(EPS * (k + 1) / ne)
            n_, d_ = fN(U, S, E), fD(U, S, E)
            assert not d_.contains(0), (i, j, k)
            q = abs(n_ / d_)
            if q.upper() > worst.upper(): worst = q
            if mind is None or abs(d_).lower() < mind: mind = abs(d_).lower()
print("sup |D| over the box  <=", worst.upper().str(8), "   min |denominator| >=", mind.str(5))
