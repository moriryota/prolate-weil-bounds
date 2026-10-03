"""Paper III, Lemma 4.1: symbolic check (SymPy) of the Schroedinger form of the outer prolate equation, of the residual
of the two-term approximation w_0, and of the layer endpoint values (D_L, D_L', D_L'') = (1, 0, -c^4/2).
Usage: python proofs-III/cancellation_identities.py"""
import sympy as s
t,c,r=s.symbols('t c r',positive=True)
p=t*t-1
w=s.Function('w')(t)
f=w/s.sqrt(p)
lhs=s.simplify(s.sqrt(p)*(s.diff(f,t,2)+2*t/p*s.diff(f,t)+c*c*(t*t-r)/p*f))
V=c*c*(1-r)/p+1/p**2
assert s.simplify(lhs-s.diff(w,t,2)-(c*c+V)*w)==0
alpha=-c*(1-r)/2
w0=s.sin(c*t)+alpha/t*s.cos(c*t)
res=s.simplify(s.diff(w0,t,2)+(c*c+V)*w0)
expected=(V-c*c*(1-r)/t**2)*s.sin(c*t)+(2*alpha/t**3+V*alpha/t)*s.cos(c*t)
assert s.simplify(res-expected)==0
k,u=s.symbols('k u',positive=True)
D=2*s.exp(-k*(1-u*u))-s.exp(-2*k*(1-u*u))
assert s.simplify(D.subs(u,1))==1
assert s.simplify(s.diff(D,u).subs(u,1))==0
assert s.simplify(s.diff(D,u,2).subs(u,1))==-8*k*k
print('Schrodinger transformation residual = 0')
print('first asymptotic correction residual difference = 0')
print('layer endpoint (D, D_prime, D_second) = (1, 0, -8*k^2)')
print('alpha =',alpha)
print('V =',V)
print('All symbolic assertions passed.')

