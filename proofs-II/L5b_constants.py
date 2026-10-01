from fractions import Fraction as Q
import math
import sympy as S
from flint import arb, ctx
ctx.prec=320
t,s,c=S.symbols('t s c',positive=True)
p=t*t-1;r=t*t-s;q=c*c*r
ell=-(S.diff(p,t)/p+S.diff(q,t)/q)/4
delta=S.factor((p*(S.diff(ell,t)+ell**2)+S.diff(p,t)*ell)/q)
D=s*s+2*s+(3-6*s)*t*t
A=1/(4*c*c*p*r);B=D/(4*c*c*r**3)
assert S.simplify(delta-(A-B))==0
print('ell =',S.factor(ell))
print('delta = A - B; A =',A,'; B =',B)
print('delta =',delta)
N=S.factor(delta*4*c*c*p*r**3)
v=S.symbols('v',nonnegative=True)
lower=S.expand((2*r*r+N).subs(t*t,1+v))
upper=S.expand((2*r*r-N).subs(t*t,1+v))
print('2 r²+N at t²=1+v =',lower)
print('2 r²-N at t²=1+v =',upper)
assert S.expand(lower-(6*s*v*v+(3-2*s-s*s)*v+3*(1-s)**2))==0
assert S.expand(upper-((4-6*s)*v*v+(5-6*s+s*s)*v+(1-s)**2))==0
assert Q(25,36)<Q(7,10)
assert 1-(Q(1801,1800)**2)/4 > Q(37,50)
J=Q(3,10)+Q(9,196)
assert J<Q(7,20)
x=Q(7,20)
exp_upper=sum(x**k/Q(math.factorial(k)) for k in range(5))+x**5/Q(math.factorial(5))/(1-x/6)
assert exp_upper<Q(10,7)
asin_lower=Q(3,4)+Q(3,4)**3/6+3*Q(3,4)**5/40+5*Q(3,4)**7/112
assert asin_lower>Q(21,25)
assert Q(405,931*900)<Q(1,2000)
mass_constant=2/(2+Q(1,2000))*(Q(7,10)*Q(21,25)-Q(30,931))*Q(1369,3000)
assert mass_constant>Q(1,4)
print('EXACT J upper =',J,'; float =',float(J))
print('EXACT exp(7/20) series upper =',exp_upper,'<',Q(10,7))
print('EXACT arcsin(3/4) series lower =',asin_lower,'>',Q(21,25))
print('EXACT exterior mass coefficient lower =',mass_constant,'; float =',float(mass_constant),' > 1/4')
print('EXACT allowed K =',1/mass_constant,'; float =',float(1/mass_constant),' < 4')
print('Arb exp(7/20) =',(arb(7)/20).exp())
print('Arb arcsin(3/4) =',(arb(3)/4).asin())
assert (arb(7)/20).exp()<arb(10)/7
assert (arb(3)/4).asin()>arb(21)/25
print('All exact constant comparisons and Arb directed inequalities PASS')
