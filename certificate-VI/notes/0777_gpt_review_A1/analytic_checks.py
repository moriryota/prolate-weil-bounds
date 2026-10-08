"""Reviewer-owned rational checks; no candidate code imported."""
from fractions import Fraction as F
from math import factorial
import json
from pathlib import Path

def exp_bounds(x, n=40):
    s=sum((x**k/F(factorial(k)) for k in range(n+1)),F(0))
    next_term=x**(n+1)/factorial(n+1)
    return s,s+next_term/(1-x/F(n+2))

assert exp_bounds(F(1))[1]<F(68,25)<3
assert exp_bounds(F(17,8))[1]<10<exp_bounds(F(7,3))[0]
assert 3*49152<64*4095
assert F(3*768*7,6*4096)<F(2,3)
assert F(4,9)<F(1,2)
pi_lo=16*(F(1,5)-F(1,375))-F(4,239)
pi_hi=16*(F(1,5)-F(1,375)+F(1,15625))-4*(F(1,239)-F(1,3*239**3))
assert pi_lo>F(157,50) and pi_hi<F(22,7)
assert F(448,11)**5>F(68,25)**18
C=F(18,5)-F(1,256)-F(7,2)
assert C==F(123,1280) and C>F(1,16)
c=F(272)
u0=50*(c+1)/(157*c*c)
u1=50*(3*c*c+6*c+2)/(157*c**3)
assert 81*u0<C and 27*u1<C
c=F(896,3)
l0=7*(c-1)/(22*c*c)
l1=21*(c-2)/(22*c*c)
assert min(l0,l1)>F(1,2048)
assert F(1,15)+F(7,2*256)<1
assert F(97,4)+136<162
assert 325+25*16==725 and 725+2+F(7,2)<1000
eb=F(320*1024*3*10**48,2**383)
ek=F(320*1024*4000*10**48,2**383)
dp=2*F(3,5)**721/factorial(721)
ej=ek+81*(2*eb+eb*eb)+160*dp
assert 96*dp+64*dp*dp<160*dp
assert 639+720==2*680-1 and 1278<=2*642-1
assert ej<F(1,2**190)
r,p,m=F(1,2**190),F(1,2**440),F(1,2)
e=1500*r+p; n=896*r*r+p; h=12000*r*r+p
alpha=e+2*e*e/(m-n)
beta=(e+n)/(m*m)+2*h*h/(m*m*(m-n))
assert n<m and alpha+beta*F(13,2)**2<F(1,2**170)
out={'status':'PASS','C':str(C),'81_e0_upper':str(81*u0),'27_e1_upper':str(27*u1),
     'e0_lower':str(l0),'e1_lower':str(l1),'J_quad_decimal':float(ej),
     'J_quad_over_2^-190':float(ej*2**190),'alpha':float(alpha),'beta':float(beta)}
Path(__file__).with_name('analytic_results.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
