"""Small independent exact/Arb checks for the final alpha,beta,delta,tau connection."""
from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
import json
ctx.prec=768;ctx.threads=1
OUT=Path(__file__).resolve().parent
e,n,h=F(1,2**706),F(1,2**1424),F(1,2**1419)
m=F(1,2);b=F(35,4)
alpha=e+2*e*e/(m-n)
beta=(e+n)/(m*m)+2*h*h/(m*m*(m-n))
payment=alpha+beta*b*b
assert n<m and payment<F(1,2**697)
lift=lambda q:arb(q.numerator)/q.denominator
L=arb(15).log()/2;c=736*L
C=(736/(2*arb.pi())).log()-arb(1)/736-arb(37)/8
d=(-2*(c+1)).exp()/(4*(c+1)*(c+4))
delta=arb(2)**-2901;tau=arb(2)**-2904
assert C>arb(1)/8 and d>delta and C*delta>tau
result={'status':'PASS_FINAL_CONNECTION_SCALARS','reviewer':'sol',
 'precision_bits':ctx.prec,'alpha_exact':str(alpha),'beta_exact':str(beta),
 'payment_ratio_to_2^-697':lift(payment*2**697).str(40),
 'C':C.str(40),'d_over_delta':(d/delta).str(40),
 'Cdelta_over_tau':(C*delta/tau).str(40),
 'delta':'2^-2901','tau':'2^-2904'}
(OUT/'connection_scalars.json').write_text(json.dumps(result,indent=2)+'\n')
for key in ('payment_ratio_to_2^-697','C','d_over_delta','Cdelta_over_tau'):print(key,result[key])
print('COMPLETE: exact alpha,beta payment and strict Cdelta>tau')
