"""Scalar handoff only, not a finite matrix certificate."""
from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
import json
ctx.prec=512;o=Path(__file__).parent
s=json.loads((o/'budget.json').read_text())
r=next(r for r in s['rows'] if r['N']==832 and arb(r['rho'])==3)
exponents={'e':-311,'n':-634,'h':-629}
v={k:F(1,2**(-p)) for k,p in exponents.items()}
for k,p in exponents.items():assert arb(r[k])<arb(2)**p
m=F(1,2);b=F(61,8);e,n,h=(v[k] for k in ('e','n','h'));assert n<m
a=e+2*e*e/(m-n);beta=(e+n)/(m*m)+2*h*h/(m*m*(m-n));pay=a+beta*b*b
assert pay<F(1,2**303)
q=next(r for r in s['quadrature'] if r['N']==832 and r['q']==256)
assert arb(q['epsilon_J_quadrature'])<arb(2)**-263
res={'scope':'RATIONAL_SCALAR_BUDGET; source bound derivation requires review; finite certificate absent',
'N':832,'rho':3,'scalar_upper_bounds':{k:str(x) for k,x in v.items()},
'alpha':str(a),'beta':str(beta),'payment':str(pay),'payment_lt_2_minus303':True,
'J_quadrature_lt_2_minus263':True,'candidate_source_J_budget':'2^-250',
'caveat':'Rational checks do not establish actual matrix C/W positivity, X/Y enclosures, or the final paid congruence.'}
(o/'rational_budget.json').write_text(json.dumps(res,indent=2));print('PASS: exact rational payment < 2^-303; Arb analytic quadrature bound < 2^-263')
