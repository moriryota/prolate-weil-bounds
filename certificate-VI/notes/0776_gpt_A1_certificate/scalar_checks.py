"""Exact checks of the stated geometry/constants and the accepted finite margin."""
import json
from pathlib import Path
from fractions import Fraction as F
from math import factorial
from flint import arb,ctx
from bounds import analytic_bounds,fraction_text
OUT=Path(__file__).parent;ctx.prec=1280
panels=[F(0),F(1,2)]+[F(2)**k for k in range(9)]
assert len(panels)==11
geometry=[]
for a,b in zip(panels,panels[1:]):
    length=b-a;center=(a+b)/2
    immax=3*length/8;absmax=center+5*length/8
    if length<=1:distance=F(1,4)-immax/2
    else:
        assert b==2*a and a>=2
        distance=(center-5*length/8)/2
    assert immax<=48 and absmax<=272 and distance>=F(1,16)
    geometry.append({'a':str(a),'b':str(b),'imaginary_bound':str(immax),'absolute_bound':str(absmax),'distance_from_digamma_poles_lower':str(distance)})
assert F(1,4)+25-F(48,2)==F(5,4)
assert F(97,4)+F(272,2)<162
assert 1+162*2==325 and 325+25*16==725
assert F(725)+2+F(7,2)<1000
# L<6/5 and exp(3/5)<2, without numerical logarithms.
assert sum(F(12,5)**k/factorial(k) for k in range(9))>10
assert 1+F(3,5)+F(3,5)**2/(2*(1-F(1,5)))<2
# exp(17/8)<10 proves log(10)/2>17/16.
x=F(17,8)
assert sum(x**k/factorial(k) for k in range(21))+x**21/factorial(21)/(1-x/22)<10
vals=analytic_bounds()
v=json.loads((OUT/'verification.json').read_text())
assert v['status']=='PASS_EXTERNAL_CERTIFICATE' and [r['parity'] for r in v['results']]==[0,1]
for r in v['results']:
    assert arb(r['T']['positive_lower_bound'])>arb(2)**-144
    assert arb(r['epsilon_T'])<arb('1e-46')
result={'status':'PASS','geometry':geometry,'constants':{k:{'exact':fraction_text(q),'display':(arb(q.numerator)/q.denominator).str(18)} for k,q in vals.items()},
        'finite_T_lower_bound_both_parities':'2^-144','L_larger_than_17_over_16':True}
path=OUT/'scalar_results.json';assert not path.exists();path.write_text(json.dumps(result,indent=2))
print('Geometry, rational constants, and finite T > 2^-144: PASS')
