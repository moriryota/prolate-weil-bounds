"""Replay exact rational gates without rewriting the saved result."""
import json
from fractions import Fraction as F
from pathlib import Path
r=F(1,2**190);p=F(1,2**440);m=F(1,2)
e=1500*r+p;n=896*r*r+p;h=12000*r*r+p
a=e+2*e*e/(m-n);b=(e+n)/(m*m)+2*h*h/(m*m*(m-n))
assert n<m and a+b*F(13,2)**2<F(1,2**170)
saved=json.loads(Path(__file__).with_name('exact_parameters.json').read_text())
assert saved['n_lt_m'] and saved['scalar_budget_lt_2_pow_minus170']
assert saved['decimal_reference_only']=={'n':float(n),'alpha':float(a),'beta':float(b),'budget':float(a+b*F(13,2)**2)}
print('Exact rational gates and saved decimal reference: PASS')
