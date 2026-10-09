"""Reviewer scalar gates. No imports from any candidate or prior project module."""
from pathlib import Path
from fractions import Fraction
from math import factorial
import json,subprocess
from flint import arb,acb,ctx
P=Path(__file__).resolve().parent;ctx.prec=320;ctx.threads=1
L=arb(12).log()/2;O=arb(416);c=O*L;pi=arb.pi();gamma=arb(65)/16
C=(O/(2*pi)).log()-1/O-gamma
d=(-2*(c+1)).exp()/(4*(c+1)*(c+4));delta=arb(2)**-1515
pp={2:2,3:3,4:2,5:5,7:7,8:2,9:3,11:11}
# Direct x-space activation events, rather than the candidate's rational q-events.
events=sorted([arb(0),L]+[abs(L-arb(n).log()) for n in pp])
weight=lambda x:1+arb(19)/40*(x/L)**2
boxes=[];caps=[]
for left,right in zip(events,events[1:]):
 mid=(left+right)/2
 shifts=[]
 for n,p in pp.items():
  h=arb(n).log();coef=arb(p).log()/arb(n).sqrt()
  if mid+h<L:shifts.append((h,coef))
  if mid-h>-L:shifts.append((-h,coef))
 # 1024 boxes in physical x. Closed boxes enclose uncertain algebraic endpoints.
 uppers=[]
 for k in range(1024):
  x=(left+(right-left)*k/1024).union(left+(right-left)*(k+1)/1024)
  total=sum((a*weight(x+h) for h,a in shifts),arb(0))/weight(x)
  assert total<arb(57)/16
  uppers.append(total.upper())
 cap=max(uppers);caps.append(cap)
 boxes.append({'left':str(left),'right':str(right),'active_shifts':[str(h) for h,_ in shifts],'row_sum_upper':str(cap)})
H=max(abs((arb(1)/4).digamma()-pi.log()-gamma).upper(),
      abs(acb(arb(1)/4,O/2).digamma().real-pi.log()-gamma).upper())
K0=O*H/pi
E0=(c+1)/(pi*c**2);E1=(3*c*c+6*c+2)/(pi*c**3)
e0lo=(c-1)/(pi*c*c);e1lo=3*(c-2)/(pi*c*c)
rho=arb(3);N=832
# First compute the uniform Chebyshev error then multiply by sqrt(2L).
r=(2*L).sqrt()*2*(c*(rho-1/rho)/2).exp()*rho**(-N)/(1-1/rho)
q=(2*L).sqrt()*(L/2).exp()*(L/2)**N/arb(factorial(N))
Hp=(2*L).sqrt()*(L/2).cosh()
e=(K0*(2*L).sqrt()+(O/pi).sqrt()*(207*E0.sqrt()+68*E1.sqrt()))*r+4*Hp*q
n=K0*r*r+4*q*q
h=(K0+275*O/pi)*r*r+4*q*q
band=4*O*416*(arb(5)/2)*arb(150).exp()*arb(2)**-511
pole_remainder=2*(arb(5)/8)**913/arb(factorial(913))
Kerr=1100*band;Uerr=207*(2*band+band*band)
perr=96*pole_remainder+64*pole_remainder*pole_remainder
jerr=Kerr+Uerr+perr
pert=275*(2*d+d*d)
F=Fraction;ee=F(1,2**311);nn=F(1,2**634);hh=F(1,2**629);m=F(1,2)
alpha=ee+2*ee*ee/(m-nn);beta=(ee+nn)/(m*m)+2*hh*hh/(m*m*(m-nn))
pay=alpha+beta*F(61,8)**2
assert pay<F(1,2**303)
comparisons={
 'prime_norm':(max(caps),arb(57)/16),'C_lower':(C,arb(1)/8),
 'd_vs_delta':(d,delta),'d_vs_upper':(d,arb(2)**-1514),
 '207E0_vs_C':(207*E0,C),'68E1_vs_C':(68*E1,C),
 'e0_lower_vs_delta':(e0lo,delta),'e1_lower_vs_delta':(e1lo,delta),
 'e_vs_budget':(e,arb(2)**-311),'n_vs_budget':(n,arb(2)**-634),'h_vs_budget':(h,arb(2)**-629),
 'Jquad_vs_budget':(jerr,arb(2)**-263),'delta_perturbation_vs_budget':(pert,arb(2)**-1504)}
assert C>arb(1)/8 and delta<d<arb(2)**-1514
assert 207*E0<C and 68*E1<C and e0lo>delta and e1lo>delta
assert e<arb(2)**-311 and n<arb(2)**-634 and h<arb(2)**-629
assert jerr<arb(2)**-263 and pert<arb(2)**-1504
data={'status':'PASS','precision_bits':ctx.prec,'x_boxes_per_cell':1024,'prime_cells':boxes,
 'values':{k:str(v) for k,v in {'L':L,'c':c,'H0':H,'K0':K0,'r':r,'qN':q,'Hp':Hp,
 'epsilon_band':band,'epsilon_K':Kerr,'epsilon_U':Uerr,'epsilon_pole':perr,
 'minus_log2_d':-d.log()/arb(2).log(),'minus_log2_Jquad':-jerr.log()/arb(2).log(),
 'tau':delta/8}.items()},
 'comparisons':{k:{'value':str(a),'comparison':str(b),'ratio':str(a/b)} for k,(a,b) in comparisons.items()},
 'exact_budget':{'alpha':str(alpha),'beta':str(beta),'block_payment_over_2^-303':str(pay/F(1,2**303))}}
(P/'scalar_review.json').write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
for k,v in data['comparisons'].items():print(k,json.dumps(v),flush=True)
for k,v in data['values'].items():print(k,v)
print('block_payment_over_2^-303',float(pay/F(1,2**303)))
print('COMPLETE independent scalar gates')
