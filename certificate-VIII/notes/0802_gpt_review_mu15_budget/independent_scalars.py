"""Independent Arb audit; no candidate imports or large matrices."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
from flint import arb,acb,ctx
import json,time,resource
P=Path(__file__).resolve().parent;ROOT=P.parents[1]
ctx.prec=768;ctx.threads=1;start=time.monotonic()
def aq(x):return arb(x.numerator)/x.denominator if isinstance(x,Q) else arb(x)
def ss(x):return x.str(32)
def l2(x):return x.log()/arb(2).log()
def dyadic(v):
 k=0
 while v<arb(2)**(-k-1):k+=1
 assert v<arb(2)**-k and v>=arb(2)**(-k-1)
 return k
L=arb(15).log()/2;Omega=736;c=Omega*L
m=arb(1)/2;s=arb(33)/8;gamma=s+m;b=2*s+m;pi=arb.pi()
assert L>arb(5)/4 and L<arb(34)/25
prime={2:2,3:3,4:2,5:5,7:7,8:2,9:3,11:11,13:13}
# Derive thresholds in physical coordinates, independently of the candidate q-events.
groups={}
for n in prime:
 event=max(Q(15,n*n),Q(n*n,15));x=abs(L-arb(n).log())
 if event in groups:
  assert x.overlaps(groups[event])
  groups[event]=groups[event].union(x)
 else:groups[event]=x
switches=sorted(groups)
boundaries=[arb(0)]+[groups[q] for q in switches]+[L]
q_events=[Q(1)]+switches+[Q(15)]
assert q_events==sorted(set(q_events)) and len(q_events)==10
rows=[];tops=[];boxes=0
for k,(left,right) in enumerate(zip(boundaries,boundaries[1:])):
 assert right>left
 qmid=(q_events[k]+q_events[k+1])/2
 active=[(sign,arb(n).log(),arb(p).log()/arb(n).sqrt(),n)
  for n,p in prime.items() for sign in [-1,1]
  if (sign==1 and qmid<Q(15,n*n)) or (sign==-1 and qmid>Q(n*n,15))]
 top=arb(0)
 for j in range(256):
  x=(left+(right-left)*j/256).union(left+(right-left)*(j+1)/256)
  w=lambda z:1+arb(2)/3*(z/L)**2
  R=sum((coef*w(x+sign*shift) for sign,shift,coef,n in active),arb(0))/w(x)
  assert R<s and R<arb('4.0915')
  top=max(top,R.upper());boxes+=1
 rows.append({'q_left':str(q_events[k]),'q_right':str(q_events[k+1]),'upper':ss(top),
  'active':[(n,sign) for sign,shift,coef,n in active]})
 tops.append(top)
print('SCHUR',boxes,'maximum',ss(max(tops)),'s',ss(s),'ratio',ss(max(tops)/s),flush=True)
C=(arb(Omega)/(2*pi)).log()-1/arb(Omega)-gamma
C_previous=(arb(704)/(2*pi)).log()-1/arb(704)-gamma
E=[(c+1)/(pi*c*c),(3*c*c+6*c+2)/(pi*c**3)]
elo=[(c-1)/(pi*c*c),3*(c-2)/(pi*c*c)]
u=[428,142];d=(-2*(c+1)).exp()/(4*(c+1)*(c+4))
delta=arb(2)**-2901;tau=arb(2)**-2904
assert C>arb(1)/8 and C_previous<arb(1)/8 and c>2
assert d>delta and d<arb(2)**-2900 and C*delta>tau
for j in [0,1]:
 assert u[j]*E[j]<C and C<(u[j]+1)*E[j] and elo[j]>delta
 print('TAIL',j,'C/E',ss(C/E[j]),'uE/C',ss(u[j]*E[j]/C),'e_lower/delta',ss(elo[j]/delta),flush=True)
H=max(abs((arb(1)/4).digamma()-pi.log()-gamma).upper(),
      abs(acb(arb(1)/4,arb(Omega)/2).digamma().real-pi.log()-gamma).upper())
K0=arb(Omega)/pi*H;Hp=(2*L).sqrt()*(L/2).cosh()
def projection(N,rho):
 R2=8*L*rho**2/(rho-1)**2*(c*(rho-1/rho)).exp()*rho**(-2*N)
 q=(2*L).sqrt()*(L/2).exp()*(L/2)**N/arb(factorial(N))
 ee=(K0*(2*L).sqrt()+(arb(Omega)/pi).sqrt()*sum(u[j]*E[j].sqrt() for j in [0,1]))*R2.sqrt()+4*Hp*q
 nn=K0*R2+4*q*q
 hh=(K0+sum(u)*arb(Omega)/pi)*R2+4*q*q
 assert nn<m
 aa=ee+2*ee*ee/(m-nn);bb=(ee+nn)/(m*m)+2*hh*hh/(m*m*(m-nn))
 return ee,nn,hh,aa,bb,aa+bb*b*b
prows=[];candidate=json.loads((ROOT/'notes/0792_claude_lowmem_generator/scalar_mu15.json').read_text())
for N in [1408,1536,1664,1792]:
 for rho in [Q(2),Q(5,2),Q(3),Q(7,2)]:
  ee,nn,hh,aa,bb,pay=projection(N,aq(rho))
  row={'N':N,'rho':str(rho),'e':ss(ee),'n':ss(nn),'h':ss(hh),'payment':ss(pay),'log2_payment':ss(l2(pay))}
  old=next(r for r in candidate['mu15_budget'] if r['N']==N)
  if old['rho']==str(rho):
   ref=arb(old['log2_payment']);assert l2(pay).overlaps(ref)
   row['candidate_log2']=ss(ref);row['log2_difference']=ss(l2(pay)-ref)
  prows.append(row)
  print('PROJECTION',N,str(rho),ss(l2(pay)),flush=True)
ee,nn,hh,aa,bb,pay=projection(1664,arb(3))
kE,kN,kH=[dyadic(v) for v in [ee,nn,hh]]
eF,nF,hF=[Q(1,2**k) for k in [kE,kN,kH]]
mF=Q(1,2);bF=Q(35,4)
alpha=eF+2*eF**2/(mF-nF);beta=(eF+nF)/mF**2+2*hF**2/(mF**2*(mF-nF))
pf=alpha+beta*bF*bF
assert pay>arb(2)**-698 and pay<arb(2)**-697
assert pf<Q(1,2**697)
print('ADOPTED',{'e_exp':kE,'n_exp':kN,'h_exp':kH},'exact payment/2^-697',ss(aq(pf*2**697)),flush=True)
ends=[Q(0),Q(1,2),Q(1)]+[Q(2**k) for k in range(1,7)]+[Q(x) for x in range(128,736,64)]+[Q(736)]
assert len(ends)==20 and all(y>x for x,y in zip(ends,ends[1:]))
panels=[]
for left,right in zip(ends,ends[1:]):
 h=(right-left)/2;ctr=(right+left)/2;remin=left-h/4;immax=3*h/4
 assert immax<=24
 # Uniform pole-distance witness for both signs in w=1/4 +- iz/2.
 bound=Q(1,4)-immax/2 if right-left<=1 else remin/2
 assert bound>=Q(1,16)
 panels.append({'left':str(left),'right':str(right),'half_length':str(h),
  'real_min':str(remin),'imaginary_max':str(immax),'pole_distance_lower':str(bound),'abs_z_upper':str(ctr+5*h/4)})
assert max(Q(x['abs_z_upper']) for x in panels)<768
assert 1+2*397+13*16==1003
gb=arb(1003)+pi.log()+gamma;assert gb<1009 and gb<1100
assert 2*(2*L).sqrt()<4 and Hp<3
HB=2*L*(48*L).exp();qrows=[]
for N in [1536,1664,1792]:
 for q in [192,224,256]:
  dim=N//2;eb=dim*4*Omega*HB*arb(2)**(-(2*q-1))
  dp=2*(L/2)**(N+81)*(L/2).exp()/arb(factorial(N+81))
  eq=1100*eb+428*(2*eb+eb*eb)+96*dp+64*dp*dp
  if q==224:assert eq<arb(2)**-319
  if q==256:assert eq<arb(2)**-383
  qrows.append({'N':N,'q':q,'epsilon_band':ss(eb),'d_p':ss(dp),'epsilon_J_quad':ss(eq),'log2':ss(l2(eq))})
  print('QUADRATURE',N,q,ss(l2(eq)),flush=True)
delta_change=570*(2*d+d*d);assert delta_change<arb(2)**-2889
precise_change=570*(2*(d-delta)+(d-delta)**2)
print('DELTA',ss(d),'bits',ss(-l2(d)),'change_bound',ss(delta_change),'ratio_to_2^-2889',ss(delta_change/(arb(2)**-2889)),flush=True)
# Fixed matrix-error proposal: these are requirements for the future certificate.
eA=eB=Q(1,2**400);eW=Q(1,2**350);eD=Q(1,2**380);eJ=Q(1,2**300);eY=Q(1,2**310)
xn=Q(512);an=Q(9);wround=Q(1,2**480)
wreq=(xn*xn+4*xn)*eA+2*xn*xn*eB+wround
dreq=eB+2*an*eA+eA*eA
assert wreq<eW and dreq<eD
et=2*bF*bF*eW+eY+eJ+beta*eD
assert et<Q(1,2**290)
eq224=arb(next(r['epsilon_J_quad'] for r in qrows if r['N']==1664 and r['q']==224))
assert eq224+delta_change<aq(eJ)/2
scheme={'epsilon_A':'2^-400','epsilon_B':'2^-400','epsilon_W':'2^-350','epsilon_D':'2^-380','epsilon_J':'2^-300','epsilon_Y_required':'2^-310','epsilon_T_cap':'2^-290',
 'W_required_under_Xcap':ss(aq(wreq)),'W_required_ratio':ss(aq(wreq/eW)),
 'D_required_under_Acap':ss(aq(dreq)),'D_required_ratio':ss(aq(dreq/eD)),
 'T_worst':ss(aq(et)),'T_ratio_to_cap':ss(aq(et/Q(1,2**290))),
 'J_analytic_ratio_to_fixed_J':ss((eq224+delta_change)/aq(eJ)),
 'X_norm_cap':'512','A0_norm_cap':'9','W_rounding_cap':'2^-480'}
surrogates=[]
for p in [0,1]:
 path=ROOT/f'notes/0792_claude_lowmem_generator/runs/m15n1664p3328q224_p{p}_diagnostic.json'
 rec=json.loads(path.read_text())
 assert rec['block_errors_paid'] is False and rec['integral_remainders_paid'] is False
 raw=arb(rec['surrogate_lower_bound']);z2=arb(rec['frobenius_squared'])
 surrogates.append({'parity':p,'diagnostic_lower':ss(raw),'minus_log2':ss(-l2(raw)),
  'error_cap_times_diagnostic_z2':ss(arb(2)**-290*z2),'worst_T_error_to_diagnostic':ss(aq(et)/raw),
  'status':'diagnostic only; final matrix gates pending'})
data={'status':'PASS_ANALYTIC_SCALARS','precision_bits':ctx.prec,'L':ss(L),'c':ss(c),'C':ss(C),'previous_C704':ss(C_previous),
 'u':u,'C_over_E':[ss(C/z) for z in E],'uE_over_C':[ss(u[j]*E[j]/C) for j in [0,1]],'schur_boxes':boxes,'schur_q_events':[str(x) for x in q_events],'schur_cells':rows,'schur_max':ss(max(tops)),
 'd':ss(d),'minus_log2_d':ss(-l2(d)),'d_over_delta':ss(d/delta),'projection_rows':prows,
 'adopted':{'N':1664,'rho':'3','e_exponent':kE,'n_exponent':kN,'h_exponent':kH,'alpha':str(alpha),'beta':str(beta),'payment_ratio_to_2^-697':ss(aq(pf*2**697)),'actual_payment_over_2^-698':ss(pay/(arb(2)**-698))},
 'panels':panels,'HB':ss(HB),'g_bound':ss(gb),'Hp':ss(Hp),'quad_rows':qrows,'delta_change_bound':ss(delta_change),'delta_change_precise':ss(precise_change),'delta_change_power':'2^-2889','scheme':scheme,'diagnostic_forecast':surrogates,
 'seconds':time.monotonic()-start,'maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(P/'independent_scalars.json').write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
print('SCHEME',json.dumps(scheme),flush=True)
print('FORECAST',json.dumps(surrogates),flush=True)
print('COMPLETE',boxes,'Schur boxes',len(prows),'projection choices',len(qrows),'quadrature cases',flush=True)
