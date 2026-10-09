"""0784 scalar budgets. Analytic bounds, not a full finite certificate."""
from pathlib import Path
import json,math
from flint import arb,acb,ctx
ctx.prec=512;out=Path(__file__).parent;L=arb(12).log()/2;c=416*L;gam=arb(65)/16;m=arb(1)/2;b=arb(61)/8
assert L<arb(5)/4
C=(arb(416)/(2*arb.pi())).log()-1/arb(416)-gam
E=[(c+1)/(arb.pi()*c*c),(3*c*c+6*c+2)/(arb.pi()*c*c*c)]
assert 207*E[0]<C and 68*E[1]<C
H=max(abs((arb(1)/4).digamma()-arb.pi().log()-gam).upper(),abs(acb(arb(1)/4,208).digamma().real-arb.pi().log()-gam).upper())
K0=416/arb.pi()*H;Hp=(2*L).sqrt()*(L/2).cosh();rows=[]
for N,rho in ((640,arb(3)/2),(640,arb(2)),(832,arb(3)/2),(832,arb(3)),(1008,arb(3)/2),(1008,arb(7)/2)):
 r2=8*L*rho*rho/((rho-1)*(rho-1))*(c*(rho-1/rho)).exp()*rho**(-2*N);r=r2.sqrt()
 qp=(2*L).sqrt()*(L/2).exp()*(L/2)**N/arb(math.factorial(N))
 e=(K0*(2*L).sqrt()+(207*E[0].sqrt()+68*E[1].sqrt())*(416/arb.pi()).sqrt())*r+4*Hp*qp
 n=K0*r2+4*qp*qp;h=(K0+275*416/arb.pi())*r2+4*qp*qp
 assert n<m
 a=e+2*e*e/(m-n);beta=(e+n)/(m*m)+2*h*h/(m*m*(m-n));payment=a+beta*b*b
 rows.append({'N':N,'rho':str(rho),'r':str(r),'log2_r':str(r.log()/arb(2).log()),'e':str(e),'n':str(n),'h':str(h),'alpha':str(a),'beta':str(beta),'payment':str(payment),'log2_payment':str(payment.log()/arb(2).log()),'below_2_minus200':bool(payment<arb(2)**-200)})
quad=[]
for N in (640,832,1008):
 dp=2*(arb(5)/8)**(N+81)/arb(math.factorial(N+81));epole=96*dp+64*dp*dp
 for q in (224,256):
  eb=(N//2)*4*416*(arb(5)/2)*arb(150).exp()*arb(2)**(-(2*q-1))
  ej=1100*eb+207*(2*eb+eb*eb)+epole
  quad.append({'N':N,'q':q,'epsilon_band':str(eb),'epsilon_J_quadrature':str(ej),'log2_epsilon_J':str(ej.log()/arb(2).log()),'below_2_minus200':bool(ej<arb(2)**-200)})
res={'status':'ANALYTIC_ERROR_BUDGET_ONLY','rows':rows,'quadrature':quad}
(out/'budget.json').write_text(json.dumps(res,indent=2))
for row in rows:print('block',row['N'],float(arb(row['rho'])),float(arb(row['log2_payment'])),row['below_2_minus200'])
for row in quad:print('quad',row['N'],row['q'],float(arb(row['log2_epsilon_J'])),row['below_2_minus200'])
