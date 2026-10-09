"""Trusted scalar specification. Exact budgets plus independently rerun Arb gates.
No import of the proposal or of prior executable modules.
"""
from fractions import Fraction as F
from math import factorial
from flint import arb,acb

def ab(x):return arb(x.numerator)/x.denominator if isinstance(x,F) else arb(x)
def fraction_text(x):return str(x.numerator)+'/'+str(x.denominator)
def analytic_bounds():
 m=F(1,2);e=F(1,2**311);n=F(1,2**634);h=F(1,2**629)
 alpha=e+2*e*e/(m-n);beta=(e+n)/(m*m)+2*h*h/(m*m*(m-n))
 assert alpha+beta*F(61,8)**2<F(1,2**303)
 return dict(e=e,n=n,h=h,alpha=alpha,beta=beta,J_quad=F(1,2**263)+F(1,2**1504),frequency_pole_U=F(1,2**263),delta_change=F(1,2**1504))

def validate_scalars():
 L=arb(12).log()/2;c=416*L;pi=arb.pi();gam=arb(65)/16;delta=arb(2)**-1515;d=(-2*(c+1)).exp()/(4*(c+1)*(c+4));C=(416/(2*pi)).log()-arb(1)/416-gam
 assert delta<=d and d<arb(2)**-1514 and C>arb(1)/8 and arb(2)**-1518<=C*delta
 pp={2:2,3:3,4:2,5:5,7:7,8:2,9:3,11:11}
 events=sorted({F(1),F(12),*[max(F(12,n*n),F(n*n,12)) for n in pp]})
 w=lambda x:1+arb(19)/40*(x/L)*(x/L)
 cells=[];caps=[]
 for a,b in zip(events,events[1:]):
  mid=(a+b)/2;sh=[(n,1) for n in pp if mid<F(12,n*n)]+[(n,-1) for n in pp if mid>F(n*n,12)]
  lo=ab(a).log()/2;hi=ab(b).log()/2;values=[]
  for j in range(256):
   x=(lo+(hi-lo)*j/256).union(lo+(hi-lo)*(j+1)/256)
   v=sum(arb(pp[n]).log()/arb(n).sqrt()*w(x+s*arb(n).log()) for n,s in sh)/w(x)
   assert v<arb(57)/16
   values.append(v.upper())
  cap=max(values);caps.append(cap);cells.append({'a':str(a),'b':str(b),'upper':str(cap)})
 E0=(c+1)/(pi*c*c);E1=(3*c*c+6*c+2)/(pi*c*c*c)
 assert 207*E0<C and 68*E1<C and (c-1)/(pi*c*c)>delta and 3*(c-2)/(pi*c*c)>delta
 H=max(abs((arb(1)/4).digamma()-pi.log()-gam).upper(),abs(acb(arb(1)/4,208).digamma().real-pi.log()-gam).upper())
 K0=416/pi*H;Hp=(2*L).sqrt()*(L/2).cosh();rho=arb(3);N=832
 r2=8*L*rho*rho/(rho-1)**2*(c*(rho-1/rho)).exp()*rho**(-2*N);r=r2.sqrt()
 qN=(2*L).sqrt()*(L/2).exp()*(L/2)**N/arb(factorial(N))
 e=(K0*(2*L).sqrt()+(416/pi).sqrt()*(207*E0.sqrt()+68*E1.sqrt()))*r+4*Hp*qN
 n=K0*r2+4*qN*qN;h=(K0+275*416/pi)*r2+4*qN*qN
 for v,k in ((e,311),(n,634),(h,629)):assert v<arb(2)**-k
 assert L<arb(5)/4 and (arb(5)/8).exp()<2 and (2*L).sqrt()*(L/2).cosh()<3
 eb=416*4*416*(arb(5)/2)*arb(150).exp()*arb(2)**-511
 dp=2*(arb(5)/8)**913/arb(factorial(913))
 eq=1100*eb+207*(2*eb+eb*eb)+96*dp+64*dp*dp
 assert eq<arb(2)**-263
 perturb=275*(2*d+d*d);assert perturb<arb(2)**-1504
 analytic_bounds()
 return {'prime_cells':cells,'prime_norm_upper':str(max(caps)),'C':str(C),'d':str(d),'delta':'2^-1515','tau':'2^-1518','e':str(e),'n':str(n),'h':str(h),'J_quadrature':str(eq),'delta_perturbation':str(perturb),'status':'PASS_SCALAR_GATES'}
