"""Fixed mu15 scalar specification, incorporating 0802 conditional acceptance. Exact budgets plus independently rerun Arb gates.
No import of the proposal or of prior executable modules.
"""
from fractions import Fraction as F
from math import factorial
from flint import arb,acb

def ab(x):return arb(x.numerator)/x.denominator if isinstance(x,F) else arb(x)
def fraction_text(x):return str(x.numerator)+'/'+str(x.denominator)
def analytic_bounds():
 m=F(1,2);e=F(1,2**706);n=F(1,2**1424);h=F(1,2**1419)
 alpha=e+2*e*e/(m-n);beta=(e+n)/(m*m)+2*h*h/(m*m*(m-n))
 assert alpha+beta*F(35,4)**2<F(1,2**697)
 return dict(e=e,n=n,h=h,alpha=alpha,beta=beta,J_quad=F(1,2**319)+F(1,2**2889),frequency_pole_U=F(1,2**319),delta_change=F(1,2**2889))

def validate_scalars():
 L=arb(15).log()/2;c=736*L;pi=arb.pi();gam=arb(37)/8;delta=arb(2)**-2901;d=(-2*(c+1)).exp()/(4*(c+1)*(c+4));C=(736/(2*pi)).log()-arb(1)/736-gam
 assert delta<=d and d<arb(2)**-2900 and C>arb(1)/8 and arb(2)**-2904<=C*delta
 pp={2:2,3:3,4:2,5:5,7:7,8:2,9:3,11:11,13:13}
 events=sorted({F(1),F(15),*[max(F(15,n*n),F(n*n,15)) for n in pp]})
 w=lambda x:1+arb(2)/3*(x/L)*(x/L)
 cells=[];caps=[]
 for a,b in zip(events,events[1:]):
  mid=(a+b)/2;sh=[(n,1) for n in pp if mid<F(15,n*n)]+[(n,-1) for n in pp if mid>F(n*n,15)]
  lo=ab(a).log()/2;hi=ab(b).log()/2;values=[]
  for j in range(256):
   x=(lo+(hi-lo)*j/256).union(lo+(hi-lo)*(j+1)/256)
   v=sum(arb(pp[n]).log()/arb(n).sqrt()*w(x+s*arb(n).log()) for n,s in sh)/w(x)
   assert v<arb(66)/16
   values.append(v.upper())
  cap=max(values);caps.append(cap);cells.append({'a':str(a),'b':str(b),'upper':str(cap)})
 E0=(c+1)/(pi*c*c);E1=(3*c*c+6*c+2)/(pi*c*c*c)
 assert 428*E0<C and 142*E1<C and (c-1)/(pi*c*c)>delta and 3*(c-2)/(pi*c*c)>delta
 H=max(abs((arb(1)/4).digamma()-pi.log()-gam).upper(),abs(acb(arb(1)/4,368).digamma().real-pi.log()-gam).upper())
 K0=736/pi*H;Hp=(2*L).sqrt()*(L/2).cosh();rho=arb(3);N=1664
 r2=8*L*rho*rho/(rho-1)**2*(c*(rho-1/rho)).exp()*rho**(-2*N);r=r2.sqrt()
 qN=(2*L).sqrt()*(L/2).exp()*(L/2)**N/arb(factorial(N))
 e=(K0*(2*L).sqrt()+(736/pi).sqrt()*(428*E0.sqrt()+142*E1.sqrt()))*r+4*Hp*qN
 n=K0*r2+4*qN*qN;h=(K0+570*736/pi)*r2+4*qN*qN
 for v,k in ((e,706),(n,1424),(h,1419)):assert v<arb(2)**-k
 assert L<arb(34)/25 and (2*L).sqrt()<2 and (2*L).sqrt()*(L/2).cosh()<3
 eb=832*4*736*(2*L)*(48*L).exp()*arb(2)**-447
 dp=2*(L/2)**1745*(L/2).exp()/arb(factorial(1745))
 eq=1100*eb+428*(2*eb+eb*eb)+96*dp+64*dp*dp
 assert eq<arb(2)**-319
 perturb=570*(2*d+d*d);assert perturb<arb(2)**-2889
 analytic_bounds()
 return {'prime_cells':cells,'prime_norm_upper':str(max(caps)),'C':str(C),'d':str(d),'delta':'2^-2901','tau':'2^-2904','e':str(e),'n':str(n),'h':str(h),'J_quadrature':str(eq),'delta_perturbation':str(perturb),'status':'PASS_SCALAR_GATES'}
