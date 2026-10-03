"""Paper III: constants of Lemma 4.1 and Proposition 4.2 (C_f, A_f, g_A), Proposition 5.3 (Z_C, K_*),
Section 6 for T = 3e12, and Theorems F(a), G(a) (6.034e14, 237522697, 1.434e23). Arb, directed rounding asserted.
Usage: python proofs-III/theoremF_a_constants.py"""
from flint import arb,ctx
import json
ctx.prec=256
A=arb;pi=A.pi();ell=A(395)/399;c0=10*pi;c50=100*pi;lo=A(50).log()
beta=A(4);eb=A(8).exp()
j=32*eb;v=2+48*eb;z=j*(8/pi).sqrt()/ell
d0=A(32).sqrt()*eb
def up(x,n=6):
    scale=10**n;k=int((x.upper()*scale).ceil().unique_fmpz())
    s=str(k) if n==0 else ('-' if k<0 else '')+f'{abs(k)//scale}.{abs(k)%scale:0{n}d}'
    assert A(s)>=x.upper()
    return s
def down(x,n=6):
    scale=10**n;k=int((x.lower()*scale).floor().unique_fmpz())
    s=str(k) if n==0 else ('-' if k<0 else '')+f'{abs(k)//scale}.{abs(k)%scale:0{n}d}'
    assert A(s)<=x.lower()
    return s
def U(c):return 6621*c**(A(9)/2)*(-2*c).exp()
def V(c):return v*(16*pi/ell).sqrt()+2*z/ell.sqrt()*(4/c+2/A(3).sqrt())
def m(c):return ((2*(592+A('152.2')*(2*(c/3).log()+A('8.08'))**2)).sqrt()+V(c))**2
def EE(c):return V(c)+A('1.78')*A(8).sqrt()*(pi/2+1+2*(c/3).log()+A('5.49'))+A('2.0408').sqrt()/2
def D(c):return ((8*c).sqrt()+EE(c))**2
CD=A(337392159829)
def DB(c):return 16*c+CD*c.log()**2
def Zraw(c):return (m(c).sqrt()+D(c).sqrt()/c)**2
ZC=A(up(Zraw(c50)/lo**2,0))
Af=10+48*eb;Cf=20+321*eb
ga=3*(Af+2*Cf)/(2*pi).sqrt()
def G(c):return (ga*c*c+120*c**(A(3)/2))**2
def Z(c):return ZC*(c/(2*pi)).log()**2
def K(c):return (Z(c).sqrt()+G(c).sqrt()+c*DB(c).sqrt())**2
Cz=A('3.325')+A('.2')/A(1).exp()
def S(d):return (A(4)/3*d.log()+1)/(pi*d)+(A('.772')*d.log()+2*Cz+A('.193'))/d**2
T=A(3)*10**12;x0=(3*T).log()
n0=(1/(4*pi)+A('.222'))*x0+A('.55')*x0.log()+A('4.9')+A('.4')/A(1).exp()
At=2*(A(1000)/999)**2*n0*(A('11.15')*x0)**2/T**2*(A(572587)/414720)
tail=(T.log()+1)/(pi*T)
def Ltail(mu):return A('.5')+A(mu).sqrt()/2*(A('11.15')*x0/2).log()
M=10000;cM=2*pi*M
def high(mu):
    c=2*pi*mu
    return 2*At*Ltail(mu)*K(c)+4*A(mu).sqrt()*DB(c)*tail
low_b_max=4*DB(c50)*S(c50**2)
CPexact=48*pi*26*ZC+(low_b_max+high(M))/lo**4
CPderived=A(up(CPexact,0))
CP=A(603400000000000)
CE=A(143400000000000000000000)
assert CP>=CPderived.upper()
assert CE>=A(237522697)*CP
checks={}
checks['endpoint_coeff']= (8/ell).sqrt()*(1+z*U(c0).sqrt()/c0)
assert checks['endpoint_coeff']<3
assert A('1.77')**2>A('3.12')
assert A('1.1')**4*A('.75')*A('.925')>1
assert A('1.16')**2>A(4)/3
assert A('1.21')>(A(4)/3+A(40)/37)/2
assert A('1.947')*(A('1.16')+A('2.21')/(2*c0))<3
assert 3*(2*pi).sqrt()+2<10
assert 1+2*A(2).exp()/3<6
assert 9*(A(4)/3)**(A(3)/2)<14
assert 3/c0**2+14/c0**4<2
assert 256/c0**2<1
near=A('4.5')+3/c0**2+12*A('3.12')*(A(17)/14)**2*(A(2)/3)*2**(A(3)/2)+24*A('3.12')*2*A(2).sqrt()
assert near<400
checks['near_integral_coefficient']=near
ratio=K(c50).sqrt()/c50**3/Z(c50).sqrt()
assert ratio<1
checks['K_over_R2Z_at_50']=ratio**2
assert (2*c50**3).log()+1<5*lo
assert 2*(2*c50**3).log()+8<12*lo
assert A('.5')+4/lo<A('1.53')
assert 5/lo+2*A('1.53')*8<26
assert (A('.25')+5*A(2).sqrt()/4)**2<8
def eta(c):return d0*(8/(ell*c)).sqrt()+z/c
kn=(A('.109').sqrt()-10*eta(c50)*U(c50).sqrt())**2
assert kn>A('.1089')
rows=[]
for mu in [50,100,1000,M]:
    c=2*pi*mu;a=A(mu).log()/2;R=c**3;L=(2*R).log()
    Y=(Z(c).sqrt()+A(2).sqrt()*(a+4)*(K(c).sqrt()+Z(c).sqrt()/4))**2
    A0=L*Z(c)+K(c)/R**2;A1=2*(a+4)**2*L*Z(c)+Y/R**2
    qlow=48*pi*(A0+2*(A0*A1).sqrt())
    blow=4*DB(c)*S(c*c);hi=high(mu);direct=qlow+blow+hi
    envelope=CP*A(mu).log()**4
    assert direct<envelope
    rows.append(dict(mu=mu,J2_B_upper=up(DB(c),0),qnorm_B_upper=up(Z(c),0),
                     log10_K_upper=up(K(c).log()/A(10).log(),6),
                     low_Q_upper=up(qlow,0),low_boundary_upper=up(blow,0),
                     high_zeros_upper=up(hi,0),
                     direct_P_upper=up(direct,0),P_tilde_upper=up(envelope,0)))
result=dict(certified=True,arb_bits=256,beta=4,smoothing_power=2,d="c^2",R="c^3",
            mu_min=50,mu_max=M,mu_max_is_optimal=False,log_power=4,
            ZC_upper=str(ZC.unique_fmpz()),derived_CP_upper=str(CPderived.unique_fmpz()),CP_upper=str(CP.unique_fmpz()),
            CE_upper=str(CE.unique_fmpz()),frequency_A_upper=up(Af,6),frequency_C_upper=up(Cf,6),
            frequency_ga_upper=up(ga,6),checks_upper={k:up(x,9) for k,x in checks.items()},
            rows=rows,limitation="Constants are conservative. High zero sum is absorbed into the logarithmic envelope on the stated finite interval, not asserted < B.")
print(json.dumps(result,ensure_ascii=False,indent=2))
