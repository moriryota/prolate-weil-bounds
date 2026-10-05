"""Paper III, Section 3: constants of the layer (Lemma 3.1, (3.2)), the coefficients (Lemma 3.2),
Proposition 3.3 (eta, theta, V, the norm bound 0.1089) and m_0(c) in (3.3). Arb, directed rounding asserted.
Usage: python proofs-III/layer_constants.py"""
import json
from fractions import Fraction as F
from math import comb
from flint import arb, ctx
ctx.prec=256
A=arb
pi=A.pi();ell=A(395)/399;bl=A(4);c0=10*pi
d=(8*bl).sqrt()*(2*bl).exp()
j=8*bl*(2*bl).exp()
z=j*(8/pi).sqrt()/ell
v=2+(8*bl+16)*(2*bl).exp()
def up(x,digits=6):
    scale=10**digits
    n=int((x.upper()*scale).ceil().unique_fmpz())
    y=A(n)/scale
    assert y>=x.upper(), ("upper rounding",x,y)
    result=str(n) if digits==0 else ("-" if n<0 else "")+f"{abs(n)//scale}.{abs(n)%scale:0{digits}d}"
    assert A(result)>=x.upper(), ("printed upper bound", result, x)
    return result
def down(x,digits=6):
    scale=10**digits
    n=int((x.lower()*scale).floor().unique_fmpz())
    y=A(n)/scale
    assert y<=x.lower(), ("lower rounding",x,y)
    result=str(n) if digits==0 else ("-" if n<0 else "")+f"{abs(n)//scale}.{abs(n)%scale:0{digits}d}"
    assert A(result)<=x.lower(), ("printed lower bound", result, x)
    return result
def eta(c):return d*(8/(ell*c)).sqrt()+z/c
def vv(c):return v*(16*pi/ell).sqrt()+2*z/ell.sqrt()*(4/c+2/A(3).sqrt())
def m0(c):
    k1=2*(c/3).log()+A('8.08')
    m=592+A('152.2')*k1*k1
    return ((2*m).sqrt()+vv(c))**2
def theta(c):return A('2.0408').sqrt()+2*z/(ell.sqrt()*c)
vk=2-8*bl/c0**2
C1=(2*(8*bl/vk).sqrt()).exp();C2=C1/vk
r=F(1,4)
mom=[1/(1-r),r/(1-r)**2,r*(1+r)/(1-r)**3,
     r*(1+4*r+r*r)/(1-r)**4,r*(1+11*r+11*r*r+r**3)/(1-r)**5]
S=sum(F(comb(4,k),40**k)*mom[k] for k in range(5))
Sarb=A(S.numerator)/S.denominator
T0=A(3)*10**12;x0=(3*T0).log()
assert x0>40*A(2).log()
nx=(1/(4*pi)+A('.224'))*x0+A('.556')*x0.log()+A('5.02')+A('.4')/A(1).exp()
At=2*(A(1000)/999)**2*nx*(A('11.15')*x0)**2/T0**2*Sarb
rows=[]
for mu0 in [50,100,1000,10000]:
    mu=A(mu0);lam=mu.sqrt();c=2*pi*mu
    U=6621*c**(A(9)/2)*(-2*c).exp()
    q=eta(c)*U.sqrt()
    ratio=((A('.109').sqrt()-(2*mu).sqrt()*q)/(1+q))**2
    err=2*A('.109').sqrt()*(1+q)*(A('.109').sqrt()+(2*mu).sqrt())*q
    assert A('.109').sqrt()>(2*mu).sqrt()*q
    assert ratio>A('.1089')
    rows.append(dict(mu=mu0,C3_upper=up((1+eta(c))**2,0),
        C4_upper=up(theta(c)**2,0),m0_upper=up(m0(c),0),
        eta_upper=up(eta(c),3),log10_q_upper=up(q.log()/A(10).log(),3),
        log10_error_upper=up(err.log()/A(10).log(),3),
        norm_ratio_lower=down(ratio,4),
        conditional_tail_coefficient_times_1e20_upper=up(At*(A('.5')+lam/2*(A('11.15')*x0/2).log())*10**20,6)))
result=dict(interval_certified=True,precision_bits=ctx.prec,
    beta=4,power=2,K=8,C1_upper=up(C1,0),C2_upper=up(C2,0),
    d_upper=up(d,3),j_upper=up(j,3),z_upper=up(z,3),v_upper=up(v,3),
    uniform_C3_upper=up((1+eta(c0))**2,0),uniform_C4_upper=up(theta(c0)**2,0),
    S4_exact=str(S),tail_A_times_1e20_upper=up(At*10**20,6),
    rows=rows,mu_max=None,
    caution="The 10000 row is a reference point, not mu_max. Tail coefficient is conditional on unproved K. No C-prime or E-prime is certified.")
print(json.dumps(result,ensure_ascii=False,indent=2))
