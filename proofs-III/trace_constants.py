"""Paper III, Lemma 5.1: the bound J^2 <= D_J(c) B <= (16c + 337392159829 (log c)^2) B, and the zero-count
constants used in Section 6. Arb, directed rounding asserted.
Usage: python proofs-III/trace_constants.py"""
import json
from flint import arb,ctx
ctx.prec=256
A=arb;pi=A.pi();ell=A(395)/399;c0=10*pi
bl=A(4);d=(8*bl).sqrt()*(2*bl).exp();j=8*bl*(2*bl).exp()
z=j*(8/pi).sqrt()/ell;v=2+(8*bl+16)*(2*bl).exp()
def up(x,places=6):
    t=10**places;n=int((x.upper()*t).ceil().unique_fmpz())
    s=str(n) if not places else ('-' if n<0 else '')+f'{abs(n)//t}.{abs(n)%t:0{places}d}'
    assert A(s)>=x.upper()
    return s
def down(x,places=6):
    t=10**places;n=int((x.lower()*t).floor().unique_fmpz())
    s=str(n) if not places else ('-' if n<0 else '')+f'{abs(n)//t}.{abs(n)%t:0{places}d}'
    assert A(s)<=x.lower()
    return s
def V(c):return v*(16*pi/ell).sqrt()+2*z/ell.sqrt()*(4/c+2/A(3).sqrt())
def E(c):return V(c)+A('1.78')*A(8).sqrt()*(pi/2+1+2*(c/3).log()+A('5.49'))+A('2.0408').sqrt()/2
def D(c):return ((8*c).sqrt()+E(c))**2
def eta(c):return d*(8/(ell*c)).sqrt()+z/c
def m0(c):return ((2*(592+A('152.2')*(2*(c/3).log()+A('8.08'))**2)).sqrt()+V(c))**2
CE=(V(c0)+A('1.78')*A(8).sqrt()*(pi/2+1+A('5.49')-2*A(3).log())+A('2.0408').sqrt()/2)/c0.log()+2*A('1.78')*A(8).sqrt()
t=A(14);Zplus=(t.log()+1)/(pi*t)+A('.386')*(t.log()+A('.5'))/t**2+(A('3.325')+A('.2')/A(1).exp())/t**2
assert Zplus<A('.107')
c50=100*pi
pert=eta(c50)*(6621*c50**(A(9)/2)*(-2*c50).exp()).sqrt()
kH=(A('.109').sqrt()-A(100).sqrt()*pert)**2
assert kH>A('.1089')
rows=[]
for mu0 in [5,7,10,50,100,1000]:
    c=2*pi*mu0;kap=(1+15/c**2)/(1-9/((1-A('3.6e-9'))*c)).sqrt()
    rows.append(dict(mu=mu0,E_upper=up(E(c),3),J2_over_B_upper=up(D(c),0),
                     kappa_upper=up(kap,9),m0_upper=up(m0(c),0)))
print(json.dumps(dict(certified_constants=True,precision_bits=256,
    CE_upper=up(CE,6),two_CE_squared_upper=up(2*CE**2,0),
    positive_zero_reciprocal_sum_upper=up(Zplus,9),
    simple_positive_zero_sum_upper="0.107",
    norm_k_over_original_H2_lower=down(kH,4),
    conditional_E_conversion_constant_upper=up(6621*(2*pi)**(A(9)/2)/A('.1089'),0),
    rows=rows,P_tilde=None,mu_max=None,
    status="No polylog trace estimate exists uniformly at infinity; C-prime/E-prime not proved."),ensure_ascii=False,indent=2))

