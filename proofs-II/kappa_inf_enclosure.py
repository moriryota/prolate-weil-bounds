"""Hermite limit using exact polynomial algebra and Arb incomplete Gaussian moments.
Independent of the other numerical scripts.
"""
from pathlib import Path
import json
from decimal import Decimal,ROUND_FLOOR,ROUND_CEILING,localcontext
from flint import arb,ctx
import sympy as sy
ctx.prec=256
D=Path("outputs/proofs-II")
pi=arb.pi()
def bs(a): return a.str(70)

# Independent closed form of I4*g0-I0*g4.
x=sy.symbols("x",real=True)
g0=2**sy.Rational(1,4)*sy.exp(-sy.pi*x*x)
g4=g0*sy.hermite(4,sy.sqrt(2*sy.pi)*x)/sy.sqrt(2**4*sy.factorial(4))
I0=g0.subs(x,0); I4=g4.subs(x,0)
h=sy.simplify(I4*g0-I0*g4)
closed=4*sy.sqrt(3)*sy.pi*x*x*(1-sy.Rational(2,3)*sy.pi*x*x)*sy.exp(-sy.pi*x*x)
assert sy.simplify(h-closed)==0
print("HERMITE_LIMIT_CLOSED_FORM",h)
norm=11*arb(2).sqrt()/8
print("H_NORM_SQUARED",bs(norm))
print("H_ZERO",sy.simplify(h.subs(x,0)))
print("H_INTEGRAL",sy.integrate(h,(x,-sy.oo,sy.oo)))
for n in (0,2,4):
    gn=2**sy.Rational(1,4)*sy.hermite(n,sy.sqrt(2*sy.pi)*x)*sy.exp(-sy.pi*x*x)/sy.sqrt(2**n*sy.factorial(n))
    integral=sy.simplify(sy.integrate(gn,(x,-sy.oo,sy.oo)))
    value=gn.subs(x,0)
    print("HERMITE_INTEGRAL",n,integral,"VALUE_AT_ZERO",value,"RATIO",sy.simplify(integral/value))

# e_inf is even by Poisson and hhat=h, h(0)=integral h=0.
# Norm squared is 2*int_1^infty (sum_m h(mx))^2 dx.
Acoef=4*arb(3).sqrt()*pi
Bcoef=-8*arb(3).sqrt()*pi*pi/3
def moments(a,lower=arb(1)):
    # J_k=int_lower^infty x^(2k)*exp(-a*x^2)dx.
    aa=a*lower**2
    J=[pi.sqrt()/(2*a.sqrt())*(a.sqrt()*lower).erfc()]
    for k in range(1,5):
        J.append(lower**(2*k-1)*(-aa).exp()/(2*a)+(2*k-1)*J[-1]/(2*a))
    return J

def finite_sum_square(M,lower=arb(1)):
    total=arb(0)
    for m in range(1,M+1):
        for n in range(1,M+1):
            J=moments(pi*(m*m+n*n),lower)
            total+=Acoef**2*(m*n)**2*J[2]
            total+=Acoef*Bcoef*(m*m*n**4+m**4*n*n)*J[3]
            total+=Bcoef**2*(m*n)**4*J[4]
    return total

def tail_l2(M):
    j=M+1
    # For x>=1, all h(mx)<0 and |h(mx)|<=|B|m^4*x^4*exp(-pi*m^2*x^2).
    # Consecutive factors m^4 exp(-pi*m^2) have ratio <= q.
    q=(arb(j+1)/j)**4*(-pi*(2*j+1)).exp()
    assert q<1
    return abs(Bcoef)*j**4/(1-q)*moments(2*pi*j*j)[4].sqrt()

records=[]
for M in (3,4,6):
    s=finite_sum_square(M)
    t=tail_l2(M)
    lower=2*s.lower()/norm
    upper=2*(s.sqrt()+t)**2/norm
    # Lower uses same sign of every h(mx) on x>=1.
    assert upper>lower
    row=dict(M=M,finite_norm_squared=bs(2*s),tail_L2_upper=bs(t),
             kappa_lower_ball=bs(lower),kappa_upper_ball=bs(upper))
    records.append(row)
    print("RIGOROUS_KAPPA",json.dumps(row,ensure_ascii=False))
    if M==6: kl,ku=lower,upper

lo=arb("0.219247199548")
hi=arb("0.219247199550")
assert lo<kl and ku<hi
assert ku<arb("0.21925")
assert kl>arb("0.21924")
print("CERTIFIED_DECIMAL_INTERVAL",str(lo),str(hi))
print("PRINTED_0.21925_MINUS_KAPPA_INTERVAL",
      bs(arb("0.21925")-ku),bs(arb("0.21925")-kl))
original=arb("0.219247199549")
print("ORIGINAL_DISPLAY",original,"DISPLAY_MINUS_LOWER",bs(original-kl),
      "DISPLAY_MINUS_UPPER",bs(original-ku))

# Local share, exact finite moments with analytically bounded neglected series.
# For Y>=0 use [1,exp(Y)] and evenness. The discarded tail is smaller
# than the total tail_l2(M), so use a symmetric absolute enclosure.
shares=[]
M=6
sfull=finite_sum_square(M)
t=tail_l2(M)
full_upper=2*(sfull.sqrt()+t)**2
for Ytext in ("0.5","0.8","1.0","1.2"):
    Y=arb(Ytext)
    spart=sfull-finite_sum_square(M,Y.exp())
    # all summands share negative sign, so finite squared integral is a lower bound.
    share_lower=2*spart.lower()/full_upper
    share_upper=2*(spart.sqrt()+t)**2/(2*sfull.lower())
    rec=dict(Y=Ytext,share_lower=bs(share_lower),share_upper=bs(share_upper))
    shares.append(rec)
    print("LOCAL_SHARE",json.dumps(rec))
(D/"independent_limit.json").write_text(json.dumps(dict(
    precision_bits=ctx.prec,closed_form=str(h),norm_squared=bs(norm),
    truncations=records,certified_interval=["0.219247199548","0.219247199550"],
    shares=shares),ensure_ascii=False,indent=2)+"\n")
print("FINISHED")
