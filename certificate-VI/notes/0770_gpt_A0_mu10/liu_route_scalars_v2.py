"""New mu=10 scalar checks for a Liu-style bounded comparison. No Liu data reused."""
from fractions import Fraction as F
from flint import arb,ctx
import json,math
ctx.prec=384
A=lambda q:arb(q.numerator)/q.denominator if isinstance(q,F) else arb(q)
ns=[2,3,4,5,7,8,9];base={2:2,3:3,4:2,5:5,7:7,8:2,9:3};L=arb(10).log()/2
assert L<A(F(7,6)) and L>A(F(17,16))
events=sorted({F(1),F(10),*[max(F(10,n*n),F(n*n,10)) for n in ns]})
w=lambda x:1+A(F(19,40))*(x/L)**2
cells=[]
for a,b in zip(events,events[1:]):
 mid=(a+b)/2
 sh=[(n,1) for n in ns if mid<F(10,n*n)]+[(n,-1) for n in ns if mid>F(n*n,10)]
 lo=A(a).log()/2;hi=A(b).log()/2;worst=None
 for j in range(256):
  x=(lo+(hi-lo)*j/256).union(lo+(hi-lo)*(j+1)/256)
  row=sum(arb(base[n]).log()/arb(n).sqrt()*w(x+sign*arb(n).log()) for n,sign in sh)/w(x)
  assert row<3,(a,j,row)
  if worst is None or float(row.upper())>float(worst.upper()):worst=row
 cells.append({'exp_2x_interval':[str(a),str(b)],'row_upper_enclosure':str(worst)})
# Exact integer comparison for Chebyshev projection at N=640, rho=3/2, Omega=256.
# 72 L exp(5 Omega L/6) (2/3)^(2N) <= 84 (68/25)^249 (2/3)^1280 < 2^-380.
projection_lhs=84*68**249*2**(1280+380)
projection_rhs=25**249*3**1280
assert projection_lhs<projection_rhs
# Lagrange continuation proof: 3*768 L/4096 < 2/3, (2/3)^2 <1/2.
assert F(3*768,4096)*F(7,6)<F(2,3)
assert F(4,9)<F(1,2)
# Lagrange sum bounded by 64^4095 (same integer inequality, new L-independent check).
assert 3*49152<64*4095
# Rank-two tail: c in [272, 896/3], retain 81 and 27 safely.
C=F(123,1280);c=F(272)
e0_hi=50*(c+1)/(157*c*c)
e1_hi=50*(3*c*c+6*c+2)/(157*c*c*c)
assert 81*e0_hi<C and 27*e1_hi<C
assert e0_hi<F(1,29**2) and e1_hi<F(9,2500)
# Lower e0,e1 bounds throughout c range, using pi<22/7 and monotonicity for c>4.
c=F(896,3)
assert (c-1)/(F(22,7)*c*c)>F(1,2**11)
assert 3*(c-2)/(F(22,7)*c*c)>F(1,2**11)
assert C>F(1,16)
# Pole projection, norms ||a+/-||<3, remainder norm <4/N!, yields <48/N!.
assert F(48,math.factorial(640))<F(1,2**440)
r=F(1,2**190);p=F(1,2**440);e=1500*r+p;n=896*r*r+p;h=12000*r*r+p
assert F(1,2)-n>0
# theta=1, chi=1, m=1/2. Error coefficients alone, not matrix acceptance.
theta=F(1);m=F(1,2)
alpha=e*theta+2*e*e/(m-n)
beta=(e/theta+n)/(m*m)+2*h*h/(m*m*(m-n))
assert alpha+beta*F(13,2)**2<F(1,2**170)
assert arb(1).exp()<A(F(68,25))
assert F(448,11)**5>F(68,25)**18
result={'status':'PASS scalar preconditions only; finite target matrix not assembled or certified',
 'L':str(L),'prime_row_weight':'1+(19/40)(x/L)^2','prime_operator_bound':3,'cells':cells,
 'M_bounds':['1/2','13/2'],'Omega':256,'N_all_degrees':640,'per_parity_dimension':320,
 'projection_r':'2^-190','theta':'1','pole_p':'2^-440','uniform_band_defect_delta':'2^-49158',
 'uniform_weil_target_tau':'2^-49162 only if remaining matrix contract passes',
 'retained_rank_weights':[81,27], 'error_e':'1500 r+p','error_n':'896 r^2+p','error_h':'12000 r^2+p',
 'uniform_error_budget_alpha_plus_beta_b_squared':str(alpha+beta*F(13,2)**2),'uniform_error_budget_lt_2_pow_minus_170':bool(alpha+beta*F(13,2)**2<F(1,2**170)), 'alpha_float_diagnostic':float(alpha),'beta_float_diagnostic':float(beta),
 'projection_integer_slack_bits_diagnostic':projection_rhs.bit_length()-projection_lhs.bit_length()}
print(json.dumps(result,indent=2))
