"""Paper III, Sections 6-7: symbolic checks (SymPy) of S_4 = 572587/414720, of the derivatives used in the monotonicity
arguments, and of the normalisation of A_1.
Usage: python proofs-III/all_mu_identities.py"""
import sympy as s,json
r,x,l,b,Z,R,a,v,Cd,c=s.symbols('r x l b Z R a v Cd c', positive=True)
f=1/(1-r);total=f
for k in range(1,5):
 f=r*s.diff(f,r);total+=s.binomial(4,k)*f/40**k
S4=s.factor(total.subs(r,s.Rational(1,4)))
assert S4==s.Rational(572587,414720)
j=(s.log(b*x/l)+s.exp(-1))/2
jdiff=s.simplify(s.diff(j/x,x)-(1-s.log(b*x/l)-s.exp(-1))/(2*x*x))
n=(1/(4*s.pi)+s.Rational(222,1000))*x+s.Rational(55,100)*s.log(x)+s.Rational(49,10)+s.Rational(2,5)/s.E
ndiff=s.simplify(s.diff(n/x,x)-(s.Rational(55,100)*(1-s.log(x))-s.Rational(49,10)-s.Rational(2,5)/s.E)/x**2)
# Z here denotes the constant Z_C; set sqrt K/(R sqrt{Z_C}l)=v.
Y=Z*l*l*(1+s.sqrt(2)*(a+4)*(v*R+s.Rational(1,4)))**2
A1=2*(a+4)**2*s.log(2*R)*Z*l*l+Y/R**2
A1split=2*Z*((a+4)/l)**2*s.log(2*R)/l+Z/l*(1/(R*l)+s.sqrt(2)*(a+4)/l*(v+1/(4*R)))**2
adiff=s.simplify(A1/l**5-A1split)
D=16*c+Cd*s.log(c)**2
Dder=s.factor(s.diff(D/c**3,c))
assert jdiff==ndiff==adiff==0
# The two cases of max(T0,c^3)^2/(T0^(2/3)c^4), expressed using u=c^3/T0.
u=s.symbols('u',positive=True)
print(json.dumps({'S4':str(S4),'j_over_x_derivative_difference':str(jdiff),'n_over_x_derivative_difference':str(ndiff),'A1_normalization_difference':str(adiff),'D_over_c3_derivative':str(Dder),'T_interpolation_ratio_when_u_le_1':str(u**(-s.Rational(4,3))),'T_interpolation_ratio_when_u_ge_1':str(u**s.Rational(2,3))},indent=2))
