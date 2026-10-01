import sympy as S
import mpmath as mp
from pathlib import Path
import json
Y,s=S.symbols('y s',positive=True);a=1-s;r=Y*Y+a;v=Y*Y+1;g=S.sqrt(r/v)
l=1/(2*Y)+Y/(2*r)
B=S.factor((S.diff(l,Y)+l*l-l*S.diff(g,Y)/g)/g**2)
C=(Y**4+(5-8*s)*Y*Y+2*s*s-6*s+4)/(4*r**3)
assert S.simplify(B+v/(4*Y*Y*r)-C)==0
u=S.symbols('u');Cu=(u*u+(5-8*s)*u+2*s*s-6*s+4)/(4*(u+a)**3)
assert S.simplify(-4*(u+a)**4*S.diff(Cu,u)-(u*u+(8-14*s)*u+7-5*s-2*s*s))==0
assert S.simplify(C-(1/(4*r)+(3-6*s)/(4*r*r)+5*s*a/(4*r**3)))==0
W=2*Y*S.sqrt(v)/r**S.Rational(5,2)
assert S.simplify(S.diff(W,Y)*r**S.Rational(7,2)*S.sqrt(v)/2-(a+(2*a-4)*Y*Y-3*Y**4))==0
Ctheta,Stheta,h=S.symbols('C S h')
J0=Ctheta+Stheta*h/8-9*Ctheta*h*h/128
J1=Stheta+3*Ctheta*h/8+15*Stheta*h*h/128
squares=S.Poly(S.expand(J0**2+J1**2),h)
assert squares.coeff_monomial(h)==Ctheta*Stheta
assert squares.coeff_monomial(h*h)==Stheta**2/4
print('Exact Delta decomposition, C monotonicity polynomial, omega unimodality and R logarithmic term: passed.')
# Integral residual is the negative t derivative shown in PROOF.md.
z,t,alpha=S.symbols('z t alpha')
residual=(alpha-1)*z*z*t*t+(2*z+1)*t*(1+z*t)-(alpha+1)*(1+z*t)**2
boundary_derivative=(alpha+1)*(1+z*t)+(alpha-1)*z*t-t*(1+z*t)
assert S.expand(residual+boundary_derivative)==0
print('Exact integral ODE residual: passed.')
mp.mp.dps=40;rows=[]
for nu in [0,1]:
    al=mp.mpf(nu)-mp.mpf('.5')
    for xx in ['1','3.7','64']:
        x=mp.mpf(xx)
        V=mp.quad(lambda t:mp.exp(-t)*t**al*(1+1j*t/(2*x))**al,[0,1,mp.inf])/mp.gamma(al+1)
        H=mp.sqrt(2/(mp.pi*x))*mp.exp(1j*(x-nu*mp.pi/2-mp.pi/4))*V
        ref=mp.hankel1(nu,x);diff=abs(H-ref)
        # End point singular quadrature is a numerical diagnostic, not a proof.
        assert diff<mp.mpf('1e-20')
        row={'nu':nu,'x':xx,'integral_value':str(H),'mpmath_value':str(ref),'difference':str(diff),'modulus_ratio':str(abs(H)/abs(ref))}
        rows.append(row);print(json.dumps(row),flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps({'symbolic_assertions':True,'numerical_hankel_checks_certified':False,'rows':rows},indent=2)+'\n')
