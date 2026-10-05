"""Paper III: P''(mu) of (7.1), the envelope constants C_L, C_M, C_T, C_B, Theorems F(b), G(b), and the comparison
G(c) < 0.832 x [Paper II, Thm E] for all mu >= 50. Arb, directed rounding asserted.
Usage: python proofs-III/theoremF_b_constants.py"""
from flint import arb as A,ctx
import json
ctx.prec=320
pi=A.pi(); b=A('11.15'); ell=A(395)/399; T0=A(3)*10**12
l0=A(50).log(); c0=100*pi; x0=(3*T0).log(); mus=T0**(A(1)/3)/(2*pi); ls=mus.log()
ZC=A(125363711220); CD=A(337392159829); e8=A(8).exp(); z=32*e8*(8/pi).sqrt()/ell
V=lambda c:(2+48*e8)*(16*pi/ell).sqrt()+2*z/ell.sqrt()*(4/c+2/A(3).sqrt())
m0=lambda c:((2*(592+A('152.2')*(2*(c/3).log()+A('8.08'))**2)).sqrt()+V(c))**2
Dp=lambda c:16*c+CD*c.log()**2
Z=lambda c:ZC*(c/(2*pi)).log()**2
ga=3*(50+690*e8)/(2*pi).sqrt()
K=lambda c:(Z(c).sqrt()+ga*c*c+120*c**(A(3)/2)+c*Dp(c).sqrt())**2
CN=A('3.385')+A('.2')/A(1).exp()
S=lambda d:(A(4)/3*d.log()+1)/(pi*d)+(A('.780')*d.log()+2*CN+A('.195'))/d**2
n=lambda x:(1/(4*pi)+A('.224'))*x+A('.556')*x.log()+A('5.02')+A('.4')/A(1).exp()
S4=A(572587)/414720; kappa=(A(1000)/999)**2
Jfactor=lambda x,l:((b*x/l).log()+1/A(1).exp())/2
def up(x,digits=0):
 scale=10**digits;k=int((x.upper()*scale).ceil().unique_fmpz());y=A(k)/scale
 assert y>=x.upper()
 s=str(k) if digits==0 else ('-' if k<0 else '')+f'{abs(k)//scale}.{abs(k)%scale:0{digits}d}'
 assert A(s)>=x.upper();return s
def down(x,digits=0):
 scale=10**digits;k=int((x.lower()*scale).floor().unique_fmpz());y=A(k)/scale
 assert y<=x.lower()
 s=str(k) if digits==0 else ('-' if k<0 else '')+f'{abs(k)//scale}.{abs(k)%scale:0{digits}d}'
 assert A(s)<=x.lower();return s
def low(mu):
 c=2*pi*mu;l=A(mu).log();a=l/2;R=c**3;zz=Z(c);kk=K(c);LR=(2*R).log()
 Y=(zz.sqrt()+A(2).sqrt()*(a+4)*(kk.sqrt()+zz.sqrt()/4))**2
 A0=LR*zz+kk/R**2;A1=2*(a+4)**2*LR*zz+Y/R**2
 return 48*pi*(A0+2*(A0*A1).sqrt())+4*Dp(c)*S(c*c)
alpha0=x0/l0;j0=Jfactor(x0,l0)
# Only these rounded upper constants enter the global envelope.
CL=A(up(low(50)/l0**4));K4=A(up(K(c0)/c0**4))
Cm_exact=2*(n(x0)/ls)*(b*x0/ls)**2*Jfactor(x0,ls)*m0(2*pi*mus)/ls**2
Cm=A(up(Cm_exact))
Ct_exact=4*kappa*S4*(n(x0)/x0)*b*b*alpha0**3*j0*K4/T0**(A(2)/3)
Ct=A(up(Ct_exact))
Cb_exact=4*alpha0/pi*Dp(c0)/c0**3
Cb=A(up(Cb_exact))
CP_exact=CL/(A(50).sqrt()*l0)+Cm+Ct/l0**2+Cb/l0**4
CP=A(up(CP_exact));CPsimple=A(10)**13*A(up(CP/A(10)**13))
conversion=6621*(2*pi)**(A(9)/2)/A('.1089');Econv=A(237522697)
assert conversion<Econv
CEexact=Econv*CPsimple;CEsimple=A(10)**21*A(up(CEexact/A(10)**21))
advantage=Econv/A('5.046e17')*(CL*l0/A(50)**(A(7)/2)+Ct/A(50)**3+Cb/(A(50)**3*l0**2)+Cm*ls**2/mus**3)
assert advantage<A('.85')
# Assertions behind the analytic monotonicity reductions.
assert x0>40*A(2).log();assert x0>1;assert c0>3
assert mus>50;assert mus<10000;assert b*3>1;assert alpha0>x0/ls
assert K(c0)/c0**6/Z(c0)<A('.000032283')
assert (2*c0**3).log()/l0>3

def calc(mu):
 mu=A(mu);l=mu.log();lam=mu.sqrt();c=2*pi*mu
 if c**3>T0:T=c**3;middle=True
 else:T=T0;middle=False
 x=(3*T).log();jf=Jfactor(x,l);J=lam*jf
 mid=2*n(x)*(b*x)**2*J*m0(c) if middle else A(0)
 At=2*kappa*n(x)*(b*x)**2/T**2*S4
 tq=2*At*J*K(c);tb=4*lam*Dp(c)*(T.log()+1)/(pi*T)
 lo=low(mu);P=lo+mid+tq+tb
 env=CPsimple*lam*l**5
 # Original paper II P, with its independent exponentially small additive term omitted.
 xx=20*mu+1
 R=A('.112')*xx+A('.278')*xx.log()+A('2.510')+A('.2')/A(1).exp()
 N=(xx-(2*pi).log())/(4*pi)+2*R
 r=1/(b*xx)
 Pold=2*N/r**2*(A('.5')+lam/2*(1/(2*r)).log())*(592+A('152.2')*(2*(c/3).log()+A('8.08'))**2)
 Cprime=A('6.034e14')*l**4
 oldE=A('5.046e17')*mu**8*l**3
 newE=Econv*mu**(A(9)/2)*P
 assert P<env;assert newE<oldE
 return {'mu':str(int(mu.unique_fmpz())),'T_upper':up(T),'middle_nonempty':middle,'low_upper':up(lo),'middle_upper':up(mid),'tail_Q_upper':up(tq),'tail_b_upper':up(tb),'P_double_upper':up(P),'P_double_lower':down(P),'global_envelope_upper':up(env),'P_thmFa_formal_upper':up(Cprime),'P_thmFa_valid':bool(mu<=10000),'P_paperII_upper':up(Pold),'P_double_over_thmFa_upper':up(P/Cprime,12),'P_double_over_paperII_interval':str(P/Pold),'E_refined_over_paperII_interval':str(newE/oldE),'E_simple_over_paperII_interval':str(CEsimple*mu**5*l**5/oldE),'P_double_over_sqrt_mu_log5_interval':str(P/(lam*l**5))}
result={'arb_bits':ctx.prec,'mu_star_interval':str(mus),'q_double':5,'CP_simple':str(CPsimple.unique_fmpz()),'CE_simple':str(CEsimple.unique_fmpz()),'CP_raw_upper':str(CP.unique_fmpz()),'E_conversion_upper':str(Econv.unique_fmpz()),'global_improvement_factor_interval':str(advantage),'global_improvement_factor_upper':up(advantage,9),'constants':{'CL_upper':str(CL.unique_fmpz()),'Cm_upper':str(Cm.unique_fmpz()),'Ct3_upper':str(Ct.unique_fmpz()),'Cb1_upper':str(Cb.unique_fmpz()),'K4_upper':str(K4.unique_fmpz()),'alpha0_upper':up(alpha0,9),'j0_upper':up(j0,9)},'rounding_comparisons':{},'rows':[calc(mu) for mu in [50,1000,10**6,10**12]],'status':'New analytic proof candidate; no claim of independent review or formal verification.'}
for key,exact,bound in [('CL',low(50)/l0**4,CL),('Cm',Cm_exact,Cm),('Ct3',Ct_exact,Ct),('Cb1',Cb_exact,Cb),('K4',K(c0)/c0**4,K4),('CP',CP_exact,CP),('CP_simple',CP,CPsimple),('CE_simple',CEexact,CEsimple)]:
 assert bound>=exact.upper()
 result['rounding_comparisons'][key]={'exact_interval':str(exact),'upper':str(bound),'gap_interval':str(bound-exact),'ratio_interval':str(exact/bound)}
print(json.dumps(result,ensure_ascii=False,indent=2))
