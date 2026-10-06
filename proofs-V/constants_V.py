"""Paper V: explicit constants of Sections 5-7 (Proposition, Theorems J and K), with the zero-counting constants of the published version of Trudgian (2014): 0.112, 0.278, 2.510. Arb 256-bit.
Run from repository root.
"""
from flint import arb, ctx
ctx.prec=256
A=arb
pi=A.pi(); e=A(1).exp(); cmin=100*pi; cmax=20000*pi
l0=A(50).log(); T0=A(3)*10**12; ell=A(395)/399
z=32*e**8*(8/pi).sqrt()/ell

def V(c): return (2+48*e**8)*(16*pi/ell).sqrt()+2*z/ell.sqrt()*(4/c+2/A(3).sqrt())
def EJ(c): return V(c)+A('1.78')*A(8).sqrt()*(pi/2+1+2*(c/3).log()+A('5.49'))+A('2.0408').sqrt()/2
def D(c): return 16*c+A(337392159829)*c.log()**2
def S(d): return ((A(4)/3)*d.log()+1)/(pi*d)+(A('.780')*d.log()+2*(A('3.385')+A('.2')/e)+A('.195'))/d**2
mbar_raw=(A(1458).sqrt()+V(cmin))**2
Zbar_raw=(A(1458).sqrt()+V(cmin)+(8/cmin).sqrt()+EJ(cmin)/cmin)**2
Z=A(2)*10**12; mbar=Z
assert mbar_raw < mbar and Zbar_raw < Z
Af=10+48*e**8; Cf=20+321*e**8; gA=3*(Af+2*Cf)/(2*pi).sqrt()
Kbar=(Z.sqrt()/cmin**2+gA+120/cmin.sqrt()+D(cmin).sqrt()/cmin)**2
kY=A(32)/e**2
Ybar=(Z.sqrt()/cmin**2+kY.sqrt()*(Kbar.sqrt()+Z.sqrt()/(4*cmin**2)))**2
A0=5*Z+Kbar/(cmin**2*l0)
A1=5*kY*Z+Ybar/(cmin**2*l0)
CL=48*pi*(A0+2*(A0*A1).sqrt())+4*D(cmin)*S(cmin**2)/l0
S4=A(572587)/414720
x0=(3*T0).log()
def n(x): return (1/(4*pi)+A('.224'))*x+A('.556')*x.log()+A('5.02')+A('.4')/e
def j(x,l): return ((A('11.15')*x/l).log()+1/e)/2
At0=2*(A(1000)/999)**2*S4*n(x0)*(A('11.15')*x0)**2/T0**2
LT100=A(1)/2+50*(A('11.15')*x0/2).log()
HT=(T0.log()+1)/(pi*T0)
CFA=CL+(2*At0*LT100*Kbar*cmax**4+400*D(cmax)*HT)/l0
mustar=T0**(A(1)/3)/(2*pi); ls=mustar.log()
CM=2*(n(x0)/ls)*(A('11.15')*x0/ls)**2*j(x0,ls)*mbar
CT=4*(A(1000)/999)**2*S4*(n(x0)/l0)*(A('11.15')*x0/l0)**2*j(x0,l0)*Kbar/(T0**(A(2)/3))
CB=4/pi*(D(cmin)/cmin**3)*(x0/l0)
CFB=CL/(A(50).sqrt()*l0**2)+CM+CT+CB/l0**2
# Bounds used in the elementary exact-phase argument.
delta_derivative=(4/(A('.75')**2*A('.925')**2)+6/A('.925')**3+6*(3+A('.69')/4)/A('.925')**4)/4
edge=8*A('7.58')*(2+A(3).sqrt()).log()
assert delta_derivative < 12
assert edge < 81
assert A('2.0408')/4 < A('.72')**2
assert (4*A('1.78')*A(7)/3+9+A('.72'))**2 < 729
assert A('0.00464')+A('0.3829')+A('0.00464')*A('0.3829')/4 < A('.4')
assert pi**2/6 < A(5)/3
for name,val in locals().copy().items():
    if name in ['mbar_raw','Zbar_raw','Kbar','Ybar','CL','CFA','CM','CT','CB','CFB','delta_derivative','edge']:
        print(name, val)
print('G_a',A(237522697)*CFA)
print('G_b',A(237522697)*CFB)
# Deliberately rounded-up published-candidate constants, asserted as balls.
assert CFA < A('1.47')*10**16
assert CFB < A('3.19')*10**16
assert A(237522697)*A('1.47')*10**16 < A('3.50')*10**24
assert A(237522697)*A('3.19')*10**16 < A('7.58')*10**24
# Thresholds printed in the paper (section 7)
assert mbar_raw < A('1.906e12') and Zbar_raw < A('1.918e12')
assert CL < A('7.785e15')
assert CM < A('3.175e16')
assert CT < A('1.165e10')
assert CB < A('3.493e6')
assert (CFA-CL) < A('6.9e15') and (CFA-CL) > A('6.7e15')   # 'about 6.8e15' from the zeros above T0
print('ALL ASSERTIONS PASS')
