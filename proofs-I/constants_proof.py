"""Arb certificates for the analytic bounds in PROOF.md; no sampled extrema."""
from flint import arb, ctx
from pathlib import Path
import json
ctx.prec=256
pi=arb.pi(); a=arb(7)/10; s=arb(3)/10
D0=(3-s)/(3*a*a)
D1=pi*(20-17*s)/(64*a**arb('1.5'))+pi*s/(16*a*a)
D2=D0
# BL3 uses E(x)=R(x)+cos(2x)/(2pi), |E|<=2 (0<x<1),
# |E|<=1/(4x) (x>=1), and a single maximum of omega.
omega_max=9/(8*arb(3).sqrt()*a*a)
r1=omega_max/(2*pi)+pi/(8*a*a)+2/(a**3*10*pi)
values={'D0':D0,'D1':D1,'D2':D2,'r1':r1}
caps={'D0':arb('1.837'),'D1':arb('1.370'),'D2':arb('1.837'),'r1':arb('1.2')}
for key,value in values.items():
    assert value<caps[key], (key,value,caps[key])
    print(key,'Arb=',value,'safe upper=',caps[key])
# Check all coefficients in -dC/d(y²) positive throughout the s interval.
S=arb('0.15','0.15')
assert (8-14*S)>0
assert (7-5*S-2*S*S)>0
# dF0/ds=(5-s)/(3(1-s)^3)>0; dIC/ds=(26-17s)/(128(1-s)^(5/2))*pi>0.
assert (5-S)>0 and (26-17*S)>0
# Elementary upper bounds used in E on x>=1: |P0|<=9/8, |P1|<=11/8;
# combined Hankel remainders e0<=9/(128x²), e1<=15/(128x²).
A=(arb(9)/64+2*(arb(9)/8*arb(9)/128+arb(11)/8*arb(15)/128)
   +(arb(9)**2+arb(15)**2)/arb(128)**2)/pi
assert A<arb(1)/4
print('BL3 E coefficient via Hankel <=',A,'<1/4')
# Near 0: J0,J1 <=1 and x<1, giving |R|<=1+1/pi, |E|<=1+3/(2pi)<2.
assert arb(5)/8*arb('.5').exp()+3/(2*pi)<2
print('BL3 origin E bound',arb(5)/8*arb('.5').exp()+3/(2*pi),'<2')
# BL2 tail: Pnu=1+ia1/x, |epsnu|<=|a2|/x², x>=64.
X=arb(64)
b0=1+arb(1)/(8*X); b1=1+arb(3)/(8*X)
tailH=(1+arb(1)/(8*X)+2*(b0*arb(9)/128+b1*arb(15)/128)/X
       +(arb(9)**2+arb(15)**2)/(arb(128)**2*X**3))/(2*pi)
assert tailH<arb('.162')
print('BL2 tail x>=64:',tailH,'<0.162')
# BL2 origin x<=1/64; inequality derived directly from convergent Frobenius series.
eps=arb(1)/64; z=eps**2/4
originH=eps**2/pi*(2*z).exp()*((1+2*z)*((2/eps).log()+1)+arb(1)/2+z*(1+z))
assert originH<arb('.001')
print('BL2 origin x<=1/64:',originH,'<0.001')
result={key:{'arb':str(value),'safe_upper':str(caps[key])} for key,value in values.items()}
result.update({'BL2_tail':str(tailH),'BL2_origin':str(originH),'BL3_E_coefficient':str(A),'precision_bits':ctx.prec})
Path(__file__).with_suffix('.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print('All analytic Arb assertions passed.')
