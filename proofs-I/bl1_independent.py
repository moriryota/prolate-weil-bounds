"""Independent BL1 certificate: no Nicholson integral used by the computation."""
from flint import arb,ctx
from pathlib import Path
import json,time
ctx.prec=128
pi=arb.pi();eps=arb(1)/64;z=eps*eps/4
# x*log(2/x)^j is increasing for j=1,2 on this small interval.
L=(2/eps).log()+1
origin=pi*eps/2*(2*z).exp()*(1+4/pi**2*(L+z)**2)
assert origin<arb('.4')
print('0<x<=1/64: bound=',origin,'<0.4',flush=True)
# H0^(1)=sqrt(2/(pi*x))*e^(i(x-pi/4))*(1-i/(8x)-9/(128x²)+e3)
# |e3|<=75/(1024x³), by DLMF 10.17(iii), ell=3, nu=0.
u=arb(1)/128+arb(0,arb(1)/128)
d2=arb(9)/128;d3=arb(75)/1024
bracket=-arb(1)/8+d2*d2*u*u+2*(1+u/8+d2*u*u)*d3*u+d3*d3*u**4
assert bracket<0
print('x>=64: |normalized Hankel|²<=1+x^-2*bracket; bracket=',bracket,'<0',flush=True)
maximum=arb(0);maxat=None;started=time.monotonic()
for k in range(16,64*1024):
    mid=arb(2*k+1)/2048;x=arb(mid,arb(1)/2048)
    jm=mid.bessel_j(0);ym=mid.bessel_y(0)
    Nm=pi*mid/2*(jm*jm+ym*ym)
    j1m=mid.bessel_j(1);y1m=mid.bessel_y(1)
    derivative=pi/2*(jm*jm+ym*ym-2*mid*(jm*j1m+ym*y1m))
    lo=arb(k)/1024;hi=arb(k+1)/1024
    j0l=lo.bessel_j(0);y0l=lo.bessel_y(0);j1l=lo.bessel_j(1);y1l=lo.bessel_y(1)
    energy=j0l*j0l+y0l*y0l+j1l*j1l+y1l*y1l
    # For either solution C: (C0²+C1²)'=-2C1²/x, so energy decreases.
    # |N''| <= (pi/2)*(1+2hi)*energy(lo); second order Taylor encloses the whole cell.
    bound=Nm.upper()+abs(derivative).upper()/2048+pi*(1+2*hi)*energy/(4*2048**2)
    assert bound<1,(k,Nm,derivative,bound)
    if bound>maximum: maximum=bound;maxat=str(mid)
    if k%8192==0: print('certified BL1 through x=',k/1024,flush=True)
print('compact max certified bound=',maximum,'at',maxat,flush=True)
print('BL1 for every x>0 passed independently of Nicholson integral.')
Path(__file__).with_suffix('.json').write_text(json.dumps({'origin_bound':str(origin),'tail_bracket':str(bracket),'compact_bound':str(maximum),'midpoint':maxat,'cells':65520,'bits':ctx.prec,'seconds':time.monotonic()-started,'Nicholson_used':False},indent=2)+'\n')
