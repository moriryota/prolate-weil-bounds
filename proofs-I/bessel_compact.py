"""Certified whole compact interval: Arb values + |H'|<=1/pi (BL1)."""
from flint import arb, ctx
from pathlib import Path
import json,time
ctx.prec=128
pi=arb.pi(); eps=arb(1)/64; X=64; N=(X*1024-16)
maximum=arb(0); at=None
start=time.monotonic()
for k in range(16,64*1024):
    x=arb(2*k+1)/2048
    H=x*x/2*(x.bessel_j(0)*x.bessel_y(0)+x.bessel_j(1)*x.bessel_y(1))
    cap=abs(H).upper()+1/(2048*pi)
    assert cap<arb('.160'),(k,H,cap)
    if cap>maximum:
        maximum=cap; at=str(x)
    if k%8192==0: print('certified through x=',k/1024,flush=True)
print('interval [1/64,64], cells=',N,'each width=1/1024')
print('maximum certified cell bound',maximum,'at midpoint',at)
print('BL2 compact <0.160; global h0=0.162 after analytic endpoints.')
result={'precision_bits':ctx.prec,'cells':N,'width':'1/1024','compact_upper':'0.160','max_cell_bound':str(maximum),'midpoint':at,'global_h0':'0.162','seconds':time.monotonic()-start}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
