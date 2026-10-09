# 0786 (Claude): independent spot check of J[0][0] (even) for mu=12: band (a-gamma) part + poles + U.
import gzip, json
from mpmath import mp, mpf, log, sqrt, quad, sin, sinh, pi, psi, mpc, re, exp, linspace
mp.dps = 40
L = log(12)/2; g = mpf(65)/16; Om = 416; c = Om*L
b0sq = lambda t: 2*L*(sin(L*t)/(L*t))**2 if t != 0 else 2*L
a = lambda t: re(psi(0, mpc(mpf(1)/4, t/2))) - log(pi)
pts = linspace(0, Om, 833)
K = quad(lambda t: (a(t) - g)*b0sq(t), pts)/pi
inband = quad(b0sq, pts)/pi
e0 = 1 - inband
d = exp(-2*(c+1))/(4*(c+1)*(c+4))
pole = 2*(4*sinh(L/2)/sqrt(2*L))**2
J = K + pole + 207*(e0 - d)**2
src = mpf(json.load(gzip.open('notes/0785_gpt_A4_mu12_certificate/sources/n832p1920q256_p0_J.json.gz', 'rt'))[0][0].strip('[').split(' ')[0])
print('J00 indep', mp.nstr(J, 25), ' source', mp.nstr(src, 25), ' diff', mp.nstr(abs(J - src), 3))
print('e0', mp.nstr(e0, 15), ' bounds', mp.nstr((c-1)/(pi*c*c), 15), mp.nstr((c+1)/(pi*c*c), 15))
assert abs(J - src) < mpf(10)**-20
print("J SPOT CHECK PASS")
