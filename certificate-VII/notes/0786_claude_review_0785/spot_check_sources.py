# 0786 (Claude): independent spot check of source entries A, B (j=0 of each parity) for mu=12.
import gzip, json
from mpmath import mp, mpf, log, sqrt, quad
mp.dps = 50
L = log(12)/2; g = mpf(65)/16
pp = {2:2,3:3,4:2,5:5,7:7,8:2,9:3,11:11}
v = {0: lambda x: 1/sqrt(2*L), 1: lambda x: sqrt(3/(2*L))*x/L}
ind = lambda x: 1 if -L <= x <= L else 0
def Cp(f, x):
    return sum(log(p)/sqrt(n)*(ind(x-log(n))*f(x-log(n)) + ind(x+log(n))*f(x+log(n))) for n, p in pp.items())
brk = sorted({-L, L, *[s*(L - log(n)) for n in pp for s in (1, -1) if log(n) < 2*L]})
for par in (0, 1):
    f = v[par]
    A = quad(lambda x: f(x)*(g*f(x) - Cp(f, x)), brk)
    B = quad(lambda x: (g*f(x) - Cp(f, x))**2, brk)
    src = {k: json.load(gzip.open(f'notes/0785_gpt_A4_mu12_certificate/sources/n832p1920q256_p{par}_{k}.json.gz', 'rt'))[0][0] for k in 'AB'}
    for k, val in (('A', A), ('B', B)):
        s = mpf(src[k].strip('[').split(' ')[0])
        print(par, k, mp.nstr(val, 30), mp.nstr(s, 30), 'diff', mp.nstr(abs(val - s), 3))
        assert abs(val - s) < mpf(10)**-25
print("SPOT CHECK PASS")
