"""Paper V, Section 8: the same computation with more quadrature nodes (n_tri = N + extra_tri, n_ext = N + extra_ext)
and, typically, more bits, to check the precision of galerkin_lambda_min.py.
Usage: python numerics-V/precision_check.py mu N prec nq extra_tri extra_ext"""
import sys, time
sys.path.insert(0, 'common'); sys.path.insert(0, 'numerics-V')
import mpmath as mp
import prolate_kvector
from weil_lowmem import WeilSpectralLowMem
mu, N, prec, nq, xt, xe = (int(x) for x in sys.argv[1:7])
class W2(WeilSpectralLowMem):
    def __init__(self, L, N, prec=256, **kw):
        super().__init__(L, N, prec=prec, n_tri=N + xt, **kw)
        self.n_ext = N + xe
prolate_kvector.WeilSpectral = W2
t0 = time.time()
r = prolate_kvector.run(mu, N=N, prec=prec, nq=nq)
print(f"mu={mu} N={N} prec={prec} n_tri=N+{xt} n_ext=N+{xe}  W={mp.nstr(r['rayleigh'], 10)}  lambda_min={mp.nstr(r['lam_min'], 10)}  min={(time.time()-t0)/60:.1f}", flush=True)
