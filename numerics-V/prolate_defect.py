"""Paper V, Section 8: the prolate defect 1 - lambda_4(2 pi mu) from the Legendre-tridiagonal eigenvectors (centre relation).
Usage: python numerics-V/prolate_defect.py mu [mu ...]"""
import sys
sys.path.insert(0, 'common')
import mpmath as mp
from prolate_kvector import pswf_even, legendre_norm_vals
for mu in [int(x) for x in sys.argv[1:]]:
    mp.mp.dps = int(1.3 * 2 * mp.pi * mu) + 120
    c = 2 * mp.pi * mu; K = int(float(c) / 2) + 120
    ks, ps = pswf_even(c, K, which=(0, 2))
    P0 = legendre_norm_vals(mp.mpf(0), ks[-1])
    co = ps[2][1]; psi0 = sum(co[i] * P0[k] for i, k in enumerate(ks)); ipsi = co[0] * mp.sqrt(2)
    print(mu, mp.nstr(1 - c / (2 * mp.pi) * (ipsi / psi0) ** 2, 8), flush=True)
