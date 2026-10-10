"""Galerkin values of W(k_lambda)/||k_lambda||^2 used in the numerical table of paper II (not part of any proof).
The Weil form is projected onto N normalized Legendre polynomials on the window, in Arb ball arithmetic
(common/weil_spectral.py); k_lambda is projected by piecewise Gauss quadrature with nq nodes per piece
(common/prolate_kvector.py). The quadrature must be fine enough (nq = N + 40 for mu >= 13). With the parameters below the values of W(k) are
not fully converged: they may differ from the limits by a few percent (at mu = 5, W(k) increases by 3.5% from N = 140
to N = 320), and lambda_min(even) is not converged. Smaller N (e.g. N = 140 at mu = 11) or lower precision is worse.
Runtimes on a laptop: about 10 min (mu <= 11), 9-12 min (mu = 13, 15), 47 min (mu = 17).
Usage: python numerics-II/galerkin_W.py [mu ...]"""
import sys
sys.path.insert(0, 'common')
import mpmath as mp
from prolate_kvector import run
PARAMS = {5: dict(N=140), 7: dict(N=140), 11: dict(N=180, nq=200), 13: dict(N=260, prec=480, nq=300),
          15: dict(N=260, prec=480, nq=300), 17: dict(N=300, prec=560, nq=340)}
for mu in [int(a) for a in sys.argv[1:]] or sorted(PARAMS):
    r = run(mu, **PARAMS[mu])
    print(f"mu={mu:>3} {PARAMS[mu]}  W(k)/||k||^2={mp.nstr(r['rayleigh'], 6)}  lambda_min(even)={mp.nstr(r['lam_min'], 6)}", flush=True)
