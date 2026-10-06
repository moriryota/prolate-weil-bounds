"""Paper V, Section 8 (numerical, not part of any proof): Galerkin bottom of the even part of W_lambda on N Legendre
modes, and W(k)/||k||^2 for the prolate vector of paper II. Uses common/prolate_kvector.py with the memory-light
Gram accumulation of numerics-V/weil_lowmem.py (same quadrature; only the summation order differs).
Usage (from the repository root): python numerics-V/galerkin_lambda_min.py mu N prec nq
Convergence needs N/(a c) >~ 2.3 (a = log(mu)/2, c = 2 pi mu), nq = N + 40, and a precision that grows with N
(mu = 20, N = 570 needs more than 840 bits; see precision_check.py). Runtimes: from minutes (mu = 7) to many hours."""
import sys, time
sys.path.insert(0, 'common'); sys.path.insert(0, 'numerics-V')
import mpmath as mp
import prolate_kvector
from weil_lowmem import WeilSpectralLowMem
prolate_kvector.WeilSpectral = WeilSpectralLowMem
from prolate_kvector import run
mu, N, prec, nq = (int(x) for x in sys.argv[1:5])
t0 = time.time()
r = run(mu, N=N, prec=prec, nq=nq)
print(f"mu={mu} N={N} prec={prec} nq={nq}  W(k)/||k||^2={mp.nstr(r['rayleigh'], 8)}  lambda_min(even)={mp.nstr(r['lam_min'], 8)}  "
      f"proj_loss={mp.nstr(r['proj_loss'], 3)}  minutes={(time.time()-t0)/60:.1f}", flush=True)
