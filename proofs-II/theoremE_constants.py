"""Theorem E: for all mu >= MU0, lambda_min <= [P(mu) C4 (2 pi mu)^{9/2} e^{-4 pi mu} + e^{-16 mu}] / KAPPA0
 <= K mu^8 (log mu)^3 e^{-4 pi mu}.  KAPPA0 = F(MU0, theta) from stage3_constants (F increasing in mu for fixed theta)."""
import io, contextlib, importlib.util, math
from flint import arb
def load(path):
    spec = importlib.util.spec_from_file_location(path, path); m = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(m)
    return m
s3 = load('proofs-II/L6_effective_constants.py'); tc = load('proofs-II/theoremC_constants.py')
PI = arb.pi()
C4 = arb('6620.3392')                 # paper I, Theorem 1 (n = 4), with the certified starting value U_4 <= 2.5967761e-18 at c0 = 10 pi
for MU0, TH in [(50, '0.145'), (80, '0.095')]:
    k, _ = s3.F(MU0, arb(TH))
    KAPPA0 = arb(str(math.floor(float(k.lower()) * 1000) / 1000)); assert KAPPA0 <= k
    P, *_ = tc.P_of(MU0)
    mu = arb(MU0)
    ratio = P / (mu**arb('3.5') * mu.log()**3)          # decreasing in mu (each factor monotone, see note)
    extra = (-16 * mu).exp() / (mu**8 * mu.log()**3 * (-4 * PI * mu).exp())   # decreasing; tiny
    K = (ratio * C4 * (2 * PI)**(arb(9) / 2) + extra) / KAPPA0
    print(f"MU0={MU0} (theta={TH}): kappa_lambda >= {KAPPA0} for all mu >= {MU0};  K <= {float(K.upper()):.4e}  [P/(mu^3.5 log^3) at MU0 = {float(ratio.upper()):.4e}]")
    for m in [MU0, 100, 1000]:
        mm = arb(m); Pm, *_ = tc.P_of(m)
        thmE = (Pm * C4 * (2 * PI * mm)**(arb(9) / 2) + (-16 * mm + 4 * PI * mm).exp()) / KAPPA0
        thmB = arb('1.362e9') * mm**23
        print(f"   mu={m}: Thm E coefficient (times e^-4pi mu) <= {float(thmE.upper()):.3e}  vs Thm B {float(thmB.mid()):.3e}")
