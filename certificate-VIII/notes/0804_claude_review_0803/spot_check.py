"""0804 (Claude): independent mpmath spot checks of the mu=15 source matrices (0792 runs, reused by 0803).

Not derived from generate_lowmem.py / weil_spectral.py: direct integration of the definitions.
Index convention: parity p, row k <-> Legendre degree j = 2k+p, v_j = sqrt((2j+1)/(2L)) P_j(x/L).
  A = <v, M v>, B = <M v, M v>, M = gamma - C_p, C_p f(x) = sum Lambda(n)/sqrt(n) [f(x-log n)+f(x+log n)] (zero-extended)
  J = (1/pi) int_0^Omega (a(t)-gamma) b_i b_j dt + 2<v_i,cosh(x/2)><v_j,cosh(x/2)> - 2<v_i,sinh(x/2)><v_j,sinh(x/2)>
      + u_p h_i h_j,  h_i = <v_i, (I - B_band - delta_old) v_{p}>,  delta_old = d(c), (B_band)_{ij} = (1/pi) int_0^Omega b_i b_j
  a(t) = Re psi(1/4+it/2) - log pi, b_j(t) = |F_j(t)| up to the parity phase: F_j(t) = 2L i^j sqrt((2j+1)/(2L)) j_j(Lt).
"""
import gzip, json, sys
from mpmath import mp, mpf, log, sqrt, quad, pi, psi, mpc, re, exp, cosh, sinh, legendre, besselj, linspace
mp.dps = 40
ROOT = sys.argv[1]
L = log(15)/2; g = mpf(37)/8; Om = 736; c = Om*L; u = {0: 428, 1: 142}
pp = {2: 2, 3: 3, 4: 2, 5: 5, 7: 7, 8: 2, 9: 3, 11: 11, 13: 13}
def v(j): return lambda x: sqrt(mpf(2*j+1)/(2*L))*legendre(j, x/L) if -L <= x <= L else mpf(0)
def Cp(f, x): return sum(log(p)/sqrt(n)*(f(x-log(n)) + f(x+log(n))) for n, p in pp.items())
brk = sorted({-L, L, *[s*(L-log(n)) for n in pp for s in (1, -1) if log(n) < 2*L]})
def src(p, name, k, l):
    path = f"{ROOT}/notes/0792_claude_lowmem_generator/runs/m15n1664p3328q224_p{p}_{name}.json.gz"
    rows = json.load(gzip.open(path, 'rt'))
    s = rows[k][l]; return mpf(s.strip('[').split(' ')[0]), s
NEED = {(p, 'A'): [(0, 0), (3, 5), (40, 41)] for p in (0, 1)}
NEED.update({(p, 'B'): NEED[(p, 'A')] for p in (0, 1)})
NEED.update({(0, 'J'): [(0, 0), (0, 2), (1, 2)], (1, 'J'): [(0, 0), (1, 2), (2, 2)]})
VALS = {}
import gc
for (p, name), idx in NEED.items():
    if name != 'J' and '--all' not in sys.argv: continue   # load one matrix at a time, keep only the needed entries
    rows = json.load(gzip.open(f"{ROOT}/notes/0792_claude_lowmem_generator/runs/m15n1664p3328q224_p{p}_{name}.json.gz", 'rt'))
    assert len(rows) == 832 and all(len(r) == 832 for r in rows)
    for (k, l) in idx: VALS[(p, name, k, l)] = rows[k][l]
    del rows; gc.collect()
    print('loaded', p, name, flush=True)
def sval(p, name, k, l):
    s = VALS[(p, name, k, l)]; return mpf(s.strip('[').split(' ')[0])
worst = mpf(0); out = []
def report(tag, mine, theirs, tol):
    global worst
    diff = abs(mine - theirs); worst = max(worst, diff)
    line = f"{tag}: indep {mp.nstr(mine, 25)}  source {mp.nstr(theirs, 25)}  diff {mp.nstr(diff, 3)}"
    print(line, flush=True); out.append(line)
    assert diff < tol, line
# ---- A, B
for p in ((0, 1) if '--all' in sys.argv else ()):
    for (k, l) in [(0, 0), (3, 5), (40, 41)]:
        fi, fj = v(2*k+p), v(2*l+p)
        Mi = lambda x: g*fi(x) - Cp(fi, x); Mj = lambda x: g*fj(x) - Cp(fj, x)
        # subdivide finely for high degree
        pts = sorted(set(brk) | set(linspace(-L, L, 41)))
        A = quad(lambda x: fi(x)*Mj(x), pts)
        B = quad(lambda x: Mi(x)*Mj(x), pts)
        report(f"p{p} A[{k}][{l}]", A, sval(p, 'A', k, l), mpf(10)**-25)
        report(f"p{p} B[{k}][{l}]", B, sval(p, 'B', k, l), mpf(10)**-25)
# ---- J
def bj(j):
    # b_j = i^{-p} F_j = (-1)^k 2L nrm j_j(Lt), j = 2k+p  (run 1 omitted the phase (-1)^k; see spot_check_run1.log)
    nrm = sqrt(mpf(2*j+1)/(2*L))
    def f(t):
        if t == 0: return 2*L*nrm if j == 0 else mpf(0)
        return (-1)**(j//2)*2*L*nrm*sqrt(pi/(2*L*t))*besselj(j+mpf(1)/2, L*t)
    return f
a = lambda t: re(psi(0, mpc(mpf(1)/4, t/2))) - log(pi)
tp = linspace(0, Om, 1473)
d_old = exp(-2*(c+1))/(4*(c+1)*(c+4))
def band(i, j): return quad(lambda t: bj(i)(t)*bj(j)(t), tp)/pi
def ip(f, w): return quad(lambda x: f(x)*w(x), [-L, 0, L])
for p, pairs in ((0, [(0, 0), (0, 2), (1, 2)]), (1, [(0, 0), (1, 2), (2, 2)])):
    j0 = p
    for (k, l) in pairs:
        i, j = 2*k+p, 2*l+p
        K = quad(lambda t: (a(t)-g)*bj(i)(t)*bj(j)(t), tp)/pi
        pole = 2*ip(v(i), lambda x: cosh(x/2))*ip(v(j), lambda x: cosh(x/2)) \
             - 2*ip(v(i), lambda x: sinh(x/2))*ip(v(j), lambda x: sinh(x/2))
        hi = (1 - d_old if i == j0 else 0) - band(i, j0)
        hj = (1 - d_old if j == j0 else 0) - band(j, j0)
        J = K + pole + u[p]*hi*hj
        report(f"p{p} J[{k}][{l}]", J, sval(p, 'J', k, l), mpf(10)**-20)
print("max diff", mp.nstr(worst, 3)); print("SPOT CHECK PASS")
open(sys.argv[2], 'w').write("\n".join(out) + f"\nmax diff {mp.nstr(worst, 3)}\nSPOT CHECK PASS\n")
