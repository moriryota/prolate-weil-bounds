"""Paper III, Section 8 (first table): out-of-band mass of f = psi_n * (1 - exp(-c^p (1-x^2)/beta))^s on [-1,1], bandwidth c,
computed as ||f||^2 minus the in-band mass (60 digits). p = 2 is the layer of the paper; p = 1 shows why width 1/c fails.
Not part of any proof.  Usage: python numerics-III/smoothed_oob.py c_over_pi beta 0 p"""
import sys
sys.path.insert(0, 'common')
import mpmath as mp
from prolate_kvector import pswf_even, legendre_norm_vals
cp = int(sys.argv[1]); beta = mp.mpf(sys.argv[2]); Wf = int(sys.argv[3])
mp.mp.dps = 60
c = cp * mp.pi
K = int(float(c) / 2) + 40
ks, ps = pswf_even(c, K, which=(0, 2))
L = 2 * K + 80                      # Legendre degree for the smoothed function
xg, wg = mp.gauss_quadrature(2 * L, 'legendre')
Pg = [legendre_norm_vals(x, L) for x in xg]
def psi(m, Pv): return sum(ps[m][1][i] * Pv[k] for i, k in enumerate(ks))
X, Wq = mp.gauss_quadrature(24, 'legendre')
def oob(co):
    """B = ||f||^2 - (1/2pi) int_{|w|<c} |f^|^2 (no tail truncation; dps 60 absorbs the cancellation)."""
    kk = [k for k in range(len(co)) if k % 2 == 0 and abs(co[k]) > mp.mpf(10)**-58]
    def fh(w):
        return sum(co[k] * mp.sqrt((2 * k + 1) / mp.mpf(2)) * 2 * (-1)**(k // 2) * mp.sqrt(mp.pi / (2 * w)) * mp.besselj(k + mp.mpf(1) / 2, w) for k in kk)
    tot = mp.mpf(0); nseg = int(mp.ceil(c / (mp.pi / 4)))
    for j in range(nseg):
        lo, hi = c * j / nseg, c * (j + 1) / nseg
        tot += sum(wq * fh((hi + lo) / 2 + (hi - lo) / 2 * x) ** 2 for x, wq in zip(X, Wq)) * (hi - lo) / 2
    return sum(x**2 for x in co) - 2 * tot / (2 * mp.pi), max(abs(co[k]) for k in range(L - 20, L + 1))
for n, m in ((0, 0), (4, 2)):
    base = None
    for s in (0, 1, 2):
        vals = [psi(m, Pv) * (1 - mp.e**(-c**int(sys.argv[4]) * (1 - x**2) / beta))**s for x, Pv in zip(xg, Pg)]
        co = [sum(w * v * Pv[k] for w, v, Pv in zip(wg, vals, Pg)) for k in range(L + 1)]
        nrm = sum(x**2 for x in co)
        B, last = oob(co)
        r = B / nrm
        if s == 0: base = r
        print(f"c={cp}pi p={sys.argv[4]} beta={beta} n={n} s={s}: B/||f||^2 = {mp.nstr(r, 6)}  ratio to s=0: {mp.nstr(r / base, 6)}  max|coef| deg>L-20: {mp.nstr(last, 3)}", flush=True)
