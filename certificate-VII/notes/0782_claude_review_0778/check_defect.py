# Independent check of 0778 DEFECT_LOWER_BOUND (Claude, 0782). Legendre-Galerkin for D_c, mpmath.
from mpmath import mp, mpf, matrix, eigsy, sqrt, exp, pi, legendre, besseli, diff, log
mp.dps = 80
def ground(c, K):
    A = matrix(K, K)
    for k in range(K):
        A[k, k] = k*(k+1) + c**2*mpf(2*k*k+2*k-1)/((2*k-1)*(2*k+3))
        if k+2 < K:
            v = c**2*mpf((k+1)*(k+2))/((2*k+3)*sqrt(mpf((2*k+1)*(2*k+5))))
            A[k, k+2] = A[k+2, k] = v
    E, Q = eigsy(A)
    i = min(range(K), key=lambda j: E[j])
    b = [Q[j, i] for j in range(K)]
    if b[0] < 0: b = [-x for x in b]
    psi = lambda x: sum(b[k]*sqrt(k+mpf(1)/2)*legendre(k, x) for k in range(K))
    lam = c/(2*pi)*(sqrt(2)*b[0]/psi(0))**2
    return E[i], psi, lam
def d(c): return exp(-2*(c+1))/(4*(c+1)*(c+4))
print(" c     1-lam0          d(c)        ratio   psi(1)^2/[e^-2c/(4c+12)]  chi0<=c^2/3")
for c in [mpf(x) for x in ('0.5', '1', '2', '4', '8', '12', '16', '20')]:
    K = int(2*c) + 50
    chi, psi, lam = ground(c, K)
    _, _, lam2 = ground(c, K+20)
    assert abs(lam-lam2) < mpf(10)**-6*(1-lam), (lam-lam2)
    defect = 1-lam
    r1 = psi(1)**2/(exp(-2*c)/(4*c+12))
    assert defect > d(c) and r1 > 1 and chi <= c**2/3
    print(f"{float(c):4.1f}  {float(defect):.4e}  {float(d(c)):.4e}  {float(defect/d(c)):8.1f}  {float(r1):10.2f}  {bool(chi<=c**2/3)}")
# Slepian identity lam0'(c) = 2 lam0 psi(1)^2 / c at c = 3
mp.dps = 40
c0 = mpf(3)
dl = diff(lambda c: ground(c, 50)[2], c0)
_, psi, lam = ground(c0, 50)
print("Slepian:", float(dl), float(2*lam*psi(1)**2/c0)); assert abs(dl-2*lam*psi(1)**2/c0) < mpf(10)**-20
# Comparison-function identity (pY')' = [c^2x^2 + 3/(4x^2) + 1/4] Y
c = mpf(7)
Y = lambda x: x**mpf(-0.5)*besseli(0, c*sqrt(1-x**2))
for x in [mpf('0.1'), mpf('0.5'), mpf('0.93')]:
    lhs = diff(lambda t: (1-t**2)*diff(Y, t), x)
    rhs = (c**2*x**2 + 3/(4*x**2) + mpf(1)/4)*Y(x)
    assert abs(lhs/rhs-1) < mpf(10)**-15
print("Y identity OK; psi <= psi(1) Y check:")
_, psi, _ = ground(c, 70)
print(all(psi(x) <= psi(1)*Y(x) for x in [mpf(j)/40 for j in range(1, 40)]))
# mu10 corollary: c = 256 * (1/2) log 10, delta_new = 2^-872
c = 256*log(10)/2
print("c=", float(c), " log2 d(c) =", float(log(d(c), 2)), " 2^-872 < d:", mpf(2)**-872 < d(c))
