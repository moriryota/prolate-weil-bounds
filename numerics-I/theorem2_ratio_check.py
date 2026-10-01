"""Numerical check of Theorem 2 (paper I) for all n <= 4 (incl. odd n), mu = 5, 7:
  T = (1-lambda_n)/lambda_n  vs  RHS = (psi_n(1)^2/c)(2/pi)K(s)(1 + 15/c^2),  s = chi_n/c^2.
PSWFs from the Legendre tridiagonal matrix of L_c in the parity sector; lambda from the finite Fourier relation
 even n: mu_n = int psi / psi(0);  odd n: mu_n = i c int y psi / psi'(0);  lambda = (c/2pi)|mu_n|^2."""
import mpmath as mp
mp.mp.dps = 200
def Pbar(k, x): return mp.sqrt((2 * k + 1) / mp.mpf(2)) * mp.legendre(k, x)
def pswf(c, n, K):
    par = n % 2
    ks = [2 * i + par for i in range(K)]
    a = lambda k: mp.mpf(k) / mp.sqrt((2 * k - 1) * (2 * k + 1)) if k > 0 else mp.mpf(0)
    A = mp.matrix(K, K)
    for i, k in enumerate(ks):
        A[i, i] = k * (k + 1) + c**2 * (a(k + 1)**2 + a(k)**2)
        if i + 1 < K: A[i, i + 1] = A[i + 1, i] = c**2 * a(k + 1) * a(k + 2)
    E, V = mp.eigsy(A)
    order = sorted(range(K), key=lambda i: E[i])
    j = order[n // 2]
    return E[j], ks, [V[i, j] for i in range(K)]
for mu in [5, 7]:
    c = 2 * mp.pi * mu
    K = int(float(c) / 2) + 60
    for n in range(5):
        chi, ks, co = pswf(c, n, K)
        psi1 = sum(co[i] * Pbar(k, 1) for i, k in enumerate(ks))
        if n % 2 == 0:
            ipsi = co[0] * mp.sqrt(2) if ks[0] == 0 else 0      # int Pbar_0 = sqrt 2
            psi0 = sum(co[i] * Pbar(k, 0) for i, k in enumerate(ks))
            m_n = ipsi / psi0
        else:
            # int y psi = coefficient of Pbar_1 times int y Pbar_1 = sqrt(2/3); psi'(0) = sum co Pbar_k'(0)
            iy = co[0] * mp.sqrt(mp.mpf(2) / 3)
            dpsi0 = sum(co[i] * mp.sqrt((2 * k + 1) / mp.mpf(2)) * mp.diff(lambda x: mp.legendre(k, x), 0) for i, k in enumerate(ks))
            m_n = c * iy / dpsi0
        lam = c / (2 * mp.pi) * m_n**2
        T = (1 - lam) / lam
        s = chi / c**2
        rhs = psi1**2 / c * (2 / mp.pi) * mp.ellipk(s) * (1 + 15 / c**2)
        print(f"mu={mu} n={n}: 1-lambda={mp.nstr(1 - lam, 6)}  s={mp.nstr(s, 5)}  T/RHS={mp.nstr(T / rhs, 8)}  T/main={mp.nstr(T / (rhs / (1 + 15 / c**2)), 8)}", flush=True)
