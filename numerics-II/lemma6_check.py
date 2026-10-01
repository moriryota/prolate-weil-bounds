"""Theorem C, numerical check of Lemma 6 and S1a' (not part of the proof):
  E(t)/A0^2 <= 3.12 on [t0, 60], t0 = 1 + 1/(2c^2)     (E = w^2 + w_z^2, w = (pq)^{1/4} psi~)
  psi~(t)^2 c sqrt(t^2-1)/A0^2 <= 3.79 on (1, 2]
  |rho(t)| <= min(2, 2c/(3t)) + 0.4/t^2 on [2, 300]  (A_inf, Phi fitted on [400, 405])"""
import sys
sys.path.insert(0, 'common')
import mpmath as mp
from prolate_kvector import pswf_even, legendre_norm_vals
def sph_j(k, w): return mp.sqrt(mp.pi / (2 * w)) * mp.besselj(k + mp.mpf(1) / 2, w)
mu = int(sys.argv[1]) if len(sys.argv) > 1 else 5
mp.mp.dps = int(2.2 * mu) + 45
c = 2 * mp.pi * mu
K = int(float(c) / 2) + 70
ks, ps = pswf_even(c, K)
for m in (0, 2):
    coef = ps[m][1]
    def val(s):
        P = legendre_norm_vals(s, ks[-1]); return sum(coef[i] * P[k] for i, k in enumerate(ks))
    F = lambda t: sum(coef[i] * mp.sqrt((2 * k + 1) / mp.mpf(2)) * 2 * (1j)**k * sph_j(k, c * t) for i, k in enumerate(ks))
    s0 = mp.mpf('0.3'); mun = F(s0) / val(s0)
    psit = lambda t: (F(t) / mun).real
    A0 = abs(val(mp.mpf(1)))
    chi = None
    # chi^SL from the ODE at a sample point: (p psi')' + q psi = 0  ->  chi = c^2 t^2 + ((p psi')')/psi
    tt = mp.mpf('1.7')
    d1 = lambda t: mp.diff(psit, t); pp = lambda t: (t**2 - 1) * d1(t)
    chi = c**2 * tt**2 + mp.diff(pp, tt) / psit(tt)
    s = chi / c**2
    def E(t):
        p = t**2 - 1; q = c**2 * t**2 - chi
        ps_, dps_ = psit(t), d1(t)
        w = (p * q) ** mp.mpf(0.25) * ps_
        dw = (p * q) ** mp.mpf(0.25) * dps_ + mp.mpf(0.25) * (p * q) ** mp.mpf(-0.75) * (2 * t * q + p * 2 * c**2 * t) * ps_
        wz = dw * mp.sqrt(p / q)
        return w**2 + wz**2
    t0 = 1 + 1 / (2 * c**2)
    grid = [t0 * (1 + mp.mpf(j) / 4000) for j in range(0, 40)] + [1 + mp.mpf(j) / 50 for j in range(1, 50)] + [mp.mpf(2) + mp.mpf(j) / 3 for j in range(0, 175)]
    Emax = max(E(t) / A0**2 for t in grid if t >= t0)
    near = max(psit(t)**2 * c * mp.sqrt(t**2 - 1) / A0**2 for t in [1 + mp.mpf(10)**(-j) for j in range(2, 8)] + [1 + mp.mpf(j) / 200 for j in range(1, 201)])
    tf = [400 + mp.mpf(j) / mp.mpf("7.3") for j in range(36)]   # spacing not commensurate with 2pi/c
    M = mp.matrix([[mp.sin(c * t), mp.cos(c * t)] for t in tf]); y = mp.matrix([psit(t) * mp.sqrt(c) * t for t in tf])
    sol = mp.lu_solve(M.T * M, M.T * y); a_, b_ = sol[0], sol[1]
    A = mp.sqrt(a_**2 + b_**2); Phi = mp.atan2(b_, a_)
    worst = 0
    for t in [2 + mp.mpf(j) / 9 for j in range(0, 2683)]:
        rho = abs(psit(t) * mp.sqrt(c) * t / A - mp.sin(c * t + Phi))
        bd = min(2, 2 * c / (3 * t)) + mp.mpf('0.4') / t**2
        worst = max(worst, rho / bd)
    print(f"mu={mu} n={2*m}: s=chi/c^2={mp.nstr(s, 5)}  E(t0)/A0^2={mp.nstr(E(t0) / A0**2, 5)} (bound 2.1984)  max_[t0,60] E/A0^2={mp.nstr(Emax, 5)} (bound 3.12)  "
          f"A_inf/A0={mp.nstr(A / A0, 5)} (bound 1.767)  near-field max={mp.nstr(near, 5)} (bound 3.79)  max |rho|/bound on [2,300]={mp.nstr(worst, 4)}", flush=True)
