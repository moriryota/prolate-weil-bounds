"""Spot check of (BL3): I = int_0^inf R(xi) omega(xi) dxi vs r1/c^3 = 1.2/c^3 (non-rigorous, double precision ODE).
R = (x^2/2)(J0^2+J1^2) - x/pi, omega = -(1/q)_xi, q = c^2 (t^2 - s), dt/dxi = sqrt(p/q)."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import j0, j1
for mu, s in [(5, 0.27512), (7, 0.19891), (5, 0.031065)]:
    c = 2 * np.pi * mu
    kap = c * np.sqrt(1 - s)
    x0 = 1e-6
    t0 = 1 + x0**2 / (2 * kap**2)          # xi ~ kappa sqrt(2(t-1))
    def rhs(x, y):
        t = y[0]
        p = t * t - 1; q = c * c * (t * t - s)
        dt = np.sqrt(p / q)
        om = (2 * c * c * t * dt) / q**2
        R = x * x / 2 * (j0(x)**2 + j1(x)**2) - x / np.pi
        return [dt, R * om, om]
    sol = solve_ivp(rhs, [x0, 4000], [t0, 0.0, 0.0], rtol=1e-11, atol=1e-16, max_step=0.05)
    I, Om = sol.y[1, -1], sol.y[2, -1]
    print(f"mu={mu} s={s}: int R omega = {I:.4e};  c^3 * |I| = {abs(I) * c**3:.4f} (proved bound r1 = 1.2);  int omega = {Om:.4e} vs 1/q(0) = {1 / kap**2:.4e}")
