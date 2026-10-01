"""Bessel-type Liouville transformation of the radial prolate equation (p psi')' + q psi = 0, t > 1,
p = t^2 - 1, q = c^2 (t^2 - s).  xi(t) = int_1^t sqrt(q/p),  P = sqrt(p q),  psi = (xi/P)^{1/2} W  gives
(xi W')' + xi W = xi Delta W,   Delta = (P^{1/2})''/P^{1/2} + 1/(4 xi^2)   (derivatives in xi).
Report c^2 sup|Delta|, c * int |Delta| dxi, and the Olver-type control integral (pi/2) int xi |Delta| (J0^2+Y0^2) dxi (xi >= 1 part)."""
import mpmath as mp
mp.mp.dps = 30
def run(c, s):
    def xi(t):
        V = mp.sqrt(t - 1)
        return mp.quad(lambda v: 2 * c * mp.sqrt((1 + v**2)**2 - s) / mp.sqrt(2 + v**2), [0, V])
    def Delta(t):
        p = t**2 - 1; r = t**2 - s
        z = c * mp.sqrt(r / p)
        lt = (2 * t / p + 2 * t / r) / 4
        ltt = (2 / p - 4 * t**2 / p**2 + 2 / r - 4 * t**2 / r**2) / 4
        zt_over_z = (2 * t / r - 2 * t / p) / 2
        lxixi = (ltt - lt * zt_over_z) / z**2
        lxi = lt / z
        X = xi(t)
        return lxixi + lxi**2 + 1 / (4 * X**2), X
    ts = [1 + mp.mpf(10)**(-k) for k in range(10, 1, -1)] + [1 + mp.mpf(j) / 20 for j in range(1, 21)] + [2 + mp.mpf(j) / 2 for j in range(0, 60)] + [40, 100, 400]
    vals = [(t,) + Delta(t) for t in ts]
    supc2 = max(abs(d) * c**2 for _, d, _ in vals)
    vals.sort(key=lambda v: v[2])
    I = sum((abs(vals[i][1]) + abs(vals[i + 1][1])) / 2 * (vals[i + 1][2] - vals[i][2]) for i in range(len(vals) - 1))
    tail = abs(vals[-1][1]) * vals[-1][2]
    return supc2, c * (I + tail), [(mp.nstr(v[2], 4), mp.nstr(v[1] * c**2, 4)) for v in vals[::10]]
for mu, s in [(5, 0.27512), (7, 0.19891), (5, 0.031065)]:
    c = 2 * mp.pi * mu
    sc, ic, sample = run(c, mp.mpf(s))
    print(f"mu={mu} s={s}: c^2 sup|Delta| = {mp.nstr(sc, 5)};  c * int|Delta| dxi ~ {mp.nstr(ic, 5)};  sample (xi, c^2 Delta): {sample}", flush=True)
