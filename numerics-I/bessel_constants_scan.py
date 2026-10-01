"""Bessel constants (numerical scan; the rigorous bounds are in proofs-I/):
  G(x) = (x^2/2)(J0^2+J1^2) = int_0^x eta J0^2,   C1p = sup_x (G(x) - x/pi)
  H(x) = (x^2/2)(J0Y0+J1Y1) = int_0^x eta J0 Y0,  h0 = sup_x |H(x)|
  (pi/2) x (J0^2+Y0^2) <= 1 (Nicholson; DLMF 10.18 monotonicity)"""
import mpmath as mp
mp.mp.dps = 25
G = lambda x: x**2 / 2 * (mp.besselj(0, x)**2 + mp.besselj(1, x)**2)
H = lambda x: x**2 / 2 * (mp.besselj(0, x) * mp.bessely(0, x) + mp.besselj(1, x) * mp.bessely(1, x))
Nm = lambda x: mp.pi / 2 * x * (mp.besselj(0, x)**2 + mp.bessely(0, x)**2)
xs = [mp.mpf(j) / 50 for j in range(1, 50 * 60)] + [60 + mp.mpf(j) / 5 for j in range(0, 2000)]
c1 = max((G(x) - x / mp.pi, x) for x in xs)
h0 = max((abs(H(x)), x) for x in xs)
nm = max((Nm(x), x) for x in xs)
print("C1p = sup(G - x/pi) ~", mp.nstr(c1[0], 6), "at x =", mp.nstr(c1[1], 5))
print("h0 = sup|H| ~", mp.nstr(h0[0], 6), "at x =", mp.nstr(h0[1], 5))
print("sup (pi/2) x (J0^2+Y0^2) ~", mp.nstr(nm[0], 8), "at x =", mp.nstr(nm[1], 5))
print("check identity G' = x J0^2 at x=3.7:", mp.nstr(mp.diff(G, 3.7) - 3.7 * mp.besselj(0, 3.7)**2, 3), "; H' = x J0 Y0:", mp.nstr(mp.diff(H, 3.7) - 3.7 * mp.besselj(0, 3.7) * mp.bessely(0, 3.7), 3))
