# 0782: spot checks of GENERAL_TAIL.md (Claude). e_j = out-of-band mass of v_0, v_1; c = Omega*L.
# Closed forms: int_c^inf sin^2u/u^2 = sin^2c/c + pi/2 - Si(2c);  j_1^2 integral by quadrature with Si-based check.
from mpmath import mp, mpf, sin, cos, pi, si, quad, inf, psi, mpc, re
mp.dps = 40
def e0(c): return 2/pi*(sin(c)**2/c + pi/2 - si(2*c))
def e1(c):
    f = lambda u: (sin(u)/u**2 - cos(u)/u)**2
    # split: oscillatory part integrated over periods up to c+20000*pi, asymptotic tail 1/(2u^2)-avg added analytically
    T = c + 20000*pi
    body = quad(f, [c + k*pi for k in range(0, 20001)])
    tail = 1/(2*T)  # average of cos^2/u^2 beyond T, error O(T^-2)
    return 3*2/pi*(body + tail)
for c in [mpf(10), mpf('294.73'), mpf('636.14'), mpf(1500)]:
    a, b = e0(c), e1(c)
    lo0, hi0 = (c-1)/(pi*c**2), (c+1)/(pi*c**2)
    lo1, hi1 = 3*(c-2)/(pi*c**2), (3*c**2+6*c+2)/(pi*c**3)
    print(float(c), bool(lo0 <= a <= hi0), bool(lo1 <= b <= hi1), float(a*pi*c), float(b*pi*c/3))
vals = [re(psi(0, mpc(0.25, t/2))) for t in [0, 0.5, 1, 2, 5, 10, 50, 100, 500]]
print("a(t) increasing:", all(x < y for x, y in zip(vals, vals[1:])))
