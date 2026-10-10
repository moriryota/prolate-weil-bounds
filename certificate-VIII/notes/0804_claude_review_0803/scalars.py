# 0804 (Claude): independent scalar gates (mpmath, 60 digits) and exact alpha+beta*b^2 (Fraction).
from fractions import Fraction as F
from mpmath import mp, mpf, log, exp, pi, log as ln
mp.dps = 60
L = log(15)/2; c = 736*L; g = mpf(37)/8
d = exp(-2*(c+1))/(4*(c+1)*(c+4)); C = log(736/(2*pi)) - mpf(1)/736 - g
print('log2 d(c)', mp.nstr(ln(d, 2), 12), ' C', mp.nstr(C, 12))
assert mpf(2)**-2901 <= d < mpf(2)**-2900 and C > mpf(1)/8 and mpf(2)**-2904 < C*mpf(2)**-2901
E0 = (c+1)/(pi*c*c); E1 = (3*c*c+6*c+2)/(pi*c**3)
print('u0*E0, u1*E1 vs C', mp.nstr(428*E0, 8), mp.nstr(142*E1, 8)); assert 428*E0 < C and 142*E1 < C
print('u from floor(C/E):', int(C/E0), int(C/E1))
m = F(1, 2); e = F(1, 2**706); n = F(1, 2**1424); h = F(1, 2**1419)
alpha = e + 2*e*e/(m-n); beta = (e+n)/(m*m) + 2*h*h/(m*m*(m-n)); tot = alpha + beta*F(35, 4)**2
import math
print('log2(alpha+beta b^2) ~', math.log2(tot.numerator) - math.log2(tot.denominator)); assert tot < F(1, 2**697)
print('2^-698 bound holds?', tot < F(1, 2**698))
print('perturb 570(2d+d^2) log2', mp.nstr(ln(570*(2*d+d*d), 2), 10)); assert 570*(2*d+d*d) < mpf(2)**-2889
print('SCALARS PASS')
