"""Paper IV, Section 8: exact lower frame bound A = 1/delta - lambda_max(K) of delta Z (delta = 1/(1+eps)) in PW_pi with M consecutive
points deleted; -ln A is linear in M. Not part of any proof.
Usage: python numerics-IV/lattice_deletion.py"""
import mpmath as mp
mp.mp.dps = 320
def sinc(x): return mp.mpf(1) if x == 0 else mp.sin(mp.pi * x) / (mp.pi * x)
print(" eps     M   -ln A        -ln A / M   -ln A /(M ln(1/eps)+M)")
for eps in (mp.mpf(1)/2, mp.mpf(1)/4, mp.mpf(1)/8):
    d = 1 / (1 + eps)
    for M in (5, 10, 20, 40, 80):
        K = mp.matrix(M, M)
        for i in range(M):
            for j in range(M): K[i, j] = sinc((i - j) * d)
        lmax = max(mp.eigsy(K)[0])
        A = 1 / d - lmax
        la = -mp.log(A)
        print(f"{float(eps):5.3f} {M:4d} {float(la):10.3f} {float(la/M):10.4f} {float(la/(M*(1+mp.log(1/eps)))):10.4f}", flush=True)
