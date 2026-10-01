"""Certify χ_4^{SL}(c) < c²/2 for c = 10π (μ = 5); by monotonicity of χ_n(c)/c² in c (Hellmann–Feynman:
d/dc(χ/c²) = (2/c)(⟨s²ψ,ψ⟩ − χ/c²) ≤ 0 since χ = ⟨(1−s²)ψ'²⟩ + c²⟨s²ψ²⟩ ≥ c²⟨s²ψ²⟩) this gives χ_0, χ_4 < c²/2 for all c ≥ 10π.
Rayleigh–Ritz: the j-th eigenvalue of the truncated matrix of L_c on even normalised Legendre polynomials P̄_0..P̄_{2K-2}
is ≥ the j-th even eigenvalue of L_c (χ_0, χ_2, χ_4 for j = 1, 2, 3). We count, in Arb, the eigenvalues of the truncated
matrix below x = c²/2 with the Sturm sequence of the symmetric tridiagonal matrix (inertia via LDLᵀ); if ≥ 3, then χ_4 < c²/2."""
from flint import arb, ctx
ctx.prec = 200
def a(k):  # z P̄_k = a(k+1) P̄_{k+1} + a(k) P̄_{k-1}
    return arb(k) / ((arb(2 * k - 1) * arb(2 * k + 1)).sqrt()) if k > 0 else arb(0)
def count_below(c, x, K):
    ks = [2 * i for i in range(K)]
    diag = [arb(k * (k + 1)) + c**2 * (a(k + 1)**2 + a(k)**2) for k in ks]
    off = [c**2 * a(k + 1) * a(k + 2) for k in ks[:-1]]
    # Sturm: d_1 = diag_1 - x ; d_i = (diag_i - x) - off_{i-1}^2 / d_{i-1}; # negative d_i = # eigenvalues < x
    cnt = 0; d = None
    for i in range(K):
        di = diag[i] - x - (off[i - 1]**2 / d if i > 0 else 0)
        assert not di.contains(0), ("sign undetermined", i)
        if di < 0: cnt += 1
        d = di
    return cnt
c = 10 * arb.pi()
for K in (12, 20, 30):
    print(f"K={K}: # eigenvalues of truncated even matrix below c^2/2 at c=10pi: {count_below(c, c**2 / 2, K)};  below c^2/4: {count_below(c, c**2 / 4, K)};  below 0.3c^2: {count_below(c, 3 * c**2 / 10, K)}")
