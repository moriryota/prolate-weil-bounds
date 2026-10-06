"""Paper V numerics: memory-light variant of WeilSpectral.hardy_interior. The original stores all (N+40)^2 x N rows and forms
two (N+40)^2 x N arb matrices (OOM at N=390, prec=720 on 32 GB). Here the Gram sum is accumulated over chunks of
n_tri rows (one outer Gauss node at a time). Same quadrature and same formula; only the summation order differs."""
import sys
sys.path.insert(0, 'common')
from flint import arb, arb_mat
import weil_spectral
from weil_spectral import WeilSpectral, gauss_legendre, legendre_dd, K_of, gram


class WeilSpectralLowMem(WeilSpectral):
    def hardy_interior(self):
        L, N = self.L, self.N
        xs, ws = gauss_legendre(self.n_tri)
        acc = arb_mat(N, N)
        for a in range(self.n_tri):
            s = (xs[a] + 1) / 2
            ws_ = ws[a] / 2
            rows, wts = [], []
            for b in range(self.n_tri):
                r = (xs[b] + 1) / 2
                wr = ws[b] / 2
                x = -L + 2 * L * s
                y = -L + 2 * L * s * r
                u = x - y
                k1 = u * K_of(u)
                weight = ws_ * wr * (2 * L) ** 2 * s * k1 * u / 2
                D = legendre_dd(x / L, y / L, N)
                rows.append([D[j] * self.norm[j] / L for j in range(N)])
                wts.append(weight)
            acc += gram(rows, wts, N)
        return acc
