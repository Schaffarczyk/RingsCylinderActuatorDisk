"""Ring-vortex (Bessel-Laplace) kernels in closed form, and checks vs. quadrature.

I(lam,mu,nu)(R,rho,z) = int_0^inf s^lam J_mu(sR) J_nu(s rho) exp(-s|z|) ds
"""
import numpy as np
from scipy.special import ellipkm1, ellipe, j0, j1
from scipy.integrate import quad


def _ke(R, rho, z):
    D2 = (R + rho) ** 2 + z * z
    d2 = (R - rho) ** 2 + z * z
    m1 = d2 / D2           # k'^2
    m = 1.0 - m1           # k^2
    K = ellipkm1(m1)
    E = ellipe(m)
    return K, E, m, np.sqrt(D2), d2


def I011(R, rho, z):
    """stream-function kernel  (Psi_ring = Gamma*rho*R/2 * I011)"""
    K, E, m, D, d2 = _ke(R, rho, z)
    k = np.sqrt(m)
    return ((2.0 / k - k) * K - 2.0 / k * E) / (np.pi * np.sqrt(R * rho))


def I110(R, rho, z):
    """axial-velocity kernel  (Vz_ring = Gamma*R/2 * I110)"""
    K, E, m, D, d2 = _ke(R, rho, z)
    return (K + (R * R - rho * rho - z * z) / d2 * E) / (np.pi * R * D)


def I111(R, rho, z):
    """radial-velocity kernel (Vr_ring = +sign(z) * Gamma*R/2 * I111)"""
    K, E, m, D, d2 = _ke(R, rho, z)
    return abs(z) / (np.pi * R * rho * D) * (-K + (R * R + rho * rho + z * z) / d2 * E)


# ---- semi-infinite cylinder kernels (rings from z'=0 .. inf, field point at axial distance z>=0 upstream of the start
# i.e. int_0^inf I(lam,1,nu)(R,rho, z + t) dt  = I(lam-1,1,nu)(R,rho,z)
def Im111(R, rho, z):
    return quad(lambda t: I011(R, rho, z + t), 0, np.inf, limit=200)[0]


def I010(R, rho, z):
    return quad(lambda t: I110(R, rho, z + t), 0, np.inf, limit=200)[0]


def I011_tail(R, rho, z):
    return quad(lambda t: I111(R, rho, z + t), 0, np.inf, limit=200)[0]


def bessel_quad(lam, mu, nu, R, rho, z):
    J = {0: j0, 1: j1}
    f = lambda s: s ** lam * J[mu](s * R) * J[nu](s * rho) * np.exp(-s * abs(z))
    return quad(f, 0, np.inf, limit=500)[0]


if __name__ == "__main__":
    for (R, rho, z) in [(1.0, 0.3, 0.5), (1.0, 1.3, 0.2), (0.8, 0.79, 0.05), (1.2, 0.2, 2.0)]:
        print(R, rho, z)
        print("  I011", I011(R, rho, z), bessel_quad(0, 1, 1, R, rho, z))
        print("  I110", I110(R, rho, z), bessel_quad(1, 1, 0, R, rho, z))
        print("  I111", I111(R, rho, z), bessel_quad(1, 1, 1, R, rho, z))
        print("  I-111", Im111(R, rho, z), bessel_quad(-1, 1, 1, R, rho, z))
        print("  I010", I010(R, rho, z), bessel_quad(0, 1, 0, R, rho, z))
        print("  I011t", I011_tail(R, rho, z), bessel_quad(0, 1, 1, R, rho, z))
