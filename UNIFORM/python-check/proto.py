"""Prototype: uniformly loaded actuator disk (h - h0 = b inside slipstream) in Conway's framework.
Vortex sheet on the slipstream boundary r = R(z), z>=0, ring strength per unit axial length
    gamma_z(z) = b / Vzm(z),  Vzm = mean axial velocity on the sheet.
Nondimensional: U_inf = 1, R_disk = 1.
"""
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from kernels import I011, I110, I111, Im111, I010, I011_tail

X2 = 40.0


class Sheet:
    def __init__(self, b, N=40, z1=0.01):
        self.b = b
        self.zn = z1 * (X2 / z1) ** (np.arange(N) / (N - 1))   # geometric nodes z1 .. X2
        self.w = np.sqrt(1.0 + 2.0 * b) - 1.0          # far-wake velocity increment (exact)
        self.fitR(np.ones(N))                          # R(z) = 1: cylinder
        self.gn = np.full(N, self.w)                   # gamma at nodes; start: far-wake value
        self.psie = 0.5 * (1.0 + self.w / 2.0)         # start: momentum theory

    # --- slipstream radius representation: natural cubic spline through (0,1),(zn,Rn) ---
    def R(self, z):
        return self.spl(np.clip(z, 0.0, X2))

    def Rz(self, z):
        return self.dspl(np.clip(z, 0.0, X2))

    def fitR(self, Rn):
        self.spl = CubicSpline(np.r_[0.0, self.zn], np.r_[1.0, Rn], bc_type='natural')
        self.dspl = self.spl.derivative()

    def gam(self, z):
        return np.interp(z, self.zn, self.gn)          # const. extrapolation outside

    # --- field quantities ---
    def _tailpars(self):
        return self.gn[-1], self.R(X2)

    def Psi(self, rho, z, sub=None):
        """Stokes stream function. sub=(gam0,) subtract log singularity at z'=z (point on sheet)."""
        f = lambda zp: self.gam(zp) * self.R(zp) * I011(self.R(zp), rho, z - zp)
        ginf, Rinf = self._tailpars()
        if sub is None:
            if 0.0 < z < X2:
                v = quad(f, 0, z, limit=200)[0] + quad(f, z, X2, limit=200)[0]
            else:
                v = quad(f, 0, X2, limit=200, points=[z] if 0 <= z <= X2 else None)[0]
        else:
            g0 = sub
            S = lambda zp: -(g0 / np.pi) * np.log(abs(z - zp))
            fr = lambda zp: 0.0 if zp == z else f(zp) - S(zp)
            v = quad(fr, 0, z, limit=200)[0] + quad(fr, z, X2, limit=200)[0]
            v += -(g0 / np.pi) * (xlnx(z) - z + xlnx(X2 - z) - (X2 - z))
        v += ginf * Rinf * Im111(Rinf, rho, X2 - z)
        return 0.5 * rho * rho + 0.5 * rho * v

    def Vz(self, rho, z, sub=None):
        """axial velocity; sub=(gam0,R0,Rz0) => point on sheet, returns mean velocity (PV)"""
        f = lambda zp: self.gam(zp) * self.R(zp) * I110(self.R(zp), rho, z - zp)
        ginf, Rinf = self._tailpars()
        if sub is None:
            if rho < 1e-12:
                f = lambda zp: self.gam(zp) * self.R(zp) ** 2 / (self.R(zp) ** 2 + (z - zp) ** 2) ** 1.5
            if 0.0 < z < X2:
                v = quad(f, 0, z, limit=200)[0] + quad(f, z, X2, limit=200)[0]
            else:
                v = quad(f, 0, X2, limit=200)[0]
        else:
            g0, R0, Rz0 = sub
            c1 = -g0 / (2.0 * np.pi * R0)
            c2 = -g0 * Rz0 / (np.pi * (1.0 + Rz0 * Rz0))
            fr = lambda zp: 0.0 if zp == z else f(zp) - c1 * np.log(abs(z - zp)) - c2 / (z - zp)
            v = quad(fr, 0, z, limit=200)[0] + quad(fr, z, X2, limit=200)[0]
            v += c1 * (xlnx(z) - z + xlnx(X2 - z) - (X2 - z))
            if c2 != 0.0:
                v += c2 * np.log(z / (X2 - z))
        if rho < 1e-12:
            v += ginf * (1.0 - (X2 - z) / np.sqrt(Rinf ** 2 + (X2 - z) ** 2))
        else:
            v += ginf * Rinf * I010(Rinf, rho, X2 - z)
        return 1.0 + 0.5 * v

    def Vr(self, rho, z):
        f = lambda zp: np.sign(z - zp) * self.gam(zp) * self.R(zp) * I111(self.R(zp), rho, z - zp)
        ginf, Rinf = self._tailpars()
        if 0.0 < z < X2:
            v = quad(f, 0, z, limit=200)[0] + quad(f, z, X2, limit=200)[0]
        else:
            v = quad(f, 0, X2, limit=200)[0]
        v -= ginf * Rinf * I011_tail(Rinf, rho, X2 - z)
        return 0.5 * v

    # --- one Conway-type iteration ---
    def iterate(self, relax=0.7, relaxg=0.7):
        b = self.b
        # stream function at the disk edge (start of the sheet)
        self.psie = self.Psi(1.0, 0.0, sub=self.gam(0.0))
        Rn = self.R(self.zn)
        Rnew = np.empty_like(Rn)
        for i, z in enumerate(self.zn):
            G = (self.Psi(Rn[i], z, sub=self.gn[i]) - 0.5 * Rn[i] ** 2) * 2.0 / Rn[i]
            Rq = -0.5 * G + np.sqrt(0.25 * G * G + 2.0 * self.psie)
            Rnew[i] = Rn[i] + relax * (Rq - Rn[i])
        dR = np.sqrt(np.mean((Rnew - Rn) ** 2))
        self.fitR(Rnew)
        # sheet strength from Bernoulli jump
        Rn = self.R(self.zn); Rzn = self.Rz(self.zn)
        gnew = np.empty_like(self.gn)
        self.vzm = np.empty_like(self.gn)
        for i, z in enumerate(self.zn):
            if z < X2:
                self.vzm[i] = self.Vz(Rn[i], z, sub=(self.gn[i], Rn[i], Rzn[i]))
            else:   # last node: slope ~0, omit 1/zeta subtraction
                self.vzm[i] = self.Vz(Rn[i], z, sub=(self.gn[i], Rn[i], 0.0))
            gnew[i] = b / self.vzm[i]
        dg = np.sqrt(np.mean((gnew - self.gn) ** 2))
        self.gn = self.gn + relaxg * (gnew - self.gn)
        return dR, dg


def xlnx(x):
    return x * np.log(x) if x > 0 else 0.0


if __name__ == "__main__":
    import sys
    b = float(sys.argv[1]) if len(sys.argv) > 1 else -4.0 / 9.0
    s = Sheet(b)
    print("b =", b, " CT(prop sign) = 2b =", 2 * b, " w_inf =", s.w)
    for L in range(1, 13):
        dR, dg = s.iterate()
        print(f"it {L:2d}  psie={s.psie:.6f}  R(X2)={s.R(X2):.6f}  dR={dR:.2e} dg={dg:.2e}")
    w = s.w
    print("psie            =", s.psie)
    print("psie (mom.th.)  =", 0.5 * (1 + w / 2))
    print("R_inf exact from psie =", np.sqrt(2 * s.psie / (1 + w)), "  R(X2) =", s.R(X2),
          "  mom.th. =", np.sqrt((1 + w / 2) / (1 + w)))
    print("CP = 4 b psie   =", 4 * b * s.psie, "  mom.th. CT(1+a):", 2 * b * (1 + w / 2))
    np.set_printoptions(precision=5, suppress=True, linewidth=150)
    print(" z      R       Rz      gamma    Vzm")
    for z, R, Rz, g, v in zip(s.zn, s.R(s.zn), s.Rz(s.zn), s.gn, s.vzm):
        print(f"{z:7.4f} {R:8.5f} {Rz:8.5f} {g:8.5f} {v:8.5f}")
    # disk-plane axial velocity
    print("Vz(r,0):", [(r, round(s.Vz(r, 0.0), 5)) for r in (0.0, 0.3, 0.6, 0.9, 0.99, 1.01, 1.2, 1.5)])
    print("Vz axis z:", [(z, round(s.Vz(0.0, z), 5)) for z in (-2, -1, 0, 1, 2, 5, 9)])
    # consistency checks on the sheet: streamline condition and velocity jump
    print("checks (z, Vr/Vz-Rz in/out, jump (Vzin-Vzout)(1+Rz^2)-gamma):")
    for z in (0.1, 0.5, 1.0, 2.0, 5.0):
        R0 = s.R(z); Rz0 = s.Rz(z); e = 2e-3
        vi, vo = s.Vz(R0*(1-e), z), s.Vz(R0*(1+e), z)
        ri, ro = s.Vr(R0*(1-e), z), s.Vr(R0*(1+e), z)
        print(f"  {z:5.2f}  {ri/vi-Rz0:9.5f} {ro/vo-Rz0:9.5f}   {(vi-vo)*(1+Rz0**2)-s.gam(z):9.5f}")
    np.save(f"proto_b{b:.4f}.npy", dict(zn=s.zn, gn=s.gn, Rn=s.R(s.zn), psie=s.psie, vzm=s.vzm), allow_pickle=True)
