"""V4 prototype: vK&L-type discretization (WE 2016, App. A+B).

- N vortex rings with circulations Gam_i AT NODES (x_i, r_i);
  ring 1 fixed at the disc edge (0, 1).
- vortex kernel delta: for evaluation distances rho < delta/2 the
  Marshall (smoothed) expressions are used -> no midpoint trick needed.
- semi-infinite tube at x0 = Ltube, radius Rw, strength gw (momentum th.).
- two-phase convergence scheme:
    A: r_i += dA*(Psi_wake - Psi_i)/(r_i*vx_i)      (i >= 2)
    B: r_i += dB*sgn*vn_i                            (fine tuning)
  Gam_i += dG*(Gam*_i - Gam_i),  Gam*_i = -ct/(2*v_i)*ds_i (force free)
"""
import numpy as np
from scipy.special import ellipk, ellipe, elliprf, elliprj

def ellpi(n, m):
    n = np.asarray(n, float); m = np.asarray(m, float)
    y = 1.0 - m
    return elliprf(0.0, y, 1.0) + (n/3.0)*elliprj(0.0, y, 1.0, 1.0-n)

class ADV4:
    def __init__(self, ct=-8/9, N=1000, Ltube=None, delta=0.002,
                 c1=2.72, c2=0.7, beta=None):
        # beta: if set, kernel radius per ring = beta*ds_i (consistent
        # with the analytic strip integral for beta = exp(-1/2) ~ 0.607);
        # if None, fixed delta as in vK&L.
        self.ct, self.N, self.delta, self.beta = ct, N, delta, beta
        vzinf = np.sqrt(1.0+ct)
        vz0 = 0.5*(1.0+vzinf)
        self.Rw = np.sqrt(vz0/vzinf)
        self.gw = 1.0-np.sqrt(1.0+ct)          # tube strength (gamma)
        self.psiwake = 0.25*(1.0+np.sqrt(1.0+ct))
        a = 0.5*(np.sqrt(1.0+ct)-1.0)
        self.cpm = 4.0*a*(1.0+a)**2
        self.Ltube = 30.0*self.Rw if Ltube is None else Ltube
        # node grid (vK&L-like stretched cosine), ring 1 at x=0
        i = np.arange(0, N)
        phi = i*np.pi/(c1*(N-1))
        g = (1.0-np.cos(phi))/(1.0-np.cos(np.pi/c1))
        self.x = self.Ltube*g**c2
        self.x[0] = 0.0
        # initial shape: van Kuik D.10-like (as in adcode nwt=3)
        r0l = 1.0 + i/(N-1.0)*(self.Rw-1.0)
        u0 = 1.0-0.5*self.gw
        uw = 1.0-0.5*self.gw*(1.0+self.x/np.sqrt(r0l**2+self.x**2))
        self.r = np.sqrt(u0/uw)
        self.r[0] = 1.0
        # circulations from far-wake strength * local spacing
        self.Gam = self.gw*self.ds_centered()

    # centered arclength weights (ds_i)
    def ds_centered(self):
        s = self.arclength()
        ds = np.zeros(self.N)
        ds[1:-1] = 0.5*(s[2:]-s[:-2])
        ds[0] = 0.5*(s[1]-s[0])
        # last ring: half spacing to previous + half gap to tube start
        gap = np.hypot(self.Ltube-self.x[-1], self.Rw-self.r[-1])
        ds[-1] = 0.5*(s[-1]-s[-2]) + 0.5*gap
        return ds

    def arclength(self):
        seg = np.hypot(np.diff(self.x), np.diff(self.r))
        return np.concatenate([[0.0], np.cumsum(seg)])

    # ---- ring-to-point induction matrices with kernel ----
    def ring_fields(self, xp, rp):
        """vz, vr, psi at points (xp, rp) from all N rings (kernel-smoothed)."""
        xj = self.x[None, :]; rj = self.r[None, :]
        G = self.Gam[None, :]
        XP = np.asarray(xp, float)[:, None]; RP = np.asarray(rp, float)[:, None]
        ks = 4.0*RP*rj/((rj+RP)**2+(XP-xj)**2)
        ks = np.clip(ks, 1e-14, 1.0-1e-14)
        K = ellipk(ks); E = ellipe(ks)
        f1 = rj/(2.0*RP)*ks/(1.0-ks)-(2.0-ks)/(2.0*(1.0-ks))
        vz = (-G*np.sqrt(ks)*np.sqrt(rj/RP)/(4.0*np.pi*rj))*(f1*E+K)
        g1 = (2.0-ks)/(2.0*(1.0-ks))*E - K
        vr = (-G*np.sqrt(ks)*(XP-xj)*(rj/RP)**1.5/(4.0*np.pi*rj*rj))*g1
        k = np.sqrt(ks)
        psi = (-G*np.sqrt(RP*rj)/(2.0*np.pi))*((2.0/k-k)*K-(2.0/k)*E)
        # kernel smoothing where rho < delta_j/2 (delta_j per ring)
        if self.beta is not None:
            dloc = (self.beta*self.ds_centered())[None, :]
        else:
            dloc = self.delta*np.ones((1, self.N))
        rho = np.hypot(XP-xj, RP-rj)
        m = rho < 0.5*dloc
        if m.any():
            # Marshall, sign-matched to the exact ring formulas above
            vzk = (-G/(4.0*np.pi*rj)*(np.log(16.0*rj/dloc)-0.25)
                   )*np.ones_like(rho)
            psik = (-G*rj/(2.0*np.pi)*(np.log(16.0*rj/dloc)-1.5)
                    )*np.ones_like(rho)
            vz = np.where(m, vzk, vz)
            vr = np.where(m, 0.0, vr)
            psi = np.where(m, psik, psi)
        return vz.sum(1), vr.sum(1), psi.sum(1)

    # ---- tube (cylinder) fields, CORRECT vr sign ----
    def tube_fields(self, xp, rp):
        R, gam = self.Rw, self.gw
        xc = np.asarray(xp, float)-self.Ltube
        r = np.asarray(rp, float)
        ks = 4.0*r*R/((R+r)**2+xc**2)
        n = 4.0*r*R/(R+r)**2
        f1 = ellipk(ks)+((R-r)/(R+r))*ellpi(n, ks)
        f2 = xc*np.sqrt(ks)/(2.0*np.pi*np.sqrt(r*R))
        f3 = np.where(np.abs(R-r) < 1e-12, 0.5, (np.sign(R-r)+1.0)/2.0)
        vz = -0.5*gam*(f3+f2*f1)
        k = np.sqrt(ks)
        f1r = ((2.0-ks)*ellipk(ks)-2.0*ellipe(ks))/k
        vr = +gam*f1r*np.sqrt(R/r)/(2.0*np.pi)      # corrected sign
        m = 4.0*r*R/((r+R)**2+xc**2)
        m0 = 4.0*r*R/((r+R)**2)
        aa = -gam*xc*np.sqrt(m*R*r)/(2.0*np.pi*m0)
        g1 = (1.0-m0+m0/m)*ellipk(m)
        g2 = -(m0/m)*ellipe(m)
        g3 = (m0-1.0)*ellpi(m0, m)
        g4 = np.where(r > R, 0.25*(-gam)*R**2, 0.25*(-gam)*r**2)
        psi = aa*(g1+g2+g3)+g4
        return vz, vr, psi

    def fields_at_rings(self):
        vz, vr, psi = self.ring_fields(self.x, self.r)
        tz, tr, tpsi = self.tube_fields(self.x, self.r)
        return vz+tz, vr+tr, psi+tpsi+0.5*self.r**2   # + inflow psi

    def gam_update(self, vz, vr, dG=0.05):
        vm = np.sqrt((1.0+vz)**2+vr**2)
        Gstar = -self.ct/(2.0*vm)*self.ds_centered()
        self.Gam += dG*(Gstar-self.Gam)

    def vn(self, vz, vr):
        """normal velocity at rings (tangent by centered differences)."""
        tx = np.gradient(self.x); tr = np.gradient(self.r)
        tn = np.hypot(tx, tr)
        return (-tr*(1.0+vz)+tx*vr)/tn

    def cp(self, nd=3000):
        rr = np.linspace(0.0, 1.0, nd+1)[1:]
        vz, vr, psi = self.ring_fields(np.zeros(nd), rr)
        tz, trr, tpsi = self.tube_fields(np.zeros(nd), rr)
        vzt = vz+tz
        a = np.trapezoid(2.0*rr*vzt, rr)
        return self.ct*(1.0+a)

    # ---- one iteration ----
    def step(self, phase="A", dA=0.05, dB=0.0025, dG=0.05, sgnB=+1.0,
             drmax=0.01):
        vz, vr, psi = self.fields_at_rings()
        self.gam_update(vz, vr, dG)
        if phase == "A":
            dr = dA*(self.psiwake-psi)/(self.r*(1.0+vz))
        else:
            dr = sgnB*dB*self.vn(vz, vr)
        dr = np.clip(dr, -drmax, drmax)
        dr[0] = 0.0                      # ring 1 pinned to the disc edge
        self.r += dr
        dpsi = self.psiwake-psi
        return dpsi, self.vn(vz, vr)
