uniform-load: Conway actuator disk with UNIFORM load (h - h0 = b = const, J = 0)
================================================================================
Companion to ../parabolic-load (Conway-main.f).  Theory and description:
../../uniform-load-derivation.pdf (.tex).

Build:   gfortran -O3 uniform-main.f -o uniform.exe      (self-contained,
         math.f is NOT needed: the three cylinder integrals of GENINT are
         evaluated in closed form in subroutine CYLINT)
Run:     ./uniform.exe   and enter
            b               (h - h0 inside the slipstream, cT = 2b;
                             wind turbine b < 0, e.g. -0.444444 = Betz cT = 8/9;
                             propeller b > 0; b > -0.5 required)
            X2 NN NITER     (0 0 0 = defaults 40 40 30; for cT -> 1 use X2 = 100)

Output (as Conway-main.f): Stream.DAT, Solution.DAT, Psi.DAT, Vz1..10.DAT,
Vr1..10.DAT (perturbation velocities Vz-1, Vr on the profiles z = -2.5 .. 2.0),
AxVz0.DAT (axis, exact), plus Gamma.DAT: z, R, dR/dz, gamma_z, Vzm on the sheet.

Results in this directory: run for b = -4/9 (Betz optimum), defaults.
   psi_e = 0.333830 (momentum theory 1/3), cP = 0.59348 (Betz 0.59259),
   R_inf = 1.4153 (sqrt 2 = 1.4142)
python-check/: independent SciPy prototype (kernels.py, proto.py) used for
   verification (agrees to 5-6 digits).
