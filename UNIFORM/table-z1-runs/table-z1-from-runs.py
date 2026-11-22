#!/usr/bin/env python3
"""table-z1-from-runs.py -- derive the columns of Table 3 (tab:z1) of
'The edge of the uniformly loaded actuator disc' from the raw runs in this
directory.

For every run directory z1_*_NN* it reads
  Gamma.DAT : z, R, dR/dz, gamma_z (strength per unit axial length), Vzm
  convergence.log : iteration history (psie per iteration)
and prints
  z1, Psi_e, Psi_e - 1/3, cP = 4|b| Psi_e, R'(z1), theta_s = arctan R',
  gamma_s(z1) = gamma_z / sqrt(1+R'^2)  (strength per unit arc length),
  V_t,m(z1)   = Vzm * sqrt(1+R'^2)      (mean tangential velocity),
  gamma_s * V_t,m  (= b identically by the Bernoulli jump condition).
For the runs z1 = 1e-5, 1e-6 the innermost nodes do not converge (the sheet
turns past the vertical, Sect. 4.6 of the paper); Psi_e is then reported as
the mean and half-range over the last 50 iterations, and the geometry
columns are marked '--'.
"""
import glob, os, re, sys
import numpy as np

b = -4.0/9.0
rows = []
for d in sorted(glob.glob('z1_*_NN*'), key=lambda s: -float(s.split('_')[1])):
    z1 = float(d.split('_')[1])
    g = np.loadtxt(os.path.join(d, 'Gamma.DAT'), comments='#')
    hist = [float(m.group(1)) for m in re.finditer(r'psie =\s*([0-9.]+)',
            open(os.path.join(d, 'convergence.log')).read())]
    nit = len(hist)
    niter = int(open(os.path.join(d, 'input.txt')).read().split()[3])
    tail = np.array(hist[-50:])
    psie, dpsie = tail.mean(), (tail.max()-tail.min())/2
    converged = nit < niter or dpsie < 2e-6      # stopped on criterion, or stationary
    if converged:
        psie, dpsie = hist[-1], 0.0
    z, R, dR, gam, vzm = g[0, :5]
    s = np.sqrt(1+dR**2)
    if converged:
        geo = (f"{dR:6.2f}  {np.degrees(np.arctan(dR)):6.1f}  {gam/s:7.3f}  "
               f"{vzm*s:6.3f}  {gam/s*vzm*s:7.3f}")
        ps = f"{psie:.5f}    "
    else:
        geo = "   --      --       --      --       --"
        ps = f"{psie:.5f}({int(round(dpsie*1e5))})"
    rows.append((z1, ps, psie-1/3, 4*abs(b)*psie, geo, nit))

print("  z1     Psi_e        Psi_e-1/3   cP       R'(z1) theta_s gamma_s  Vt,m   gam*Vt  iterations")
for z1, ps, d3, cp, geo, nit in rows:
    print(f"{z1:7.0e}  {ps}  {d3:9.2e}  {cp:.5f}  {geo}  {nit}")
