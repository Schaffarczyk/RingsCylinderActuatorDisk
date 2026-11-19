# FWRV-AD runs 010–014 — free-wake ring-vortex actuator-disc output

Paper: Bontempo, Conway, Schaffarczyk, *The edge of the uniformly loaded
actuator disc: edge-flow structure, regularity, and the exactness of
one-dimensional momentum theory* (submitted to Wind Energy Science).

The files in this directory are the output of the free-wake ring-vortex
actuator-disc code FWRV-AD (`rngVxAD_UniLoad`) of R. Bontempo and M. Manna
for the uniformly loaded disc with the Betz load (ΔH = −8/9, i.e.
b = −4/9, c_T = 8/9). The runs were computed by **Rodolfo Bontempo**
(Università degli Studi di Napoli Federico II) in October and November 2026
and are included here with his permission. The FWRV-AD code itself is not
part of this repository.

They are the data behind the last line of Table `fwrv`, Table `fwrvedge`,
Fig. `fwrvedge`, the FWRV-AD entries of Table `adiabatic` and the paragraphs
"Edge-refined FWRV-AD runs (geometric panelling)" and "Beyond the vertical"
of the paper.

## Runs

| run | name in the paper | code version | panelling law | panels | first panel (converged) | c_P | R_∞ | iterations |
|---|---|---|---|---|---|---|---|---|
| 010 | – (predecessor of 014) | 3.4 | cosine | 10 000 | 1.5e-7 | 0.592644 | 1.41274 | 220 |
| 011 | G6 | 3.4.1 | geometric (first length) | 6 000 | 5.0e-10 | 0.592827 | 1.40261 | 222 |
| 012 | G10 | 3.4.1 | geometric (first length) | 10 000 | 5.0e-10 | 0.592732 | 1.40172 | 220 |
| 013 | G10⁺ | 3.4.1 | geometric (first length) | 10 000 | 3.8e-12 | 0.592769 | 1.39830 | 219 |
| 014 | C10 | 3.4.1 | cosine | 10 000 | 4.3e-7 | 0.592645 | 1.41274 | 220 |

Common settings: quadruple precision, wake-radius residual 1e-12,
relaxation factor 0.2, semi-infinite vortex cylinder from z = 50 with the
momentum-theory strength, complete elliptic integrals by AGM / Carlson.
Momentum theory: c_P = 16/27 = 0.592593, R_∞ = √2 = 1.41421.
Run 010 was made with code version 3.4 and is superseded by run 014
(version 3.4.1, same panelling); it is kept for comparison.

## Files (per run `schaffarczykNNNquad`)

| file | content |
|---|---|
| `.ech` | echo of the input file (all settings of the run) |
| `.log` | run log with the iteration history |
| `Res.dat` | residual per iteration |
| `Wake.dat` | converged wake sheet, one line per panel, from the disc edge downstream: `zWk rWk dnstyStrgthWk pnlLngWk cosineWk sineWk slopeWk vzWk vrWk vmWk` (axial and radial coordinate, sheet strength γ, panel length, cosine and sine of the panel angle, slope, axial, radial and mean velocity on the sheet); disc radius 1, free-stream speed 1, edge at (z, r) = (0, 1) |
| `velDisk.dat` | velocities on the disc plane, 1000 radial points |
| `velBdOnDisk.dat` | body-induced velocities on the disc (zero, open configuration) |
| `perfCoeff.dat` | thrust and power coefficients, mean axial velocity at the disc, wake radius |
| `CPUTime.dat` | computational cost |

The wake-convergence histories (`WkConv.dat`, 140–230 MB per run) are not
included because of their size; they are available from the authors on
request.

## Evaluation scripts

* `fwrv-edge-analysis.py` — reads the five `Wake.dat` files and prints, per
  decade of the distance ρ from the edge, the number of panels, the sheet
  angle θ_s, its slope per decade, β_geo, γ, V_in, V_out, V_out/V_in and
  β_vel = ln(V_out/V_in)/2π. Output: `fwrv-edge-analysis.txt`; it also
  writes `runs.npy` for the second script.
* `fwrv-edge-blocks.py` — from run 013 (G10⁺): the rows of Table `fwrvedge`
  and Table `adiabatic` and the plot coordinates of Fig. `fwrvedge`,
  including the integration of the adiabatic law
  dβ/d ln ρ = −(β²/π) tanh(2πβ). Output: `fwrv-edge-blocks.txt`.

Requirements: Python 3 with NumPy. Run both scripts in this directory,
`fwrv-edge-analysis.py` first.

## Licence

Data and documentation: CC BY 4.0; scripts: MIT (see `LICENSE` in the
repository root). Please credit R. Bontempo for the FWRV-AD runs and cite
the paper above and Bontempo & Manna, *Wind Energy* (2017), for the method.
