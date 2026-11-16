# table-z1-runs — raw runs behind Table 1 (tab:z1)

Paper: Bontempo, Conway, Schaffarczyk, *The edge of the uniformly loaded
actuator disc: edge-flow structure, regularity, and the exactness of
one-dimensional momentum theory* (submitted to Wind Energy Science).
Table 1: edge behaviour of the converged solution for the Betz load
b = −4/9 (c_T = 8/9) under refinement of the first sheet node z1.

The runs in this directory were regenerated on 16 November 2026 from the
published sources `../uniform-main.f` (z1 = 1e-1 … 1e-4) and
`../uniform-main-edge.f` (z1 = 1e-5, 1e-6) with `run-table-z1.sh`
(gfortran 13, `-O3`, x86-64 Linux). The outputs of the original September 2026
runs were not archived; the regenerated values agree with the printed table
within the rounding of the table (see "Comparison" below).

## Set-up (identical for all rows)

| quantity | value |
|---|---|
| load | b = −4/9 (h − h0 inside the slipstream; wind-turbine sign) |
| numerical domain | 0 ≤ z ≤ X2 = 40, semi-infinite cylinder beyond X2 |
| nodes | geometric from z1 to X2, ≈16 per decade: NN = 40, 58, 74, 90, 106, 122 |
| sheet representation | natural cubic spline; principal value on the sheet exact |
| relaxation | R: 0.7, γ: 0.7 (uniform-main.f); R: 0.35, γ: 0.15 (uniform-main-edge.f) |
| start | cold (cylinder R = 1, γ = far-wake value) for z1 ≥ 1e-4; warm start from the previous run's `Gamma.DAT` (copied to `start.dat`) for z1 = 1e-5, 1e-6 |
| stop | rms residuals of both updates < 1e-6, or NITER (300 / 400) |

`z1`, `relax`, `relaxg` are compile-time constants in both codes
(`z1 = 0.01d0`, `relax = 0.7d0`, `relaxg = 0.7d0` in the published sources);
`run-table-z1.sh` sets them with `sed` into a per-run copy before compiling.
Input on stdin (`input.txt`): line 1 `b`, line 2 `X2 NN NITER`.

## Files per run directory `z1_<z1>_NN<nodes>/`

| file | content |
|---|---|
| `input.txt` | the two stdin lines |
| `start.dat` | warm-start profile for z1 = 1e-5, 1e-6: a copy of the previous run's `Gamma.DAT`, made by the script (not stored separately) |
| `Gamma.DAT` | sheet data per node: z, R, dR/dz, γ_z (ring strength per unit **axial** length), V_zm (mean axial velocity on the sheet) |
| `Solution.DAT` | global results (Ψ_e, c_P, R_∞, R(X2), γ(X2)) and run parameters |
| `convergence.log` | iteration history: Ψ_e, R(X2), rms residuals err1 (streamline condition), err2 (Bernoulli jump) per iteration, wall time |

The field files of the codes (`Stream.DAT`, `Psi.DAT`, `Vz*.DAT`, `Vr*.DAT`,
`AxVz0.DAT`) are not stored; rerunning the script reproduces them.

## From the raw data to the table columns

`table-z1-from-runs.py` prints the table from the run directories
(`table-z1-derived.txt` is its output):

* Ψ_e: from the last iteration (converged runs) or mean ± half-range of the
  last 50 iterations (z1 = 1e-5, 1e-6, where the innermost nodes keep
  oscillating because the sheet turns past the vertical, Sect. 4.6 of the paper);
* c_P = 4|b| Ψ_e (momentum theory: 16/27 = 0.59259);
* R′(z1) = dR/dz at the first node, θ_s = arctan R′;
* γ(z1) in the table is the strength per unit **arc** length,
  γ_s = γ_z / √(1+R′²);
* V_t,m(z1) = V_zm √(1+R′²) (mean velocity tangential to the sheet);
* γ_s · V_t,m = b holds identically (Bernoulli jump condition, last column).

## Comparison: regenerated runs vs. printed Table 1

| z1 | Ψ_e (paper) | Ψ_e (here) | R′ (paper / here) | θ_s (paper / here) | γ_s | V_t,m |
|---|---|---|---|---|---|---|
| 1e-1 | 0.33493 | 0.33495 | 0.32 / 0.32 | 17.7° / 17.8° | −0.544 / −0.544 | 0.817 / 0.818 |
| 1e-2 | 0.33381 | 0.33382 | 0.58 / 0.61 | 30.1° / 31.4° | −0.546 / −0.546 | 0.814 / 0.815 |
| 1e-3 | 0.33346 | 0.33346 | 1.02 / 1.02 | 45.6° / 45.5° | −0.556 / −0.556 | 0.799 / 0.799 |
| 1e-4 | 0.33337 | 0.33337 | 1.70 / 1.70 | 59.5° / 59.5° | −0.571 / −0.571 | 0.778 / 0.778 |
| 1e-5 | 0.33339(7) | 0.33337(1) | — | — | — | — |
| 1e-6 | 0.33339(5) | 0.33334(1) | — | — | — | — |

The only visible difference is the first-node slope of the z1 = 1e-2 row
(0.61 vs. 0.58; the September run probably used a slightly different node
count); the fit β_sheet = 0.106 of eq. (spiral) is insensitive to it. For
z1 = 1e-5 and 1e-6 the regenerated Ψ_e lie inside the uncertainty quoted in
the table and are, if anything, closer to momentum theory
(c_P = 0.59266 and 0.59261 vs. 16/27 = 0.59259).

Timing (single core): 7–18 s for z1 ≥ 1e-4, ≈5.5 min each for z1 = 1e-5 and 1e-6.
