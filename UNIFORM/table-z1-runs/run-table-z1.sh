#!/bin/bash
# run-table-z1.sh -- regenerate the raw runs behind Table 3 (tab:z1) of
# "The edge of the uniformly loaded actuator disc" (Bontempo/Conway/Schaffarczyk)
# b = -4/9 (Betz load, cT = 8/9), X2 = 40, ~16 nodes per decade.
#   z1 = 1e-1 .. 1e-4 : uniform-main.f       (relax = relaxg = 0.7, cold start)
#   z1 = 1e-5, 1e-6   : uniform-main-edge.f  (relax 0.35 / relaxg 0.15, warm start
#                        from the previous run's Gamma.DAT, see header of that file)
# z1 is a compile-time constant in both codes; it is set with sed into a copy.
set -e
SRC=$(cd "$(dirname "$0")"/.. && pwd)      # UNIFORM/ (this script lives in UNIFORM/table-z1-runs/)
B=-0.4444444444444444
run () {   # run <code.f> <z1> <NN> <NITER> <relax> <relaxg> <tag> [startfile]
  local code=$1 z1=$2 nn=$3 nit=$4 rx=$5 rg=$6 tag=$7 start=$8
  mkdir -p "$tag"; cd "$tag"
  sed -e "s/^      z1     = 0.01d0/      z1     = ${z1}/" \
      -e "s/^      relax  = 0.7d0/      relax  = ${rx}/" \
      -e "s/^      relaxg = 0.7d0/      relaxg = ${rg}/" "$SRC/$code" > prog.f
  grep -q "z1     = ${z1}" prog.f
  gfortran -O3 prog.f -o prog.exe
  [ -n "$start" ] && cp "../$start" start.dat
  printf '%s\n%s %s %s\n' "$B" 40 "$nn" "$nit" > input.txt
  ( time ./prog.exe < input.txt ) > run.log 2>&1
  rm -f prog.exe Vz*.DAT Vr*.DAT AxVz0.DAT Stream.DAT Psi.DAT
  tail -3 run.log; cd ..
}
run uniform-main.f      1.d-1  40 300 0.7d0  0.7d0  z1_1e-1_NN040
run uniform-main.f      1.d-2  58 300 0.7d0  0.7d0  z1_1e-2_NN058
run uniform-main.f      1.d-3  74 300 0.7d0  0.7d0  z1_1e-3_NN074
run uniform-main.f      1.d-4  90 300 0.7d0  0.7d0  z1_1e-4_NN090
run uniform-main-edge.f 1.d-5 106 400 0.35d0 0.15d0 z1_1e-5_NN106 z1_1e-4_NN090/Gamma.DAT
run uniform-main-edge.f 1.d-6 122 400 0.35d0 0.15d0 z1_1e-6_NN122 z1_1e-5_NN106/Gamma.DAT
