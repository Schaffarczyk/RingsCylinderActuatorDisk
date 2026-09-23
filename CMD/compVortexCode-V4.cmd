echo 'V4 needs mem.f and SubA-V3.f in this directory (copies included)'
gfortran -c -O3 -fno-automatic -fdefault-real-8 mem.f
gfortran -c -O3 -fno-automatic -fdefault-real-8 VortexCode-V4.f SubA-V3.f
gfortran mem.o VortexCode-V4.o SubA-V3.o -o VortexCode-V4.exe
