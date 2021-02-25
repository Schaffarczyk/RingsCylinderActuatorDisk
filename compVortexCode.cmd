echo 'compile module mem'
gfortran -c -O3 mem.f
echo 'compile rest'
gfortran -c -O3 -fno-automatic -fbounds-check -fdefault-real-8 VortexCode.f  SubA.f 
echo 'link all together'
gfortran mem.o VortexCode.o  SubA.o -o VortexCode.exe
