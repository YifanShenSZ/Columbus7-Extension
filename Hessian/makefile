libdir = /usr/local/lib
incdir = /usr/local/include

libflg = -I$(incdir) -L$(libdir) -lFL
cplflg = -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -O3

vibration.exe: vibration.f90
	ifort $(libflg) $(cplflg) $^ -o $@
