libdir = /usr/local/lib
incdir = /usr/local/include

libflg = -I$(incdir) -L$(libdir) -lFL
cplflg = -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -O3

loop.exe: loop.f90
	ifort $(libflg) $(cplflg) $^ -o $@

path.exe: path.f90
	ifort $(libflg) $(cplflg) $^ -o $@
