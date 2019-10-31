libdir = /home-4/yshen57@jhu.edu/lib
incdir = /home-4/yshen57@jhu.edu/include

GenGeom.exe: GenGeom.f90
	ifort -I$(incdir) -L$(libdir) -lFL -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -O3 $^ -o GenGeom.exe
