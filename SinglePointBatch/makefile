libdir = /usr/local/lib
incdir = /usr/local/include

loop.exe: loop.f90
	ifort -I$(incdir) -L$(libdir) -lFL -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -O3 $^ -o loop.exe

path.exe: path.f90
	ifort -I$(incdir) -L$(libdir) -lFL -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -O3 $^ -o path.exe
