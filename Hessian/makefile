libdir = /usr/local/lib
incdir = /usr/local/include

.PHONY: make
make: GenLoop.exe FiniteDifference.exe NormalModeLoop.exe

GenLoop.exe: GenLoop.f90
	ifort -I$(incdir) -L$(libdir) -lFL -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -O3 $^ -o GenLoop.exe

FiniteDifference.exe: FiniteDifference.f90
	ifort -I$(incdir) -L$(libdir) -lFL -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -O3 $^ -o FiniteDifference.exe

NormalModeLoop.exe:  NormalModeLoop.f90
	ifort -I$(incdir) -L$(libdir) -lFL -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -O3 $^ -o NormalModeLoop.exe
