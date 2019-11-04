libdir = /home-4/yshen57@jhu.edu/lib
incdir = /home-4/yshen57@jhu.edu/include

GenLoop.exe FiniteDifference.exe NormalModeLoop.exe: GenLoop.f90 FiniteDifference.f90 NormalModeLoop.f90
	ifort -I$(incdir) -L$(libdir) -lFL -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -O3 GenLoop.f90 -o GenLoop.exe
	ifort -I$(incdir) -L$(libdir) -lFL -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -O3 FiniteDifference.f90 -o FiniteDifference.exe
	ifort -I$(incdir) -L$(libdir) -lFL -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -O3 NormalModeLoop.f90 -o NormalModeLoop.exe
