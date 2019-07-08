MyLibDir = ~/Library/Fortran-Library
MyLib = $(MyLibDir)/General.f90 $(MyLibDir)/Mathematics.f90 $(MyLibDir)/LinearAlgebra.f90 \
$(MyLibDir)/mkl_rci.f90 $(MyLibDir)/NonlinearOptimization.f90 $(MyLibDir)/Chemistry.f90 \
$(MyLibDir)/GeometryTransformation.f90

GenLoop: $(MyLib) GenLoop.f90
	ifort -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -ipo -O3 $^ -o GenLoop.exe

FiniteDifference: $(MyLib) FiniteDifference.f90
	ifort -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -ipo -O3 $^ -o FiniteDifference.exe

NormalModeLoop: $(MyLib) NormalModeLoop.f90
	ifort -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -ipo -O3 $^ -o NormalModeLoop.exe

HarmonicityCheck: $(MyLib) HarmonicityCheck.f90
	ifort -m64 -xCORE-AVX2 -mtune=core-avx2 -mkl -ipo -O3 $^ -o HarmonicityCheck.exe

clean:
	rm *.mod
