MyLibDir = /home-4/yshen57@jhu.edu/Library/Fortran-Library
MyLib = $(MyLibDir)/General.f90 $(MyLibDir)/Mathematics.f90 $(MyLibDir)/LinearAlgebra.f90 \
$(MyLibDir)/mkl_rci.f90 $(MyLibDir)/NonlinearOptimization.f90 $(MyLibDir)/Chemistry.f90 \
$(MyLibDir)/GeometryTransformation.f90

GenLoop: $(MyLib) GenLoop.f90
	ifort -mkl $^ -o GenLoop.exe

FiniteDifference: $(MyLib) FiniteDifference.f90
	ifort -mkl $^ -o FiniteDifference.exe

NormalModeLoop: $(MyLib) NormalModeLoop.f90
	ifort -mkl $^ -o NormalModeLoop.exe

clean:
	rm *.mod