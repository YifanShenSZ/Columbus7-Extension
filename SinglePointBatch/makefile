MyLibDir = /home-4/yshen57@jhu.edu/Library/Fortran-Library
MyLib = $(MyLibDir)/General.f90 $(MyLibDir)/Mathematics.f90 $(MyLibDir)/LinearAlgebra.f90 \
$(MyLibDir)/mkl_rci.f90 $(MyLibDir)/NonlinearOptimization.f90 $(MyLibDir)/Chemistry.f90 \
$(MyLibDir)/GeometryTransformation.f90

GenColumbusInput.exe: $(MyLib) GenColumbusInput.f90
	ifort -mkl $^ -o GenColumbusInput.exe

clean:
	rm *.mod