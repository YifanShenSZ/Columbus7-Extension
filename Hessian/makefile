flag = -m64 -xCORE-AVX2 -mtune=core-avx2 -no-prec-div -fp-model fast=2 -mkl -parallel -O3

vibration.exe: vibration.f90
	ifort $(flag) -lFL $^ -o $@
