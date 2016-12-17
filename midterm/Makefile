#########################################################
#   make: compiles fortran source and makes plots
#       'stress.ps' and 'pre_cor.ps' by executing codes
#       and calling gnuplot plotting scripts
#   make clean: removes all files aside from source code
#########################################################
all: prob1 prob2 plots

prob1: prob1.f90 spline.o
	gfortran -o prob1 prob1.f90 spline.o
spline.o: spline.f90
	gfortran -c spline.f90
prob2: prob2.f90
	gfortran -o prob2 prob2.f90

plots: stress.dat pre_cor.dat stress.ps pre_cor.ps

stress.dat:
	./prob1
pre_cor.dat:
	./prob2
stress.ps: stress.gps
	gnuplot stress.gps
pre_cor.ps: pre_cor.gps
	gnuplot pre_cor.gps

clean:
	rm -f *.o *.mod *.dat *.ps prob1 prob2
