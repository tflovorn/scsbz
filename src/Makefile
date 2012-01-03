FC = gfortran
FFLAGS =

LDFLAGS =

all: environment.mod

double.mod double.o: double.f90
	$(FC) $(FFLAGS) -c -o double.o double.f90

environment.mod environment.o: double.mod double.o environment.f90
	$(FC) $(FFLAGS) -c -o environment.o environment.f90
