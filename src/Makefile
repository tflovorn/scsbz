FC = gfortran
FFLAGS =

LDFLAGS =

all: rootfinder.mod brillouin.mod

double.mod double.o: double.f90
	$(FC) $(FFLAGS) -c -o double.o double.f90

rootfinder.mod rootfinder.o: double.mod double.o rootfinder.f90
	$(FC) $(FFLAGS) -c -o rootfinder.o rootfinder.f90

environment.mod environment.o: double.mod double.o environment.f90
	$(FC) $(FFLAGS) -c -o environment.o environment.f90

brillouin.mod brillouin.o: environment.mod environment.o brillouin.f90
	$(FC) $(FFLAGS) -c -o brillouin.o brillouin.f90
