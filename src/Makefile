FC = gfortran
FFLAGS =

LDFLAGS =

all: scsolve.mod brillouin.mod

double.mod double.o: double.f90
	$(FC) $(FFLAGS) -c -o double.o double.f90

environment.mod environment.o: double.mod double.o environment.f90
	$(FC) $(FFLAGS) -c -o environment.o environment.f90

brillouin.mod brillouin.o: environment.mod environment.o brillouin.f90
	$(FC) $(FFLAGS) -c -o brillouin.o brillouin.f90

scsolve.mod scsolve.o: environment.mod environment.o scsolve.f90
	$(FC) $(FFLAGS) -c -o scsolve.o scsolve.f90
