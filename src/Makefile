FC = gfortran
FFLAGS =

LDFLAGS =

all: main.out

double.mod double.o: double.f90
	$(FC) $(FFLAGS) -c -o double.o double.f90

environment.mod environment.o: double.mod double.o environment.f90
	$(FC) $(FFLAGS) -c -o environment.o environment.f90

brillouin.mod brillouin.o: environment.mod environment.o brillouin.f90
	$(FC) $(FFLAGS) -c -o brillouin.o brillouin.f90

scsolve.mod scsolve.o: environment.mod environment.o scsolve.f90
	$(FC) $(FFLAGS) -c -o scsolve.o scsolve.f90

scsystem.mod scsystem.o: scsolve.mod scsolve.o scsystem.f90
	$(FC) $(FFLAGS) -c -o scsystem.o scsystem.f90

sbzequations.mod sbzequations.o: brillouin.mod brillouin.o scsystem.mod scsystem.o sbzequations.f90
	$(FC) $(FFLAGS) -c -o sbzequations.o sbzequations.f90

driver.mod driver.o: sbzequations.mod sbzequations.o driver.f90
	$(FC) $(FFLAGS) -c -o driver.o driver.f90

main.out: driver.mod driver.o main.f90
	$(FC) $(FFLAGS) $(LDFLAGS) -o main.out main.f90 driver.o sbzequations.o scsystem.o scsolve.o brillouin.o environment.o double.o
