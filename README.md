Solution of mean-field self-consistent equations in the pseudogap state for the slave-boson formalism of high-temperature superconductivity in the cuprates.

Based on [Kotliar and Liu, PRB **38**, 7 (1988)](http://prb.aps.org/abstract/PRB/v38/i7/p5142_1). Includes hopping in the c direction from the start.

Unit testing is done with [fUnit](http://nasarb.rubyforge.org/funit/).

Considerable advice taken from the Cambridge [course on Modern Fortran](http://www-uxsup.csx.cam.ac.uk/courses/Fortran/).

## Modules

* double.f90
    * Define double-precision data type.
* environment.f90
    * Define data type holding relevant environment variables needed for evaluating error in self-consistent equations.
* brillouin.f90
    * Evaluate a sum over the first Brillouin zone (square lattice).
* scsolve.f90
    * Define data type for a self-consistent equation. Solve individual equations.
* scsystem.f90
    * Define data type for a system of self-consistent equations. Solve that system to given tolerances.
* sbzequations.f90
    * Define the self-consistent equations for the slave-boson system.
* driver.f90
    * Read environment(s) from a specified input file, solve the corresponding self-consistent system(s), and output the results to a specified output file.
* main.f90
    * Call driver using command-line arguments for input and output files.
