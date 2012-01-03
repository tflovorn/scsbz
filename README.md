Solution of mean-field self-consistent equations for the slave-boson formalism of high-temperature superconductivity in the cuprates.

Includes hopping in the c direction.

Unit testing is done with [fUnit](http://nasarb.rubyforge.org/funit/).

## Modules

* double.f90
    * Define double-precision data type.
* rootfinder.f90
    * Find a root of the given function. Assume only one root exists in the specified range.
* environment.f90
    * Define data type holding relevant environment variables needed for evaluating error in self-consistent equations.
* brillouin.f90
    * Evaluate a sum over the first Brillouin zone.
* scsolve.f90
    * Define data type for a self-consistent equation. Solve individual equations and systems of equations.
* scequations.f90
    * Define the self-consistent equations for the system.
* driver.f90
    * Read environment(s) from a specified input file, solve the corresponding self-consistent system(s), and output the results to a specified output file.
* main.f90
    * Call driver using command-line arguments for input and output files.
