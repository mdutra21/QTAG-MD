This repository contains the Fortran source files for QTAG calculations.

QTAG.f90 - Main program file.
constants.f90 - contains some public variables and subroutine for reading input file.
initialize.f90 - contains subroutines for initializing basis description of wavefunction.
tier1.f90 - contains subroutines for various small tasks.
tier2.f90 - contains subroutines for larger tasks, depends on tier1.
prop.f90 - contains subroutines for propagation of basis functions.
reexpand.f90 - contains subroutines for types of reexpansions.
finalize.f90 - currently unused
sobol_dataset.f90 - subroutine for generating a list of Sobol numbers for basis coordinates.
rxsobl.f90 - reexpansion subroutine specifically for using coordinates from sobol_dataset.
simplex.f90 - currently unused
quartic.f90 - contains functions describing a quartic potential and its derivatives.
harmonic.f90 - contains functions describing a harmonic potential and its derivatives.

