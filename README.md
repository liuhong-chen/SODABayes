# SODABayes
subsetted orthogonal data augmentation for fast parallel implementation of Bayesian linear regression models for whole-genome analyses

supportted models: BayesA, BayesB, BayesC, BayesCpi, BayesR, and BayesRc
# Source code
SODABayes.f90 – main program;
globals.f90 – definition of global variables;
parallelrng.f90 – parallel random number and random variable generator;
utilities.f90 – functions and subroutines that can be used by other modules;
useroptions.f90 – a module that handles various user defined options;
soda.f90 – a module that generates the orthogonal augmenting matrix;
models.f90 – Bayesian models using either regular or SODA algorithms.
# Compiling
ifort -o sodaBayes -qopenmp globals.f90 parallelrng.f90 utilities.f90 useroptions.f90 soda.f90 models.f90 SODABayes.f90

or

gfortran -o sodaBayes -fopenmp globals.f90 parallelrng.f90 utilities.f90 useroptions.f90 soda.f90 models.f90 SODABayes.f90
# Usage
The program can be run with command line options or a parameter file. Refer to the manual for instructions.

# License
This program is free for research use. Please contact the authors for non-research use. 
You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
Please report to the author(s) any bugs you found or any suggestions that you think could 
potentially improve the performance of the program.
