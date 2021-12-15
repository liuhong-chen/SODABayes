# SODABayes
subsetted orthogonal data augmentation for Bayesian genomic prediction

supportted models: BayesA, BayesB, BayesC, BayesCpi, BayesR, and BayesRc
# source code
SODABayes.f90
globals.f90
parallelrng.f90
useroptions.f90
utilities.f90
soda.f90
models.f90
# compiling
ifort -o soda -qopenmp globals.f90 parallelrng.f90 utilities.f90 useroptions.f90 soda.f90 models.f90 SODABayes

or

gfortran -o soda -fopenmp globals.f90 parallelrng.f90 utilities.f90 useroptions.f90 soda.f90 models.f90 SODABayes
# Usage
read the manual for instructions

# license
This program is free for academic uses. You can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
Please report to the author(s) any bugs you found or any suggestions that you think could 
potentially improve the performance of the program.
