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
# read the manual for instructions
