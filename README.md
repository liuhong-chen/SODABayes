# SODABayes
subsetted orthogonal data augmentation for Bayesian genomic prediction
# source code
SODABayes.f90
globals.f90
parallelrng.f90
useroptions.f90
utilities.f90
soda.f90
models.f90
# compile
ifort -o soda10 -qopenmp globals.f90 parallelrng.f90 utilities.f90 useroptions.f90 soda.f90 models.f90 SODABayes

or

gfortran -o soda10 -fopenmp globals.f90 parallelrng.f90 utilities.f90 useroptions.f90 soda.f90 models.f90 SODABayes
# read manual for instructions
