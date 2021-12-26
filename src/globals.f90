 module globals
  use iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  !Global variables
  character (len=8)  :: cdate
  character (len=10) :: ctime
  integer(int64) :: seed1
  integer, parameter :: maxdist=20
  integer :: nloci, nind, nmix, numit, burnin, thin
  logical :: mcmc
  !Input
  logical :: useparfile, alterpheno
  character(len=200) :: job_name, parfil, genfil, phenfil, alterphenfil, grpfil, bimfil, inprefix, outprefix
  !Data
  character(len=80), dimension(:), allocatable :: animalID, snpID
  character(len=80) :: pheno_name
  integer :: pheno_ncol
  real(real64), dimension(:), allocatable :: why, pred, freqstore
  real(real64), target, dimension(:,:), allocatable :: X
  !Models
  integer BayesMethod
  integer, parameter :: BayesA=1, BayesB=2, BayesC=3, BayesCpi=4, BayesR=5, BayesRc=6
  real(real64) :: pizero, vara, vare, mu
  real(real64), dimension(:), allocatable :: gp, vare_gp, mix, alpha, g, yadj
  real(real64) :: scale_va, scale_ve, df_va, df_ve
  real(real64), dimension(:), allocatable :: xpx
  integer, dimension(:), allocatable :: permvec
  !Results and Outputs
  character(len=200) :: logfil, freqfil, mbvfil, hypfil, locfil, modfil, efffil, seedfil
  real(real64) :: mu_store, vare_store, vara_store, included_store 
  real(real64), dimension(:), allocatable :: gstore
  real(real64), dimension(:,:), allocatable :: indiststore
  !BayesRc
  integer, parameter :: maxgrp=100
  character(len=80) :: grp(maxgrp)
  integer ngroup
  integer, dimension(:,:), allocatable :: snpindist, snptracker
  real(real64), dimension(:,:), allocatable :: p, log_p, dirx, varindist
  real(real64), dimension(:,:), allocatable :: pstore, snpstore, varstore   
  !BayesA and BayesB
  real(real64), dimension(:), allocatable :: vara_s, vara_sstore
  !Soda related
  integer*1, parameter :: FIVEMASK8=z'55'
  integer(int64), parameter :: FIVEMASK64=z'5555555555555555'  
  logical :: SODAOFF, rblocks
  logical :: setblocks, setblocksize
  integer :: maxthreads, nthreads, nprocs, nblocks, blocksize, blockt, naind, nbatch, nskipped, naless
  logical, dimension(:), allocatable :: skip_j
  integer, dimension(:), allocatable :: i_skipped    
  integer(int64), allocatable :: genbits(:,:)
  real(real64), dimension(:), allocatable :: xxd, yA, yAadj
  real(real64), target, dimension(:,:), allocatable :: XA
end module globals
