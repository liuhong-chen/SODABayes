
module ParallelRNG
!------------------------------------------------------------------------------------ 
!This module contains a parallel random number generator and random variables from
!various distributions
! 
! Each thread uses its own pseudo-random sequence.
! The method is described in Marsaglia, G. & Tsang, W.W. (2000)
! `The ziggurat method for generating random variables', J. Statist. Software, v5(8).
!------------------------------------------------------------------------------------
  use globals, only: int32, int64, real32, real64
  implicit none

  private :: rng_xoroshiro128plus, rng_jump

  public :: random, rand_set_seed, rand_uniform, rand_normal, rand_exponential, rand_gamma, rand_chi_square
  public :: rand_scaled_inverse_chi_square, rand_inverse_gamma, rand_weibull, rand_cauchy, rand_student_t
  public :: rand_laplace, rand_log_normal, rand_beta, rand_dirichlet
  
  integer :: rng_n
  real(real64), parameter :: PI=3.141592653589793238462
  type state
      integer(int64) :: s1
      integer(int64) :: s2
  end type state
  type(state), save, allocatable :: ss(:)
  
  interface random
    module procedure rng_uni
    module procedure rng_uni_array
    module procedure rng_uni_32
    module procedure rng_uni_ser
    module procedure rng_uni_32_ser
  end interface

  interface rand_set_seed
    module procedure rng_zigset
  end interface
  
  interface rng_xoroshiro128plus
    module procedure rng_xoroshiro128plus
    module procedure rng_xoroshiro128plus_kpar
    module procedure rng_xoroshiro128plus_kpar_32
  end interface
  
  interface rng_jump
    module procedure rng_jump
    module procedure rng_jump_n
  end interface
   

contains


  subroutine rng_zigset(npar, rng_jsrseed)

    integer, intent(in)  :: npar
    integer(int64), intent(in)  :: rng_jsrseed(2)

    type(state) :: s
    integer :: i, kpar
    real(real64) dn, tn, de, te

    rng_n = npar

    allocate(ss(0:npar-1))
    do i = 0, npar-1
        ss(i)%s1 = 0.0d0
        ss(i)%s2 = 0.0d0
    end do
    s%s1 = rng_jsrseed(1)
    s%s2 = rng_jsrseed(2)

    do kpar = 0,npar-1
        call rng_jump(s)
        ss(kpar)%s1 = s%s1
        ss(kpar)%s2 = s%s2        
    enddo
  end subroutine rng_zigset

  ! Generate random 64-bit integers
  subroutine rng_xoroshiro128plus_kpar(ival, kpar)
    integer(int64), intent(out)   :: ival
    integer, intent(in) :: kpar
        
    call rng_xoroshiro128plus(ival, ss(kpar))
  end subroutine

  ! Generate random 32-bit integers by taking the *higher* 32-bits of the result
  subroutine rng_xoroshiro128plus_kpar_32(ival, kpar)
    integer(int32), intent(out)   :: ival
    integer, intent(in) :: kpar
    integer(int64) :: tmp

    call rng_xoroshiro128plus(tmp, ss(kpar))
    
    ival = transfer(shiftr(tmp,32), ival)
  end subroutine

  ! Generate random 64-bit integers
  subroutine rng_xoroshiro128plus(ival, s) 
    integer(int64), intent(out)   :: ival
    type(state), intent(inout) :: s
    integer(int64) :: s1, s2
    
    s1 = s%s1
    s2 = s%s2
        
    ival = s1 + s2
    s2 = ieor(s2, s1)
    s1 = ieor( ieor(rotl(s1, 55), s2), shiftl(s2, 14))
    s2 = rotl(s2, 36)
 
    s%s1 = s1
    s%s2 = s2

  contains
    function rotl(x, k)
      integer(int64) :: rotl, x
      integer :: k

      rotl = ior(shiftl(x,k), shiftr(x, 64-k))
    end function
  end subroutine

  ! Jump by 2^64 steps
  subroutine rng_jump(s)
    type(state), intent(inout) :: s
    ! The first constant is df900294d8f554a5 using not() to avoid overflow
    integer(int64), parameter :: jmp(2) = (/-4707382666127344949_int64, &
                                            -2852180941702784734_int64/)
    integer(int64) :: s1, s2, dummy
    
    integer :: i, b
    
    s1 = 0
    s2 = 0
    
    do i = 1, 2
      do b = 0, 63
        if (iand(jmp(i), shiftl(1_int64,b))/=0) then
          s1 = ieor(s1, s%s1)
          s2 = ieor(s2, s%s2)
        end if  
        call rng_xoroshiro128plus(dummy, s)
      end do
    end do
    s%s1 = s1
    s%s2 = s2
  end subroutine
  
  ! Do n jumps defined above
  subroutine rng_jump_n(s, n)
    type(state), intent(inout) :: s
    integer, intent(in) :: n
    integer :: i
    
    do i = 1, n
      call rng_jump(s)
    end do
  end subroutine

  !  Generate uniformly distributed random numbers [0,1), sequence kpar
  subroutine rng_uni(fn_val, kpar)
    integer :: kpar
    real(real64), intent(out) ::  fn_val
    integer(int64) :: x

    if (kpar >= rng_n) then
        write(*,*) 'RNG Error: thread number',kpar, 'exceeds initialized max: ',rng_n-1
        write(*,*) '           or RNG not initialized.'
        stop
    endif
    
    call rng_xoroshiro128plus(x, kpar)
    
    fn_val = rng_int64_to_real64(x)
  end subroutine rng_uni
  
  subroutine rng_uni_array(n, fn_val, kpar)
    integer :: kpar, n, i
    real(real64), intent(out) ::  fn_val(n)
    integer(int64) :: x

    if (kpar >= rng_n) then
        write(21,*) 'RNG Error: thread number',kpar, 'exceeds initialized max: ',rng_n-1
        write(21,*) '           or RNG not initialized.'
        stop
    endif
    do i=1, n
        call rng_xoroshiro128plus(x, kpar)   
        fn_val(i) = rng_int64_to_real64(x)
    end do
    
  end subroutine rng_uni_array
  subroutine rng_uni_32(fn_val, kpar)
    integer :: kpar
    real(real32), intent(out) ::  fn_val
    integer(int64) :: x

    if (kpar >= rng_n) then
        write(*,*) 'RNG Error: thread number',kpar, 'exceeds initialized max: ',rng_n-1
        write(*,*) '           or RNG not initialized.'
        stop
    endif
    
    call rng_xoroshiro128plus(x, kpar)
    
    fn_val = rng_int64_to_real32(x)
  end subroutine rng_uni_32

  subroutine rng_uni_ser(x)
    real(real64), intent(out) :: x

    call rng_uni(x, 0)
  end subroutine
 
  subroutine rng_uni_32_ser(x)
    real(real32), intent(out) :: x
    real(real64) :: x64
    call rng_uni(x64, 0)
    x = real(x64, real32)
  end subroutine

  function rng_int32_to_real32(i) result(res)
    real(real32) :: res
    integer(int32), value :: i

    i  = ior(int(Z'3F800000',int32), shiftr(i, 9))
    res = transfer(i, 1.0_real32) - 1;
  end function
 
  function rng_int64_to_real32(i) result(res)
    real(real32) :: res
    integer(int64), value :: i
    integer(int32) :: tmp

    tmp = transfer(shiftr(i,32), tmp)
    tmp  = ior(int(Z'3F800000',int32), shiftr(tmp, 9))
    res = transfer(tmp, 1.0_real32) - 1;
  end function
 
  function rng_int64_to_real64(i) result(res)
    real(real64) :: res
    integer(int64), value :: i
    integer(int32) :: tmp
    
    i  = ior(int(Z'3FF0000000000000',int64), shiftr(i, 12))
    res = transfer(i, 1.0_real64) - 1;
  end function 
  
!Generate random variables from specific distributions
!
! Random Sample from a uniform distribution
!
  function rand_uniform(a,b,kpar) result(c)
    real(real64) :: a,b,c,temp
    integer :: kpar
    call random(temp, kpar)
    c= a+temp*(b-a)
  end function rand_uniform
!
! Random Sample from normal (Gaussian) distribution
!
  function rand_normal(mean,stdev,kpar) result(c)
    real(real64) :: mean,stdev,c,r,theta,temp(2)
    integer :: kpar   
    if(stdev <= 0.0d0) then
       !Write(*,*) "Standard Deviation must be +ve"
       c=mean
    else
       call random(temp(1),kpar)
       call random(temp(2),kpar)
       r=(-2.0d0*log(temp(1)))**0.5d0
       theta = 2.0d0*PI*temp(2)
       c= mean+stdev*r*sin(theta)
    end if
  end function rand_normal
!
!  Random smaple from an exponential distribution
!
  function rand_exponential(mean, kpar) result(c)
    real(real64) :: mean,c,temp
    integer :: kpar
    if (mean <= 0.0d0) then
       write(*,*) "mean must be positive"
    else
       call random(temp, kpar)
       c=-mean*log(temp)
    end if
  end function rand_exponential
!
!  Return a random sample from a gamma distribution
!
  recursive function rand_gamma(shape, scale, kpar) result(ans)
    real(real64) shape,scale,ans,u,w,d,c,x,xsq,g,v
    integer :: kpar
    if (shape <= 0.0d0) then
       write(*,*) "Shape parameter must be positive"
    end if
    if (scale <= 0.0d0) then
       write(*,*) "Scale parameter must be positive"
    end if
!
!    ## Implementation based on "A Simple Method for Generating Gamma Variables"
!    ## by George Marsaglia and Wai Wan Tsang.
!    ## ACM Transactions on Mathematical Software
!    ## Vol 26, No 3, September 2000, pages 363-372.
!
    if (shape >= 1.0d0) then
       d = shape - 1.0d0/3.0d0
       c = 1.0d0/(9.0d0*d)**0.5d0
       do while (.true.)
          x = rand_normal(0.0d0, 1.0d0, kpar)
          v = 1.0d0 + c*x
          do while (v <= 0.0d0)
             x = rand_normal(0.0d0, 1.0d0, kpar)
             v = 1.0d0 + c*x
          end do
          v = v*v*v
          call random(u, kpar)
          xsq = x*x
          if ((u < 1.0d0 -.0331d0*xsq*xsq) .or. (dlog(u) < 0.5d0*xsq + d*(1.0d0 - v + dlog(v))) )then
             ans=scale*d*v
             return
          end if
       end do
    else
        g = rand_gamma(shape+1.0d0, 1.0d0, kpar)
        call random(w, kpar)
        ans=scale*g*(w)**(1.0d0/shape)
        return
     end if
   end function rand_gamma
!
! ## return a random sample from a chi square distribution
! ## with the specified degrees of freedom
!
   function rand_chi_square(dof, kpar) result(ans)
     real(real64) ans,dof
     integer :: kpar
     ans=rand_gamma(0.5d0*dof, 2.0d0, kpar)
   end function rand_chi_square
!
! ## return a random sample from a scaled 
!    inverse chi square distribution with
!    df and scale parameter
!
   function rand_scaled_inverse_chi_square(dof,scale, kpar) result(ans)
     real(real64) ans,dof,scale
     integer :: kpar
     ans=rand_inverse_gamma(dble(0.5)*dof, dble(0.5)*dof*scale, kpar)
   end function rand_scaled_inverse_chi_square

! ## return a random sample from an inverse gamma random variable
!
   function rand_inverse_gamma(shape, scale, kpar) result(ans)
     real(real64) shape,scale,ans
     integer :: kpar
!    ## If X is gamma(shape, scale) then
!    ## 1/Y is inverse gamma(shape, 1/scale)
     ans= 1.0d0 / rand_gamma(shape, 1.0d0 / scale, kpar)
   end function rand_inverse_gamma
!
!## return a sample from a Weibull distribution
!
   function rand_weibull(shape, scale, kpar) result(ans)
     real(real64) shape,scale,temp,ans
     integer :: kpar
     if (shape <= 0.0d0) then
        write(*,*) "Shape parameter must be positive"
     end if
     if (scale <= 0.0d0) then
        write(*,*) "Scale parameter must be positive"
     end if
     call random(temp, kpar)
     ans= scale * (-log(temp))**(1.0d0 / shape)
   end function rand_weibull
!
!## return a random sample from a Cauchy distribution
!
   function rand_cauchy(median, scale, kpar) result(ans)
     real(real64) ans,median,scale,p
     integer :: kpar
     if (scale <= 0.0d0) then
        write(*,*) "Scale parameter must be positive"
     end if
     call random(p, kpar)
     ans = median + scale*tan(PI*(p - 0.5d0))
   end function rand_cauchy
!
!## return a random sample from a Student t distribution
!
   function rand_student_t(dof, kpar) result(ans)
     real(real64) ans,dof,y1,y2
     integer :: kpar
     if (dof <= 0.d0) then
        write(*,*) "Degrees of freedom must be positive"
     end if
!
! ## See Seminumerical Algorithms by Knuth
      y1 = rand_normal(0.0d0, 1.0d0, kpar)
      y2 = rand_chi_square(dof, kpar)
      ans= y1 / (y2 / dof)**0.50d0
!
    end function rand_student_t
!
!## return a random sample from a Laplace distribution
!## The Laplace distribution is also known as the double exponential distribution.
!
    function rand_laplace(mean, scale, kpar)  result(ans)
      real(real64) ans,mean,scale,u
      integer :: kpar
      if (scale <= 0.0d0) then
        write(*,*) "Scale parameter must be positive"
     end if
     call random(u, kpar)
     if (u < 0.5d0) then
        ans = mean + scale*log(2.0d0*u)
     else
        ans = mean - scale*log(2.0d0*(1.0d0-u))
     end if
   end function rand_laplace
!
! ## return a random sample from a log-normal distribution
!
   function rand_log_normal(mu, sigma, kpar) result(ans)
     real(real64) ans,mu,sigma
     integer :: kpar
     ans= exp(rand_normal(mu, sigma, kpar))
   end function rand_log_normal
!
! ## return a random sample from a beta distribution
!
   function rand_beta(a, b, kpar) result(ans)
     real(real64) a,b,ans,u,v
     integer :: kpar
     if ((a <= 0.0d0) .or. (b <= 0.0d0)) then
        write(*,*) "Beta parameters must be positive"
     end if
     u = rand_gamma(a, 1.0d0, kpar)
     v = rand_gamma(b, 1.0d0, kpar)
     ans = u / (u + v)
   end function rand_beta

   function rand_dirichlet(n,irx,kpar) result(x)
     integer :: n, i, kpar
     real(real64) :: sx
     real(real64), dimension(n) :: irx, x
     
     do i=1,n
        x(i)=rand_gamma(irx(i),1.0d0, kpar)
     enddo
     sx=sum(x)
     x=x/sx
   end function rand_dirichlet
   
end module ParallelRNG  
