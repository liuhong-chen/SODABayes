module soda
    use globals
    use ParallelRNG
    use utilities
    use omp_lib
    
contains
    
  subroutine choleskyPD(A, U, n, status)
  !choleskyPD computes the Cholesky factorization of a real symmetric
  !positive definite matrix A. Output is an upper triangular matrix U
  integer,      intent(in)    :: n
  real(real64),     intent(in)    :: A(n,n)
  real(real64),     intent(out)   :: U(n,n)
  integer,      intent(out)   :: status
  real(real64) :: tmp
  integer :: i,j

  !quick return if n<=0
  if(n<=0) then 
      status = -1
      return
  end if
  !check the positiveness of the diagonals
  do i=1,n
      if(A(i,i) <= 0.0d0) then
        status = -1
        return    
      end if
  end do
  
  status = 0
  U(:,:)=0.0d0
  do i = 1, n
      
     tmp = A(i,i) - dot_product(U(1:i-1,i),U(1:i-1,i))
     if(tmp <= 0.0d0 .or. isnan(tmp)) then
          status = i
          return
     end if
     U(i,i) = sqrt(tmp)
     do j = i+1, n
        U(i,j)  = ( A(j,i) - dot_product(U(1:i-1,j),U(1:i-1,i)) ) / U(i,i)
     end do
  end do
  
  return
  end subroutine choleskyPD

  subroutine choleskyPSD(A, piv, n, rank, status)
  !choleskyPSD computes the Cholesky factorization with complete
  !pivoting of a real symmetric positive semi-definite matrix A.
  !Output is an upper triangular matrix, rank, and a pivoting vector.
  !Reference: http://eprints.maths.manchester.ac.uk/689/1/para_paper.pdf
  integer, intent(in) :: n
  integer, intent(out) :: piv(n)
  integer, intent(out) :: rank, status
  real(real64), intent(inout) :: A(n,n)

  real(real64) :: Ajj, DSTOP, EPS, dtemp
  real(real64), allocatable :: work(:)
  integer :: itemp, i, j, jj, p, pvt

  allocate(work(2*n))
  status = 0
  if( n.LE.0 ) then
    status = -1
    return
  end if
  
  do p = 1, n
    piv(p) = p
  end do

  EPS = epsilon(1.0d0)
  
  call IDMAX( n, A(1, 1), n+1, pvt, dtemp)
  Ajj = A(pvt, pvt)
  if( Ajj.LE.EPS ) then
    rank = 0
    status = -1
    return
  end if

  DSTOP = n*EPS*Ajj

  do p = 1, n
    work(p) = 0
  end do

  do j = 1, n
    do p = j, n
      if( j.GT.1 ) then
        work(p) = work(p) + A(j-1, p)**2
      end if
      work(n+p) = A(p,p) - work(p)
    end do
    if(j.GT.1) then
      call IDMAX( n-j+1, work(n+j), 1, itemp, dtemp)
      pvt = itemp + j - 1
      Ajj = work(n+pvt)
      if(Ajj.LE.DSTOP) then
        A(j,j) = Ajj
        rank = j - 1
        return
      end if
    end if
    if(j.NE.pvt) then
      A(pvt, pvt) = A(j,j)
      call dswap( j-1, A(1,j), 1, A(1,pvt), 1)
      CALL dswap( pvt-j-1, A(j, j+1), n, A(j+1, pvt), 1)
      if(pvt<n) then
            CALL dswap( n-pvt, A(j, pvt+1), n, A(pvt, pvt+1), n)
      end if
         
      dtemp = work(j)
      work(j) = work(pvt)
      work(pvt) = dtemp
      itemp = piv(pvt)
      piv(pvt) = piv(j)
      piv(j) = itemp
    end if
    
    Ajj = sqrt( Ajj )
    A(j,j) = Ajj

    if( j.LT.n ) then
      call DGEMV( 'Trans', j-1, n-j, -1.0d0, A( 1, j+1 ), N, A( 1, j ), 1, 1.0d0, A( j, j+1 ), n )
      call DSCAL( n-j, 1.0 / Ajj, A( j, j+1 ), n )
    end if
    
  end do
  rank = n

  deallocate(work)
  return
end subroutine choleskyPSD
!subroutines for cholesky decomposition
!adapted from lapack packages      
subroutine dswap(n,dx,incx,dy,incy)
       integer incx,incy,n
       real(real64) dx(*), dy(*)

       real(real64) DTEMP
       integer I,IX,IY,M,MP1

       if (n.LE.0) return
       if (incx.EQ.1 .AND. incy.EQ.1) then
 
          m = mod(n,3)
          if (m.NE.0) then
             do i = 1,m
                dtemp = dx(i)
                dx(i) = dy(i)
                dy(i) = dtemp
             end do
             if (n.LT.3) return
          end if
          mp1 = m + 1
          do i = mp1,n,3
             dtemp = dx(i)
             dx(i) = dy(i)
             dy(i) = dtemp
             dtemp = dx(i+1)
             dx(i+1) = dy(i+1)
             dy(i+1) = dtemp
             dtemp = dx(i+2)
             dx(i+2) = dy(i+2)
             dy(i+2) = dtemp
          end do
       else
          ix = 1
          iy = 1
          if (incx.LT.0) ix = (-n+1)*incx + 1
          if (incy.LT.0) iy = (-n+1)*incy + 1
          do i = 1,n
             dtemp = dx(ix)
             dx(ix) = dy(iy)
             dy(iy) = dtemp
             ix = ix + incx
             iy = iy + incy
          end do
       end if
       return

end subroutine dswap

subroutine IDMAX( n, X, incx, k, r )
!IDMAX finds the largest component of x, r, and
!determines the smallest index, k, such that x(k) = r.

  real(real64) :: r
  integer :: incx, k, n

  real(real64) X(*)
  integer i, ix

  k = 0
  if( n.LT.1 .OR. incx.LE.0 ) return
  k = 1
  if( n.EQ.1 ) return

  if( incx.EQ.1 ) then
    r = X( 1 )
    do i = 2, N
      if( X( i ).GT.r ) then     
          k = i
          r = X( i )
      end if
    end do
  else
    r = X( 1 )
    ix = 1
    do i = 2, n
        ix = ix + incx
        if( X( ix ).GT.r ) then
          k = i
          r = X( ix )
        end if
    end do
  end if
  
  return
end subroutine IDMAX

subroutine dgemv(trans, m, n, alpha, A, LDA, x, incx, beta, y, incy)
       real(real64) :: alpha, beta
       integer :: incx, incy, LDA, m, n
       character :: trans
       real(real64) :: A(LDA,*), x(*), y(*)
       
       real(real64) :: temp
       integer :: i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
 
       info = 0
       if (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND. .NOT.lsame(trans,'C')) return
       if (m.LE.0 .OR. n.LE.0 .OR. LDA.LT.max(1,m) .OR. incx.EQ.0 .OR. incy.EQ.0) return
       if ((alpha.EQ.0.0d0) .AND. (beta.EQ.1.0d0)) return

       if (lsame(trans,'N')) then
           lenx = n
           leny = m
       else
           lenx = m
           leny = n
       end if
       if (incx.GT.0) then
           kx = 1
       else
           kx = 1 - (lenx-1)*incx
       end if
       if (incy.GT.0) then
           ky = 1
       else
           ky = 1 - (leny-1)*incy
       end if

 !     First form  y := beta*y.
 
       if (beta.NE.1.0d0) then
           if (incy.EQ.1) then
               if (beta.EQ.0.0d0) then
                   do i = 1,leny
                       y(i) = 0.0d0
                   end do
               else
                   do i = 1,leny
                       y(i) = beta*y(i)
                   end do
               end if
           else
               iy = ky
               if (beta.EQ.0.0d0) then
                   do  i = 1,leny
                       y(iy) = 0.0d0
                       iy = iy + incy
                   end do
               else
                   do  i = 1,leny
                       y(iy) = beta*y(iy)
                       iy = iy + incy
                   end do
               end if
           end if
       end if
       if (alpha.EQ.0.0d0) return
       if (lsame(trans,'N')) then
 
 !        Form  y := alpha*A*x + y.
 
           jx = kx
           if (incy.EQ.1) then
               do  j = 1,n
                   temp = alpha*x(jx)
                   do  i = 1,m
                       y(i) = y(i) + temp*a(i,j)
                   end do
                   jx = jx + incx
               end do
           else
               do  j = 1,n
                   temp = alpha*x(jx)
                   iy = ky
                   do  i = 1,m
                       y(iy) = y(iy) + temp*a(i,j)
                       iy = iy + incy
                   end do
                   jx = jx + incx
               end do
           end if
       else
 
 !        Form  y := alpha*A**T*x + y.
 
           jy = ky
           if (incx.EQ.1) then
               do  j = 1,n
                   temp = 0.0d0
                   do  i = 1,m
                       temp = temp + a(i,j)*x(i)
                   end do
                   y(jy) = y(jy) + alpha*temp
                   jy = jy + incy
               end do
           else
               do  j = 1,n
                   temp = 0.0d0
                   ix = kx
                   do  i = 1,m
                       temp = temp + a(i,j)*x(ix)
                       ix = ix + incx
                   end do
                   y(jy) = y(jy) + alpha*temp
                   jy = jy + incy
               end do
           end if
       end if
 
       return
end subroutine dgemv

subroutine dscal(n,da,dx,incx)
       real(real64) DA
       integer INCX,N
 
       real(real64) DX(*)
 
       integer i,m,mp1,nincx
 
       if (n.LE.0 .OR. incx.LE.0) return
       if (incx.EQ.1) then
          m = mod(n,5)
          if (m.NE.0) then
             do i = 1,m
                dx(i) = da*dx(i)
             end do
             if (n.LT.5) return
          end if
          mp1 = m + 1
          do i = mp1,n,5
             dx(i) = da*dx(i)
             dx(i+1) = da*dx(i+1)
             dx(i+2) = da*dx(i+2)
             dx(i+3) = da*dx(i+3)
             dx(i+4) = da*dx(i+4)
          end do
       else
          nincx = n*incx
          do i = 1,nincx,incx
             dx(i) = da*dx(i)
          end do
       end if
       return
end subroutine DSCAL

logical function lsame(CA,CB)

       character ca,cb 
       integer inta,intb,zcode
       lsame = ca .EQ. cb
       
       if (lsame) return
 
       zcode = ichar('Z')
       inta = ichar(ca)
       intb = ichar(cb)

       if (zcode.EQ.90 .OR. zcode.EQ.122) then
           if (inta.GE.97 .AND. inta.LE.122) inta = inta - 32
           if (intb.GE.97 .AND. intb.LE.122) intb = intb - 32
       else if (zcode.EQ.233 .OR. zcode.EQ.169) then
           if (inta.GE.129 .AND. inta.LE.137 .OR. &
               inta.GE.145 .AND. inta.LE.153 .OR. &
               inta.GE.162 .AND. inta.LE.169) inta = inta + 64
           if (intb.GE.129 .AND. intb.LE.137 .OR. &
               intb.GE.145 .AND. intb.LE.153 .OR. &
               intb.GE.162 .AND. intb.LE.169) intb = intb + 64
       else if (zcode.EQ.218 .OR. zcode.EQ.250) then
           if (inta.GE.225 .AND. inta.LE.250) inta = inta - 32
           if (intb.GE.225 .AND. intb.LE.250) intb = intb - 32
       end if
       lsame = inta .EQ. intb
end function lsame
    
subroutine power_method ( n, A, lambda, p)
!subroutine for finding the dominant eigenvalue of matrix A

  integer n, p
  real (real64) ::  A(:,:)
  real (real64) lambda, delta

  real (real64), allocatable :: y(:), v(:)
  integer :: i, maxit
  real(real64), parameter :: tol = 1.0e-14
  maxit = 1000*n
  allocate(y(n),v(n))
  !initialize
  call random(n, y, p)
  do i = 1, maxit
    v = y / sqrt ( sum ( y**2 ) )
    y = matmul ( A, v )
    lambda = dot_product ( y, v )
    delta = sum((y-lambda*v)**2) / abs(lambda)
    if ( delta <= tol ) exit
  end do

  lambda = abs(lambda)
  deallocate(y,v)
  return
end subroutine power_method

subroutine setup_soda
  integer k, ios
  if(SODAOFF) return

  if(setblocksize) then
      nblocks = ceiling(dble(nloci) / blocksize)
      blockt = nblocks
  else if(setblocks) then
      blocksize = ceiling(dble(nloci) / nblocks)
      blockt = nloci - (blocksize-1) * nblocks
  else
      blocksize = nthreads
      nblocks = ceiling(dble(nloci) / blocksize)  
      blockt = nblocks      
  end if
  if(blocksize <= 1) then
      SODAOFF = .true.
      return
  end if
  naind = blocksize
  allocate(xxd(nblocks), XA(naind,nloci), yA(naind), yAadj(naind), skip_j(nloci), i_skipped(nloci),stat=ios)
  if( ios /= 0 )then
     write(21,'(a50)') 'ERROR :: Unable to allocate required storage for pseudo data'
     call flush(21)
     stop 'Unable to allocate required storage for pseudo data'
  end if
  if(rblocks) then
    call permutate(permvec,nloci,blocksize)
  else
     do k=1,nloci
      permvec(k)=k  
     enddo
  endif  

  nbatch = ceiling(dble(nind)/(int64*4))
  allocate(genbits(nbatch,nloci))
  call load_geno_bits
  call get_aug_matrix
  naind = naind-naless
  deallocate(genbits)
  
  return
end subroutine setup_soda  

subroutine get_aug_matrix
  integer i, j, k, l, jj, kk, ls, le, n, rank, status, ios, thread
  integer, allocatable :: piv(:)
  real(real64) d, temp, J0
  real(real64), allocatable :: xxtemp(:,:)
  
  allocate(xxtemp(blocksize,blocksize),piv(blocksize),stat=ios)
  if( ios /= 0 )then
      write(21,'(a50)') 'ERROR :: Unable to allocate enough memory for cholesky decomposition'
      call flush(21)
      stop 'Unable to allocate enough memory for cholesky decomposition'
  end if
  skip_j=.false.
  i_skipped=0
  nskipped=0  
  naless=nloci
  XA = 0.0d0
  do i = 1, nblocks      
    if(i<=blockt) then  
        ls = (i-1)*blocksize+1
        le = min(ls+blocksize-1,nloci)
    else 
        ls = blocksize*blockt+(i-blockt-1)*(blocksize-1)+1
        le = min(ls+blocksize-2,nloci)
    end if 
    n = le - ls + 1
    xxtemp = 0.0d0

    !$omp parallel private(jj,kk) 
    !$omp do
    do j=1, n
        jj = permvec(ls + j - 1)
        do k=1, j
            kk = permvec(ls + k - 1)
            xxtemp(k,j) = xxproduct(jj, kk)
            xxtemp(j,k) = xxtemp(k,j)
        end do
    end do
    !$omp end do
    !$omp end parallel
    call power_method(n, xxtemp(1:n,1:n), d, 0)
    d = d + 0.001
    xxd(i) = d
    xxtemp(:,:) = -xxtemp(:,:)
    do j=1, n
        xxtemp(j,j) = d + xxtemp(j,j)
    end do
    call choleskyPD(xxtemp(1:n,1:n), XA(1:n,ls:le),n,status)
    if(status /=0) then 
        write(21,'(a80)') 'Warnimg: D-XX is not positive definite'
        write(21,'(a80)') '         switching to pivoted cholesky decomposition'
        call choleskyPSD(xxtemp(1:n,1:n), piv, n, rank, status)   
        if(status /= 0) then
            write(21,'(a80)',advance='no') 'Error: unable to construct an augmented matrix for X'
            write(21,'(a80)',advance='no') '       (D-XX is not positive semidefinite)'
            stop 'Unable to construct an augmented matrix for X'
        end if     
        do j=1,n
            jj = piv(j)
            if(j<=rank) then    
                XA(1:j,ls+jj-1) = xxtemp(1:j,j)
            else
                skip_j(ls+jj-1) = .true.
                nskipped = nskipped + 1             
                i_skipped(nskipped) = permvec(ls+jj-1)
            end if            
        end do
     else
         naless = 0
    end if
     naless = min(naless,blocksize-n+nskipped)
  end do 

  deallocate(xxtemp)  
  return
end subroutine get_aug_matrix


function xxproduct(jj,kk) result(res)
    integer jj,kk
    real(real64) res
    integer(int64) :: a, b, c, d, e, f, g, v
    integer i, nc, nd, ng
    real(real64) :: pj, pk, qj, qk, ep

    pj = freqstore(jj)
    pk = freqstore(kk)
    qj = 1.0d0 - pj
    qk = 1.0d0 - pk
    ep = epsilon(1.0d0)
    if(pj.le.ep.or.qj.le.ep.or.pk.le.ep.or.qk.le.ep) then
        res= 0.0d0
        return
    end if
    pj = pj * 2.0d0
    pk = pk * 2.0d0
    nc = 0
    nd = 0
    ng = 0
    do i=1,nbatch
        a = genbits(i,jj)
        b = genbits(i,kk)
        c = iand(a,b)
        d = iand(shiftr(c,1),FIVEMASK64)
        e = iand(iand(shiftr(a,1),b),FIVEMASK64)
        f = iand(iand(a,shiftr(b,1)),FIVEMASK64)
        g = ior(e,f)    
        nc = nc + popcnt(c)
        nd = nd + popcnt(d)
        ng = ng + popcnt(g)
    end do
    res = nd*3 + ng*2 + nc
    res = (res - nind*pj*pk) / sqrt(pj*pk*qj*qk)
    return
end function xxproduct

end module soda