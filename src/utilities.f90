
module utilities
use globals
use ParallelRNG
implicit none

    contains
subroutine to_upper(c)
    character(len=*), intent(inout) :: c
    integer i, l, intc, zcode
    zcode = ichar('Z')
    l = len(c)
    do i=1,l
       intc = ichar(c(i:i))

       if (zcode.EQ.90 .OR. zcode.EQ.122) then
 !        ASCII is assumed - ZCODE is the ASCII code of either lower or
 !        upper case 'Z'.           
           if (intc.GE.97 .AND. intc.LE.122) intc = intc - 32
       else if (zcode.EQ.233 .OR. zcode.EQ.169) then
 !        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
 !        upper case 'Z'.           
           if (intc.GE.129 .AND. intc.LE.137 .OR. &
               intc.GE.145 .AND. intc.LE.153 .OR. &
               intc.GE.162 .AND. intc.LE.169) intc = intc + 64
       else if (zcode.EQ.218 .OR. zcode.EQ.250) then
 !        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
 !        plus 128 of either lower or upper case 'Z'           
           if (intc.GE.225 .AND. intc.LE.250) intc = intc - 32
       end if      
       c(i:i) = achar(intc)
    end do
    return
end subroutine to_upper

subroutine to_lower(c)
    character(len=*), intent(inout) :: c
    integer i, l, intc, zcode
    zcode = ichar('Z')
    l = len(c)
    do i=1,l
       intc = ichar(c(i:i))

       if (zcode.EQ.90 .OR. zcode.EQ.122) then
 !        ASCII is assumed - ZCODE is the ASCII code of either lower or
 !        upper case 'Z'.           
           if (intc.GE.97 .AND. intc.LE.122) intc = intc + 32
       else if (zcode.EQ.233 .OR. zcode.EQ.169) then
 !        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
 !        upper case 'Z'.           
           if (intc.GE.129 .AND. intc.LE.137 .OR. &
               intc.GE.145 .AND. intc.LE.153 .OR. &
               intc.GE.162 .AND. intc.LE.169) intc = intc - 64
       else if (zcode.EQ.218 .OR. zcode.EQ.250) then
 !        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
 !        plus 128 of either lower or upper case 'Z'           
           if (intc.GE.225 .AND. intc.LE.250) intc = intc + 32
       end if      
       c(i:i) = achar(intc)
    end do
    return
end subroutine to_lower

logical function ssmatch(sa, sb)
    character(len=*), intent(in) :: sa, sb
    character(len=:), allocatable :: s1, s2
    integer n1,n2, comp

    n1=len_trim(sa)
    n2=len_trim(sb)
    allocate(character(len=n1) :: s1)
    allocate(character(len=n2) :: s2)
    s1 = trim(sa)
    s2 = trim(sb)

    call to_upper(s1)
    call to_upper(s2)

    comp = index(trim(s1),trim(s2))*index(trim(s2),trim(s1))

    if(comp /=0) then
        ssmatch = .true.
    else 
        ssmatch = .false.
    end if
    deallocate(s1, s2)
    
end function ssmatch

subroutine left_trim(line)
    character(len=*), intent(inout) :: line
    integer i
    i=1
    do while(i<=len_trim(line))
        if(line(i:i) /= ' ') then
            line = line(i:len_trim(line))
            exit
        end if
        i = i+1
    end do
end subroutine left_trim

  
  logical function str_match(str1,str2)
    implicit none
    character(*) :: str1, str2
    integer :: comp
    !exact match
    comp = index(trim(str1),trim(str2))*index(trim(str2),trim(str1))
    str_match = .false.
    if (comp /= 0) str_match = .true.
  end function str_match

  integer function nfields(line, sep)
    character(len=*), intent(in) :: line
    character(len=*), intent(in) :: sep
    integer i, n

    i = 1;
    n = len_trim(line)
    nfields = 0
    do while(i <= n)
       do while(scan(sep,line(i:i)) /= 0)
          i = i + 1
          if (n < i) return
       enddo
       nfields = nfields + 1
       do
          i = i + 1
          if (n < i) return
          if (scan(sep,line(i:i)) /= 0) exit
       enddo
    enddo
  end function nfields
  
  subroutine get_fields(str,sep,n,val)
    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: sep
    character(len=*),dimension(*), intent(out) :: val
    integer, intent(out) :: n
    integer :: i, pos1, pos2, ln

    ln=len_trim(str)
    i=1
    n=0
    do
       do while(scan(sep,str(i:i)) /= 0)
          i = i + 1
          if (ln < i) return
       enddo
       pos1=i
       pos2 = scan(trim(str(pos1:)), sep)
       if (pos2 == 0) then
          n=n+1
          val(n) = str(pos1:)
          exit
       endif
       n=n+1
       val(n) = str(pos1:pos1+pos2-2)
       i = pos2+pos1
    enddo
  end subroutine get_fields

  integer function cast_int(value)
    character(len=*) :: value
    read(value,*)cast_int
  end function cast_int

  real function cast_float(value)
    character(len=*) :: value
    read(value,*)cast_float
  end function cast_float

  logical function cast_logical(value)
    character(len=*) :: value
    cast_logical=.false.
    if(str_match(trim(value),'t')) cast_logical=.true.
  end function cast_logical

  logical function is_key(key)
    character(len=*) :: key
    character(len=2),parameter :: str='-'
    is_key=.false.
    if(str_match(trim(key(1:1)),str)) then
        if(.not. is_numeric(trim(key(2:2)))) then
            is_key=.true.
        end if
    end if
  end function is_key

  logical function is_numeric(key)
    character(len=*) :: key
    integer :: value,ios
    is_numeric=.true.
    read(key,*,iostat=ios) value
    if(ios /=0)is_numeric=.false.
  end function is_numeric

subroutine write_log
   character(len=8) :: bmethod(6)=(/'BayesA  ', 'BayesB  ', 'BayesC  ','BayesCpi', 'BayesR  ','BayesRc '/)
   if(mcmc) then
       write(21,902) 'Prefix for input files',trim(inprefix)
       write(21,903) 'No. of loci',nloci
       write(21,903) 'No. of individuals',nind
       if(alterpheno) then
           write(21,902) 'Alternative phenotype file', alterphenfil
           if(pheno_name /= "") then
               write(21,902) 'Phenotype name', pheno_name
           else
               write(21,903) 'Phenotype colunm', pheno_ncol
           end if
       end if     
       write(21,906) 'Prior Vara(Scale DF)', scale_va, df_va
       write(21,906) 'Prior Vare(Scale DF)', scale_ve, df_ve
       write(21,903) 'No. of MCMC cycles',numit
       write(21,903) 'Burnin ',burnin
       write(21,903) 'Thinning rate',thin
       write(21,902) 'Bayesian method ', trim(bmethod(BayesMethod))
       if(BayesMethod==BayesB .or. BayesMethod==BayesC) then
           write(21,905) 'pi (%SNP with zero effects)', pizero
       else if(BayesMethod == BayesR .or. BayesMethod==BayesRc) then
           write(21,903) 'No. of mixtures',nmix         
           write(21,905) 'Variance of dist ', mix
           write(21,905) 'Dirichlet prior', alpha
       end if     
       if(SODAOFF .eqv. .false.) then
           write(21,908) 'SODA on', .not.SODAOFF                      
           write(21,903) 'No. of threads', nthreads
           write(21,903) 'No. of blocks', nblocks
           write(21,903) 'Block size', blocksize
           write(21,908) 'Random blocks', rblocks 
       else 
            write(21,908) 'SODA off', SODAOFF
       end if
       write(21,901) 'Seed ', seed1  
       write(21,902) 'Prefix for output files',trim(outprefix)        
   else
       write(21,902) 'Prefix for input files',trim(inprefix)
       write(21,903) 'No. of loci',nloci
       write(21,903) 'No. of individuals',nind
       write(21,902) 'model file',trim(modfil)
       write(21,902) 'freq file',trim(freqfil)
       write(21,902) 'effects file',trim(efffil)       
       write(21,902) 'Prefix for output files',trim(outprefix)       
   end if
   call flush(21)
901 format(a,t30,'= ',i12)   
902 format(a,t30,': ',a)
903 format(a,t30,'= ',i8)
904 format(a,t30,'= ',f20.6)
905 format(a,t30,'= ',10f10.5)
906 format(a,t30,'= ',2f10.6)
907 format(a,t30,'= ',f10.2,a)
908 format(a,t30,'= ',l)
end subroutine write_log

subroutine get_size
  character(len=80) :: dum,grpinfo
  integer i, ngrp, ios
  logical newgrp
  nind=0
  open(35,file=trim(phenfil),status='old',form='formatted')
  do
     read(35,*,iostat=ios) dum
     if (ios.ne.0) exit
     nind=nind+1
  enddo
  close(35,status='keep')
  nloci=0
  open(36,file=trim(bimfil),status='old',form='formatted')
  do
     read(36,*,iostat=ios) dum
     if (ios.ne.0) exit
     nloci=nloci+1          
  enddo
  close(36,status='keep')
  ngroup=1
  if(BayesMethod==BayesRc .and. mcmc) then
      open(42,file=trim(grpfil),status='old',form='formatted')
      ngrp=0
    do
      read(42,*,iostat=ios) dum, grpinfo
      if(ios.ne.0) exit
      newgrp=.true.
      do i=1,ngrp
          if(str_match(trim(grpinfo),grp(i))) then
              newgrp=.false.
              exit
          end if
      end do
      if(newgrp) then
          ngrp = ngrp+1
          if(ngrp>maxgrp) then
              print *, 'error: max no. of SNP groups: ', maxgrp
              stop
          end if
          grp(ngrp)=trim(grpinfo)
      end if
    end do
    ngroup=ngrp
    close(42, status='keep')
  end if             
end subroutine get_size

subroutine load_snp_grp
    integer i,j,ios
    character(len=80) :: id, grpinfo
    
    open(42,file=grpfil,status='old',form='formatted')
    do i=1,nloci
        read(42,*,iostat=ios) id, grpinfo
        if(ios/=0) then
            print *, 'error: reading grpfile:', grpfil 
            stop
        end if
        if(str_match(trim(id),snpid(i)) .neqv. .true.) then
            print *, 'error: SNP IDs in grpfile and bimfile not matching'
            stop
        end if
        do j=1,ngroup
            if(str_match(trim(grpinfo),grp(j))) then
                snptracker(i,1)=j
                exit
            end if
        end do
    end do
    close(42,status='keep')
end subroutine load_snp_grp

subroutine load_pheno
  character(len=1024) :: str
  character(len=80) :: tmp, tID
  integer i, j, ncolumns, nfields,ios
  integer, parameter :: maxfields=120, maxfieldlength=80
  character(len=1), parameter :: headersymbol='$'
  character(len=3) :: sep=' , '
  character(len=maxfieldlength),dimension(maxfields) :: fields

  sep(3:3) = achar(9) 
  open(31,file=trim(phenfil),status='old',form='formatted')
  do i=1,nind
      read(31,'(a)') str
      read(str,*) tmp,animalID(i),tmp,tmp,tmp, why(i)
  enddo    
  close(unit=31,status='keep')
  if(alterpheno .eqv. .true.)  then
      open(unit=30,file=trim(alterphenfil),status='old',form='formatted')
      ncolumns = 0
      if(pheno_name /= "") then
          read(30,'(a)') str
          call left_trim(str)
          if(index(str,headersymbol) /=1) then
              write(*,*) 'error: header line must start with',headersymbol
              stop
          end if
          call get_fields(str(2:),sep,nfields,fields)
          ncolumns = nfields
          do j=1,nfields
              if(str_match(fields(j),pheno_name)) then
                  pheno_ncol = j
                  exit
              end if
          end do
          if(pheno_ncol <= 0) then
              write(*,*) 'error: phenotype column name cannot be found.'
              stop
          end if
      end if
      if(pheno_ncol <=0) then
          write(*,'(a)') 'error: phenotype location is not specified.'
          stop
      end if 
      i=0
      do while(i<nind) 
          read(30,'(a)',iostat=ios) str
          if(ios /=0) return
          call left_trim(str)
          if(index(str,headersymbol) ==1) cycle
          i=i+1
          read(str,*) tmp,tID
          if(str_match(trim(tID), trim(animalID(i))) .neqv. .true.) then
              write(*,'(a)') 'error: animal ID '//trim(tID)//' in phenotype file does not match '&
                             //trim(animalID(i))//' in the plink file.'
          end if
          call get_fields(str,sep,nfields,fields)
          if(ncolumns == 0) then
              ncolumns = nfields
          else
              if(nfields /= ncolumns) then
                  write(*,*) 'line:',i,'no. of fields',nfields,'does not match',ncolumns
                  stop
              end if
          end if
          read(fields(pheno_ncol),*) why(i)
      end do
      if(i < nind) then
          write(*,*) 'error: number of records less than',nind
      end if     
      close(unit=30,status='keep')
  end if     
end subroutine load_pheno

subroutine load_geno
integer*1, parameter :: m1=b'01101100'
integer*1, parameter :: m2=b'00011011'
integer*1, parameter :: m3=b'00000001'
integer*1 :: magic1, magic2, mode
character(len=:), allocatable :: line
character(len=80) :: dum, id
integer*1 :: b1
integer :: i, j, k, pos, ios, nbytes

nbytes = ceiling(dble(nind)/4)
allocate(character(len=nbytes) :: line)
open(unit=40, file=trim(bimfil),status='old',form='formatted')
do j=1,nloci
    read(40,*,iostat=ios) dum, id
    if(ios/=0) then
        write(*,*) 'error: reading snpID from: '//trim(bimfil)
    end if
    snpID(j) = trim(id)
end do

open (unit=41,file=trim(genfil),status='old',access='stream',form='unformatted')
read(41) magic1, magic2, mode
if (magic1.ne.m1.or.magic2.ne.m2) then
   write(*,*)  'Binary genotype file may not be a plink file'
   write(21,'(a)') 'Binary genotype file may not be a plink file'
   call flush(21)
   stop
end if
if (mode /= m3) then
  write(*,'(a)') 'should not be individual mode - check >>>'
  write(21,'(a)') 'SNP file should not be in individual mode'
  call flush(21)
  stop
end if

do j=1,nloci
   read(41) line
   k = 0
   pos = 1
   do i=1,nind
      if(k==0) then 
          b1 = ichar(line(pos:pos))
          pos = pos + 1
          !The following shiftr+and equation gives same results 
          !from left or right of sequences coded 00, 01/10, 11
          !Missing values were treated the same as heterozygoes
          b1 = b1 - iand(shiftr(b1,1), FIVEMASK8)
      end if
      X(i,j) = ibits(b1,k,2)
      k=mod(k+2,8)
   enddo
enddo
deallocate(line)
close(41,status='keep')
if(BayesMethod==BayesRc .and. mcmc) then
        call load_snp_grp
end if
end subroutine load_geno

subroutine load_geno_bits
    implicit none
    integer*1, parameter :: m1=b'01101100'
    integer*1, parameter :: m2=b'00011011'
    integer*1, parameter :: m3=b'00000001'    
    integer*1 :: magic1, magic2, mode
    integer :: i,j,k
    integer :: ibatch, nbytes
    integer*1, allocatable :: lastbatch(:)

    ibatch = nind/(int64*4)
    nbytes = ceiling(dble(mod(nind,int64*4))/4)
    allocate(lastbatch(nbytes))
    genbits = 0
    lastbatch = 0
    
    open (unit=41,file=trim(genfil),status='old',access='stream',form='unformatted')
    read(41) magic1, magic2, mode
    if (magic1.ne.m1.or.magic2.ne.m2.or.mode.ne.m3) then
        write(*,*)  'Binary genotype file is broken'
        write(21,'(a)') 'Binary genotype file is broken'
        call flush(21)
        stop
    end if
    do i=1,nloci
        do j=1,ibatch
            read(41) genbits(j,i)
            genbits(j,i) = genbits(j,i) - iand(shiftr(genbits(j,i),1), FIVEMASK64)
        end do
        if(nbytes>0) then
            read(41) lastbatch
            do k=1,nbytes
              call mvbits(int(lastbatch(k),8),0,8,genbits(nbatch,i),(k-1)*8)              
            end do
            lastbatch = 0 
            genbits(nbatch,i) = genbits(nbatch,i) - iand(shiftr(genbits(nbatch,i),1), FIVEMASK64)
        endif
    end do
    deallocate(lastbatch)
    close(41,status='keep')
end subroutine load_geno_bits

subroutine set_random_seed
  integer :: i, n, clock
  integer(int64) :: seed(2)

  if (seed1 /= 0) then
     seed(1)=abs(seed1)
     seed(2)=seed(1) + 1
  else
     call system_clock(count=clock)
     seed = clock + 37 * (/ (i - 1, i = 1, 2) /)
     seed1 = seed(1)
  endif
  call rand_set_seed(maxthreads,seed)  
end subroutine set_random_seed

subroutine xscale
  !RESCALE GENOTYPE COVARIATES TO HAVE ZERO MEANS AND UNIT VARIANCES
  !NO MISSING GENOTYPES WERE ASSUMED OR MISSING VALUES TREATED AS HETEROZYGOTES
  integer :: j
  real(real64) :: p
  character(len=20) :: dum
  if(mcmc) then
     do j=1,nloci
        p=sum(X(:,j)) / (2.0d0*nind)
        if( p < epsilon(1.0d0) .or. (1.0d0-p) < epsilon(1.0d0) ) then
           X(:,j)=0.0d0
        else
           X(:,j)=(X(:,j)-2.0d0*p)/sqrt(2.0d0*p*(1.0d0-p))
        endif
        freqstore(j)=p
     enddo
     open(45,file=trim(freqfil),status='unknown')
     do j=1,nloci
        write(45,'(A10,1X,F10.6)') snpID(j), freqstore(j)
     enddo
     close(45,status='keep')
  else
     open(45,file=trim(freqfil),status='unknown')
      do j=1,nloci
        read(45,*) dum,freqstore(j)
     enddo
     do j=1,nloci
        p=freqstore(j)
        if( p < epsilon(1.0d0) .or. (1.0d0-p) < epsilon(1.0d0) ) then
           X(:,j)=0.0d0
        else
           X(:,j)=(X(:,j)-2.0d0*p)/sqrt(2.0d0*p*(1.0d0-p))
        endif
     enddo
     close(45,status='keep')
  endif

end subroutine xscale
    
subroutine allocate_data
  integer :: ios, ios2
  
  ios=0
  ios2=0
  if(mcmc) then
      allocate(animalID(nind), snpID(nloci), why(nind), pred(nind), X(nind,nloci), g(nloci), yadj(nind), xpx(nloci), &
          permvec(nloci), varindist(ngroup, nmix), varstore(ngroup, nmix), gstore(nloci), freqstore(nloci), stat=ios)
      if(BayesMethod .eq.BayesA) then
          allocate(vara_s(nloci), vara_sstore(nloci), stat=ios2) 
      else if(BayesMethod .eq. BayesB) then
          allocate(vara_s(nloci), vara_sstore(nloci), p(ngroup, nmix), log_p(ngroup,nmix), snpindist(ngroup, nmix), &
              snptracker(nloci,2), snpstore(ngroup,nmix), indiststore(nloci,nmix), stat=ios2)
      else if(BayesMethod .eq. BayesC) then
          allocate(p(ngroup, nmix), log_p(ngroup,nmix), snpindist(ngroup, nmix), snptracker(nloci,2), &
              snpstore(ngroup,nmix),indiststore(nloci,nmix), vare_gp(nmix), gp(nmix), stat=ios2)
      else
          allocate(p(ngroup, nmix), log_p(ngroup,nmix), snpindist(ngroup, nmix), snptracker(nloci,2), snpstore(ngroup,nmix), &
              indiststore(nloci,nmix), vare_gp(nmix), gp(nmix), pstore(ngroup, nmix), dirx(ngroup,nmix), stat=ios2)
      end if
  else
      allocate(animalID(nind), snpID(nloci), X(nind,nloci), freqstore(nloci), gstore(nloci), pred(nind), stat=ios)
  end if
  
  if( ios /= 0 .or. ios2 /= 0) then
      print *, 'ERROR :: Unable to allocate required storage'
      stop 'Unable to allocate required storage for data'
  endif
end subroutine allocate_data

subroutine permutate(p,n,s)
  integer s,n,p(n)
  integer i,j,k,ipj,itemp,m
  real(real64),allocatable :: u(:)
  allocate(u(s))
  do i=1,n
     p(i)=i
  enddo
  do i=1,n,s
     m=min(n-i+1,s)
     call random(s,u,0)
     do j=1,m
        ipj=i+j-1
        k=int(u(j)*(n-ipj+1))+ipj
        itemp=p(ipj)
        p(ipj)=p(k)
        p(k)=itemp
     enddo
  enddo
end subroutine permutate

subroutine compute_residuals
  integer :: i
  do i=1,nind
     yadj(i)=why(i)-dot_product(X(i,1:nloci),g)-mu
  enddo
end subroutine compute_residuals

subroutine compute_dgv
  integer :: i
  pred=-9999.0d0
  do i=1,nind
    pred(i)=mu+dot_product(X(i,1:nloci),gstore(1:nloci))
  end do
end subroutine compute_dgv

subroutine write_dgv()
 integer :: i

 open(unit=61,file=mbvfil,status='unknown',form='formatted')
 do i=1,nind
    write(61,'(A20,1x,E15.7)') animalID(i), pred(i)      
 enddo
 close(unit=61,status='keep')
end subroutine write_dgv

subroutine setup_prediction
  character(len=100):: str
  character(len=100):: dum
  integer :: i, ios
  real(real64), dimension(20):: gtemp

  open(31,file=trim(phenfil),status='old',form='formatted')
  do i=1,nind
      read(31,'(a)') str
      read(str,*) dum,animalID(i)
  enddo    
  close(31)
  
  open(51,file=efffil,status='old',form='formatted')
  read(51,'(a)') str
  do i=1,nloci
     read(51,*,iostat=ios) dum, gstore(i)
     if(ios/=0) then
        write(21,'(a)') 'Error opening file ',adjustl(efffil)
     call flush(21)
     end if
  enddo
  close(51)
  open(12,file=modfil,status='old',form='formatted')
  read(12,*) dum, mu
  close(12)
end subroutine setup_prediction

subroutine output_model
  integer :: i, j, ios
  character(len=10) :: ci
  character(len=15) :: ca

  open(unit=51,file=efffil,status='unknown',iostat=ios)
  if(ios/=0) then
     write(21,*) 'Error opening ',trim(efffil)
     stop
  end if
  ci='snpID'
  if(BayesMethod .eq. BayesA) then
      write(51,'(A10,1x,A15,1x,A15)') adjustl(ci), 'snpEffect', 'snpVariance'    
      do i=1,nloci
          write(51,'(A10,2(1X,E15.6))') snpID(i),gstore(i),vara_sstore(i)
      enddo
  else if(BayesMethod .eq. BayesB) then
      write(51,'(A10,1x,A15,1x,A15)',advance='no') adjustl(ci), 'snpEffect', 'snpVariance'          
      do i=1,nmix
          write(ca,'(I15)') i
          ca=adjustl(ca)
          ca="PIP"//trim(ca)
          ca=adjustr(ca)
          write(51,'(1x,A15)',advance="no") ca
      end do
      write(51,*)
      do i=1,nloci
          write(51,'(A10,2(1X,E15.6),50(1X,E15.6))') snpID(i),gstore(i),vara_sstore(i),indiststore(i,:) 
      enddo
  else
      write(51,'(A10,1x,A15)',advance='no') adjustl(ci), 'snpEffect'         
      do i=1,nmix
          write(ca,'(I15)') i
          ca=adjustl(ca)
          ca="PIP"//trim(ca)
          ca=adjustr(ca)
          write(51,'(1x,A15)',advance="no") ca
      end do
      write(51,*)
      do i=1,nloci
          write(51,'(A10,1X,E15.6,50(1X,E15.6))') snpID(i),gstore(i),indiststore(i,:) 
      enddo      
  end if
  close(51,status='keep')

      
  open(unit=12,file=modfil,status='unknown',iostat=ios)
  if(ios/=0) then
     write(21,*) 'Error opening ',trim(modfil)
     stop
  end if
  write(12,800) 'Mean',mu_store
  write(12,800) 'Ve', vare_store
  if(BayesMethod .ne. BayesA .and. BayesMethod .ne. BayesB) then  
      write(12,800) 'Vara',vara_store
  end if
  if(BayesMethod .ne. BayesA) then
      write(12,800) 'NSNP',included_store
      do i=1,ngroup
          write(ca,'(I8)') i
          ca=adjustl(ca)
          ca="Group "//trim(adjustl(ca))
          do j=1,nmix
              write(ci,'(I8)') j
              ci=adjustl(ci)
              ci="N"//trim(adjustl(ci))
              if(ngroup==1) then
                  write(12,800) ci,snpstore(i,j)
              else          
                  write(12,801) ca,ci,snpstore(i,j)
              end if
          end do
      end do
      if(BayesMethod .ne. BayesB .and. BayesMethod .ne. BayesC) then
          do i=1,ngroup
              write(ca,'(I8)') i
              ca=adjustl(ca)
              ca="Group "//trim(adjustl(ca))       
              do j=1,nmix   
                  write(ci,'(I8)') j
                  ci=adjustl(ci)
                  ci="P"//trim(adjustl(ci))
                  if(ngroup==1) then
                      write(12,800) ci,pstore(i,j)
                  else          
                      write(12,801) ca,ci,pstore(i,j)
                  end if
              end do
          end do
      end if
  end if
  do i=1,ngroup
      write(ca,'(I8)') i
      ca=adjustl(ca)
      ca="Group "//trim(adjustl(ca))
      do j=1,nmix   
          write(ci,'(I8)') j
          ci=adjustl(ci)
          ci="Va"//trim(adjustl(ci))
          if(BayesMethod .eq. BayesA) ci="Va"
          if(ngroup==1) then
              write(12,800) ci,varstore(i,j)
          else
              write(12,801) ca,ci,varstore(i,j)
          end if
      end do
  end do
  close(12,status='keep')

800 format(a,t10,E15.6)
801 format(a,t10,a,t15,E15.6)    
end subroutine output_model


subroutine clean_up
    if(mcmc) then
        deallocate(ss, animalID, snpID, why, pred, X, g, yadj, xpx, permvec, varindist, varstore, gstore, freqstore)
        if(BayesMethod .eq. BayesA) then
            deallocate(vara_s, vara_sstore)
        else if(BayesMethod .eq. BayesB) then
            deallocate(vara_s, vara_sstore, p, log_p, snpindist, snptracker, snpstore, indiststore) 
        else if(BayesMethod .eq. BayesC) then
            deallocate(p, log_p, snpindist, snptracker, snpstore, indiststore, gp, mix, vare_gp)
        else
            deallocate(p, log_p, snpindist, snptracker, snpstore, indiststore, gp, mix, vare_gp, dirx, alpha, pstore)
        end if           
	    if(.not.SODAOFF) then
            deallocate(xxd, XA, yA, yAadj, skip_j, i_skipped)
        end if
    else
        deallocate(pred, animalID, snpID, gstore, X, freqstore)
	end if
end subroutine clean_up

end module utilities
 