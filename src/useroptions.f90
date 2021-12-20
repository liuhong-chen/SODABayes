    
module useroptions
    use globals
    use utilities
    implicit none

    integer :: narg
    integer, parameter:: nopt=31
    character(len=1024) :: arg
    character(len=1024), allocatable, dimension(:) :: cmd_line, tmp_array
    integer, parameter :: a_int = 1, a_float = 2, a_char = 3, a_flag = 4
    type option
        integer:: pos
        character(len=20):: key
        character(len=20):: argtype
        character(len=100):: desc
        integer:: kind
        character(len=200):: default
    end type option
    type(option) :: opt(nopt)

    contains 
  subroutine get_default_options
    opt = (/&
    option(1, '-parfile'   ,'[filename]','parameter file'                          ,a_char,         ''),&
    option(2, '-job_name'  ,'[char]',    'job name'                                ,a_char,         ''),&        
    option(3, '-method'    ,'[char]',    'Bayesian method'                         ,a_char, 'BayesCpi'),&     
    option(4, '-bfile'     ,'[prefix]',  'prefix PLINK binary files'               ,a_char,         ''),&
    option(5, '-pheno'     ,'[filename]','alternative phenotype file'              ,a_char,         ''),&
    option(6, '-pheno_name','[char]',    'phenotype column name'                   ,a_char,         ''),&
    option(7, '-pheno_ncol','[num]',     'phenotype column no.'                    ,a_int,         '3'),&
    option(8, '-grpfile'   ,'[filename]','SNP group info file for BayesRc'         ,a_char,         ''),&
    option(9, '-pi'        ,'[num]',     '% SNP with zero effects (BayesB/BayesC)' ,a_float,    '0.99'),&
    option(10,'-scale_va'  ,'[num]',     'scale for SNP variance prior'            ,a_float,    '0.01'),&
    option(11,'-scale_ve'  ,'[num]',     'scale for residual variance prior'       ,a_float,    '0.01'),&
    option(12,'-df_va'     ,'[num]',     'degrees of freedom Va'                   ,a_float,    '0.01'),&
    option(13,'-df_ve'     ,'[num]',     'degrees of freedom Ve'                   ,a_float,    '0.01'),&
    option(14,'-alpha'     ,'[num]',     'prior for Dirichlet/Beta distribution'   ,a_float,     '1.0'),&
    option(15,'-numit'     ,'[num]',     'length of MCMC chain'                    ,a_int,     '50000'),&
    option(16,'-burnin'    ,'[num]',     'burnin steps'                            ,a_int,      '5000'),&
    option(17,'-thin'      ,'[num]',     'thinning rate'                           ,a_int,         '1'),&
    option(18,'-nmix'      ,'[num]',     'number of mixture distributions'         ,a_int,         '2'),&
    option(19,'-mix'       ,'[num]',     'effect sizes of mixtures (% of Va)'      ,a_float, '0.0,1.0'),&
    option(20,'-seed'      ,'[num]',     'seed value for random number generator'  ,a_int,         '0'),&
    option(21,'-SODAOn'    ,'[flag]',    'turn soda on'                            ,a_flag,        't'),&
    option(22,'-SODAOff'   ,'[flag]',    'turn soda off'                           ,a_flag,        'f'),&
    option(23,'-nthreads'  ,'[num]',     'number of threads'                       ,a_int,         '1'),&
    option(24,'-blocksize' ,'[num]',     'number of SNPs in a block'               ,a_int,         '0'),&
    option(25,'-nblocks'   ,'[num]',     'number of blocks'                        ,a_int,         '0'),&
    option(26,'-rblocks'   ,'[flag]',    'randomly assign SNPs into blocks'        ,a_flag,        'f'),&
    option(27,'-predict'   ,'[flag]',    'perform prediction'                      ,a_flag,        'f'),&
    option(28,'-model'     ,'[filename]','model summary file (for prediction) '    ,a_char,         ''),&
    option(29,'-freq'      ,'[filename]','SNP frequency file (for prediction)'     ,a_char,         ''),&
    option(30,'-effects'   ,'[filename]','SNP effect file (for prediction)'        ,a_char,         ''),&
    option(31,'-out'       ,'[prefix]',  'prefix for output'                       ,a_char,         '')&    
    /)    
  end subroutine get_default_options
  
  subroutine show_options
   integer :: k
   do k=1,nopt
      write(*,111) adjustl(opt(k)%key),adjustl(opt(k)%argtype), &
           adjustl(opt(k)%desc), adjustl(opt(k)%default)
   end do
111 format(a12,a12,a45,a25)
  end subroutine show_options    
    
  subroutine parse_parfile
    integer ::  i, j
    character(len=1024) :: key
    
    useparfile = .false.  
    i=1
    do while (i.le.narg)
       key=trim(cmd_line(i))
       if(str_match(key,'-parfile')) then
           parfil=trim(cmd_line(i+1))
           useparfile=.true.
           exit
       end if
       j=index(key,'.par')
       if(j /=0 .and. j==len_trim(key)-3) then     
           parfil = trim(cmd_line(i))
           useparfile = .true.
           exit
       end if
       i=i+1
    end do
  end subroutine parse_parfile
  
  subroutine parse_parfile_options(narg, parameters)
    integer, intent(in) :: narg
    character(len=1024), intent(in) :: parameters(:)
    integer ::  i, k, kind
    character(len=1024) :: key
    character(len=100) :: err1, err2
    logical :: valid

    err1='Unknown command line option : '
    err2='Missing argument for :'
    i=1
    do while(i.le.nopt)
        if(str_match(opt(i)%key,'-parfile')) then
            opt(i)%default=parfil
            exit         
        end if 
        i=i+1
    end do    
    i=1
    do while (i.le.narg)
       key=trim(parameters(i))
       k=1
       valid=.false.        
       do while(k.le.nopt)
          if(str_match(key,trim(opt(k)%key))) then
             kind=opt(k)%kind
             valid=.true.
             if(kind==4) then
                opt(k)%default='t'
                i=i+1
             else if(i==narg) then
                print *, trim(err2),trim(key),i
                stop 'ERROR: Problem parsing the option in the .par file'
                valid=.false.
             else if(is_key(trim(parameters(i+1)))) then
                print *, trim(err2),trim(key),i
                stop 'ERROR: Problem parsing the option in the .par file.'
             else
                opt(k)%default=trim(parameters(i+1))
                valid=.true.
                i=i+2
             endif
             exit           
          end if
          k=k+1
          valid=.false.
       enddo
       if(.not.valid) then
          print *, trim(err1),trim(key),i
          stop 'ERROR: Problem parsing the option in the .par file.'
       end if
    enddo
    
  end subroutine parse_parfile_options
  
  subroutine get_cmdline
  
    integer i
    
    narg = command_argument_count()
    allocate(cmd_line(narg))
    do i=1,narg
       call get_command_argument(i,arg)
       cmd_line(i)=arg
    end do
  end subroutine get_cmdline
  
  subroutine parse_cmdline
    integer ::  i, k, kind
    character(len=1024) :: key
    character(len=100) :: err1, err2
    logical :: valid

    err1='Unknown command line option : '
    err2='Missing argument for :'
    
    i=1
    do while (i.le.narg)
       key=trim(cmd_line(i))
       k=1
       valid=.false.
       if(index(key,'.par')==(len(trim(key))-3)) then                 
           valid = .true.
           i=i+1   
           cycle
       end if        
       do while(k.le.nopt)
          if(str_match(key,trim(opt(k)%key))) then
             kind=opt(k)%kind
             valid=.true.
             if(kind==4) then
                opt(k)%default='t'
                i=i+1
             else if(i==narg) then
                print *, trim(err2),trim(key),i
                stop 'ERROR: Problem parsing the command line arguments'
                valid=.false.
             else if(is_key(trim(cmd_line(i+1)))) then
                print *, trim(err2),trim(key),i
                stop 'ERROR: Problem parsing the command line arguments'
             else
                opt(k)%default=trim(cmd_line(i+1))
                valid=.true.
                i=i+2
             endif
             exit           
          end if
          k=k+1
          valid=.false.
       enddo
       if(.not.valid) then
          print *, trim(err1),trim(key),i
          stop 'ERROR: Problem parsing the command line arguments'
       end if
    enddo
  end subroutine parse_cmdline

  subroutine parse_help
   integer :: i,k
   character(len=1024) :: key
   if(narg==0) then
      write(*,'(a)')adjustl('argument    type        description &
           &                                default')
      do k=1,nopt
         write(*,111) adjustl(opt(k)%key),adjustl(opt(k)%argtype), &
              adjustl(opt(k)%desc), adjustl(opt(k)%default)
      end do
      stop
   else
      do i=1,narg
         key=trim(cmd_line(i))
         if(str_match(trim(key(1:2)),'-h')) then
            write(*,'(a)')adjustl('argument    type        description &
                 &                                default')
            do k=1,nopt
               write(*,111) adjustl(opt(k)%key),adjustl(opt(k)%argtype), &
                    adjustl(opt(k)%desc), adjustl(opt(k)%default)
            end do
            stop
         end if
111      format(a12,a12,a45,a25)
      end do
   end if
 end subroutine parse_help

subroutine parse_out
  integer::i
  do i=1,nopt
    if(str_match(trim(opt(i)%key),'-out')) then
        outprefix = trim(opt(i)%default)
        logfil = trim(outprefix)//'.log'
        freqfil= trim(outprefix)//'.frq'
        mbvfil= trim(outprefix)//'.gv'
        hypfil= trim(outprefix)//'.hyp'
        locfil= trim(outprefix)//'.snp'
        modfil=trim(outprefix)//'.mod'
        efffil=trim(outprefix)//'.eff'
        exit
     end if
  end do
end subroutine parse_out

subroutine parse_plink
  integer::i
  logical:: fileExist

  do i=1,nopt
     if(str_match(trim(opt(i)%key),'-bfile')) then
        inprefix=trim(opt(i)%default)
        genfil=trim(inprefix)//'.bed'
        inquire(file=genfil,exist=fileExist)
        if(.not.fileExist)then
           print *,'plink file ',trim(genfil),' not found'
           stop
        end if

        bimfil=trim(inprefix)//'.bim'
        inquire(file=bimfil,exist=fileExist)
        if(.not.fileExist)then
           print *,'plink file ',trim(bimfil),' not found'
           stop
        end if
        phenfil=trim(inprefix)//'.fam'
        inquire(file=phenfil,exist=fileExist)
        if(.not.fileExist)then
           print *,'plink file ',trim(phenfil),' not found'
           stop
        end if        
        exit
     end if
  enddo
end subroutine parse_plink

subroutine parse_pheno
  integer:: i, j
  logical:: fileExist
  character(len=1024) :: key

  alterpheno=.false.
  pheno_name=""
  pheno_ncol=0

  do i=1,nopt
     if(str_match(trim(opt(i)%key),'-pheno')) then          
        alterphenfil=trim(opt(i)%default)
        if(len_trim(alterphenfil)==0) return
        alterpheno = .true.        
        inquire(file=alterphenfil,exist=fileExist)
        if(.not.fileExist)then
           print *,'phenotype file ',trim(alterphenfil),' not found'
           stop
        end if  
        do j=1,nopt
            if(str_match(trim(opt(j)%key), '-pheno_name')) then
                pheno_name=trim(opt(j)%default)
                exit
            end if
        end do
        do j=1,nopt
            if(str_match(trim(opt(j)%key), '-pheno_ncol')) then
                pheno_ncol = cast_int(opt(j)%default)
                exit
            end if
        end do
        alterpheno=.true.
        exit
     end if
  end do 
end subroutine parse_pheno

     
subroutine parse_predict
  integer::i
  logical :: flag, fileExist
  character(len=200) :: filein
  mcmc=.true.
  do i=1,nopt
     if(str_match(trim(opt(i)%key),'-predict')) then
        flag = cast_logical(opt(i)%default)
        if(flag) mcmc=.false.
     end if
  end do
  if(.not.mcmc) then
     do i=1,nopt
        if(str_match(trim(opt(i)%key),'-model')) then
           filein=trim(opt(i)%default)
           if(str_match(trim(filein),' ')) then
               print *,'No model file specified'
               stop 'Error: No model file specified'
           end if
           modfil=trim(filein)
           inquire(file=modfil,exist=fileExist)
           if(.not.fileExist)then
              print *,'model file ',trim(modfil),' not found'
              stop
           end if
        else if(str_match(trim(opt(i)%key),'-freq')) then
           filein=trim(opt(i)%default)
           if(str_match(trim(filein),' ')) then 
               print *, 'No SNP frequency file specified'
               stop 'Error: No SNP frequency file specified'
           end if    
           freqfil=trim(filein)
           inquire(file=freqfil,exist=fileExist)
           if(.not.fileExist)then
              write(*,*)'freq file ',trim(freqfil),' not found'
              stop
           end if
        else if(str_match(trim(opt(i)%key),'-effects')) then
           filein=trim(opt(i)%default)
           if(str_match(trim(filein),' ')) then
               print *, 'No SNP effect file specified'
               stop 'Error: No SNP effect file specified'
           end if    
           efffil=trim(filein)
           inquire(file=freqfil,exist=fileExist)
           if(.not.fileExist)then
              print *,'effect file ',trim(efffil),' not found'
              stop
           end if
        end if
     end do
  end if
end subroutine parse_predict

subroutine parse_nmix
  integer::i
  do i=1,nopt
     if(str_match(trim(opt(i)%key),'-nmix')) then
        nmix = cast_int(opt(i)%default)
        if(nmix>maxdist) then
            print *, 'Error: nmix is too big. maximum allowed: ', maxdist
            stop
        end if        
        exit
     end if
  end do
  if(BayesMethod .eq. BayesA) then
      nmix=1
      return
  else if(BayesMethod .eq. BayesB .or. BayesMethod .eq. BayesC .or. BayesMethod .eq. BayesCpi) then
      nmix=2
      allocate(mix(nmix))
      if(BayesMethod .eq. BayesCpi) then
          allocate(alpha(nmix))
      end if
  else
      allocate(mix(nmix),alpha(nmix))
  end if
  
end subroutine parse_nmix

subroutine parse_method
  integer i
  character(len=20) :: method
  do i=1,nopt
      if(str_match(trim(opt(i)%key),'-method')) then
          method = opt(i)%default
          if(ssmatch(method,'BayesA')) then
              BayesMethod=BayesA
          else if(ssmatch(method,'BayesB')) then
              BayesMethod=BayesB
          else if(ssmatch(method,'BayesC')) then
              BayesMethod=BayesC
          else if(ssmatch(method,'BayesCpi')) then
              BayesMethod=BayesCpi
          else if(ssmatch(method,'BayesR')) then
              BayesMethod=BayesR
          else if(ssmatch(method,'BayesRc')) then
              BayesMethod=BayesRc
          else
              write(*,*) 'error: selecting Bayesian method'
              write(*,*) 'supported methods: BayesA, BayesB, BayesC, BayesCpi, BayesR, and BayesRc'
              stop
          end if
      end if
  end do
end subroutine parse_method

subroutine parse_pi
    integer i    
    do i=1,nopt
        if(str_match(trim(opt(i)%key),'-pi')) then
            pizero = cast_float(opt(i)%default)
            if(pizero<=0.0d0 .or. pizero>=1.0d0) then
                print *, 'Error: pi must be a value between 0 and 1'
                stop
            end if
        end if    
    end do
end subroutine parse_pi
    
subroutine parse_snpgrp_file
    integer i   
    logical fileExist
    do i=1,nopt
        if(str_match(trim(opt(i)%key),'-grpfile')) then
            grpfil = trim(opt(i)%default)
            inquire(file=grpfil,exist=fileExist)
            if(.not.fileExist)then
                print *,'snp group file ',trim(grpfil),' not found'
                stop
            end if   
        end if
    end do
end subroutine parse_snpgrp_file
                
subroutine parse_soda
  integer :: i
  logical :: flag
  SODAOFF = .false.
  rblocks = .false.
  setblocksize = .false.
  setblocks = .false.
  do i=1,nopt
     if(ssmatch(trim(opt(i)%key),'-nthreads')) then
        nthreads = cast_int(opt(i)%default)
        nthreads = min(nthreads,maxthreads)
     else if(ssmatch(trim(opt(i)%key),'-SODAOFF')) then
        flag = cast_logical(opt(i)%default)
        if(flag) SODAOFF = .true. 
     else if(ssmatch(trim(opt(i)%key),'-SODAON')) then
        flag = cast_logical(opt(i)%default)
        if(flag) SODAOFF = .false.        
     else if(ssmatch(trim(opt(i)%key),'-rblocks')) then
        flag = cast_logical(opt(i)%default)
        if(flag) rblocks = .true.
     else if(ssmatch(trim(opt(i)%key),'-blocksize')) then
        blocksize = cast_int(opt(i)%default)
        if(blocksize > 0) setblocksize = .true.
     else if(ssmatch(trim(opt(i)%key),'-nblocks')) then
        nblocks = cast_int(opt(i)%default)
        if(nblocks > 0) setblocks = .true.
     end if     
  end do
end subroutine parse_soda

subroutine parse_priors
  integer::i,k,n,nitem
  character(len=128), dimension(nmix) :: c_string

  do i=1,nopt
     if(str_match(trim(opt(i)%key),'-scale_va')) then
        scale_va = cast_float(opt(i)%default)
     else if(str_match(trim(opt(i)%key),'-scale_ve')) then
        scale_ve = cast_float(opt(i)%default)
     else if(str_match(trim(opt(i)%key),'-df_va')) then
        df_va = cast_float(opt(i)%default)
     else if(str_match(trim(opt(i)%key),'-df_ve')) then
        df_ve = cast_float(opt(i)%default)
     else if(str_match(trim(opt(i)%key),'-mix')) then
         if(BayesMethod .eq. BayesA .or. BayesMethod .eq. BayesB) cycle
         if(BayesMethod .eq. BayesC .or. BayesMethod .eq. BayesCpi) then
             mix(1)=0.0d0
             mix(2)=1.0d0
             cycle
         end if
        nitem=nfields(opt(i)%default, ',')
        if(nitem /= nmix) then
           print *, 'Error: Number of mixtue classes ',nmix
           print *,'        but ',nitem,'effect sizes specified (mix)'
           stop
        else
           call get_fields(trim(opt(i)%default),',',n,c_string)
           do k=1,nmix
              mix(k)=cast_float(trim(c_string(k)))
           enddo
        endif
     else if(str_match(trim(opt(i)%key),'-alpha')) then
        if(BayesMethod .eq. BayesA .or. BayesMethod .eq. BayesB .or. BayesMethod .eq. BayesC) cycle       
        nitem=nfields(opt(i)%default, ',')
        if(nitem==1) then
           call get_fields(trim(opt(i)%default),',',n,c_string)
           alpha=cast_float(trim(c_string(1)))
        else if(nitem /= nmix) then
           print *, 'Error: Number of mixtue classes',nmix
           print *,'        but',nitem,'prior values for Dirichlet pecified (alpha)'
           stop
        else
           call get_fields(trim(opt(i)%default),',',n,c_string)
           do k=1,nmix
              alpha(k)=cast_float(trim(c_string(k)))
           enddo
        endif        
     end if
  end do
end subroutine parse_priors

subroutine parse_initialise
  integer::i
  do i=1,nopt
     if(str_match(trim(opt(i)%key),'-numit')) then
        numit = cast_int(opt(i)%default)
     else if(str_match(trim(opt(i)%key),'-burnin')) then
        burnin= cast_int(opt(i)%default)
     else if(str_match(trim(opt(i)%key),'-thin')) then
        thin= cast_int(opt(i)%default)
     else if(str_match(trim(opt(i)%key),'-seed')) then
        seed1= cast_int(opt(i)%default)
     endif
  end do
end subroutine parse_initialise

subroutine parse
  call parse_help
  call parse_out
  call parse_plink
  call parse_predict
  if(mcmc) then
      call parse_method 
      if(BayesMethod==BayesB .or. BayesMethod==BayesC) then
          call parse_pi
      else if(BayesMethod==BayesRc) then
          call parse_snpgrp_file
      end if
      call parse_nmix      
      call parse_pheno  
      call parse_soda
      call parse_initialise
      call parse_priors
  end if
end subroutine parse

end module useroptions
