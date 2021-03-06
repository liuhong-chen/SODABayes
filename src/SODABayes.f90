program SODABayes
    use useroptions
    use utilities
    use models
    use omp_lib
    use soda
    
    implicit none
    logical :: fileExist
    integer i, ios, npar, status, pos1, pos2
    integer nitems, start_job, end_job
    integer, parameter :: maxpar=100
    character(len=2), parameter :: commentchar='#!'
    character(len=1024) :: line, item, tmp
    character(len=1024), allocatable :: parameters(:)
    character(len=200), dimension(100) :: values
    
    call get_cmdline
    call parse_parfile
    if(useparfile) then
        inquire(file=parfil,exist=fileExist)
        if(.not.fileExist)then
           print *,'par file ',trim(parfil),' not found'
           stop 
        end if 
        open(unit=10, file=trim(parfil), status='old', form='formatted', iostat=ios)
        if(ios/=0) then
            print *, 'cannot open parfile:', parfil
            stop
        end if
        allocate(parameters(maxpar))
        start_job = 0
        end_job = 1
        job_name=""
        do while(.true.)
            read(10,'(a)',iostat=ios) line
            if(ios /=0) then
                exit
            end if
            call left_trim(line)
            if(len_trim(line)==0) cycle            
            if(scan(line,commentchar)==1) cycle
            pos2 = scan(line,commentchar)-1
            if(pos2==-1) pos2=len_trim(line) 
            line=trim(line(1:pos2))
            read(line,*) item
            if(ssmatch(item,'begin')) then
                if(end_job /=1) then
                    write(*,*) 'error: Need an END statement for job: ', job_name
                    stop
                end if
                end_job = 0
                tmp = line
                call to_upper(tmp)
                pos1 = verify(tmp,'BEGIN')
                job_name = line(pos1:pos2)
                call left_trim(job_name)
                npar = 0
                do i=1,maxpar
                   parameters(i)=""
                end do
                start_job = 1
            else if(ssmatch(item,'end')) then
                if(start_job /= 1) then
                    write(*,*) 'error: need a BEGIN statement to start a job.'
                    stop
                end if
                call get_default_options               
                call parse_parfile_options(npar, parameters)
                call parse_cmdLine     
                write(*,'(a)') 'job '//trim(job_name) //' started...'              
                call run_program_sodaBayes(status)
                if(status /= 0) then
                    write(*,'(a)') 'job '//trim(job_name)//' did not complete.'
                    stop
                else
                    write(*,'(a)') 'job '//trim(job_name)//' completed.'                    
                end if
                end_job = 1
                start_job = 0
            else
                nitems = nfields(line,' ')
                if(nitems > 2) then         
                    write(*,'(a)') line//' has more than two items'
                    stop
                end if
                call get_fields(line,' ',nitems,values)
                if(nitems==1) then
                    npar = npar + 1
                    parameters(npar) = values(1)
                else
                    parameters(npar+1) = values(1)
                    parameters(npar+2) = values(2)
                    npar = npar + 2
                end if              
            end if
        end do
        close(unit=10,status='keep')
    else
        call get_default_options
        call parse_cmdLine  
        write(*,'(a)') 'job '//trim(job_name)//' started...'                       
        call run_program_sodaBayes(status)
        if(status /= 0) then
            write(*,'(a)') 'job '//trim(job_name)//' did not complete.'
            stop
        else
            write(*,'(a)') 'job '//trim(job_name)//' completed.'                    
        end if        
    end if
    contains
  subroutine run_program_sodaBayes(status)
    implicit none
    integer, intent(inout) :: status
    real(real64) :: start, setup, end
    character(len=100) :: s1, s2, s3
    
    status = 1
    call date_and_time(date=cdate,time=ctime)
    start = omp_get_wtime()
    call omp_set_dynamic(.false.)
    maxthreads = omp_get_max_threads()
    call parse    
    open(unit=21,file=logfil,status='unknown',form='formatted')
    write(21,'(a)') 'Running program SODABayes: '//trim(job_name)
    write(21,101) 'Run started at',cdate(1:4),cdate(5:6),cdate(7:8),ctime(1:2),ctime(3:4),ctime(5:6)
    call get_size
    call allocate_data
    call load_geno
    call xscale
    if(mcmc) then
       call load_pheno
       call set_random_seed    
       if(SODAOFF) then
           setup = omp_get_wtime()
           call write_log           
           call regular_bayes          
       else
           call setup_soda
           setup = omp_get_wtime()
           call write_log           
           call soda_bayes
       end if
    else
       call setup_prediction
       call compute_dgv
       call write_dgv
       call write_log       
    end if
    call clean_up
    end = omp_get_wtime()
    call date_and_time(date=cdate,time=ctime)
    write(21,101) 'Run ended at',cdate(1:4),cdate(5:6),cdate(7:8),ctime(1:2),ctime(3:4),ctime(5:6)
    if(mcmc) then
        write(s1,'(f16.2)') setup-start
        write(s2,'(f16.2)') end-setup
        write(s3,'(f16.2)') end-start
        write(21,102) 'Walltime used (seconds)'
        write(21,'(a)') 'Setup: '//adjustl(trim(s1))//' MCMC: '//adjustl(trim(s2))//' Total: '//adjustl(trim(s3))
        close(21)
    end if
    status = 0
101 format(a20,1x,a4,'-',a2,'-',a2,' ',a2,':',a2':',a2)
102 format(a31)
  end subroutine run_program_sodaBayes
end program SODABayes
    
