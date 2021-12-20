 
module models
    use globals
    use utilities
    use parallelRNG
    use omp_lib
    
    implicit none
    contains
subroutine soda_bayes
  logical :: overflow
  character (len=12) :: c12
  character (len=15) :: c15
  integer :: i, j, k, jk, kk, jj, blocki, ls, le, ns, ne, nsize
  integer :: rep, snploc, included, thread, indistflag, counter, sg
  integer, dimension(:), allocatable :: gflag
  real(real64), dimension(:), allocatable :: gtemp  
  real(real64) :: logdetV, uhat, zz, zz_vare, rhs, v1, gk, gt, gdiff, scale
  real(real64) :: yhat, vary, sk, skk, r, ssculm, nnind, nnaind, clike
  real(real64), dimension(maxdist) :: s, stemp
  real(real64), pointer :: z(:), za(:)
  
  if(omp_get_dynamic()) call omp_set_dynamic(.false.)
  call omp_set_num_threads(nthreads)  
  open(unit=25,file=hypfil,status='unknown',form='formatted')
  if(BayesMethod .eq. BayesA) then
      write(*,'(A12,3(1x,A15))') 'Iteration', 'mu', 'Ve', 'Va'
      write(25,'(A12,3(1x,A15))') 'Iteration', 'mu', 'Ve', 'Va'
  else if(BayesMethod .eq. BayesB) then
      write(*,'((A12,1x,A12),4(1x,A15))') 'Iteration', 'ModelSize','mu', 'Ve', 'Va1', 'Va2'
      write(25,'((A12,1x,A12),4(1x,A15))') 'Iteration', 'ModelSize','mu', 'Ve', 'Va1', 'Va2'
  else  
      write(*,'((A12,1x,A12),3(1x,A15))',advance='no') 'Iteration', 'ModelSize','mu', 'Vara','Ve'
      write(25,'((A12,1x,A12),3(1x,A15))',advance='no') 'Iteration', 'ModelSize','mu', 'Vara','Ve'
      if(ngroup==1) then
          do i=1,nmix
              write(c12,'(I12)') i
              c12=adjustl(c12)
              c12="N"//trim(c12)
              c12=adjustr(c12)
              write(*,'(1x,A12)',advance="no") c12
              write(25,'(1x,A12)',advance="no") c12
          end do
          do i=1,nmix
              write(c15,'(I15)') i
              c15=adjustl(c15)
              c15="Va"//trim(c15)
              c15=adjustr(c15)
              write(*,'(1x,A15)',advance="no") c15
              write(25,'(1x,A15)',advance="no") c15
           end do
      end if
      write(*,*)
      write(25,*) 
  end if

  nnind=dble(nind)
  nnaind=dble(naind)+nnind

  if(BayesMethod .ne. BayesA) then 
      included_store=0d0
      snpstore=0d0
      indiststore=0d0  
      varstore=0.0d0
      if(BayesMethod .eq. BayesB .or. BayesMethod .eq. BayesC) then 
          do i=1,ngroup
              p(i,1)=pizero
              p(i,2)=1.0d0-pizero
          end do
      else
          do i=1,ngroup
              p(i,1)=0.5d0
              p(i,2:nmix)=1.0d0/mix(2:nmix)
              p(i,2:nmix)=0.5*p(i,2:nmix)/sum(p(i,2:nmix)) 
          end do
          pstore=0d0
      end if
  end if  
  if(BayesMethod .ne. BayesRc) then   
      snptracker(:,1)=1
  end if  
  snptracker(:,2)=2
  counter=0
  do i=1,nskipped
      j = i_skipped(i)
      xpx(j) = dot_product(X(:,j),X(:,j))
  enddo
  mu=sum(why)/nnind
  g=0.0d0 
  gstore=0d0
  mu_store=0d0  
  call compute_residuals
  if(df_va <=0.0d0) then
      yhat=sum(why)/nnind
      vary= sum((why-yhat)*(why-yhat))/(nnind-1.0d0)
      if(BayesMethod .eq. BayesA) then
          vara=2.0d0*scale_va*vary/nloci   !scale_va as heritability
      else if(BayesMethod .eq. BayesB) then
          vara=2.0d0*scale_va*vary/nloci/p(1,2) 
      else
          do i=1,ngroup
              vara = vara + dot_product(p(i,:),mix)
          end do
          vara=scale_va*vary/nloci/vara
      end if
      scale_va=0.0d0                   !set scale_va to 0
  else 
      vara=scale_va*df_va/(df_va+2.0d0) !mode
  end if
  if(BayesMethod .eq. BayesA .or. BayesMethod .eq. BayesB) then
      vara_s = vara      
      vara_sstore=0.0d0
  else
      gp=mix*vara      
      vara_store=0d0
  end if  
  scale=(dot_product(yadj,yadj)+scale_ve*df_ve)/ (nnind+df_ve)
  vare=rand_scaled_inverse_chi_square(nnind+df_ve,scale,0)
  vare_store=0d0
  do i=1,naind
     yAadj(i)=rand_normal(0.0d0,dsqrt(vare),0)
  end do 
  if(BayesMethod .eq. BayesA) then
       allocate(gtemp(blocksize))
       gtemp=0.0d0
  else
       allocate(gtemp(blocksize), gflag(blocksize))
       gtemp=0.0d0
       gflag=0
  end if
  do rep=1,numit
      yadj=yadj+mu
      mu=rand_normal(sum(yadj)/nnind, dsqrt(vare/nnind),0)
      yadj=yadj-mu
      if(BayesMethod .eq. BayesA) then            
          do blocki=1, nblocks
          if(blocki<=blockt) then  
              ls = (blocki-1)*blocksize+1
              le = min(ls+blocksize-1,nloci)
          else 
              ls = blocksize*blockt+(blocki-blockt-1)*(blocksize-1)+1
              le = min(ls+blocksize-2,nloci)
          end if 
          zz=xxd(blocki)
          zz_vare=zz/vare
          !$omp parallel private(thread,jk,snploc,z,za,rhs,v1,gk,gdiff,scale)
          !$omp do
          do k=ls,le
              if(skip_j(k)) cycle 
              jk=k-ls+1
              thread=omp_get_thread_num() 
              snploc=permvec(k)
              z => X(:,snploc)
              za => XA(:,k)
              gk=g(snploc)       
              rhs= dot_product(yadj,z) + dot_product(yAadj,za)
              rhs = rhs + zz*gk
              gtemp(jk)=gk
              v1=zz+vare/vara_s(snploc)
              gk=rand_normal(rhs/v1, dsqrt(vare/v1), thread)          
              g(snploc)=gk
              scale=(gk*gk + scale_va*df_va)/(df_va+1.0d0)
              vara_s(snploc)=rand_scaled_inverse_chi_square(df_va+1.0d0,scale,thread)                 
          enddo
          !$omp end do    
          !$omp do reduction(+:yadj,yAadj)
          do k=ls,le 
              if(skip_j(k)) cycle
              jk = k-ls+1
              snploc=permvec(k)
              z => X(:,snploc)
              za => XA(:,k)
              gdiff=gtemp(jk)-g(snploc)
              yadj=yadj+z*gdiff
              yAadj=yAadj+za*gdiff
          end do
          !$omp end do
          !$omp end parallel    
          enddo    
          do k=1,nskipped
              snploc = i_skipped(k)
              z => X(:,snploc)
              zz=xpx(snploc)
              zz_vare=zz/vare
              gk=g(snploc)
              rhs= dot_product(yadj,z)
              rhs=rhs+zz*gk
              v1=zz+vare/vara_s(snploc)
              gt=rand_normal(rhs/v1, dsqrt(vare/v1),0)
              gk = gk - gt
              yadj=yadj+z*gk  
              g(snploc)=gt
              scale=(gt*gt + scale_va*df_va)/(df_va+1.0d0)
              vara_s(snploc)=rand_scaled_inverse_chi_square(df_va+1.0d0,scale,0)                      
          end do  
          varindist(1,1)=sum(g*g)
      else if(BayesMethod .eq. BayesB) then      
          log_p(1,1)=dlog(p(1,1))          
          log_p(1,2)=dlog(p(1,2))
          do blocki=1, nblocks
              if(blocki<=blockt) then  
                  ls = (blocki-1)*blocksize+1
                  le = min(ls+blocksize-1,nloci)
              else 
                  ls = blocksize*blockt+(blocki-blockt-1)*(blocksize-1)+1
                  le = min(ls+blocksize-2,nloci)
              end if 
              zz=xxd(blocki)
              zz_vare=zz/vare
              !$omp parallel private(jk,thread,snploc,z,za,gk,rhs,s,sk,logdetV,uhat,v1,gdiff,r,clike,indistflag,scale)
              !$omp do
              do k=ls,le
                  if(skip_j(k)) cycle 
                  jk=k-ls+1
                  thread=omp_get_thread_num() 
                  snploc=permvec(k)
                  z => X(:,snploc)
                  za => XA(:,k)
                  gk=g(snploc)       
                  rhs= dot_product(yadj,z) + dot_product(yAadj,za)
                  if(snptracker(snploc,2) > 1) then
                     rhs = rhs + zz*gk
                  endif
                  s(1)=log_p(1,1)
                  logdetV=dlog(vara_s(snploc)*zz_vare+1.0d0)
                  uhat=rhs/(zz+vare/vara_s(snploc))
                  s(2)=-0.5d0*(logdetV-(rhs*uhat/vare))+log_p(1,2)
                  clike=s(1)-s(2)
                  if (clike .gt. 700) then 
                      sk = 0.0d0
                  else if(clike .lt. -700) then
                      sk = 1.0d0
                  else
                      sk=1.0d0/(1.0d0+dexp(clike))
                  end if          
                  call random(r,thread)
                  if (r<sk) then
                      indistflag=2
                  else
                      indistflag=1
                  endif
                  gtemp(jk)=gk
                  if(indistflag==1) then
                      gk=rand_normal(0.0d0, dsqrt(vara_s(snploc)),thread)
                      g(snploc)=0.0d0
                      if(snptracker(snploc,2)==1) then
                          gflag(jk) = 1 
                      end if
                  else
                      v1=zz+vare/vara_s(snploc)
                      gk=rand_normal(rhs/v1, dsqrt(vare/v1), thread) 
                      g(snploc)=gk
                  endif
                  snptracker(snploc,2)=indistflag
                  scale=(gk*gk + scale_va*df_va)/(df_va+1.0d0)
                  vara_s(snploc)=rand_scaled_inverse_chi_square(df_va+1.0d0,scale,thread)         
              enddo
              !$omp end do    
              !$omp do reduction(+:yadj,yAadj)
              do k=ls,le 
                  if(skip_j(k)) cycle
                  jk = k-ls+1
                  if(gflag(jk) ==1) then
                      gflag(jk) = 0
                      cycle
                  end if
                  snploc=permvec(k)
                  z => X(:,snploc)
                  za => XA(:,k)
                  gdiff=gtemp(jk)-g(snploc)
                  yadj=yadj+z*gdiff
                  yAadj=yAadj+za*gdiff        
              end do
              !$omp end do
              !$omp end parallel     
          enddo
          do k=1,nskipped
              snploc = i_skipped(k)
              z => X(:,snploc)
              zz=xpx(snploc)
              zz_vare=zz/vare
              gk=g(snploc)
              rhs= dot_product(yadj,z)
              if(snptracker(snploc,2) > 1) then
                  rhs=rhs+zz*gk
              endif          
              s(1)=log_p(1,1)
              logdetV=dlog(vara_s(snploc)*zz_vare+1.0d0)
              uhat=rhs/(zz+vare/vara_s(snploc))
              s(2)=-0.5d0*(logdetV-(rhs*uhat/vare))+log_p(1,2)
              clike=s(1)-s(2)
              if (clike .gt. 700) then 
                  sk = 0.0d0
              else if(clike .lt. -700) then
                  sk = 1.0d0
              else
                  sk=1.0d0/(1.0d0+dexp(clike))
              end if          
              call random(r,0)
              if (r<sk) then
                  indistflag=2
              else
                  indistflag=1
              endif
              if(indistflag==1) then
                  gt=rand_normal(0.0d0, dsqrt(vara_s(snploc)),0)
                  if(snptracker(snploc,2)>1) then
                      yadj=yadj+z*gk
                  endif
                  g(snploc)=0.0d0
              else
                  v1=zz+vare/vara_s(snploc)
                  gt=rand_normal(rhs/v1, dsqrt(vare/v1),0)
                  gk = gk - gt
                  yadj=yadj+z*gk 
                  g(snploc)=gt
              endif
              snptracker(snploc,2)=indistflag  
              scale=(gt*gt + scale_va*df_va)/(df_va+1.0d0)
              vara_s(snploc)=rand_scaled_inverse_chi_square(df_va+1.0d0,scale,0)         
          end do
          included = count(snptracker(:,2)==2)
          snpindist(1,1) = nloci - included
          snpindist(1,2) = included
          varindist(1,1) = 0.0d0
          varindist(1,2) = sum(g*g, snptracker(:,2)==2)
      else 
          do i=1,ngroup
              log_p(i,1)=dlog(p(i,1))          
              do j=2,nmix
                  log_p(i,j)=dlog(p(i,j))
              end do
          enddo
          do i=2,nmix
              vare_gp(i)=vare/gp(i)
          end do
          do blocki=1, nblocks
              if(blocki<=blockt) then  
                  ls = (blocki-1)*blocksize+1
                  le = min(ls+blocksize-1,nloci)
              else 
                  ls = blocksize*blockt+(blocki-blockt-1)*(blocksize-1)+1
                  le = min(ls+blocksize-2,nloci)
              end if 
              zz=xxd(blocki)
              zz_vare=zz/vare
              !$omp parallel &
              !$omp private(thread,snploc,logdetV,uhat,z,za,rhs,v1,gk,gdiff,s,r) &
              !$omp private(stemp,sk,skk,jk,clike,overflow,ssculm,indistflag,sg)
              !$omp do
              do k=ls,le
                  if(skip_j(k)) cycle 
                  jk=k-ls+1
                  thread=omp_get_thread_num() 
                  snploc=permvec(k)
                  z => X(:,snploc)
                  za => XA(:,k)
                  gk=g(snploc)       
                  rhs= dot_product(yadj,z) + dot_product(yAadj,za)
                  if(snptracker(snploc,2) > 1) then
                      rhs = rhs + zz*gk
                  endif
                  sg=snptracker(snploc,1)
                  s(1)=log_p(sg,1)
                  do kk=2,nmix
                      logdetV=dlog(gp(kk)*zz_vare+1.0d0)
                      uhat=rhs/(zz+vare_gp(kk))
                      s(kk)=-0.5d0*(logdetV-(rhs*uhat/vare))+log_p(sg,kk)
                  enddo
                  stemp=0.0d0
                  do kk=1,nmix
                      skk=s(kk)
                      sk=0.0d0
                      overflow=.false.
                      do j=1,nmix
                          if(j==kk) cycle
                          clike=s(j)-skk
                          if(clike .lt. -700) then
                              cycle
                          else if (clike .gt. 700) then
                              overflow=.true.
                              exit
                          endif
                          sk=sk+dexp(clike)
                      enddo
                      if(overflow .eqv. .true.) then
                          stemp(kk) = 0.0d0
                      else
                          stemp(kk)=1.0d0/(1.0d0+sk)
                      endif
                  enddo   
                  ssculm=0.0d0  
                  call random(r,thread)
                  indistflag=1
                  do kk=1,nmix
                      ssculm=ssculm+stemp(kk)
                      if(r<ssculm) then
                          indistflag=kk
                          exit
                      endif
                  enddo
                  gtemp(jk)=gk
                  if(indistflag==1) then
                      if(snptracker(snploc,2)==1) then
                          gflag(jk) = 1 
                      end if
                      gk=0.0d0            
                  else
                      v1=zz+vare/gp(indistflag)
                      gk=rand_normal(rhs/v1, dsqrt(vare/v1), thread)          
                  endif
                  snptracker(snploc,2)=indistflag
                  g(snploc)=gk
              enddo
              !$omp end do    
              !$omp do reduction(+:yadj,yAadj)
              do k=ls,le 
                  if(skip_j(k)) cycle
                  jk = k-ls+1
                  if(gflag(jk) ==1) then
                      gflag(jk) = 0
                      cycle
                  end if
                  snploc=permvec(k)
                  z => X(:,snploc)
                  za => XA(:,k)
                  gdiff=gtemp(jk)-g(snploc)
                  yadj=yadj+z*gdiff
                  yAadj=yAadj+za*gdiff        
              end do
              !$omp end do
              !$omp end parallel 
          enddo
          do k=1,nskipped
              snploc = i_skipped(k)
              z => X(:,snploc)
              zz=xpx(snploc)
              zz_vare=zz/vare
              gk=g(snploc)
              rhs= dot_product(yadj,z)
              if(snptracker(snploc,2) > 1) then
                  rhs=rhs+zz*gk
              endif         
              sg=snptracker(snploc,1)
              s(1)=log_p(sg,1)
              do kk=2,nmix
                  logdetV=dlog(gp(kk)*zz_vare+1.0d0)
                  uhat=rhs/(zz+vare_gp(kk))
                  s(kk)=-0.5d0*(logdetV-(rhs*uhat/vare))+log_p(sg,kk)
              enddo
              stemp=0.0d0
              do kk=1,nmix
                  skk=s(kk)
                  sk=0.0d0
                  overflow=.false.
                  do j=1,nmix
                      if(j==kk) cycle
                      clike=s(j)-skk
                      if(clike .lt. -700) then
                          cycle
                      else if (clike .gt. 700) then 
                          overflow=.true.
                          exit
                      endif
                      sk=sk+dexp(clike)    
                  enddo
                  if (overflow .eqv. .true.) then
                      stemp(kk) = 0.0
                  else
                      stemp(kk)=1.0d0/(1.0d0+sk)
                  endif
              enddo
              ssculm=0.0d0
              call random(r,0)
              indistflag=1
              do kk=1,nmix
                  ssculm=ssculm+stemp(kk)
                  if (r<ssculm) then
                      indistflag=kk
                      exit
                  endif
              enddo
              if(indistflag==1) then
                  gt=0.0d0
                  if(snptracker(snploc,2)>1) then
                      yadj=yadj+z*gk
                  endif
              else
                  v1=zz+vare/gp(indistflag)
                  gt=rand_normal(rhs/v1, dsqrt(vare/v1),0)
                  gk = gk - gt
                  yadj=yadj+z*gk  
              endif
              g(snploc)=gt
              snptracker(snploc,2)=indistflag          
          end do
          included=0
          do i=1,ngroup
              do j=1,nmix
                  snpindist(i,j)=count(snptracker(:,1)==i .and. snptracker(:,2)==j)
                  varindist(i,j)=sum(g*g, mask= snptracker(:,1)==i .and. snptracker(:,2)==j)
              enddo
              included=included+snpindist(i,1)          
          end do
          included = nloci - included
          scale = 0.0d0
          do i=1,ngroup
              do j=2,nmix
                  scale = scale + varindist(i,j) / mix(j)
              end do 
          end do
          scale=(scale + scale_va*df_va)/(df_va+dble(included))
          vara=rand_scaled_inverse_chi_square(dble(included)+df_va,scale,0)
          gp=mix*vara
      end if
      scale=(dot_product(yadj,yadj)+scale_ve*df_ve)/ (nnind+df_ve)
      vare=rand_scaled_inverse_chi_square(nnind+df_ve,scale,0)
      if(BayesMethod .ne. BayesA .and. BayesMethod .ne. BayesB .and. BayesMethod .ne. BayesC) then
          do i=1,ngroup
                  dirx(i,:)=dble(snpindist(i,:))+alpha
                  p(i,:)=rand_dirichlet(nmix,dirx(i,:),0)
          end do          
      end if    
      do i=1,naind
          yAadj(i)=rand_normal(0.0d0,dsqrt(vare),0)
      end do    
      if(rep>burnin .and. mod(rep,thin)==0) then
         counter=counter+1
         gstore=gstore+g
         mu_store=mu_store+mu
         varstore=varstore+varindist
         vare_store=vare_store+vare
         if(BayesMethod .ne. BayesA) then               
             included_store=included_store+included
             snpstore=snpstore+snpindist
             do i=1,nloci
                 jj=snptracker(i,2)
                 indiststore(i,jj)=indiststore(i,jj)+1
             enddo 
             if(BayesMethod .ne. BayesB .and. BayesMethod .ne. BayesC) then
                 pstore=pstore+p
             end if
         end if
         if(BayesMethod .eq. BayesA .or. BayesMethod .eq. BayesB) then
             vara_sstore=vara_sstore+vara_s
         else       
             vara_store=vara_store+vara    
         end if           
      end if
      if(mod(rep,thin)==0) then
         if(BayesMethod .eq. BayesA) then
             write(*,'(I12,3(1x,E15.6))')  rep, mu, vare, varindist
             write(25,'(I12,3(1x,E15.6))')  rep, mu, vare, varindist
         else if(BayesMethod .eq. BayesB) then
             write(*,'((I12,1x,I12),4(1x,E15.6))')  rep, included, mu, vare, varindist
             write(25,'((I12,1x,I12),4(1x,E15.6))')  rep, included, mu, vare, varindist
         else 
             write(*,'((I12,1x,I12),3(1x,E15.6))',advance='no')  rep, included, mu, vara, vare
             write(25,'((I12,1x,I12),3(1x,E15.6))',advance='no')  rep, included, mu, vara, vare
             if(ngroup==1) then
                 write(*,'(20(1x,I12))',advance='no')  snpindist
                 write(*,'(20(1x,E15.6))',advance='no') varindist
                 write(25,'(20(1x,I12))',advance='no')  snpindist
                 write(25,'(20(1x,E15.6))',advance='no') varindist
             end if
             write(*,*)
             write(25,*)
         end if
         call flush(25)
      end if      
      ! re-calibrate residuals
      if(mod(rep,1000)==0) then
         call compute_residuals
      endif
  end do !end do mcmc
   !posterior means
   gstore=gstore/counter
   mu_store=mu_store/counter
   varstore=varstore/counter
   vare_store=vare_store/counter
   if(BayesMethod .ne. BayesA) then
       included_store=included_store/counter
       snpstore=snpstore/counter
       do i=1,nloci
           indiststore(i,:)=indiststore(i,:)/counter 
       enddo
       if(BayesMethod .ne. BayesB .and. BayesMethod .ne. BayesC) then
           pstore=pstore/counter
       end if
   end if
   if(BayesMethod .eq. BayesA .or. BayesMethod .eq. BayesB) then  
       vara_sstore=vara_sstore/counter 
   else     
       vara_store=vara_store/counter
   end if
   
   call output_model
   mu=mu_store
   call compute_dgv
   call write_dgv
   if(BayesMethod .eq. BayesA) then
       deallocate(gtemp)
   else 
       deallocate(gtemp, gflag)  
   end if
end subroutine soda_bayes   

subroutine regular_bayes()
  character (len=12) :: c12
  character (len=15) :: c15
  logical :: overflow
  integer :: i, j, k, kk, jj, b, rep, snploc, included, indistflag, counter, sg
  real(real64) :: logdetV, uhat, zz, zz_vare, rhs, v1, gk, gt, scale
  real(real64) :: yhat, vary, sk, skk, r, ssculm, nnind, clike
  real(real64), dimension(maxdist) :: s, stemp
  real(real64), pointer :: z(:)
  
   nnind=dble(nind)
   open(unit=25,file=hypfil,status='unknown',form='formatted')
   if(BayesMethod .eq. BayesA) then
       write(*,'(A12,3(1x,A15))') 'Iteration', 'mu', 'Ve', 'Va'
       write(25,'(A12,3(1x,A15))') 'Iteration', 'mu', 'Ve', 'Va'
   else if(BayesMethod .eq. BayesB) then
       write(*,'((A12,1x,A12),4(1x,A15))') 'Iteration', 'ModelSize','mu', 'Ve', 'Va1', 'Va2'
       write(25,'((A12,1x,A12),4(1x,A15))') 'Iteration', 'ModelSize','mu', 'Ve', 'Va1', 'Va2'
   else  
       write(*,'((A12,1x,A12),3(1x,A15))',advance='no') 'Iteration', 'ModelSize','mu', 'Vara','Ve'
       write(25,'((A12,1x,A12),3(1x,A15))',advance='no') 'Iteration', 'ModelSize','mu', 'Vara','Ve'
       if(ngroup==1) then
           do i=1,nmix
               write(c12,'(I12)') i
               c12=adjustl(c12)
               c12="N"//trim(c12)
               c12=adjustr(c12)
               write(*,'(1x,A12)',advance="no") c12
               write(25,'(1x,A12)',advance="no") c12
           end do
           do i=1,nmix
               write(c15,'(I15)') i
               c15=adjustl(c15)
               c15="Va"//trim(c15)
               c15=adjustr(c15)
               write(*,'(1x,A15)',advance="no") c15
               write(25,'(1x,A15)',advance="no") c15
            end do
       end if
       write(*,*)
       write(25,*) 
   end if
   
   if(BayesMethod .ne. BayesA) then 
      included_store=0d0
      snpstore=0d0
      indiststore=0d0  
      varstore=0.0d0
      if(BayesMethod .eq. BayesB .or. BayesMethod .eq. BayesC) then 
          do i=1,ngroup
              p(i,1)=pizero
              p(i,2)=1.0d0-pizero
          end do
      else
          do i=1,ngroup
              p(i,1)=0.5d0
              p(i,2:nmix)=1.0d0/mix(2:nmix)
              p(i,2:nmix)=0.5*p(i,2:nmix)/sum(p(i,2:nmix)) 
          end do
          pstore=0d0
      end if
   end if  
   if(BayesMethod .ne. BayesRc) then   
      snptracker(:,1)=1
   end if  
   snptracker(:,2)=2
   counter=0
   do i=1,nloci
      xpx(i)=dot_product(X(:,i),X(:,i))
      permvec(i)=i
   enddo
   mu=sum(why)/nnind
   g=0.0d0 
   gstore=0d0
   mu_store=0d0  
   call compute_residuals
   if(df_va <=0.0d0) then
       yhat=sum(why)/nnind
       vary= sum((why-yhat)*(why-yhat))/(nnind-1.0d0)
       if(BayesMethod .eq. BayesA) then
           vara=2.0d0*scale_va*vary/nloci   !scale_va as heritability
       else if(BayesMethod .eq. BayesB) then
           vara=2.0d0*scale_va*vary/nloci/p(1,2)
       else
           do i=1,ngroup
               vara = vara + dot_product(p(i,:),mix)
           end do
           vara=scale_va*vary/nloci/vara
       end if
       scale_va=0.0d0                   !set scale_va to 0
   else 
       vara=scale_va*df_va/(df_va+2.0d0) !mode
   end if
   if(BayesMethod .eq. BayesA .or. BayesMethod .eq. BayesB) then
       vara_s = vara      
       vara_sstore=0.0d0
   else       
       gp=mix*vara      
       vara_store=0d0
   end if  
   scale=(dot_product(yadj,yadj)+scale_ve*df_ve)/ (nnind+df_ve)
   vare=rand_scaled_inverse_chi_square(nnind+df_ve,scale,0)
   vare_store=0d0

   each_cycle : do rep=1,numit
      yadj=yadj+mu
      mu=rand_normal(sum(yadj)/nnind, dsqrt(vare/nnind),0)
      yadj=yadj-mu
      if(BayesMethod .eq. BayesA) then       
          do k=1,nloci
              snploc=permvec(k)
              z => X(:,snploc)
              zz=xpx(snploc)
              zz_vare=zz/vare
              gk=g(snploc)
              rhs= dot_product(yadj,z)
              rhs=rhs+zz*gk
              v1=zz+vare/vara_s(snploc)
              gt=rand_normal(rhs/v1, dsqrt(vare/v1),0)
              gk = gk - gt
              yadj=yadj+z*gk  
              g(snploc)=gt
              scale=(gt*gt + scale_va*df_va)/(df_va+1.0d0)
              vara_s(snploc)=rand_scaled_inverse_chi_square(df_va+1.0d0,scale,0)        
          end do
      else if(BayesMethod .eq. BayesB) then 
          log_p(1,1)=dlog(p(1,1))
          log_p(1,2)=dlog(p(1,2))
          do k=1,nloci
              snploc=permvec(k)
              z => X(:,snploc)
              zz=xpx(snploc)
              zz_vare=zz/vare
              gk=g(snploc)
              rhs= dot_product(yadj,z)
              if(snptracker(snploc,2) > 1) then
                  rhs=rhs+zz*gk
              endif         
              s(1)=log_p(1,1)
              logdetV=dlog(vara_s(snploc)*zz_vare+1.0d0)
              uhat=rhs/(zz+vare/vara_s(snploc))
              s(2)=-0.5d0*(logdetV-(rhs*uhat/vare))+log_p(1,2)
              clike=s(1)-s(2)
              if (clike .gt. 700) then 
                 sk = 0.0d0
              else if(clike .lt. -700) then
                 sk = 1.0d0
              else
                 sk=1.0d0/(1.0d0+dexp(clike))
              end if          
              call random(r,0)
              if (r<sk) then
                 indistflag=2
              else
                 indistflag=1
              endif
              if(indistflag==1) then
                 gt=rand_normal(0.0d0, dsqrt(vara_s(snploc)),0)
                 if(snptracker(snploc,2)>1) then
                     yadj=yadj+z*gk
                 endif
                 g(snploc)=0.0d0
              else
                 v1=zz+vare/vara_s(snploc)
                 gt=rand_normal(rhs/v1, dsqrt(vare/v1),0)
                 gk = gk - gt
                 yadj=yadj+z*gk            
                 g(snploc)=gt
              endif
              snptracker(snploc,2)=indistflag
              scale=(gt*gt + scale_va*df_va)/(df_va+1.0d0)
              vara_s(snploc)=rand_scaled_inverse_chi_square(df_va+1.0d0,scale,0)         
          enddo
          included = count(snptracker(:,2)==2)
          snpindist(1,1) = nloci - included
          snpindist(1,2) = included
          varindist(1,1) = 0.0d0
          varindist(1,2) = sum(g*g, snptracker(:,2)==2)
      else
          do i=1,ngroup
              log_p(i,1)=dlog(p(i,1))
              do j=2,nmix
                  log_p(i,j)=dlog(p(i,j))
              enddo
          end do
          do i=2,nmix
              vare_gp(i)=vare/gp(i)
          end do
          do k=1,nloci
              snploc=permvec(k)
              z => X(:,snploc)
              zz=xpx(snploc)
              zz_vare=zz/vare
              gk=g(snploc)
              rhs= dot_product(yadj,z)
              if(snptracker(snploc,2) > 1) then
                  rhs=rhs+zz*gk
              endif         
              sg=snptracker(snploc,1)
              s(1)=log_p(sg,1)
              do kk=2,nmix
                 logdetV=dlog(gp(kk)*zz_vare+1.0d0)
                 uhat=rhs/(zz+vare_gp(kk))
                 s(kk)=-0.5d0*(logdetV-(rhs*uhat/vare))+log_p(sg,kk)
              enddo
              stemp=0.0d0
              do kk=1,nmix
                 skk=s(kk)
                 sk=0.0d0
                 overflow=.false.
                 do j=1,nmix
                    if(j==kk) cycle
                    clike=s(j)-skk
                    if(clike .lt. -700) then
                       cycle
                    else if (clike .gt. 700) then 
                       overflow=.true.
                       exit
                    endif
                    sk=sk+dexp(clike)
                 enddo
                 if (overflow .eqv. .true.) then
                    stemp(kk) = 0.0
                 else
                    stemp(kk)=1.0d0/(1.0d0+sk)
                 endif
              enddo
              ssculm=0.0d0
              call random(r,0)
              indistflag=1
              do kk=1,nmix
                 ssculm=ssculm+stemp(kk)
                 if (r<ssculm) then
                    indistflag=kk
                    exit
                 endif
              enddo
              if(indistflag==1) then                
                 gt=0.0d0
                 if(snptracker(snploc,2)>1) then
                     yadj=yadj+z*gk
                 endif
              else
                 v1=zz+vare/gp(indistflag)
                 gt=rand_normal(rhs/v1, dsqrt(vare/v1),0)
                 gk = gk - gt
                 yadj=yadj+z*gk  
              endif
              g(snploc)=gt
              snptracker(snploc,2)=indistflag
           enddo
           included=0
           do i=1,ngroup
               do j=1,nmix
                  snpindist(i,j)=count(snptracker(:,1)==i .and. snptracker(:,2)==j)
                  varindist(i,j)=sum(g*g, mask= snptracker(:,1)==i .and. snptracker(:,2)==j)
               enddo
               included=included+snpindist(i,1)          
           end do
           included = nloci - included      
           scale = 0.0d0
           do i=1,ngroup
               do j=2,nmix
                   scale = scale + varindist(i,j) / mix(j)
               end do 
           end do
           scale=(scale + scale_va*df_va)/(df_va+dble(included))
           vara=rand_scaled_inverse_chi_square(dble(included)+df_va,scale,0)
           gp=mix*vara
      end if
      scale=(dot_product(yadj,yadj)+scale_ve*df_ve)/ (nnind+df_ve)
      vare=rand_scaled_inverse_chi_square(nnind+df_ve,scale,0)
      if(BayesMethod .ne. BayesA .and. BayesMethod .ne. BayesB .and. BayesMethod .ne. BayesC) then
          do i=1,ngroup
              dirx(i,:)=dble(snpindist(i,:))+alpha
              p(i,:)=rand_dirichlet(nmix,dirx(i,:),0)
          end do
      end if
      
      if(rep>burnin .and. mod(rep,thin)==0) then
         counter=counter+1
         gstore=gstore+g
         mu_store=mu_store+mu
         vare_store=vare_store+vare
         if(BayesMethod .ne. BayesA) then               
             included_store=included_store+included
             snpstore=snpstore+snpindist
             varstore=varstore+varindist
             do i=1,nloci
                 jj=snptracker(i,2)
                 indiststore(i,jj)=indiststore(i,jj)+1
             enddo 
             if(BayesMethod .ne. BayesB .and. BayesMethod .ne. BayesC) then
                  pstore=pstore+p
             end if
         end if
         if(BayesMethod .eq. BayesA .or. BayesMethod .eq. BayesB) then
             vara_sstore=vara_sstore+vara_s
         else       
             vara_store=vara_store+vara    
         end if           
      end if
      if(mod(rep,thin)==0) then
         if(BayesMethod .eq. BayesA) then
             write(*,'(I12,3(1x,E15.6))')  rep, mu, vare, varindist
             write(25,'(I12,3(1x,E15.6))')  rep, mu, vare, varindist
         else if(BayesMethod .eq. BayesB) then
             write(*,'((I12,1x,I12),4(1x,E15.6))')  rep, included, mu, vare, varindist
             write(25,'((I12,1x,I12),4(1x,E15.6))')  rep, included, mu, vare, varindist
         else 
             write(*,'((I12,1x,I12),3(1x,E15.6))',advance='no')  rep, included, mu, vara, vare
             write(25,'((I12,1x,I12),3(1x,E15.6))',advance='no')  rep, included, mu, vara, vare
             if(ngroup==1) then
                 write(*,'(20(1x,I12))',advance='no')  snpindist
                 write(*,'(20(1x,E15.6))',advance='no') varindist
                 write(25,'(20(1x,I12))',advance='no')  snpindist
                 write(25,'(20(1x,E15.6))',advance='no') varindist
             end if
             write(*,*)
             write(25,*)
         end if
         call flush(25)
      end if       
      if(mod(rep,1000)==0) then
         call compute_residuals
      endif
   enddo each_cycle

   !posterior means
   gstore=gstore/counter
   mu_store=mu_store/counter
   vare_store=vare_store/counter
   if(BayesMethod .ne. BayesA) then
       included_store=included_store/counter
       varstore=varstore/counter
       snpstore=snpstore/counter
       do i=1,nloci
           indiststore(i,:)=indiststore(i,:)/counter 
       enddo
       if(BayesMethod .ne. BayesB .and. BayesMethod .ne. BayesC) then
           pstore=pstore/counter
       end if
   end if
   if(BayesMethod .eq. BayesA .or. BayesMethod .eq. BayesB) then  
       vara_sstore=vara_sstore/counter 
   else     
       vara_store=vara_store/counter
   end if
   
   call output_model
   mu=mu_store
   call compute_dgv
   call write_dgv
end subroutine regular_bayes
end module models
