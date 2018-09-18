  subroutine dyn_susc_expansion
! Dynamic KS susceptibility
! Using density basis
! No fit of GF
  use global

  implicit none

! -----------------------------------------------------------------
!   i 2pi
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0)
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi)
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0), cminus = (-1.d0,0.d0)
  complex(kind=c8b), parameter :: e0up = (0.481995d0,-0.05d0), e0dn = (0.711995,-0.007d0)
! -----------------------------------------------------------------
  complex(kind=c8b), allocatable :: susc0ylm(:,:,:,:,:), suscylm(:,:,:,:,:), susc0inv(:,:,:,:,:), suscinv(:,:,:,:,:)
  complex(kind=c8b) :: zw(nomega)
  integer(kind=i4b) :: ib, jb, is, js, i4(4), ia, ja, ia2, ja2, jlm, ilm, j, i, iw, iq, jq, a
  integer(kind=i4b) :: ipiv(ndensum), info, nevals, lwork
  real(kind=r8b)    :: maxelem, start, finish
  complex(kind=c8b) :: block(4,4)
  character*100     :: filename
  complex(kind=c8b), allocatable :: work(:), tmp(:,:), evals(:), evecl(:,:), evecr(:,:,:)
  real(kind=r8b),    allocatable :: rwork(:)!, evals(:)
!  real(kind=r8b)    :: r(nalmsb), c(nalmsb), rwork(2*nalmsb), rcond, ferr(nalmsb), berr(nalmsb)
!  complex(kind=c8b) :: work(2*nalmsb), temp(nalmsb,nalmsb)
!  complex(kind=c8b) :: x(nalmsb,nalmsb)
!  character*1       :: equed
  complex(kind=c8b), allocatable :: curr_corr_lm(:,:,:,:,:,:,:,:),curr_corr_0_lm(:,:,:,:,:,:,:,:)
  complex(kind=c8b), allocatable :: curr_corr_div_lm(:,:,:,:,:,:,:),curr_corr_div_0_lm(:,:,:,:,:,:,:)

  call cpu_time(start)
! Equidistant frequencies
  do iw=1,nomega
    zw(iw) = omegamin + (omegamax-omegamin)*(iw-1.d0)/(nomega-1.d0)
  end do
  open(file='evals.dat',unit=iofile,status='replace')
! eigenvalues
  nevals = 2*nasusc2 !4*lmmax0*nasusc2
  lwork  = 2*nevals
  allocate(tmp(nevals,nevals),work(lwork),rwork(2*lwork))
  allocate(evals(nevals),evecl(nevals,nevals),evecr(nevals,nevals,nomega))
  allocate(susc0ylm(4,4,nasusc2,nasusc2,nomega),susc0inv(2,2,nasusc2,nasusc2,nomega))
  allocate(suscylm(4,4,nasusc2,nasusc2,nomega),suscinv(2,2,nasusc2,nasusc2,nomega))

! Added by Sascha
  if(lcurrcorr) allocate(curr_corr_0_lm(1:3,1:4,1:4,1:lmmax2,1:lmmax2,1:nasusc2,1:nasusc2,1:nomega))
  if(lcurrcorr) allocate(curr_corr_lm(1:3,1:4,1:4,1:lmmax2,1:lmmax2,1:nasusc2,1:nasusc2,1:nomega))
  if(lcurrcorr) allocate(curr_corr_div_lm(1:4,1:4,1:lmmax2,1:lmmax2,1:nasusc2,1:nasusc2,1:nomega))
  if(lcurrcorr) allocate(curr_corr_div_0_lm(1:4,1:4,1:lmmax2,1:lmmax2,1:nasusc2,1:nasusc2,1:nomega))

! ********************
! Frequency loop
  do iw=1,nomega
! ********************
    write(*,'("iw=",i4,", zw=",2es16.8)') iw, zw(iw)
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!             Kohn-Sham susceptibility for current frequency
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   Put together KS susc from Taylor expansion
    ! insert kssusc0, kssusc1 and kssusc 2 into kssusc
    call zcopy(ndensum*ndensum,kssusc0,1,kssusc,1)
    call zaxpy(ndensum*ndensum,zw(iw),kssusc1,1,kssusc,1)
    call zaxpy(ndensum*ndensum,0.5d0*zw(iw)**2,kssusc2,1,kssusc,1)
!   multipoles
    susc0ylm(:,:,:,:,iw) = czero; susc0inv(:,:,:,:,iw) = czero
    ! iq and jq loop on kssusc, which so far contains transverse chi_KS with two entries (+- and -+)
    do jq=1,ndensum
      i4 = i2almsbden(:,jq)
      jb = i4(1); jlm = i4(2); js = i4(3); ja2 = i4(4)
      do iq=1,ndensum
        i4 = i2almsbden(:,iq)
        ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4)
        if (ilm == 1 .and. jlm == 1) susc0ylm(is,js,ia2,ja2,iw) = susc0ylm(is,js,ia2,ja2,iw) + suscnorm(iq)*kssusc(iq,jq)*suscnorm(jq)
      end do
    end do
    if (lcartesian) then
      do ja2=1,nasusc2
      do ia2=1,nasusc2
        do jlm=1,1!lmmax0
        do ilm=1,1!lmmax0
          block(:,:) = czero
          do j=1,4
          do i=1,4
            do js=1,4
            do is=1,4
              block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*susc0ylm(is,js,ia2,ja2,iw)*pc2s(i2is(1,js),i2is(2,js),j)
            end do
            end do
          end do
          end do
          susc0ylm(:,:,ia2,ja2,iw) = block
        end do
        end do
      end do
      end do
    end if
!   final filtering
    maxelem = maxval(abs(susc0ylm(:,:,:,:,iw)))
    where (abs(susc0ylm(:,:,:,:,iw)) < susctol*maxelem) susc0ylm(:,:,:,:,iw) = 0.d0
!   --------------------------------------------------------------------
!   inverse of spherical susceptibility
    jq = 1
    do ja2=1,nasusc2
    do jlm=1,1!lmmax0
    do j=1,2
      iq = 1
      do ia2=1,nasusc2
      do ilm=1,1!lmmax0
      do i=1,2
!        tmp(iq,jq) = (suscylm(i,j,ilm,jlm,ia2,ja2) - conjg(suscylm(j,i,jlm,ilm,ja2,ia2)))/i2pi
        tmp(iq,jq) = susc0ylm(i,j,ia2,ja2,iw)
        iq = iq + 1
      end do
      end do
      end do
      jq = jq + 1
    end do
    end do
    end do
    ! the below zgetrf and zgetri are used to invert the matrix tmp
    call zgetrf(nevals,nevals,tmp,nevals,ipiv,info)
    if (info /= 0) stop 'dyn_susc_expansion: failure in zgetrf'
    call zgetri(nevals,tmp,nevals,ipiv,work,lwork,info)
    if (info /= 0) stop 'dyn_susc_expansion: failure in zgetri'
    jq = 1
    do ja2=1,nasusc2
    do jlm=1,1!lmmax0
    do j=1,2
      iq = 1
      do ia2=1,nasusc2
      do ilm=1,1!lmmax0
      do i=1,2
        susc0inv(i,j,ia2,ja2,iw) = tmp(iq,jq)
        iq = iq + 1
      end do
      end do
      end do
      jq = jq + 1
    end do
    end do
    end do
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               Full susceptibility for current frequency
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   Solve Dyson equation
    if (lenhanced) then
      ! the below basically is a matrix multiplication between kssusc and kernel, output written/ to denominator
      call zgemm('N','N',ndensum,ndensum,ndensum,cminus,kssusc,ndensum,kernel,ndensum,czero,denominator,ndensum)
      do iq=1,ndensum
        denominator(iq,iq) = denominator(iq,iq) + cone
      end do
      ! this solves A*X=B linearly coupled equation, where 
      ! on input     A=denominator, B=kssusc
      ! on output    denominator= weird non-useful matrix, kssusc=X (solution)
      ! note that on output, kssusc contains the ENHANCED susceptibility, not KS (misleading...)
      call zgesv(ndensum,ndensum,denominator,ndensum,ipiv,kssusc,ndensum,info)
      if (info /= 0) stop 'dyn_susc_expansion: failure in zgesv'
!      call zgesvx('N','N',nalmsb,nalmsb,denominator,nalmsb,temp,nalmsb,ipiv,equed,r,c,kssusc0,nalmsb,x,nalmsb,rcond,ferr,berr,work,rwork,info)
!      write(iodb,'("condition number of enhancement factor=",es16.3)') rcond
!      write(iodb,'("fwd & bkwd error in solution=",2es16.3)') maxval(abs(ferr)), maxval(abs(berr))
!      if (info /= 0) stop 'dyn_susc_real2: failure in zgesvx'
!      kssusc0 = x
!
!      Added by Sascha:
      if(lcurrcorr) then
        call current_correlation(kssusc,zw,iw,curr_corr_lm,curr_corr_0_lm,curr_corr_div_lm,curr_corr_div_0_lm)
        if(lcurrcorrint) call susc_interpolation(kssusc,iw)
      end if
    end if
!   multipoles
    suscylm(:,:,:,:,iw) = czero; suscinv(:,:,:,:,iw) = czero
    do jq=1,ndensum
      i4 = i2almsbden(:,jq)
      jb = i4(1); jlm = i4(2); js = i4(3); ja2 = i4(4)
      do iq=1,ndensum
        i4 = i2almsbden(:,iq)
        ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4)
         if (ilm == 1 .and. jlm == 1) suscylm(is,js,ia2,ja2,iw) = suscylm(is,js,ia2,ja2,iw) + suscnorm(iq)*kssusc(iq,jq)*suscnorm(jq)
      end do
    end do
    if (lcartesian) then
      do ja2=1,nasusc2
      do ia2=1,nasusc2
        do jlm=1,1!lmmax0
        do ilm=1,1!lmmax0
          block(:,:) = czero
          do j=1,4
          do i=1,4
            do js=1,4
            do is=1,4
              block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*suscylm(is,js,ia2,ja2,iw)*pc2s(i2is(1,js),i2is(2,js),j)
            end do
            end do
          end do
          end do
          suscylm(:,:,ia2,ja2,iw) = block
        end do
        end do
      end do
      end do
    end if
!   final filtering
    maxelem = maxval(abs(suscylm(:,:,:,:,iw)))
    where (abs(suscylm(:,:,:,:,iw)) < susctol*maxelem) suscylm(:,:,:,:,iw) = 0.d0
!   --------------------------------------------------------------------
!   eigenvalues
    jq = 1
    do ja2=1,nasusc2
    do jlm=1,1!lmmax0
    do j=1,2
      iq = 1
      do ia2=1,nasusc2
      do ilm=1,1!lmmax0
      do i=1,2
!        tmp(iq,jq) = (suscylm(i,j,ilm,jlm,ia2,ja2) - conjg(suscylm(j,i,jlm,ilm,ja2,ia2)))/i2pi
        tmp(iq,jq) = suscylm(i,j,ia2,ja2,iw)
        iq = iq + 1
      end do
      end do
      end do
      jq = jq + 1
    end do
    end do
    end do
!    call zheev('V','U',nevals,tmp,nevals,evals,work,lwork,rwork,info)
    call zgeev('V','V',nevals,tmp,nevals,evals(:),evecl,nevals,evecr(:,:,iw),nevals,work,lwork,rwork,info)
    if (info /= 0) stop 'dyn_susc_expansion: failure in zgeev'
    write(iofile,'(1000es16.8)') zw(iw), evals(:), 1.d0/evals(:)
!   --------------------------------------------------------------------
!   inverse of spherical susceptibility
    jq = 1
    do ja2=1,nasusc2
    do jlm=1,1!lmmax0
    do j=1,2
      iq = 1
      do ia2=1,nasusc2
      do ilm=1,1!lmmax0
      do i=1,2
!        tmp(iq,jq) = (suscylm(i,j,ilm,jlm,ia2,ja2) - conjg(suscylm(j,i,jlm,ilm,ja2,ia2)))/i2pi
        tmp(iq,jq) = suscylm(i,j,ia2,ja2,iw)
        iq = iq + 1
      end do
      end do
      end do
      jq = jq + 1
    end do
    end do
    end do
    call zgetrf(nevals,nevals,tmp,nevals,ipiv,info)
    if (info /= 0) stop 'dyn_susc_expansion: failure in zgetrf'
    call zgetri(nevals,tmp,nevals,ipiv,work,lwork,info)
    if (info /= 0) stop 'dyn_susc_expansion: failure in zgetri'
    jq = 1
    do ja2=1,nasusc2
    do jlm=1,1!lmmax0
    do j=1,2
      iq = 1
      do ia2=1,nasusc2
      do ilm=1,1!lmmax0
      do i=1,2
        suscinv(i,j,ia2,ja2,iw) = tmp(iq,jq)
        iq = iq + 1
      end do
      end do
      end do
      jq = jq + 1
    end do
    end do
    end do
! ******
  end do
! ******
  close(iofile)
! output blocks of susceptibility
  do ja2=1,nasusc2
    do ia2=1,nasusc2
      write(filename,'("ia",i4.4,"ja",i4.4,".chi0")') ia2, ja2
      open(file=filename,unit=iofile,status='replace')
      write(filename,'("ia",i4.4,"ja",i4.4,".chif")') ia2, ja2
      open(file=filename,unit=iofile2,status='replace')
      write(filename,'("ia",i4.4,"ja",i4.4,".chid")') ia2, ja2
      open(file=filename,unit=iofile3,status='replace')
      write(iofile,'("# real(energy) ((susc0(i,j),j=1,4),i=1,4)")')
      write(iofile2,'("# real(energy) ((suscylm(i,j),j=1,4),i=1,4)")')
      write(iofile3,'("# real(energy) (suscylm(i,i),i=1,4)")')
      do iw=1,nomega
        write(iofile,'(1000es16.8)') real(zw(iw)), ((susc0ylm(i,j,ia2,ja2,iw),j=1,4),i=1,4)!, ((susc0inv(i,j,ia2,ja2,iw),j=1,2),i=1,2)
        write(iofile2,'(1000es16.8)') real(zw(iw)), ((suscylm(i,j,ia2,ja2,iw),j=1,4),i=1,4)!, ((suscinv(i,j,ia2,ja2,iw),j=1,2),i=1,2)
        write(iofile3,'(1000es16.8)') real(zw(iw)), (suscylm(i,i,ia2,ja2,iw),i=1,4)
      end do
      close(iofile)
      close(iofile2)
      close(iofile3)
    end do
  end do
! Added by Sascha:
! output blocks of spincurrent spin correlation function
  if(lcurrcorr) then
  do ja2=1,nasusc2
    do ia2=1,nasusc2
      write(filename,'("ia",i4.4,"ja",i4.4,".curf")') ia2, ja2
      open(file=filename,unit=iofile,status='replace')
      write(filename,'("ia",i4.4,"ja",i4.4,".cur0")') ia2, ja2
      open(file=filename,unit=iofile2,status='replace')
      write(filename,'("ia",i4.4,"ja",i4.4,".curdiv")') ia2, ja2
      open(file=filename,unit=iofile3,status='replace')
      write(filename,'("ia",i4.4,"ja",i4.4,".curdiv0")') ia2, ja2
      open(file=filename,unit=iofile4,status='replace')
      write(iofile,'("# omega chi_xx^x chi_xy^x chi_yx^x chi_yy^x chi_xx^y ...")')
      write(iofile2,'("# omega chi_xx^x chi_xy^x chi_yx^x chi_yy^x chi_xx^y ...")')
      write(iofile3,'("# omega div( chi_xx chi_xy chi_yx chi_yy )")')
      write(iofile4,'("# omega div( chi_xx chi_xy chi_yx chi_yy )")')
      do iw=1,nomega
        do ilm=1,1!lmmax0
          write(iofile,'(es16.8,2i5,1000es16.8)') real(zw(iw)),i2lm(:,ilm), (((curr_corr_lm(a,i,j,ilm,1,ia2,ja2,iw),j=1,2),i=1,2),a=1,3)
          write(iofile2,'(es16.8,2i5,1000es16.8)') real(zw(iw)),i2lm(:,ilm), (((curr_corr_0_lm(a,i,j,ilm,1,ia2,ja2,iw),j=1,2),i=1,2),a=1,3)
          write(iofile3,'(es16.8,2i5,1000es16.8)') real(zw(iw)),i2lm(:,ilm), ((curr_corr_div_lm(i,j,ilm,1,ia2,ja2,iw),j=1,2),i=1,2)
          write(iofile4,'(es16.8,2i5,1000es16.8)') real(zw(iw)),i2lm(:,ilm), ((curr_corr_div_0_lm(i,j,ilm,1,ia2,ja2,iw),j=1,2),i=1,2)
        end do
      end do
      close(iofile)
      close(iofile2)
      close(iofile3)
      close(iofile4)
    end do
  end do
  end if
! output eigenvectors
  do ia2=1,nasusc2
    write(filename,'("evecs",i4.4,"ns",i2.2,".dat")') ia2, 2
    open(file=filename,unit=iofile,status='replace')
    do iw=1,nomega
      write(iofile,'(1000es16.8)') zw(iw), evecr(:,2*ia2-1,iw), evecr(:,2*ia2,iw)
    end do
    close(iofile)
  end do
  deallocate(tmp,work,rwork,evals,evecl,evecr,susc0ylm,suscylm,susc0inv,suscinv)
  if(lcurrcorr) deallocate(curr_corr_0_lm,curr_corr_lm,curr_corr_div_lm,curr_corr_div_0_lm)
  call cpu_time(finish)
  write(*,'(/," Dynamic KS susc time=",f10.3," s")') finish - start
! All done!
  end subroutine dyn_susc_expansion
