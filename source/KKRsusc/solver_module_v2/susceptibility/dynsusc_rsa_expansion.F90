  subroutine dynsusc_rsa_expansion
! Dynamic KS susceptibility
! Rigid spin approximation
! Transverse susceptibility only
  use global

  implicit none

! -----------------------------------------------------------------
!   i 2pi
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0)
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi)
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0), cminus = (-1.d0,0.d0)
  complex(kind=c8b), parameter :: e0up = (0.481995d0,-0.05d0), e0dn = (0.711995,-0.007d0)
! -----------------------------------------------------------------
  complex(kind=c8b), allocatable :: susc0ylm(:,:,:,:,:), suscylm(:,:,:,:,:)
  complex(kind=c8b) :: zw(nomega)
  integer(kind=i4b) :: ib, jb, is, js, i4(4), ia, ja, ia2, ja2, jlm, ilm, j, i, iw, iq, jq
  integer(kind=i4b) :: ipiv(ndensum), info, nevals, lwork
  real(kind=r8b)    :: maxelem, start, finish
  complex(kind=c8b) :: block(4,4)
  character*17      :: filename
  complex(kind=c8b), allocatable :: work(:), tmp(:,:), evals(:), evecl(:,:), evecr(:,:,:)
  real(kind=r8b),    allocatable :: rwork(:)!, evals(:)
!  real(kind=r8b)    :: r(nalmsb), c(nalmsb), rwork(2*nalmsb), rcond, ferr(nalmsb), berr(nalmsb)
!  complex(kind=c8b) :: work(2*nalmsb), temp(nalmsb,nalmsb)
!  complex(kind=c8b) :: x(nalmsb,nalmsb)
!  character*1       :: equed

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
  allocate(susc0ylm(4,4,nasusc2,nasusc2,nomega))
  allocate(suscylm(4,4,nasusc2,nasusc2,nomega))
! ********************
! Frequency loop
  do iw=1,nomega
! ********************
    write(*,'("iw=",i4,", zw=",2es16.8)') iw, zw(iw)
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!             Kohn-Sham susceptibility for current frequency
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   Put together KS susc from Taylor expansion
    call zcopy(ndensum*ndensum,kssusc0,1,kssusc,1)
    call zaxpy(ndensum*ndensum,zw(iw),kssusc1,1,kssusc,1)
    call zaxpy(ndensum*ndensum,0.5d0*zw(iw)**2,kssusc2,1,kssusc,1)
!   multipoles
    susc0ylm(:,:,:,:,iw) = czero
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
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               Full susceptibility for current frequency
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   Solve Dyson equation
    if (lenhanced) then
      call zgemm('N','N',ndensum,ndensum,ndensum,cminus,kssusc,ndensum,kernel,ndensum,czero,denominator,ndensum)
      do iq=1,ndensum
        denominator(iq,iq) = denominator(iq,iq) + cone
      end do
      call zgesv(ndensum,ndensum,denominator,ndensum,ipiv,kssusc,ndensum,info)
      if (info /= 0) stop 'dyn_susc_expansion: failure in zgesv'
!      call zgesvx('N','N',nalmsb,nalmsb,denominator,nalmsb,temp,nalmsb,ipiv,equed,r,c,kssusc0,nalmsb,x,nalmsb,rcond,ferr,berr,work,rwork,info)
!      write(iodb,'("condition number of enhancement factor=",es16.3)') rcond
!      write(iodb,'("fwd & bkwd error in solution=",2es16.3)') maxval(abs(ferr)), maxval(abs(berr))
!      if (info /= 0) stop 'dyn_susc_real2: failure in zgesvx'
!      kssusc0 = x
    end if
!   multipoles
    suscylm(:,:,:,:,iw) = czero
    do jq=1,ndensum
      i4 = i2almsbden(:,jq)
      jb = i4(1); jlm = i4(2); js = i4(3); ja2 = i4(4)
      do iq=1,ndensum
        i4 = i2almsbden(:,iq)
        ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4)
        suscylm(is,js,ia2,ja2,iw) = suscylm(is,js,ia2,ja2,iw) + suscnorm(iq)*kssusc(iq,jq)*suscnorm(jq)
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
    if (info /= 0) stop 'dyn_susc_expansion: failure in zheev'
    write(iofile,'(1000es16.8)') zw(iw), evals(:), 1.d0/evals(:)
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
      do iw=1,nomega
        write(iofile,'(1000es16.8)') real(zw(iw)), ((susc0ylm(i,j,ia2,ja2,iw),j=1,4),i=1,4)
        write(iofile2,'(1000es16.8)') real(zw(iw)), ((suscylm(i,j,ia2,ja2,iw),j=1,4),i=1,4)
      end do
      close(iofile)
      close(iofile2)
    end do
  end do
! output eigenvectors
  do ia2=1,nasusc2
    write(filename,'("evecs",i4.4,"ns",i2.2,".dat")') ia2, 2
    open(file=filename,unit=iofile,status='replace')
    do iw=1,nomega
      write(iofile,'(1000es16.8)') zw(iw), evecr(:,2*ia2-1,iw), evecr(:,2*ia2,iw)
    end do
    close(iofile)
  end do
  deallocate(tmp,work,rwork,evals,evecl,evecr,susc0ylm,suscylm)
  call cpu_time(finish)
  write(*,'(/," Dynamic KS susc time=",f10.3," s")') finish - start
! All done!
  end subroutine dynsusc_rsa_expansion
