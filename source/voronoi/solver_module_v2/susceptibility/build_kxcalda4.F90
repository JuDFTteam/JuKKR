  subroutine build_kxcalda4(rhoden,kxclm,kxcylm,kxcy00)
! assembles the xc ALDA kernel in the density basis from the sum rule

  implicit none

  complex(kind=c8b), intent(in)  :: rhoden(ndenmax,nasusc)
  real(kind=r8b),    intent(in)  :: kxclm(nrmax,lmmax4,4,4,nasusc)
  complex(kind=c8b), intent(out) :: kxcylm(4,4,lmmax0,lmmax0,nasusc)
  complex(kind=c8b), intent(out) :: kxcy00(4,4,lmmax,lmmax,nasusc)
! -----------------------------------------------------------------
  real(kind=r8b),    parameter :: fourpi = 16.d0*atan(1.d0)
  real(kind=r8b),    parameter :: tol = 1.d-8
  complex(kind=c8b), parameter :: cone = (1.d0,0.d0), czero = (0.d0,0.d0)
  integer(kind=i4b) :: s1, s2, s3, s4, is, js, nden0, nden1
  integer(kind=i4b) :: i2(2), i4(4), ie, ia, ja, jlm, ilm, klm, ib, jb, j, i, nr, iq, jq
  real(kind=r8b)    :: dr(nrmax), maxnorm, dummy(5)
  real(kind=r8b)    :: kxc(nrmax), maxelem, re, im
  complex(kind=c8b) :: suscylm(4,4,lmmax0,lmmax0,nasusc,nasusc)
  complex(kind=c8b) :: suscy00(4,4,lmmax,lmmax,nasusc,nasusc)
  complex(kind=c8b) :: suscden(4,lmmax,nasusc)
  complex(kind=c8b) :: work(nrmax), norm, norm2, bxceff(ngfsum), magbasis(ngfsum)
  integer(kind=i4b) :: iqmax, jqmax, ipiv(ngfsum), info
  logical           :: exists


  stop 'fix bug with ia/nasusc being used as ia2/nasusc2!'

  kxcylm = 0.d0; kxcy00 = 0.d0; kernel = 0.d0
! if the kernel file exists just read in values
  inquire(file='kernel.dat',exist=exists)
  if (exists) then
    open(file='kernel.dat',unit=iofile,status='old')
    do jq=1,ndensum
      do iq=1,ndensum
        read(iofile,*) re, im
        kernel(iq,jq) = cmplx(re,im)
      end do
    end do
    close(iofile)
    return
  end if
! copy components of magnetization
  do ia=1,nasusc
    nden0 = sum(nalmsbden(1:ia-1))
    nden1 = sum(nalmsbden(1:ia))
    do iq=1,nden1-nden0
      magbasis(iq+nden0) = rhoden(iq,ia)
      write(iodb,'("magbasis=",4i4,2es16.8)') i2almsbden(:,iq), magbasis(iq+nden0)
    end do
  end do
! static susceptibility
  call dyn_susc_real2(0.d0,suscylm,suscy00,suscden,.true.,.false.,.false.)
! inverse of static susceptibility
  denominator = 0.d0
  do iq=1,ndensum
    denominator(iq,iq) = 1.d0 - tol
  end do
  call zgesv(ndensum,ndensum,kssusc,ngfsum,ipiv,denominator,ngfsum,info)
!  call zgesv(ndensum,1,kssusc,ngfsum,ipiv,magbasis,ngfsum,info)
  if (info /= 0) stop 'build_kxcalda4: failure in zgesv'
! multiply magnetization by inverse susceptibility and put it in bxc
  call zgemv('N',ndensum,ndensum,cone,denominator,ngfsum,magbasis,1,czero,bxceff,1)
! divide bxc by the magnetization
  do ia=1,nasusc    ! atom i
    nden0 = sum(nalmsbden(1:ia-1))
    nden1 = sum(nalmsbden(1:ia))
    do jq=1+nden0,nden1
      i4 = i2almsbden(:,jq)
      jb = i4(1); jlm = i4(2); js = i4(3)!; ia = i4(4)
      do iq=1+nden0,nden1
        i4 = i2almsbden(:,iq)
        ib = i4(1); ilm = i4(2); is = i4(3)!; ia = i4(4)
!       ***********************************
        if (ilm == jlm .and. ib == jb .and. is /= js) then
!       ***********************************
        write(iodb,'("bxceff=",4i4,2es16.8)') i2almsbden(:,iq), bxceff(iq)
        if (abs(magbasis(jq)) > tol) then
          bxceff(iq) = bxceff(iq)/magbasis(jq)
        else
          bxceff(iq) = 0.d0
        end if
        if (is > 2) bxceff(iq) = 0.d0
        write(iodb,'("kernel=",4i4,2es16.8)') i2almsbden(:,iq), bxceff(iq)
        if (abs(bxceff(iq)) > tol) kernel(iq,jq) = bxceff(iq)
!       ******
        end if
!       ******
      end do
    end do
  end do
! output kernel
  open(file='kernel.dat',unit=iofile,status='new')
  do jq=1,ndensum
    do iq=1,ndensum
      write(iofile,'(2es16.8)') kernel(iq,jq)
    end do
  end do
  close(iofile)
! use the ASA
! use the fact that the basis was constructed per l channel
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do ia=1,nasusc    ! atom i
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    nr = nrpts(ia)
    dr(1:nr) = drmesh(1:nr,ia)
    maxnorm = 0.d0
    nden0 = sum(nalmsbden(1:ia-1))
    nden1 = sum(nalmsbden(1:ia))
    do jq=1+nden0,nden1
      i4 = i2almsbden(:,jq)
      jb = i4(1); jlm = i4(2); js = i4(3)!; ia = i4(4)
      i2 = i2is(:,js)
      s3 = i2(1); s4 = i2(2)
      do iq=1+nden0,nden1
        i4 = i2almsbden(:,iq)
        ib = i4(1); ilm = i4(2); is = i4(3)!; ia = i4(4)
        i2 = i2is(:,is)
        s1 = i2(1); s2 = i2(2)
!       ***********************************
        if (ilm == jlm .and. ib == jb) then
!       ***********************************
!        write(iodb,'("build_kxcalda2: indices=",7i8)') ia, is, js, ilm, jlm, ib, jb
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        kxc = 0.d0
!        do klm=1,lmmax4
!        if (abs(gaunt(ilm,jlm,klm)) > ylmtol) then
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         convert from cartesian to spin
!          do j=1,4
!            do i=1,4
!              kxc(1:nr) = kxc(1:nr) + pc2s(s1,s2,i)*kxclm(1:nr,klm,i,j,ia)*ds2c(j,s3,s4)*gaunt(ilm,jlm,klm)
!            end do
!          end do
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        end if
!        end do
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       matrix elements
!        work(1:nr) = suscbasis(1:nr,ib,ilm,is,ia)*suscbasis(1:nr,jb,jlm,js,ia)*kxc(1:nr)
        norm = kernel(iq,jq)  !radint(nr,work,dr,npanat(ia),ircutat(:,ia))
        if (abs(norm) < tol) norm = 0.d0
        kxcylm(is,js,ilm,jlm,ia) = kxcylm(is,js,ilm,jlm,ia) + norm*suscnorm(iq)*suscnorm(jq)
!        kernel(iq,jq) = kernel(iq,jq) + norm
        if (abs(norm) > maxnorm) then
          maxnorm = abs(norm)
          iqmax = iq; jqmax = jq
        end if
        if (abs(norm) > tol) write(iodb,'("build_kxcalda4: ia,is,js,ilm,jlm,ib,jb,norm=",7i8,2es16.8)') ia, is, js, ilm, jlm, ib, jb, norm
!       ******
        end if
!       ******
      end do
    end do
    write(iodb,'("build_kxcalda4: norm=",7i8,2es16.8)') ia, i2almsbden(1:3,iqmax), i2almsbden(1:3,jqmax), maxnorm
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do            ! atom i
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  kxcsusc = real(kxcsusc)
! All done
  end subroutine build_kxcalda4

