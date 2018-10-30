  subroutine build_kxcalda3(kxclm,kxcylm,kxcy00)
! assembles the xc ALDA kernel in the density basis

  implicit none

  real(kind=r8b),    intent(in)  :: kxclm(nrmax,lmmax4,4,4,nasusc)
  complex(kind=c8b), intent(out) :: kxcylm(4,4,lmmax0,lmmax0,nasusc)
  complex(kind=c8b), intent(out) :: kxcy00(4,4,lmmax,lmmax,nasusc)
! -----------------------------------------------------------------
  real(kind=r8b),    parameter :: fourpi = 16.d0*atan(1.d0)
  integer(kind=i4b) :: s1, s2, s3, s4, is, js, nden0, nden1
  integer(kind=i4b) :: i2(2), i4(4), ie, ia, jlm, ilm, klm, ib, jb, j, i, nr, iq, jq
  real(kind=r8b)    :: dr(nrmax), maxnorm, dummy(5)
  real(kind=r8b)    :: kxc(nrmax), maxelem
  complex(kind=c8b) :: work(nrmax), norm, norm2
  complex(kind=c8b) :: suscylm(4,4,lmmax0,lmmax0,nasusc,nasusc)
  complex(kind=c8b) :: suscy00(4,4,lmmax,lmmax,nasusc,nasusc)
  complex(kind=c8b) :: suscden(4,lmmax,nasusc)
  integer(kind=i4b) :: iqmax, jqmax, ipiv(ngfsum), info

  kxcylm = 0.d0; kxcy00 = 0.d0; kernel = 0.d0
! static susceptibility
  call dyn_susc_real2(0.d0,suscylm,suscy00,suscden,.true.,.false.,.false.)
! kernel as inverse of static susceptibility
  do iq=1,ndensum
    kernel(iq,iq) = 1.d0 - 1.d-6
  end do
  call zgesv(ndensum,ndensum,kssusc,ngfsum,ipiv,kernel,ngfsum,info)
  if (info /= 0) stop 'build_kxcalda3: failure in zgesv'
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
!        write(iodb,'("build_kxcalda3: indices=",7i8)') ia, is, js, ilm, jlm, ib, jb
!       transverse only
        if (is > 2 .or. js > 2) kernel(iq,jq) = 0.d0
        norm = kernel(iq,jq)
        kxcylm(is,js,ilm,jlm,ia) = kxcylm(is,js,ilm,jlm,ia) + norm*suscnorm(iq)*suscnorm(jq)
        if (abs(norm) > maxnorm) then
          maxnorm = abs(norm)
          iqmax = iq; jqmax = jq
        end if
        if (abs(norm) > susctol) write(iodb,'("ia,is,js,ilm,jlm,ib,jb,norm=",7i8,2es16.8)') ia, is, js, ilm, jlm, ib, jb, norm
      end do
    end do
    write(iodb,'("build_kxcalda3: norm=",7i8,2es16.8)') ia, i2almsbden(1:3,iqmax), i2almsbden(1:3,jqmax), maxnorm
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do            ! atom i
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  kxcsusc = real(kxcsusc)
! All done
  end subroutine build_kxcalda3

