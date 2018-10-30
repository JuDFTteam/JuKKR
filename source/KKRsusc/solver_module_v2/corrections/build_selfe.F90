  subroutine build_selfe(ie,ia,senumd,sedend,efshift,se_infty,se_fit,vselfeb)
! Puts together the model self-energy

  use global

  implicit none

! Which energy
  integer(kind=i4b), intent(in)    :: ie
! Which atom
  integer(kind=i4b), intent(in)    :: ia
! Degree of numerator and denominator
  integer(kind=i4b), intent(in)    :: senumd, sedend
! Fermi energy shift
  real(kind=r8b),    intent(in)    :: efshift
! Constant part of SE
  real(kind=r8b),    intent(in)    :: se_infty(2,nlms)
! Fit parameters
  real(kind=r8b),    intent(in)    :: se_fit(2,nlms,1:1+senumd+sedend)
! Self-energy potential (added to)
  complex(kind=c8b), intent(inout) :: vselfeb(nlmsb,nlmsb)
! -----------------------------------------------------------------
  integer(kind=i4b) :: is, js, ib, jb, i3(3), i2(2)
  integer(kind=i4b) :: jlms, ilms, jlm, ilm, i, j, inum, iden
  complex(kind=c8b) :: norm, fac, ze, zen, znum, zden

! Frequency starts at Fermi energy
!  ze = esusc(ie) - real(escf(nescf)) - efshift
!  ze = esusc(ie) - real(escf(nescf))
  ze = real(esusc(ie)-escf(nescf))
  do j=1,nlmsba(ia)
    i3 = i2lmsb(:,j,ia)
    jb = i3(1); jlm = i3(2); js = i3(3)
    jlms = lms2i(jlm,js)
    do i=1,nlmsba(ia)
      i3 = i2lmsb(:,i,ia)
      ib = i3(1); ilm = i3(2); is = i3(3)
      ilms = lms2i(ilm,is)
!     ------------------------------------------------------------------
!     assuming diagonal self-energy
      if (i == j) then
!     revise the basis
        fac = cmplx(se_fit(1,ilms,1),se_fit(2,ilms,1))
        znum = fac
        zen = ze
        do inum=1,senumd
          fac = cmplx(se_fit(1,ilms,1+inum),se_fit(2,ilms,1+inum))
          znum = znum + fac*zen
          zen = zen*ze
        end do
        zden = 1.d0
        zen = ze
        do iden=1,sedend
          fac = cmplx(se_fit(1,ilms,1+senumd+iden),se_fit(2,ilms,1+senumd+iden))
          zden = zden + fac*zen
          zen = zen*ze
        end do
        fac = cmplx(se_infty(1,ilms),se_infty(2,ilms))
        norm = fac + znum/zden
!        vselfeb(i,j) = vselfeb(i,j) + overlap(i,j,ia)*norm
        vselfeb(i,j) = vselfeb(i,j) + overlap(i,j,ia)*(norm - efshift)
      end if
!     ------------------------------------------------------------------
    end do
  end do
!  write(iodb,'(" LDA+U kernel, ia=",i4)') ia
!  write(iodb,'("vldaub norm:",6i4,2es16.8)') i2lmsb(:,imax,ia), i2lmsb(:,jmax,ia), maxnorm
!  write(iodb,'("  ib ilm  is, jb jlm  js, matrix element")')
!  do j=1,nlmsba(ia)
!    do i=1,nlmsba(ia)
!      if (abs(vldaub(i,j)) > 1.d-8) then
!        write(iodb,'(6i4,2es10.2)') i2lmsb(:,i,ia), i2lmsb(:,j,ia), vldaub(i,j)
!      end if
!    end do
!  end do
! All done
  end subroutine build_selfe
