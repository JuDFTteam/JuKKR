  subroutine build_vldaub(ie,ia,vldaub)
! assembles the Dudarev LDA+U potential in the projection basis
  use global, only: i4b, r8b, c8b, nlms, nlmsb, nlmsba, i2lmsb, i2lm, lms2i, ueff, jeff, overlap, rhomat, vshift

  implicit none

! Which energy
  integer(kind=i4b), intent(in)    :: ie
! Which atom
  integer(kind=i4b), intent(in)    :: ia
! LDA+U potential (added to)
  complex(kind=c8b), intent(inout) :: vldaub(nlmsb,nlmsb)
! -----------------------------------------------------------------
  integer(kind=i4b) :: is, js, ib, jb, nbi, nbj, i3(3), i2(2)
  integer(kind=i4b) :: jm, jl, im, il, jlm, ilm, i, j
  real(kind=r8b)    :: uminusj
  complex(kind=c8b) :: norm
  integer(kind=i4b) :: imax, jmax

!  maxnorm = 0.d0
  do j=1,nlmsba(ia)
    i3 = i2lmsb(:,j,ia)
    jb = i3(1); jlm = i3(2); js = i3(3)
    i2 = i2lm(:,jlm)
    jm = i2(1); jl = i2(2)
    do i=1,nlmsba(ia)
      i3 = i2lmsb(:,i,ia)
      ib = i3(1); ilm = i3(2); is = i3(3)
      i2 = i2lm(:,ilm)
      im = i2(1); il = i2(2)
!   ----------------------------------------------------------------
!   selection rules
      if (jl == il) then
        uminusj = ueff(il,ia) - jeff(il,ia)
!     revise the basis
        norm = -rhomat(lms2i(ilm,is),lms2i(jlm,js),ia)
        if (ilm == jlm .and. is == js) norm = norm + 0.5d0 !+ vshift(il,is,ia)
!        if (abs(norm) > maxnorm) then
!          maxnorm = abs(norm)
!          imax = i; jmax = j
!        end if
        vldaub(i,j) = vldaub(i,j) + uminusj*overlap(i,j,ia)*norm
!        vldaub(i,j) = vldaub(i,j) + (ueff - jeff)*norm
!        if (is == 1) js = 2
!        if (is == 2) js = 1
!        norm = rhomat(lms2i(jlm,js),lms2i(ilm,js))
!        vldaub(i,j) = vldaub(i,j) + jeff*overlap(i,j,ia)*norm
      end if
!   ----------------------------------------------------------------
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
  end subroutine build_vldaub

