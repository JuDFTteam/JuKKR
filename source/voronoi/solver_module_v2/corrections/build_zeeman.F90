  subroutine build_zeeman(ia,magdir1,deltapot,ipick)
! assembles the Zeeman coupling in the projection basis
  use global !, only: i4b, r8b, c8b, nrmax, nrpts, drmesh, phiref, nlmsba, nlmsb, nlms, i2lmsb, i2lms, i2lm, ldots, atol, pauli, lorb

  implicit none

! Which atom
  integer(kind=i4b), intent(in)    :: ia
! Spin quantization axis
  real(kind=r8b),    intent(in)    :: magdir1(3)
! Zeeman coupling in the basis (added to)
  complex(kind=c8b), intent(inout) :: deltapot(nlmsb,nlmsb)
! Transverse, longitudinal or full Zeeman operator
  integer(kind=i4b), intent(in)    :: ipick
! -----------------------------------------------------------------
  integer(kind=i4b) :: is, js, ib, jb, i3(3), i2(2)
  integer(kind=i4b) :: jm, jl, im, il, jlm, ilm, i, j
  real(kind=r8b)    :: bfield(3), newlen

! use the structure of L+S
! also use the fact that the basis was constructed per l channel
! ----------------------------------------------------------------------
! Include a possible constraining field
  bfield(:) = blen(ia)*bdir(:,ia) + bconlen(ia)*bcondir(:,ia)
! The decomposition by ipick is used mostly for the spin sum rule
! All components
!  if (ipick == 1) do nothing, use the full vector
! Only z component
  if (ipick == 2) bfield(:) = dot_product(bfield,magdir1)*magdir1(:)
! Only xy components
  if (ipick == 3) bfield(:) = bfield(:) - dot_product(bfield,magdir1)*magdir1(:)
! ----------------------------------------------------------------------
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
!   selection rules: basis functions for each il are orthonormal
      if (il == jl .and. ib == jb) then
!       spin Zeeman
        if (im == jm .and. (ibfield(ia) == 1 .or. ibfield(ia) == 3 .or. ibfield(ia) == 4)) then
          deltapot(i,j) = deltapot(i,j) - sum(bfield*pauli(is,js,1:3))
        end if
!       orbital Zeeman
        if (is == js .and. (ibfield(ia) == 2 .or. ibfield(ia) == 3)) then
          deltapot(i,j) = deltapot(i,j) - sum(bfield*lorb(ilm,jlm,1:3))
        end if
      end if
!   ----------------------------------------------------------------
    end do
  end do
! All done
  end subroutine build_zeeman
