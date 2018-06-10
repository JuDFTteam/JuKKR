!-------------------------------------------------------------------------------
!> @brief Calculates the source terms J,H and the left solution J2, H2 for:
!> - non-relativistic
!> - scalar-relativistic
!> - full-relativistic
!> calculations
!-------------------------------------------------------------------------------
subroutine rllsllsourceterms(nsra, nvec, eryd, rmesh, nrmax, nrmaxd, lmax, &
  lmsize, use_fullgmat, jlk_index, hlk, jlk, hlk2, jlk2, gmatprefactor)

  use :: constants

  implicit none

  integer :: nsra, lmax, nrmax, nrmaxd, nvec
  double complex :: eryd
  double precision, dimension (nrmaxd) :: rmesh
  integer, dimension (2*lmsize) :: jlk_index
  integer :: l1, lm1, m1, ivec, ispinfullgmat, ir
  integer :: use_fullgmat
  integer :: lmsize
  double complex :: ek, ek2, gmatprefactor
  double complex, dimension (1:4*(lmax+1), nrmax) :: hlk, jlk
  double complex, dimension (1:4*(lmax+1), nrmax) :: hlk2, jlk2

  if (nsra==2) then
    nvec = 2
  else if (nsra==1) then
    nvec = 1
  end if

  lm1 = 1
  do ivec = 1, nvec
    do ispinfullgmat = 0, use_fullgmat
      do l1 = 0, lmax
        do m1 = -l1, l1
          jlk_index(lm1) = l1 + (ivec-1)*(lmax+1) + 1
          lm1 = lm1 + 1
        end do
      end do
    end do !ispinorbit=0,use_fullgmat
  end do !nvec

  if (nsra==1) then
    ek = sqrt(eryd)
    ek2 = sqrt(eryd)
  else if (nsra==2) then
    ek = sqrt(eryd+(eryd/cvlight)**2)
    ek2 = sqrt(eryd+(eryd/cvlight)**2)*(1.0d0+eryd/cvlight**2)
  end if

  do ir = 1, nrmax

    call beshank(hlk(:,ir), jlk(:,ir), ek*rmesh(ir), lmax)
    if (nsra==2) then
      call beshank_smallcomp(hlk(:,ir), jlk(:,ir), ek*rmesh(ir), rmesh(ir), &
        eryd, lmax)
    end if

    do l1 = 1, nvec*(lmax+1)
      hlk(l1, ir) = -ci*hlk(l1, ir)
    end do

    if (nsra==1) then
      do l1 = 1, nvec*(lmax+1)
        jlk2(l1, ir) = jlk(l1, ir)
        hlk2(l1, ir) = hlk(l1, ir)
      end do
    else if (nsra==2) then
      do l1 = 1, lmax + 1
        jlk2(l1, ir) = jlk(l1, ir)
        hlk2(l1, ir) = hlk(l1, ir)
      end do
      do l1 = lmax + 2, 2*(lmax+1)
        jlk2(l1, ir) = -jlk(l1, ir)
        hlk2(l1, ir) = -hlk(l1, ir)
      end do
    end if

  end do
  gmatprefactor = ek2
end subroutine
