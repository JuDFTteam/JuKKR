!-------------------------------------------------------------------------------
! SUBROUTINE: DECIMATE
!> @brief Decimation method
!> - Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine decimate(gllke, naez, tinvbup, tinvbdown, vacflag, factl, nlbasis, &
  nrbasis, alm, ndim, lmmaxd)

  use :: global_variables

  implicit none

!     .. Parameters ..
!
! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
! *          function, set up in the spin-independent non-relativstic *
! *          (l,m_l)-representation                                   *
! *                                                                   *
! *********************************************************************
!
  integer, intent (in) :: alm !< NAEZ*LMMAXD
  integer, intent (in) :: ndim !< NPRINCD*LMMAXD
  integer, intent (in) :: naez !< Number of atoms in unit cell
  integer, intent (in) :: lmmaxd !< (KREL+KORBIT+1)(LMAX+1)^2
  integer, intent (in) :: nlbasis !< Number of basis layers of left host (repeated units)
  integer, intent (in) :: nrbasis !< Number of basis layers of right host (repeated units)
  logical, dimension (2), intent (in) :: vacflag
  double complex, dimension (lmmaxd, lmmaxd), intent (in) :: factl
  double complex, dimension (lmmaxd, lmmaxd, *), intent (in) :: tinvbup
  double complex, dimension (lmmaxd, lmmaxd, *), intent (in) :: tinvbdown
  double complex, dimension (alm, alm), intent (inout) :: gllke
! .. Local Scalars
  integer :: ldi1, ldi1t, ldi2, ldi2t, lm1, lm2, nlayer, icall
  integer :: ichck, ihost, ii1, ii2, il1, il2, ip1, ip1t, ip2, ip2t, itermax
  double precision :: errmax
! .. Local Arrays
  double complex, dimension (ndim, ndim) :: a1, an, b1, bn, c1, cn, x1, xn
! ..
  data icall/0/
! .. External Subroutines
  external :: bofm, surfgf
! .. Save statement
  save :: icall, nlayer, itermax, errmax, ichck
! .. External Functions
  logical :: opt
  external :: opt
! .. Intrinsic Functions
  intrinsic :: mod
!
  icall = icall + 1
!----------------------------------------------------------------------------
  if (icall==1) then
    nlayer = naez/nprincd
! Parameters for the "decimation" technique.
    itermax = 300
    errmax = 1.0d-180
    ichck = 1
  end if
!----------------------------------------------------------------------------
  if (.not. vacflag(1)) then
!-------------------------------------------------------------------------
! Get the matrix B1
!-------------------------------------------------------------------------
    call bofm(1, 1, b1, ndim, gllke, alm)

! Now Subtract t-mat of left host
    do ip1 = 1, nprincd
      ihost = mod(ip1-1, nlbasis) + 1
      do lm1 = 1, lmmaxd
        do lm2 = 1, lmmaxd
          il1 = lmmaxd*(ip1-1) + lm1
          il2 = lmmaxd*(ip1-1) + lm2
          b1(il1, il2) = (b1(il1,il2)-tinvbup(lm1,lm2,ihost))
        end do
      end do
    end do

    call bofm(1, 2, c1, ndim, gllke, alm)
    call bofm(2, 1, a1, ndim, gllke, alm)

! It performs the 'space decimation' iterative procedure.
    call surfgf(ndim, a1, b1, c1, x1, itermax, errmax, ichck, lmmaxd)
! Adds to the matrix GLLKE the elements that couples the
! interface to the two half-spaces.
    do ip1 = 1, nprincd
      do ip2 = 1, nprincd
        ii1 = ip1
        ii2 = ip2
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            ldi1 = lmmaxd*(ip1-1) + lm1
            il1 = lmmaxd*(ii1-1) + lm1
            ldi2 = lmmaxd*(ip2-1) + lm2
            il2 = lmmaxd*(ii2-1) + lm2
            gllke(il1, il2) = gllke(il1, il2) - x1(ldi1, ldi2)
          end do
        end do
      end do
    end do
  end if

  if (.not. vacflag(2)) then

!  If 'ONEBULK' is activated then it calculates the xn decimated element
!  from the x1 element: this is just in the case of equal bulks on the

    if (.not. opt('ONEBULK ')) then

!----------------------------------------------------------------------
! Get the matrix BN
!----------------------------------------------------------------------
      call bofm(nlayer, nlayer, bn, ndim, gllke, alm)

! Now Substract t-mat right host
! Notes : the indexing is easier like that
      do ip1 = 1, nprincd
        ihost = nrbasis - mod(ip1, nrbasis)
        ihost = mod(ip1-1, nrbasis) + 1
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            il1 = lmmaxd*(ip1-1) + lm1
            il2 = lmmaxd*(ip1-1) + lm2
            bn(il1, il2) = (bn(il1,il2)-tinvbdown(lm1,lm2,ihost))
          end do
        end do
      end do

      call bofm(nlayer, nlayer-1, an, ndim, gllke, alm)
      call bofm(nlayer-1, nlayer, cn, ndim, gllke, alm)

! It performs the 'space decimation' iterative procedure.
      call surfgf(ndim, cn, bn, an, xn, itermax, errmax, ichck, lmmaxd)
!
    else
!
      do ip1 = 1, nprincd
        do ip2 = 1, nprincd
          ip1t = (nprincd+1) - ip2
          ip2t = (nprincd+1) - ip1
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              ldi1 = lmmaxd*(ip1-1) + lm1
              ldi2 = lmmaxd*(ip2-1) + lm2
              ldi1t = lmmaxd*(ip1t-1) + lm2
              ldi2t = lmmaxd*(ip2t-1) + lm1
              xn(ldi1t, ldi2t) = factl(lm1, lm2)*x1(ldi1, ldi2)
            end do
          end do
        end do
      end do
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             Added on 1.02.2000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Adds to the matrix GLLKE the elements that couples the
! interface to the two half-spaces.
    do ip1 = 1, nprincd
      do ip2 = 1, nprincd
        ii1 = (nlayer-1)*nprincd + ip1
        ii2 = (nlayer-1)*nprincd + ip2
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            ldi1 = lmmaxd*(ip1-1) + lm1
            il1 = lmmaxd*(ii1-1) + lm1
            ldi2 = lmmaxd*(ip2-1) + lm2
            il2 = lmmaxd*(ii2-1) + lm2
            gllke(il1, il2) = gllke(il1, il2) - xn(ldi1, ldi2)
          end do
        end do
      end do
    end do
  end if

  return

end subroutine
