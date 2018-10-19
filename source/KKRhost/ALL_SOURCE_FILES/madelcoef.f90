!------------------------------------------------------------------------------------
!> Summary: Calculation of the Madelung coefficients 
!> Author: 
!> Calculation of the structure dependent Madelung coefficients for the determination 
!> of the electrostatic potential. 
!------------------------------------------------------------------------------------
module mod_madelcoef
  use :: mod_datatypes, only: dp
  use :: constants, only : pi
  private :: dp

contains
  !-------------------------------------------------------------------------------
  !> Summary: Calculation of the Madelung coefficients 
  !> Author: 
  !> Category: electrostatics, KKRhost 
  !> Deprecated: False 
  !> Calculation of the structure dependent Madelung coefficients for the determination 
  !> of the electrostatic potential.
  !-------------------------------------------------------------------------------
  subroutine madelcoef(linterface,lpot,a,b,smat,cleb,icleb,iend,lpotd,lmpotd,lmxspd,&
    nclebd)

    implicit none
    ! ..
    ! .. Input variables
    integer, intent(in) :: iend   !! Number of nonzero gaunt coefficients
    integer, intent(in) :: lpot   !! Maximum l component in potential expansion
    integer, intent(in) :: lpotd  !! Maximum l component in potential expansion
    integer, intent(in) :: lmpotd !! (lpot+1)**2
    integer, intent(in) :: lmxspd !! (2*lpot+1)**2
    integer, intent(in) :: nclebd !! lmxspd*lmpotd
    logical, intent(in) :: linterface !! If True a matching with semi-inifinite surfaces must be performed
    integer, dimension(nclebd,3), intent(in) :: icleb !! Pointer array
    real (kind=dp), dimension(lmxspd), intent(in) :: smat !! Lattice summations coming from the Ewald method
    real (kind=dp), dimension(nclebd), intent(in) :: cleb !! !! GAUNT coefficients (GAUNT)
    ! ..
    ! .. Output variables
    real (kind=dp), dimension(lmpotd), intent(out) :: b !! Madelung coefficients
    real (kind=dp), dimension(lmpotd, lmpotd), intent(out) :: a !! Madelung coefficients

    ! ..
    ! .. Local scalars
    real (kind=dp), parameter :: fpi = 4.0_dp*pi
    integer :: i, l, l1, l2, lm1, lm2, lm3, lmpot, m
    integer, dimension(lmxspd) :: loflm
    ! INTEGER ICALL_madelcoef
    ! ..
    ! .. Local arrays
    real (kind=dp), dimension(0:lpotd,0:lpotd) :: dfac
    ! ..
    ! .. Data statements
    ! DATA ICALL_madelcoef /0/
    ! integer, save :: icall_madelcoef=0
    ! ..
    ! .. Intrinsic functions
    intrinsic :: real
    ! ..................................................................

    lmpot = (lpot+1)**2

    i = 1

    ! --> determine the l-value for given lm

    do l = 0, 2*lpot
      do m = -l, l
        loflm(i) = l
        i = i + 1
      end do
    end do

    ! --> calculate:                             (2*(l+l')-1)!!
    ! dfac(l,l') = 4pi**2 *  ----------------------
    ! (2*l+1)!! * (2*l'+1)!!

    dfac(0, 0) = fpi*fpi
    do l1 = 1, lpot
      dfac(l1, 0) = dfac(l1-1, 0)*real(2*l1-1, kind=dp)/real(2*l1+1, kind=dp)
      dfac(0, l1) = dfac(l1, 0)
      do l2 = 1, l1
        dfac(l1, l2) = dfac(l1, l2-1)*real(2*(l1+l2)-1, kind=dp)/real(2*l2+1, kind=dp)
        dfac(l2, l1) = dfac(l1, l2)
      end do
    end do

    ! --> initialize

    do lm1 = 1, lmpot
      do lm2 = 1, lmpot
        ! write(*,*) 'test',LM1,LM2,LMPOT,LMPOTD
        a(lm1, lm2) = 0.0e0_dp
      end do
    end do

    ! --> calculate a(lm1,lm2)

    do i = 1, iend
      lm1 = icleb(i, 1)
      lm2 = icleb(i, 2)
      lm3 = icleb(i, 3)
      l1 = loflm(lm1)
      l2 = loflm(lm2)

      ! --> this loop has to be calculated only for l1+l2=l3

      a(lm1, lm2) = a(lm1, lm2) + 2.0e0_dp*dfac(l1, l2)*smat(lm3)*cleb(i)
    end do

    if (linterface) return

    ! --> initialize

    do lm1 = 1, lmpot
      b(lm1) = 0.0e0_dp
    end do

    ! --> calculate b(lm1)

    do lm1 = 1, lmpot
      l1 = loflm(lm1)
      b(lm1) = b(lm1) - 2.0e0_dp*fpi/real(2*l1+1, kind=dp)*smat(lm1)
    end do
  end subroutine madelcoef

end module mod_madelcoef
