!> Module that contains subroutines that are used in the
!> KKr self consistency loop (previously done in the main program)
module KKRSelfConsistency_mod
  implicit none

  contains

  !----------------------------------------------------------------------------
  !> Extracts the DOS at the Fermi level and stores result in DENEF
  !> @param IELAST energy index of Fermi level
  !> @param LMAXD1 lmax+1
  double precision function calcDOSatFermi(DEN, IELAST, IEMXD, LMAXD1, NSPIN)
    implicit none

    double complex, intent(in) :: DEN(0:LMAXD1,IEMXD,NSPIN)
    integer, intent(in) :: IELAST
    integer, intent(in) :: IEMXD
    integer, intent(in) :: LMAXD1
    integer, intent(in) :: NSPIN

    integer :: ISPIN
    integer :: L
    double precision :: PI
    double precision :: DENEF

    PI = 4.0D0*ATAN(1.0D0)

    ! get density of states at Fermi-level
    DENEF = 0.0d0
    do ISPIN = 1,NSPIN
      do L = 0,LMAXD1
        DENEF = DENEF - 2.0D0 * &
        DIMAG(DEN(L,IELAST,ISPIN))/PI/DBLE(NSPIN)
      end do
    end do

    calcDOSatFermi = DENEF

  end function calcDOSatFermi

  !----------------------------------------------------------------------------
  !> Copy output potential to input potential for new iteration.
  subroutine resetPotentials(IRC1, IRMD, IRMIN1, IRMIND, LMPOTD, NSPIN, VINS, VISP, VONS)
    implicit none

    integer :: IRC1
    integer :: IRMIN1
    integer :: ISPIN
    integer :: J
    integer :: LM

    integer :: IRMD

    integer :: IRMIND
    integer :: LMPOTD
    integer :: NSPIN

    double precision, intent(out) :: VINS(IRMIND:IRMD,LMPOTD,2)
    double precision, intent(out) :: VISP(IRMD,2)
    double precision, intent(inout) :: VONS(IRMD,LMPOTD,2)

    !DEBUG
    if (IRMIN1 > IRMIND) then
      write(*,*) "resetPotentials: IRMIN1 > IRMIND"
      stop
    end if

    !initialise VINS
    do ISPIN = 1,2
      do LM = 1, LMPOTD
        do J = IRMIND, IRMD
          VINS(J,LM,ISPIN) = 0.0D0
        enddo
      enddo
    enddo

    !initialise VISP
    do ISPIN = 1,2
      do J = 1, IRMD
        VISP(J,ISPIN) = 0.0D0
      enddo
    enddo

    ! copy output potential to input potential for new iteration
    do ISPIN = 1,NSPIN

      call DCOPY(IRC1,VONS(1,1,ISPIN),1,VISP(1,ISPIN),1)

      if (LMPOTD>1) then

        do LM = 2,LMPOTD
          do J = IRMIN1,IRC1
            VINS(J,LM,ISPIN) = VONS(J,LM,ISPIN)
          end do
        end do
      end if
    enddo
  end subroutine resetPotentials

  !----------------------------------------------------------------------------
  !> Calculates the l-resolved charges.
  !> This is done by energy integration in the complex plane over the imaginary
  !> part of the diagonal of the structural Green's function (DEN)

  subroutine calcChargesLres(CHARGE, DEN, IELAST, LMAXD1, NSPIN, WEZ, IEMXD)
    implicit none

    double precision, intent(out) :: CHARGE(0:LMAXD1,2)
    double complex, intent(in) :: DEN(0:LMAXD1,IEMXD,NSPIN)
    doublecomplex, intent(in) :: WEZ(IEMXD)

    integer, intent(in) :: IEMXD
    integer, intent(in) :: NSPIN
    integer, intent(in) :: LMAXD1
    integer, intent(in) :: IELAST

    integer :: ISPIN
    integer :: L
    integer :: IE

    ! ---> l/m_s/atom-resolved charges

    do ISPIN = 1,NSPIN
      do L = 0,LMAXD1
        CHARGE(L,ISPIN) = 0.0D0

        do IE = 1,IELAST
          CHARGE(L,ISPIN) = CHARGE(L,ISPIN) + &
          DIMAG(WEZ(IE)*DEN(L,IE,ISPIN))/ &
          DBLE(NSPIN)
        end do

      end do
    end do
  end subroutine calcChargesLres

end module KKRSelfConsistency_mod
