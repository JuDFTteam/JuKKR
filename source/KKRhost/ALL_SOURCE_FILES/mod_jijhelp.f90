module mod_jijhelp

  implicit none
  private

  public :: set_Jijcalc_flags
  
contains

  subroutine set_Jijcalc_flags(t_dtmatJij,NATYPD,NATOMIMPD,NATOMIMP,ATOMIMP,IQAT)

    use mod_types, only: t_inc, type_dtmatJijDij

    implicit none
    type(type_dtmatJijDij), intent(inout) :: t_dtmatJij(t_inc%NATYP)
    integer, intent(in) :: NATYPD,NATOMIMPD,NATOMIMP, ATOMIMP(NATOMIMPD),IQAT(NATYPD)

    integer :: I1, IQ

    DO IQ=1,t_inc%NATYP
     DO I1=1,NATOMIMP
       IF(IQAT(IQ)==ATOMIMP(I1))then
         t_dtmatJij(IQ)%calculate = .true.
       END IF
     END DO!I1
    END DO!IQ

    !Test
!   write(*,*) '=========== TEST FOR TMAT-CALC============='
!   do IQ=1,t_inc%NATYP
!     write(*,'(A,I8,L3)') "atom",IQ, t_dtmatJij(IQ)%calculate
!   end do

  end subroutine set_Jijcalc_flags






end module mod_jijhelp
