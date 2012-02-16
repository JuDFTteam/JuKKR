!> Dummy procedure to replace mminvmod for debugging purposes.
!> Gives only some garbage results.
subroutine MMINVMOD_DUMMY(GLLH1,X2,TMATLL,NUMN0,INDN0,N2B, &
                    IAT,SCITER,ITCOUNT, &
                    GLLHBLCK,BCP,IGUESS,CNVFAC, &
                    TOL, &
                    naezd, lmmaxd, naclsd, xdim, ydim, zdim, &
                    natbld, nthrds)


  implicit none

  integer, intent(in) :: naezd
  integer, intent(in) :: lmmaxd
  integer, intent(in) :: naclsd
  integer, intent(in) :: xdim
  integer, intent(in) :: ydim
  integer, intent(in) :: zdim
  integer, intent(in) :: natbld
  integer, intent(in) :: nthrds

  double complex :: GLLH1(LMMAXD,NACLSD*LMMAXD,NAEZD) ! in?
  double complex :: X2(NAEZD*LMMAXD,LMMAXD)           ! out
  double complex :: TMATLL(LMMAXD,LMMAXD,NAEZD)       ! in
  integer:: NUMN0(NAEZD)                              ! in
  integer:: INDN0(NAEZD,NACLSD)                       ! in
  double precision::N2B(LMMAXD)                       ! in, calculated outside mminvmod (not optimal)
  integer::IAT                                        ! in
  integer::SCITER                                     ! in
  integer::ITCOUNT                                    ! out
  double complex :: GLLHBLCK(NATBLD*LMMAXD, NATBLD*XDIM*YDIM*ZDIM*LMMAXD) ! inout, work-array
  integer::BCP                                        ! in
  integer::IGUESS                                     ! in
  double precision::CNVFAC                            ! inout - useless
  double precision::TOL                               ! in

  !----------------------------------------------------------------------------

  integer :: x, y

  ! fill output with some garbage

  CNVFAC = 1000.0D0
  ITCOUNT = 1

  do y = 1, lmmaxd
    do x = 1, naezd*lmmaxd
      X2(x,y) = 0.4
    end do
  end do

end subroutine
