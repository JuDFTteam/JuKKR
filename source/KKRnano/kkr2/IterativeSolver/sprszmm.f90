subroutine SPRSZMM(IAT,GLLH,NUMN0,INDN0,X,DONE,OMEGA,DELTA, &  ! <
                   AX, &                                       ! >
                   ! new input parameters after inc.p removal
                   naez, lmmaxd, naclsd, nthrds)

  ! This routine is called very often
  ! TODO: Optimise this routine if possible

  implicit none

  integer, intent(in) :: naez
  integer, intent(in) :: lmmaxd
  integer, intent(in) :: naclsd
  integer, intent(in) :: nthrds

  !     INTEGER           LMMAXD
  !     INTEGER           NDIM,NAEZ
  !     PARAMETER        (NAEZ=NAEZD,NDIM=NAEZD*LMMAXD)
  !     INTEGER           NGTBD
  !     PARAMETER        (NGTBD = NACLSD*LMMAXD)

  double complex:: CONE
  double complex:: CZERO
  parameter        (CONE=(1.0D0,0.0D0), CZERO=(0.0D0,0.0D0))
  !     ..
  !     ... Scalars ..
  integer, intent(in) ::       IAT
  double complex, intent(in) ::OMEGA  ! scalar in Matrix-Matrix-Mult.
  double complex, intent(in) ::DELTA  ! scalar in Matrix-Matrix-Mult.

  !     ... Arrays ..
!  double complex::   X(NDIM,LMMAXD)
!  double complex::  AX(NDIM,LMMAXD)
!  double complex::GLLH(LMMAXD,NGTBD,NAEZD)

  double complex::   X(NAEZ*LMMAXD,LMMAXD)
  double complex::  AX(NAEZ*LMMAXD,LMMAXD)
  double complex::  GLLH(LMMAXD, NACLSD*LMMAXD, NAEZ)
  integer::           NUMN0(NAEZ)
  integer::           INDN0(NAEZ,NACLSD)
  logical::           DONE(LMMAXD)

  !     .. Local Arrays .. Fortran 90 automatic array
  !double complex :: SPRSX(NGTBD,LMMAXD)
  double complex :: SPRSX(NACLSD*LMMAXD, LMMAXD) ! medium size
  ! ..
  !     .. Local Scalars ..
  integer::I1
  integer::I2
  integer::I3
  integer::LM2
  integer::IL1B
  integer::I2H
  integer::I3H

  integer::MYTHRD
!$  integer::       OMP_GET_THREAD_NUM

  integer::NGTBD
  integer::NDIM

  NGTBD = NACLSD*LMMAXD
  NDIM  = NAEZ*LMMAXD

     
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!$ call OMP_SET_NUM_THREADS(NTHRDS)
!$omp parallel private (I1,LM2,I2,I3H,I2H,IL1B,SPRSX,MYTHRD)
   MYTHRD = 0
!$ MYTHRD = OMP_GET_THREAD_NUM()
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  do I1=1,NAEZ

    if (MOD(I1,NTHRDS) == MYTHRD) then

      do LM2=1,LMMAXD
        if ( .not. DONE(LM2)) then
          do I2=1,NUMN0(IAT)
            I3=INDN0(I1,I2)
            I3H = (I3-1)*LMMAXD + 1
            I2H = (I2-1)*LMMAXD + 1

            call ZCOPY(LMMAXD,X(I3H,LM2),1, &
                       SPRSX(I2H,LM2),1)

          enddo
        endif
      enddo

      IL1B=LMMAXD*(I1-1)

      call ZGEMM('N','N',LMMAXD,LMMAXD,NUMN0(IAT)*LMMAXD, &
                 OMEGA,GLLH(1,1,I1),LMMAXD, &
                 SPRSX,NGTBD, &
                 DELTA,AX(IL1B+1,1),NDIM)

    endif

  enddo

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!$omp end parallel
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end subroutine SPRSZMM
