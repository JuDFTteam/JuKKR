subroutine SPRSZMM(IAT,GLLH,NUMN0,INDN0,X,DONE,OMEGA,DELTA, &  ! <
                   AX, &                                       ! >
                   ! new input parameters after inc.p removal
                   naez, lmax, naclsd, nthrds)

  ! This routine is called very often
  ! TODO: Optimise this routine if possible

  implicit none

  integer, intent(in) :: naez
  integer, intent(in) :: lmax
  integer, intent(in) :: naclsd
  integer, intent(in) :: nthrds

  !     INTEGER           LMMAXD
  !     PARAMETER        (LMMAXD= (LMAXD+1)**2)
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
!  integer::           NUMN0(NAEZD)
!  integer::           INDN0(NAEZD,NACLSD)
!  logical::           DONE(LMMAXD)

  double complex::   X(NAEZ*(LMAX+1)**2,(LMAX+1)**2)
  double complex::  AX(NAEZ*(LMAX+1)**2,(LMAX+1)**2)
  double complex::  GLLH((LMAX+1)**2, NACLSD*(LMAX+1)**2, NAEZ)
  integer::           NUMN0(NAEZ)
  integer::           INDN0(NAEZ,NACLSD)
  logical::           DONE((LMAX+1)**2)

  !     .. Local Arrays .. Fortran 90 automatic array
  !double complex :: SPRSX(NGTBD,LMMAXD)
  double complex :: SPRSX(NACLSD*(LMAX+1)**2, (LMAX+1)**2) ! medium size
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

  integer::LMMAXD
  integer::NGTBD
  integer::NDIM

  LMMAXD= (LMAX+1)**2
  NGTBD = NACLSD*LMMAXD
  NDIM  = NAEZ*LMMAXD

     
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!$ call OMP_SET_NUM_THREADS(NTHRDS)
!$OMP PARALLEL PRIVATE (I1,LM2,I2,I3H,I2H,
!$OMP&                  IL1B,SPRSX,MYTHRD)
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
!$OMP END PARALLEL
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end subroutine SPRSZMM
