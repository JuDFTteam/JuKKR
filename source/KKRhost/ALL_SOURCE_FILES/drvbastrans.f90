SUBROUTINE drvbastrans(rc,crel,rrel,srrel,nrrel,irrel,  &
    nlmax,nkmmax,nmuemax,nkmpmax,nkmax,linmax)
!   ********************************************************************
!   *                                                                  *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

! Dummy arguments
INTEGER LINMAX,NKMAX,NKMMAX,NKMPMAX,NLMAX,NMUEMAX
COMPLEX*16 CREL(NKMMAX,NKMMAX),RC(NKMMAX,NKMMAX), &
           RREL(NKMMAX,NKMMAX),SRREL(2,2,NKMMAX)
INTEGER IRREL(2,2,NKMMAX),NRREL(2,NKMMAX)

! Local variables
REAL*8 CGC(NKMPMAX,2)
INTEGER I,IKM1LIN(LINMAX),IKM2LIN(LINMAX),IL,IMUE,IPRINT, &
        KAPTAB(NMUEMAX),LTAB(NMUEMAX),MMAX,NMUETAB(NMUEMAX), &
        NSOLLM(NLMAX,NMUEMAX)

IF (nkmmax /= 2*nlmax**2) STOP ' Check NLMAX,NKMMAX in < DRVBASTRANS > '
IF (nmuemax /= 2*nlmax) STOP ' Check NLMAX,NMUEMAX in < DRVBASTRANS > '
IF (nkmpmax /= (nkmmax+2*nlmax))  &
    STOP ' Check NLMAX,NKMMAX,NKMPMAX in < DRVBASTRANS > '
IF (nkmax /= 2*nlmax-1) STOP ' Check NLMAX,NKMAX in < DRVBASTRANS > '
IF (linmax /= (2*nlmax*(2*nlmax-1)))  &
    STOP ' Check NLMAX,LINMAX in < DRVBASTRANS > '

iprint = 0

DO i = 1,nmuemax
  ltab(i) = i/2
  IF ( 2*ltab(i) == i ) THEN
    kaptab(i) = ltab(i)
  ELSE
    kaptab(i) = -ltab(i) - 1
  END IF
  nmuetab(i) = 2*ABS(kaptab(i))
END DO

DO il = 1,nlmax
  mmax = 2*il
  DO imue = 1,mmax
    IF ( (imue == 1) .OR. (imue == mmax) ) THEN
      nsollm(il,imue) = 1
    ELSE
      nsollm(il,imue) = 2
    END IF
  END DO
END DO

CALL ikmlin(iprint,nsollm,ikm1lin,ikm2lin,nlmax,nmuemax,linmax, nlmax)

CALL calccgc(ltab,kaptab,nmuetab,cgc,nkmax,nmuemax,nkmpmax)

! ---------------------------- now calculate the transformation matrices

CALL strsmat(nlmax-1,cgc,srrel,nrrel,irrel,nkmmax,nkmpmax)

CALL bastrmat(nlmax-1,cgc,rc,crel,rrel,nkmmax,nkmpmax)

RETURN
END SUBROUTINE drvbastrans

