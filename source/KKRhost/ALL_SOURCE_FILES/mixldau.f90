SUBROUTINE mixldau(  &
        mmaxd,nspind,natypd,natyp,nspin,lopt,wldauold,  &
        wldau)
implicit none
! Input:
INTEGER NATYPD,NSPIND,MMAXD
INTEGER LOPT(NATYPD)
INTEGER NATYP,NSPIN
DOUBLE PRECISION WLDAUOLD(MMAXD,MMAXD,NSPIND,NATYPD)
! Input/Output:
DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
! Inside:
INTEGER IAT,IS,M1,M2,MMAX
INTEGER IER
DOUBLE PRECISION XMIX,XMIX2,RMSERR
CHARACTER*256 UIO ! NCOLIO=256

EXTERNAL IOINPUT


! First calculate rms error in interaction matrix

DO iat = 1,natyp
  rmserr = 0.d0
  IF (lopt(iat) >= 0) THEN
    mmax = 2*lopt(iat) + 1
    DO is = 1,nspin
      DO m2 = 1,mmax
        DO m1 = 1,mmax
          rmserr = rmserr + (wldau(m1,m2,is,iat) - wldauold(m1,m2,is,iat))**2
        END DO
      END DO
    END DO
    rmserr = DSQRT(rmserr)
    WRITE(1337,9000) iat,rmserr
    9000       FORMAT('LDA+U interaction matrix rms error for atom',  &
        i6,' = ',e10.2)
  END IF
END DO

! Now mix old/new interaction matrices
ier = 0
CALL ioinput('MIXLDAU         ',uio,1,7,ier)
IF ( ier /= 0 ) THEN
  WRITE(*,*) 'MIXLDAU not found, setting to 1.'
  RETURN
ELSE
  READ (UNIT=uio,FMT=*) xmix
  WRITE(1337,*) 'Using MIXLDAU = ',xmix
END IF

xmix2 = 1.d0 - xmix

DO iat = 1,natyp
  IF (lopt(iat) >= 0) THEN
    DO is = 1,nspin
      DO m2 = 1,mmaxd
        DO m1 = 1,mmaxd
          wldau(m1,m2,is,iat) = xmix * wldau(m1,m2,is,iat)  &
              + xmix2 * wldauold(m1,m2,is,iat)
        END DO
      END DO
    END DO
  END IF
END DO

END SUBROUTINE mixldau
