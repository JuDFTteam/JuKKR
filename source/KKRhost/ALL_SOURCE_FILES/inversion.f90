! ************************************************************************
SUBROUTINE inversion(gllke,invmod,icheck)
! ************************************************************************
! This subroutine calculates the inversion of a matrix
! in 4 different ways depending on the form of the matrix

!     INVMOD = 0  ----> total inversion scheme
!     INVMOD = 1  ----> band matrix inversion scheme
!     INVMOD = 2  ----> corner band matrix inversion scheme
!     INVMOD = 3  ----> sparse matrix inversion scheme

! ------------------------------------------------------------------------
IMPLICIT NONE

!     .. parameters ..
INCLUDE 'inc.p'

! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************

INTEGER LMMAXD
PARAMETER (LMMAXD = (KREL+KORBIT+1) * (LMAXD+1)**2)
INTEGER ALMD,NDIM
PARAMETER (ALMD= NAEZD*LMMAXD,NDIM = NPRINCD*LMMAXD)
DOUBLE COMPLEX CI,CZERO,CONE
PARAMETER (CI=(0.D0,1.D0),CZERO=(0.D0,0.D0),CONE=(1.D0,0.D0))

DOUBLE COMPLEX GLLKE(ALMD,ALMD), &
     GDI(NDIM,NDIM,NLAYERD),GUP(NDIM,NDIM,NLAYERD), &
     GDOW(NDIM,NDIM,NLAYERD)
double complex, allocatable :: GTEMP(:,:)
INTEGER I,I1,IP1,II1,IL1,LDI1,IP2,II2,IL2,LDI2,J,INVMOD
INTEGER LM1,LM2,INFO,IPVT(ALMD),NLAYER
INTEGER ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD)  ! changed 3.11.99

EXTERNAL ZGETRF,ZGETRS,ZCOPY,INVSLAB

allocate(gtemp(almd,almd))

nlayer=naezd/nprincd

IF (invmod == 0) THEN     ! total matrix inversion
  
  
  
  DO i=1,almd
    DO j=1,almd
      gtemp(i,j)=czero
      IF (i == j) THEN
        gtemp(i,j)=cone
      END IF
    END DO
  END DO
  
  
!     write (6,*) '-------full inversion calculation--------'
  
  CALL zgetrf(almd,almd,gllke,almd,ipvt,info)
  CALL zgetrs('N',almd,almd,gllke,almd,ipvt,gtemp,almd,info)
  
  CALL zcopy(almd*almd,gtemp,1,gllke,1)
  
  
  
ELSE IF ((invmod >= 1).AND.(invmod <= 2)) THEN ! slab or supercell
! inversion
  
  
  DO  i1 = 1,nlayerd
    DO  ip1 = 1,nprincd
      DO  ip2 = 1,nprincd
        ii1 = (i1-1)*nprincd+ip1
        ii2 = (i1-1)*nprincd+ip2
        DO  lm1 = 1,lmmaxd
          DO  lm2 = 1,lmmaxd
            ldi1 = lmmaxd*(ip1-1)+lm1
            il1 = lmmaxd*(ii1-1)+lm1
            ldi2 = lmmaxd*(ip2-1)+lm2
            il2 = lmmaxd*(ii2-1)+lm2
            gdi(ldi1,ldi2,i1) = gllke(il1,il2)
          END DO
        END DO
      END DO
    END DO
  END DO
  
!     this part now is correct also for    ! changes 20/10/99
!     supercell geometry : 20/10/99
!---> upper linear part
  DO  i1 = 1,nlayerd
    DO  ip1 = 1,nprincd
      DO  ip2 = 1,nprincd
        DO  lm1 = 1,lmmaxd
          DO  lm2 = 1,lmmaxd
            ldi1 = lmmaxd*(ip1-1)+lm1
            ldi2 = lmmaxd*(ip2-1)+lm2
            IF(i1 <= (nlayerd-1)) THEN
              ii1 = (i1-1)*nprincd+ip1
              ii2 = i1*nprincd+ip2
              il1 = lmmaxd*(ii1-1)+lm1
              il2 = lmmaxd*(ii2-1)+lm2
              gup(ldi1,ldi2,i1) =  gllke(il1,il2)
            ELSE
              ii1 = ip1
              ii2 = (nlayerd-1)*nprincd+ip2
              il1 = lmmaxd*(ii1-1)+lm1
              il2 = lmmaxd*(ii2-1)+lm2
              gdow(ldi1,ldi2,i1) = gllke(il1,il2)
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO
  
  
!---> lower linear part
  DO  i1 = 1,nlayerd
    DO  ip1 = 1,nprincd
      DO  ip2 = 1,nprincd
        DO  lm1 = 1,lmmaxd
          DO  lm2 = 1,lmmaxd
            ldi1 = lmmaxd*(ip1-1)+lm1
            ldi2 = lmmaxd*(ip2-1)+lm2
            IF(i1 <= (nlayerd-1)) THEN
              ii1 = i1*nprincd+ip1
              ii2 = (i1-1)*nprincd+ip2
              il1 = lmmaxd*(ii1-1)+lm1
              il2 = lmmaxd*(ii2-1)+lm2
              gdow(ldi1,ldi2,i1) =  gllke(il1,il2)
            ELSE
              ii1 = (nlayerd-1)*nprincd+ip1
              ii2 = ip2
              il1 = lmmaxd*(ii1-1)+lm1
              il2 = lmmaxd*(ii2-1)+lm2
              gup(ldi1,ldi2,i1) = gllke(il1,il2)
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO
  
!     end of the corrected part  20/10/99
  
  IF (invmod == 1) THEN
    
    CALL invslab(gdi,gup,gdow,gllke,icheck)
    
!          write (6,*) '-------slab calculation--------'
    
  ELSE IF (invmod == 2) THEN ! supercell geometry inversion
    
    CALL invsupercell(gdi,gup,gdow,gllke,icheck)
    
!          write (6,*) '-------supercell calculation--------'
    
  END IF
  
  
ELSE                      ! sparse matrix inversion
  
!     NOT YET IMPLEMENTED!!!!!!!!!
  
  
END IF


deallocate(gtemp)


RETURN

END SUBROUTINE inversion
