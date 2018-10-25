! ************************************************************************
      SUBROUTINE INVERSION(GLLKE,INVMOD,ICHECK)
! ************************************************************************
! This subroutine calculates the inversion of a matrix
! in 4 different ways depending on the form of the matrix
!
!     INVMOD = 0  ----> total inversion scheme
!     INVMOD = 1  ----> band matrix inversion scheme
!     INVMOD = 2  ----> corner band matrix inversion scheme
!     INVMOD = 3  ----> gofrin module
!
! ------------------------------------------------------------------------
      use godfrin ! GODFRIN Flaviano

      implicit none

!     .. parameters ..
      include 'inc.p'
!
! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************
!

      INTEGER LMMAXD
      PARAMETER (LMMAXD = (KREL+KORBIT+1) * (LMAXD+1)**2)
      INTEGER ALMD,NDIM
      PARAMETER (ALMD= NAEZD*LMMAXD,NDIM = NPRINCD*LMMAXD)
      DOUBLE COMPLEX CI,CZERO,CONE
      PARAMETER (CI=(0.D0,1.D0),CZERO=(0.D0,0.D0),CONE=(1.D0,0.D0))

      DOUBLE COMPLEX GLLKE(ALMD,ALMD),!GTEMP(ALMD,ALMD),
     +     GDI(NDIM,NDIM,NLAYERD),GUP(NDIM,NDIM,NLAYERD),
     +     GDOW(NDIM,NDIM,NLAYERD)
      double complex, allocatable :: GTEMP(:,:)
      INTEGER I,I1,IP1,II1,IL1,LDI1,IP2,II2,IL2,LDI2,J,INVMOD
      INTEGER LM1,LM2,INFO,IPVT(ALMD),NLAYER
      INTEGER ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD)  ! changed 3.11.99

      EXTERNAL ZGETRF,ZGETRS,ZCOPY,INVSLAB

      allocate(GTEMP(ALMD,ALMD))

!     Naive number of layers in each principal layer
      NLAYER=NAEZD/NPRINCD

!  ---------------------------------------------------------------------
!                        full matrix inversion
!  ---------------------------------------------------------------------
      IF (INVMOD==0) THEN
!       initialize unit matrix
        DO I=1,ALMD
          DO J=1,ALMD
            GTEMP(I,J)=CZERO
            IF (I.EQ.J) THEN
               GTEMP(I,J)=CONE
            END IF
          END DO
        END DO
!     write (6,*) '-------full inversion calculation--------'
        CALL ZGETRF(ALMD,ALMD,GLLKE,ALMD,IPVT,INFO)
        CALL ZGETRS('N',ALMD,ALMD,GLLKE,ALMD,IPVT,GTEMP,ALMD,INFO)
        CALL ZCOPY(ALMD*ALMD,GTEMP,1,GLLKE,1)
!  ---------------------------------------------------------------------
!            block tridiagonal inversion (slab or supercell)
!  ---------------------------------------------------------------------
      ELSE IF (INVMOD==1 .OR. INVMOD==2) THEN
!       copy blocks of assumed tridiagonal matrix
        DO I1=1,NLAYERD
          DO IP1=1,NPRINCD
          DO IP2=1,NPRINCD
            II1 = (I1-1)*NPRINCD+IP1
            II2 = (I1-1)*NPRINCD+IP2
            DO LM1=1,LMMAXD
            DO LM2=1,LMMAXD
              LDI1 = LMMAXD*(IP1-1)+LM1
              IL1 = LMMAXD*(II1-1)+LM1
              LDI2 = LMMAXD*(IP2-1)+LM2
              IL2 = LMMAXD*(II2-1)+LM2
              GDI(LDI1,LDI2,I1) = GLLKE(IL1,IL2)
            END DO
            END DO
          END DO
          END DO
        END DO
!       this part now is correct also for    ! changes 20/10/99
!       supercell geometry : 20/10/99       
!       --> upper linear part
        DO I1=1,NLAYERD
          DO IP1=1,NPRINCD
          DO IP2=1,NPRINCD
            DO LM1=1,LMMAXD
            DO LM2=1,LMMAXD
              LDI1 = LMMAXD*(IP1-1)+LM1
              LDI2 = LMMAXD*(IP2-1)+LM2
              IF (I1.LE.(NLAYERD-1)) THEN
                II1 = (I1-1)*NPRINCD+IP1
                II2 = I1*NPRINCD+IP2
                IL1 = LMMAXD*(II1-1)+LM1
                IL2 = LMMAXD*(II2-1)+LM2
                GUP(LDI1,LDI2,I1) =  GLLKE(IL1,IL2)
              ELSE
                II1 = IP1
                II2 = (NLAYERD-1)*NPRINCD+IP2
                IL1 = LMMAXD*(II1-1)+LM1
                IL2 = LMMAXD*(II2-1)+LM2
                GDOW(LDI1,LDI2,I1) = GLLKE(IL1,IL2)
              END IF
            END DO
            END DO
          END DO
          END DO
        END DO
!       --> lower linear part
        DO I1=1,NLAYERD 
          DO IP1=1,NPRINCD
          DO IP2=1,NPRINCD    
            DO LM1=1,LMMAXD
            DO LM2=1,LMMAXD
              LDI1 = LMMAXD*(IP1-1)+LM1
              LDI2 = LMMAXD*(IP2-1)+LM2
              IF (I1.LE.(NLAYERD-1)) THEN
                II1 = I1*NPRINCD+IP1
                II2 = (I1-1)*NPRINCD+IP2
                IL1 = LMMAXD*(II1-1)+LM1
                IL2 = LMMAXD*(II2-1)+LM2
                GDOW(LDI1,LDI2,I1) =  GLLKE(IL1,IL2)
              ELSE
                II1 = (NLAYERD-1)*NPRINCD+IP1
                II2 = IP2
                IL1 = LMMAXD*(II1-1)+LM1
                IL2 = LMMAXD*(II2-1)+LM2
                GUP(LDI1,LDI2,I1) = GLLKE(IL1,IL2)
              END IF
            END DO
            END DO
          END DO
          END DO
        END DO
!       end of the corrected part  20/10/99
!
!       slab: matrix is block tridiagonal
        IF (INVMOD==1) THEN
          CALL INVSLAB(GDI,GUP,GDOW,GLLKE,ICHECK)
!          write (6,*) '-------slab calculation--------'
!       supercell: matrix is tridiagonal with corner blocks
        ELSE IF (INVMOD==2) THEN
          CALL INVSUPERCELL(GDI,GUP,GDOW,GLLKE,ICHECK)
!          write (6,*) '-------supercell calculation--------'
        ENDIF
!  ---------------------------------------------------------------------
!                          godfrin module
!  ---------------------------------------------------------------------
      ELSE IF (INVMOD==3) THEN
        call sparse_inverse(gllke,t_godfrin%na,t_godfrin%nb,
     +    t_godfrin%bdims,t_godfrin%ldiag,t_godfrin%lper,
     +    t_godfrin%lpardiso)  ! GODFRIN Flaviano
!  ---------------------------------------------------------------------
      ELSE
!       if it gets here, did you have a coffee before running the code?
        STOP 'UNKNOWN INVERSION MODE !!!'
      ENDIF
!  ---------------------------------------------------------------------


      deallocate(GTEMP)

C
      RETURN
C
      END
