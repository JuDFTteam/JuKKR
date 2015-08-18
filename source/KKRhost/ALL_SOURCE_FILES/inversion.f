c ************************************************************************
      SUBROUTINE INVERSION(GLLKE,INVMOD,ICHECK)
c ************************************************************************
c This subroutine calculates the inversion of a matrix
c in 4 different ways depending on the form of the matrix
c
c     INVMOD = 0  ----> total inversion scheme
c     INVMOD = 1  ----> band matrix inversion scheme
c     INVMOD = 2  ----> corner band matrix inversion scheme
c     INVMOD = 3  ----> sparse matrix inversion scheme
c
c ------------------------------------------------------------------------

      implicit none

c     .. parameters ..
      include 'inc.p'
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  NPOTD = 2 * NATYPD                                               *
C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
C *  NSPIND = 1                                                       *
C *                                                                   *
C *********************************************************************
C

      INTEGER LMMAXD
      PARAMETER (LMMAXD = (KREL+KORBIT+1) * (LMAXD+1)**2)
      INTEGER ALMD,NDIM
      PARAMETER (ALMD= NAEZD*LMMAXD,NDIM = NPRINCD*LMMAXD)
      DOUBLE COMPLEX CI,CZERO,CONE
      PARAMETER (CI=(0.D0,1.D0),CZERO=(0.D0,0.D0),CONE=(1.D0,0.D0))

      DOUBLE COMPLEX GLLKE(ALMD,ALMD),GTEMP(ALMD,ALMD),
     +     GDI(NDIM,NDIM,NLAYERD),GUP(NDIM,NDIM,NLAYERD),
     +     GDOW(NDIM,NDIM,NLAYERD)
      INTEGER I,I1,IP1,II1,IL1,LDI1,IP2,II2,IL2,LDI2,J,INVMOD
      INTEGER LM1,LM2,INFO,IPVT(ALMD),NLAYER
      INTEGER ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD)  ! changed 3.11.99

      EXTERNAL ZGETRF,ZGETRS,ZCOPY,INVSLAB


      NLAYER=NAEZD/NPRINCD

      IF (INVMOD.EQ.0) THEN     ! total matrix inversion



         DO I=1,ALMD
            DO J=1,ALMD
               GTEMP(I,J)=CZERO
               IF (I.EQ.J) THEN
                  GTEMP(I,J)=CONE
               ENDIF
            ENDDO
         ENDDO


c     write (6,*) '-------full inversion calculation--------'

         CALL ZGETRF(ALMD,ALMD,GLLKE,ALMD,IPVT,INFO)
         CALL ZGETRS('N',ALMD,ALMD,GLLKE,ALMD,IPVT,GTEMP,ALMD,INFO)

         CALL ZCOPY(ALMD*ALMD,GTEMP,1,GLLKE,1)



      ELSE IF ((INVMOD.GE.1).AND.(INVMOD.LE.2)) THEN ! slab or supercell
                                ! inversion


         DO 337 I1 = 1,NLAYERD
         DO 337 IP1 = 1,NPRINCD
         DO 337 IP2 = 1,NPRINCD
            II1 = (I1-1)*NPRINCD+IP1
            II2 = (I1-1)*NPRINCD+IP2
            DO 347 LM1 = 1,LMMAXD
            DO 347 LM2 = 1,LMMAXD
               LDI1 = LMMAXD*(IP1-1)+LM1
               IL1 = LMMAXD*(II1-1)+LM1
               LDI2 = LMMAXD*(IP2-1)+LM2
               IL2 = LMMAXD*(II2-1)+LM2
               GDI(LDI1,LDI2,I1) = GLLKE(IL1,IL2)
 347        CONTINUE
 337     CONTINUE

c     this part now is correct also for    ! changes 20/10/99
c     supercell geometry : 20/10/99       
c---> upper linear part
       DO 357 I1 = 1,NLAYERD
       DO 357 IP1 = 1,NPRINCD
       DO 357 IP2 = 1,NPRINCD
          DO 367 LM1 = 1,LMMAXD
          DO 367 LM2 = 1,LMMAXD
             LDI1 = LMMAXD*(IP1-1)+LM1
             LDI2 = LMMAXD*(IP2-1)+LM2
             IF(I1.LE.(NLAYERD-1)) THEN
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
             ENDIF
 367      CONTINUE
 357   CONTINUE


c---> lower linear part
       DO 333 I1 = 1,NLAYERD 
       DO 333 IP1 = 1,NPRINCD
       DO 333 IP2 = 1,NPRINCD    
          DO 335 LM1 = 1,LMMAXD
          DO 335 LM2 = 1,LMMAXD
             LDI1 = LMMAXD*(IP1-1)+LM1
             LDI2 = LMMAXD*(IP2-1)+LM2
             IF(I1.LE.(NLAYERD-1)) THEN
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
             ENDIF
 335      CONTINUE
 333   CONTINUE

c     end of the corrected part  20/10/99

       IF (INVMOD.EQ.1) THEN

          CALL INVSLAB(GDI,GUP,GDOW,GLLKE,ICHECK)

c          write (6,*) '-------slab calculation--------'

       ELSE IF (INVMOD.EQ.2) THEN ! supercell geometry inversion

           CALL INVSUPERCELL(GDI,GUP,GDOW,GLLKE,ICHECK)

c          write (6,*) '-------supercell calculation--------'

       ENDIF


      ELSE                      ! sparse matrix inversion

C     NOT YET IMPLEMENTED!!!!!!!!!


      ENDIF

C
      RETURN
C
      END
