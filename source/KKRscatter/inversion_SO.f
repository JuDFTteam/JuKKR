c ************************************************************************
      SUBROUTINE INVERSION_SO(GLLKE,INVMOD,ICHECK,LMMAXSO)
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

      INTEGER LMMAXSO,LMAXSQ
      PARAMETER (LMAXSQ=(LMAXD+1)**2)
c      PARAMETER (LMMAXSO=NSPD*(LMAXD+1)**2)
      INTEGER ALMD,NDIM
      PARAMETER (ALMD= NAEZD*LMAXSQ,NDIM = NPRINCD*LMAXSQ)
      INTEGER ALMDSO,NDIMSO
      PARAMETER (ALMDSO=NAEZD*LMAXSQ*NSPOD,NDIMSO=NPRINCD*LMAXSQ*NSPOD)
      INTEGER NPRINCDSO
      PARAMETER (NPRINCDSO=NPRINCD*NSPOD)
      DOUBLE COMPLEX CI,CZERO,CONE
      PARAMETER (CI=(0.D0,1.D0),CZERO=(0.D0,0.D0),CONE=(1.D0,0.D0))

      DOUBLE COMPLEX GLLKE(ALMDSO,ALMDSO),GTEMP(ALMDSO,ALMDSO),
     +     GDI(NDIMSO,NDIMSO,NLAYERD),GUP(NDIMSO,NDIMSO,NLAYERD),
     +     GDOW(NDIMSO,NDIMSO,NLAYERD)
      INTEGER I,I1,IP1,II1,IL1,LDI1,IP2,II2,IL2,LDI2,J,INVMOD
      INTEGER LM1,LM2,INFO,IPVT(ALMDSO),NLAYER!,LMMAXSO
      INTEGER ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD)  ! changed 3.11.99
      double precision time0,time1,time2,dclock,time3 
      EXTERNAL ZGETRF,ZGETRS,ZCOPY,INVSLAB


      NLAYER=NAEZD/NPRINCD

      IF (INVMOD.EQ.0) THEN     ! total matrix inversion



         DO I=1,ALMDSO
            DO J=1,ALMDSO
               GTEMP(I,J)=CZERO
               IF (I.EQ.J) THEN
                  GTEMP(I,J)=CONE
               ENDIF
            ENDDO
         ENDDO


c     write (6,*) '-------full inversion calculation--------'
         ! time0 = dclock()
         CALL ZGETRF(ALMDSO,ALMDSO,GLLKE,ALMDSO,IPVT,INFO)
         ! time1= dclock()
         CALL ZGETRS('N',ALMDSO,ALMDSO,GLLKE,ALMDSO,IPVT,GTEMP,
     +                                               ALMDSO,INFO)
         ! time2=dclock() 
         CALL ZCOPY(ALMDSO*ALMDSO,GTEMP,1,GLLKE,1)
         ! time3=dclock()
         ! write(6,*) 'Full inv time',time1-time0,time2-time0,time3-time0

      ELSE IF ((INVMOD.GE.1).AND.(INVMOD.LE.2)) THEN ! slab or supercell
                                ! inversion


         DO 337 I1 = 1,NLAYERD
         DO 337 IP1 = 1,NPRINCD
         DO 337 IP2 = 1,NPRINCD
            II1 = (I1-1)*NPRINCD+IP1
            II2 = (I1-1)*NPRINCD+IP2
            DO 347 LM1 = 1,LMMAXSO
            DO 347 LM2 = 1,LMMAXSO
               LDI1 = LMMAXSO*(IP1-1)+LM1
               IL1 = LMMAXSO*(II1-1)+LM1
               LDI2 = LMMAXSO*(IP2-1)+LM2
               IL2 = LMMAXSO*(II2-1)+LM2
               GDI(LDI1,LDI2,I1) = GLLKE(IL1,IL2)
 347        CONTINUE
 337     CONTINUE

c     this part now is correct also for    ! changes 20/10/99
c     supercell geometry : 20/10/99       
c---> upper linear part
       DO 357 I1 = 1,NLAYERD
       DO 357 IP1 = 1,NPRINCD
       DO 357 IP2 = 1,NPRINCD
          DO 367 LM1 = 1,LMMAXSO
          DO 367 LM2 = 1,LMMAXSO
             LDI1 = LMMAXSO*(IP1-1)+LM1
             LDI2 = LMMAXSO*(IP2-1)+LM2
             IF(I1.LE.(NLAYERD-1)) THEN
                II1 = (I1-1)*NPRINCD+IP1
                II2 = I1*NPRINCD+IP2
                IL1 = LMMAXSO*(II1-1)+LM1
                IL2 = LMMAXSO*(II2-1)+LM2
                GUP(LDI1,LDI2,I1) =  GLLKE(IL1,IL2)
             ELSE
                II1 = IP1
                II2 = (NLAYERD-1)*NPRINCD+IP2
                IL1 = LMMAXSO*(II1-1)+LM1
                IL2 = LMMAXSO*(II2-1)+LM2
                GDOW(LDI1,LDI2,I1) = GLLKE(IL1,IL2)
             ENDIF
 367      CONTINUE
 357   CONTINUE


c---> lower linear part
       DO 333 I1 = 1,NLAYERD 
       DO 333 IP1 = 1,NPRINCD
       DO 333 IP2 = 1,NPRINCD    
          DO 335 LM1 = 1,LMMAXSO
          DO 335 LM2 = 1,LMMAXSO
             LDI1 = LMMAXSO*(IP1-1)+LM1
             LDI2 = LMMAXSO*(IP2-1)+LM2
             IF(I1.LE.(NLAYERD-1)) THEN
                II1 = I1*NPRINCD+IP1
                II2 = (I1-1)*NPRINCD+IP2
                IL1 = LMMAXSO*(II1-1)+LM1
                IL2 = LMMAXSO*(II2-1)+LM2
                GDOW(LDI1,LDI2,I1) =  GLLKE(IL1,IL2)
             ELSE
                II1 = (NLAYERD-1)*NPRINCD+IP1
                II2 = IP2
                IL1 = LMMAXSO*(II1-1)+LM1
                IL2 = LMMAXSO*(II2-1)+LM2
                GUP(LDI1,LDI2,I1) = GLLKE(IL1,IL2)
             ENDIF
 335      CONTINUE
 333   CONTINUE

c     end of the corrected part  20/10/99

       IF (INVMOD.EQ.1) THEN

          CALL INVSLAB_SO(GDI,GUP,GDOW,GLLKE,ICHECK)

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
