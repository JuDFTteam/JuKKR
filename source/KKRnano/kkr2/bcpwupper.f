      SUBROUTINE BCPWUPPER(
     >                     GLLH,
     <                     GLLHBLCK,
     >                     NAEZ,NUMN0,INDN0,
C                          new input after inc.p removal
     &                     lmax, nthrds, natbld, xdim, ydim, zdim,
     &                     naclsd)
C
      IMPLICIT NONE
C ------------------------------------------------------------------------
C
C
C
C ------------------------------------------------------------------------
C     .. parameters ..
C      include 'inc.p'
C      include 'inc.cls'

      INTEGER NAEZ

C     number of preconditioning blocks in each direction
      integer xdim
      integer ydim
      integer zdim

      integer lmax

C     number of atoms in the reference cluster
      integer naclsd
C     number of atoms per preconditioning block
      integer natbld
C     number of OpenMP threads
      INTEGER nthrds

C
      INTEGER          NBLCKD
C     NBLCKD = XDIM*YDIM*ZDIM

      INTEGER   NSYMAXD
      PARAMETER (NSYMAXD=48)

      INTEGER   LMMAXD
C     LMMAXD= (LMAX+1)**2

      INTEGER   ALM
C     PARAMETER (ALM = NAEZD*LMMAXD)

      INTEGER   NGTBD
C      PARAMETER (NGTBD = NACLSD*LMMAXD)
C                       = NACLSD*(LMAX+1)**2

      INTEGER   NINTACTD
      PARAMETER (NINTACTD = 19)     ! Ecken des "Interaktions-wuerfels"
C                                     nicht mit einbezogen, incl. waere
C                                     27
C

C      DOUBLE COMPLEX     GLLH(LMMAXD,NGTBD,NAEZD),
C     +                   GLLHBLCK(LMMAXD*NATBLD,LMMAXD*NATBLD*NBLCKD)

      DOUBLE COMPLEX     GLLH((LMAX+1)**2,NACLSD*(LMAX+1)**2,NAEZ),
     &                   GLLHBLCK(NATBLD*(LMAX+1)**2,
     &                            NATBLD*XDIM*YDIM*ZDIM*(LMAX+1)**2)

      INTEGER            INDN0(NAEZ,NACLSD),NUMN0(NAEZ)
C ..
C ..
C local arrays ..
C      DOUBLE COMPLEX     BLAV(NATBLD*LMMAXD,NATBLD*LMMAXD,NINTACTD),
C     +                   TMPBLCK(NATBLD*LMMAXD,NATBLD*LMMAXD)

      DOUBLE COMPLEX     BLAV(NATBLD*(LMAX+1)**2,
     &                        NATBLD*(LMAX+1)**2, NINTACTD),

     &                   TMPBLCK(NATBLD*(LMAX+1)**2,NATBLD*(LMAX+1)**2)

      INTEGER            IPIV(NATBLD*(LMAX+1)**2)
C ..
C local scalars ..
      DOUBLE COMPLEX     CONE,IMONE,CZERO
      DOUBLE PRECISION   PI
      INTEGER            I1,I2,LM1,LM2,IL1,IL2,IL2B,
     +                   IX,IY,IZ,NDIM,
     +                   IX2,IY2,IZ2,
     +                   IINT,
     +                   CEN,LEF,RIG,DOW,UPP,FOW,BAC,
     +                   LEDO,LEUP,LEBA,LEFO,RIDO,RIUP,RIBA,RIFO,
     +                   DOBA,DOFO,UPBA,UPFO,
     +                   INFO
C ..
C openMP ..
      INTEGER        MYTHRD,OMP_GET_THREAD_NUM
C ..
C intrinsic functions ..
      INTRINSIC          ZEXP,ATAN
C ..
C=======================================================================
C

      NBLCKD = XDIM*YDIM*ZDIM
      LMMAXD= (LMAX+1)**2

      PI = 4.D0*ATAN(1.D0)
C
      CONE  = ( 1.0D0,0.0D0 )
      CZERO  = ( 0.0D0,0.0D0 )
      IMONE  = ( 0.0D0,1.0D0 )
C
      DO I1=1,LMMAXD*NATBLD
        DO I2=1,LMMAXD*NATBLD*NBLCKD
          GLLHBLCK(I1,I2) = CZERO
        ENDDO
      ENDDO
      DO I1=1,LMMAXD*NATBLD
        DO I2=1,LMMAXD*NATBLD
          DO IINT = 1, NINTACTD
            BLAV(I1,I2,IINT) = CZERO
          ENDDO
        ENDDO
      ENDDO
C
C
      DO IX = 1, XDIM
        IX2 = IX - 2
        IF (IX2.LT.0) IX2 = IX2 + XDIM

        DO IY = 1, YDIM
          IY2 = IY - 2
          IF (IY2.LT.0) IY2 = IY2 + YDIM

          DO IZ = 1, ZDIM
            IZ2 = IZ - 2
            IF (IZ2.LT.0) IZ2 = IZ2 + ZDIM
C
            CEN = MOD(IX-1,XDIM) + 1 +
     +            MOD(IY-1,YDIM)*XDIM+
     +            MOD(IZ-1,ZDIM)*XDIM*YDIM
            LEF = MOD(IX2,XDIM) + 1 +
     +            MOD(IY-1,YDIM)*XDIM+
     +            MOD(IZ-1,ZDIM)*XDIM*YDIM
            RIG = MOD(IX,XDIM) + 1 +
     +            MOD(IY-1,YDIM)*XDIM+
     +            MOD(IZ-1,ZDIM)*XDIM*YDIM
            DOW = MOD(IX-1,XDIM) + 1 +
     +            MOD(IY2,YDIM)*XDIM+
     +            MOD(IZ-1,ZDIM)*XDIM*YDIM
            UPP = MOD(IX-1,XDIM) + 1 +
     +            MOD(IY,YDIM)*XDIM+
     +            MOD(IZ-1,ZDIM)*XDIM*YDIM
            BAC = MOD(IX-1,XDIM) + 1 +
     +            MOD(IY-1,YDIM)*XDIM+
     +            MOD(IZ2,ZDIM)*XDIM*YDIM
            FOW = MOD(IX-1,XDIM) + 1 +
     +            MOD(IY-1,YDIM)*XDIM+
     +            MOD(IZ,ZDIM)*XDIM*YDIM
C
C
            CALL GENBLAV(GLLH,BLAV(1,1,1),CEN,CEN,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,2),CEN,LEF,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,3),CEN,RIG,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,4),CEN,DOW,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,5),CEN,UPP,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,6),CEN,BAC,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,7),CEN,FOW,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
C
C
            LEDO = MOD(IX2,XDIM) + 1 +
     +             MOD(IY2,YDIM)*XDIM+
     +             MOD(IZ-1,ZDIM)*XDIM*YDIM
            LEUP = MOD(IX2,XDIM) + 1 +
     +             MOD(IY,YDIM)*XDIM+
     +             MOD(IZ-1,ZDIM)*XDIM*YDIM
            LEBA = MOD(IX2,XDIM) + 1 +
     +             MOD(IY-1,YDIM)*XDIM+
     +             MOD(IZ2,ZDIM)*XDIM*YDIM
            LEFO = MOD(IX2,XDIM) + 1 +
     +             MOD(IY-1,YDIM)*XDIM+
     +             MOD(IZ,ZDIM)*XDIM*YDIM
            RIDO = MOD(IX,XDIM) + 1 +
     +             MOD(IY2,YDIM)*XDIM+
     +             MOD(IZ-1,ZDIM)*XDIM*YDIM
            RIUP = MOD(IX,XDIM) + 1 +
     +             MOD(IY,YDIM)*XDIM+
     +             MOD(IZ-1,ZDIM)*XDIM*YDIM
            RIBA = MOD(IX,XDIM) + 1 +
     +             MOD(IY-1,YDIM)*XDIM+
     +             MOD(IZ2,ZDIM)*XDIM*YDIM
            RIFO = MOD(IX,XDIM) + 1 +
     +             MOD(IY-1,YDIM)*XDIM+
     +             MOD(IZ,ZDIM)*XDIM*YDIM
            DOBA = MOD(IX-1,XDIM) + 1 +
     +             MOD(IY2,YDIM)*XDIM+
     +             MOD(IZ2,ZDIM)*XDIM*YDIM
            DOFO = MOD(IX-1,XDIM) + 1 +
     +             MOD(IY2,YDIM)*XDIM+
     +             MOD(IZ,ZDIM)*XDIM*YDIM
            UPBA = MOD(IX-1,XDIM) + 1 +
     +             MOD(IY,YDIM)*XDIM+
     +             MOD(IZ2,ZDIM)*XDIM*YDIM
            UPFO = MOD(IX-1,XDIM) + 1 +
     +             MOD(IY,YDIM)*XDIM+
     +             MOD(IZ,ZDIM)*XDIM*YDIM
C  
            CALL GENBLAV(GLLH,BLAV(1,1,8),CEN,LEDO,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,9),CEN,LEUP,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,10),CEN,LEBA,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,11),CEN,LEFO,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,12),CEN,RIDO,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,13),CEN,RIUP,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,14),CEN,RIBA,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,15),CEN,RIFO,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,16),CEN,DOBA,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,17),CEN,DOFO,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,18),CEN,UPBA,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
            CALL GENBLAV(GLLH,BLAV(1,1,19),CEN,UPFO,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
C
C
          ENDDO
        ENDDO
      ENDDO
C
C normalize ...
C
      DO IINT=1,NINTACTD
        DO I1=1,NATBLD*LMMAXD
          DO I2=1,NATBLD*LMMAXD
            BLAV(I1,I2,IINT) = BLAV(I1,I2,IINT)/NBLCKD
          ENDDO
        ENDDO
      ENDDO
C
C
C
C ..
C
      NDIM = (NATBLD*LMMAXD)*(NATBLD*LMMAXD)
C
      CALL OMP_SET_NUM_THREADS(NTHRDS)
!$OMP PARALLEL PRIVATE (IX,IY,IZ,TMPBLCK,MYTHRD,INFO,IPIV)
      MYTHRD = OMP_GET_THREAD_NUM()
C
      DO IX = 1, XDIM
      IF (MOD(IX,MIN(XDIM,NTHRDS)).EQ.MYTHRD) THEN
        DO IY = 1, YDIM
          DO IZ = 1, ZDIM
C .. scales a vector by a constant >> TMPBLCK = CZERO ..
C            CALL ZSCAL(NDIM,CZERO,TMPBLCK,1)
C
            DO LM1 = 1, NATBLD*LMMAXD
              DO LM2 = 1, NATBLD*LMMAXD
                TMPBLCK(LM1,LM2) = CZERO
              ENDDO
            ENDDO
C ..
C .. constant times a vector plus a vector >> TMPBLCK = TMPBLCK + fac*BLAV ..
            CALL ZAXPY(NDIM,
     +                 CONE,BLAV(1,1,1),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(-1.0D0)*(IX-1)/XDIM),
     +                 BLAV(1,1,2),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(1.0D0)*(IX-1)/XDIM),
     +                 BLAV(1,1,3),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(-1.0D0)*(IY-1)/YDIM),
     +                 BLAV(1,1,4),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(1.0D0)*(IY-1)/YDIM),
     +                 BLAV(1,1,5),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(-1.0D0)*(IZ-1)/ZDIM),
     +                 BLAV(1,1,6),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(1.0D0)*(IZ-1)/ZDIM),
     +                 BLAV(1,1,7),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(-1.0D0)*(IX-1)/XDIM)*
     +                 ZEXP(-2*PI*IMONE*(-1.0D0)*(IY-1)/YDIM),
     +                 BLAV(1,1,8),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(-1.0D0)*(IX-1)/XDIM)*
     +                 ZEXP(-2*PI*IMONE*(1.0D0)*(IY-1)/YDIM),
     +                 BLAV(1,1,9),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(-1.0D0)*(IX-1)/XDIM)*
     +                 ZEXP(-2*PI*IMONE*(-1.0D0)*(IZ-1)/ZDIM),
     +                 BLAV(1,1,10),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(-1.0D0)*(IX-1)/XDIM)*
     +                 ZEXP(-2*PI*IMONE*(1.0D0)*(IZ-1)/ZDIM),
     +                 BLAV(1,1,11),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(1.0D0)*(IX-1)/XDIM)*
     +                 ZEXP(-2*PI*IMONE*(-1.0D0)*(IY-1)/YDIM),
     +                 BLAV(1,1,12),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(1.0D0)*(IX-1)/XDIM)*
     +                 ZEXP(-2*PI*IMONE*(1.0D0)*(IY-1)/YDIM),
     +                 BLAV(1,1,13),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(1.0D0)*(IX-1)/XDIM)*
     +                 ZEXP(-2*PI*IMONE*(-1.0D0)*(IZ-1)/ZDIM),
     +                 BLAV(1,1,14),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(1.0D0)*(IX-1)/XDIM)*
     +                 ZEXP(-2*PI*IMONE*(1.0D0)*(IZ-1)/ZDIM),
     +                 BLAV(1,1,15),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(-1.0D0)*(IY-1)/YDIM)*
     +                 ZEXP(-2*PI*IMONE*(-1.0D0)*(IZ-1)/ZDIM),
     +                 BLAV(1,1,16),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(-1.0D0)*(IY-1)/YDIM)*
     +                 ZEXP(-2*PI*IMONE*(1.0D0)*(IZ-1)/ZDIM),
     +                 BLAV(1,1,17),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(1.0D0)*(IY-1)/YDIM)*
     +                 ZEXP(-2*PI*IMONE*(-1.0D0)*(IZ-1)/ZDIM),
     +                 BLAV(1,1,18),1,TMPBLCK,1)
            CALL ZAXPY(NDIM,
     +                 ZEXP(-2*PI*IMONE*(1.0D0)*(IY-1)/YDIM)*
     +                 ZEXP(-2*PI*IMONE*(1.0D0)*(IZ-1)/ZDIM),
     +                 BLAV(1,1,19),1,TMPBLCK,1)
C ..
C .. initialize GLLHBLCK as identity-matrix
            DO I1 = 1,NATBLD*LMMAXD
              GLLHBLCK(I1,NATBLD*LMMAXD*((IX-1)+
     +                 (IY-1)*XDIM+(IZ-1)*XDIM*YDIM)+I1) = CONE
            ENDDO
C ..
C
C ..
            CALL ZGESV(NATBLD*LMMAXD,NATBLD*LMMAXD,
     +                 TMPBLCK,NATBLD*LMMAXD,
     +                 IPIV,
     +                 GLLHBLCK(1,NATBLD*LMMAXD*((IX-1)+
     +                 (IY-1)*XDIM+(IZ-1)*XDIM*YDIM)+1),
     +                 NATBLD*LMMAXD,INFO)
C
          ENDDO
        ENDDO
      ENDIF
      ENDDO
C
!$OMP END PARALLEL
C
      RETURN
C
      END
C
C
C
C
C
C
C
C
C
      SUBROUTINE GENBLAV(GLLH,BLAV,I,J,
     &                   INDN0,NUMN0,
     &                   naez, lmax, natbld, naclsd )
C
      IMPLICIT NONE
C
C      include 'inc.p'
C      include 'inc.cls'

      INTEGER naez
      INTEGER lmax
      INTEGER natbld
      INTEGER naclsd
C
      INTEGER            LMMAXD
C     PARAMETER         (LMMAXD= (LMAXD+1)**2)
C
      INTEGER            I,J,BS,II,JJ,I1,I2,LM1,LM2,IL1,IL2,IL2B,
     +                   ISRH

      INTEGER            INDN0(NAEZ,NACLSD),NUMN0(NAEZ)

C     DOUBLE COMPLEX     GLLH(LMMAXD,NACLSD*LMMAXD,NAEZ),
C     +                   BLAV(NATBLD*LMMAXD,NATBLD*LMMAXD)

      DOUBLE COMPLEX     GLLH((LMAX+1)**2,NACLSD*(LMAX+1)**2,NAEZ),
     &                   BLAV(NATBLD*(LMAX+1)**2,NATBLD*(LMAX+1)**2)


      LMMAXD = (LMAX+1)**2

      DO I1 = 1, NATBLD
        DO I2 =1, NATBLD
C
          IL2B = 0
C
          DO ISRH=1,NUMN0((I-1)*NATBLD+I1)
          IF (INDN0((I-1)*NATBLD+I1,ISRH).EQ.(J-1)*NATBLD+I2)
     +      IL2B = ISRH
          ENDDO
C
C also a pointer array POINTTO(NAEZD,NAEZD) could come up for this
C job
C          IL2B = POINTTO((J-1)*NATBLD+I2,(I-1)*NATBLD+I1)
C
          IF (IL2B.NE.0) THEN
            DO LM1 = 1, LMMAXD
              IL1=LMMAXD*(I1-1)+LM1
              DO LM2 = 1, LMMAXD
                IL2  = LMMAXD*(I2-1)+LM2
                BLAV(IL1,IL2) = BLAV(IL1,IL2) + 
     +          GLLH(LM1,LMMAXD*(IL2B-1)+LM2,(I-1)*NATBLD+I1)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO



C
      END
