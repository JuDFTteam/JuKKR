C*==calcmvec.f    processed by SPAG 6.05Rc at 18:24 on  6 Apr 2003
      SUBROUTINE CALCMVEC(TXTL,NFILCBWF,SHFTEF,SPLITSS,IEPATH,NEPATH,
     &                    IREL,IPRINT,NT,NL,MEZZ,MEZJ,TAUT,TSST,
     &                    IQAT,NKMQ,NKM,IECURR,NETAB,IGRID,WE,TXTT,FACT,
     &                    MVEVDL0,MVEVIL,BMVEVDL0,BMVEVIL,MVPHI,MVTET,
     &                    MVGAM,QMTET,QMPHI,R2DRDI,JRWS,IMT,AMEMVEC,
     &                    IKMLLIM1,IKMLLIM2,IMKMTAB,NTMAX,NLMAX,NMUEMAX,
     &                    NQMAX,NKMMAX,NMMAX,NMVECMAX,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT COMPLEX*16(A-H,O-Z)
C
C PARAMETER definitions
C
      COMPLEX*16 C0,C1,CI
      PARAMETER (C0=(0.0D0,0.0D0),C1=(1.0D0,0.0D0),CI=(0.0D0,1.0D0))
      REAL*8 PI
      PARAMETER (PI=3.141592653589793238462643D0)
      COMPLEX*16 CPRE
      PARAMETER (CPRE=-C1/PI)
C
C Dummy arguments
C
      COMPLEX*16 WE
      INTEGER IECURR,IEPATH,IGRID,IPRINT,IREL,NEPATH,NETAB,NFILCBWF,NKM,
     &        NKMMAX,NL,NLMAX,NMMAX,NMUEMAX,NMVECMAX,NQMAX,NRMAX,NT,
     &        NTMAX
      REAL*8 SHFTEF,FACT(0:100)
      LOGICAL SPLITSS
      REAL*8 AMEMVEC(NKMMAX,NKMMAX,3,NMVECMAX),MVGAM(NTMAX,NMVECMAX),
     &       MVPHI(NTMAX,NMVECMAX),MVTET(NTMAX,NMVECMAX),QMPHI(NQMAX),
     &       QMTET(NQMAX),R2DRDI(NRMAX,NMMAX)
      COMPLEX*16 BMVEVDL0(NLMAX,NTMAX,3,NMVECMAX),
     &           BMVEVIL(NLMAX,NTMAX,3,NMVECMAX),
     &           MEZJ(NKMMAX,NKMMAX,NTMAX,NMVECMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMVECMAX),
     &           MVEVDL0(NLMAX,NTMAX,3,NMVECMAX),
     &           MVEVIL(NLMAX,NTMAX,3,NMVECMAX),
     &           TAUT(NKMMAX,NKMMAX,NTMAX),TSST(NKMMAX,NKMMAX,NTMAX)
      INTEGER IKMLLIM1(NKMMAX),IKMLLIM2(NKMMAX),IMKMTAB(NKMMAX),
     &        IMT(NTMAX),IQAT(NQMAX,NTMAX),JRWS(NMMAX),NKMQ(NQMAX)
      CHARACTER*1 TXTL(0:NLMAX)
      CHARACTER*4 TXTT(NTMAX)
C
C Local variables
C
CF77--------------------------------------------------------------------
Cccc      COMPLEX*16 JF(NRMAX,2,NKMMAX),JG(NRMAX,2,NKMMAX),
Cccc     &           ZF(NRMAX,2,NKMMAX),ZG(NRMAX,2,NKMMAX),
CF77--------------------------------------------------------------------
CF90--------------------------------------------------------------------
      COMPLEX*16 JF(:,:,:),JG(:,:,:),ZF(:,:,:),ZG(:,:,:)
      ALLOCATABLE JF,JG,ZF,ZG
CF90--------------------------------------------------------------------
      COMPLEX*16 AMIN,APLS,BMVEVD(NTMAX,3,NMVECMAX),
     &           BMVEVDL(NLMAX,3,NMVECMAX),
     &           BMVEVDM(NLMAX,NMUEMAX,3,NMVECMAX),
     &           BMVEVI(NTMAX,3,NMVECMAX),CS,CSUM,CWGT,DROT4(4,4),
     &           MEIRR(NKMMAX,NKMMAX,3,NMVECMAX),
     &           MEREG(NKMMAX,NKMMAX,3,NMVECMAX),
     &           MVEVD(NTMAX,3,NMVECMAX),MVEVDL(NLMAX,3,NMVECMAX),
     &           MVEVDM(NLMAX,NMUEMAX,3,NMVECMAX),
     &           MVEVI(NTMAX,3,NMVECMAX),USC(3,3),W1(NKMMAX,NKMMAX),
     &           W2(NKMMAX,NKMMAX),W3(NKMMAX,NKMMAX),W3X3(3,3),
     &           ZFJF(2,2),ZFZF(2,2),ZGJG(2,2),ZGZG(2,2)
      LOGICAL CHECK
      INTEGER I,I0,IA_ERR,IKM,IKM1,IKM2,IKMCB(2,NKMMAX),IKMT1,IKMT2,
     &        IL,IM,IMKM1,IMKM2,IMV,IMVEC,IPOL,IQ,IT,ITI,J,J1,J2,JTOP,K,
     &        K1,K2,KAPCB,L,LI,LMAX,M,MUE,MUETOP,N,NMVEC,NOSWF,NPOL,
     &        NSOL,NSOLCB(NKMMAX)
      REAL*8 MJ,MROT(3,3),MV,MVGLO(3,NMVECMAX),MVGLOL(NLMAX,3,NMVECMAX),
     &       MVX,MVXY,MVY,MVZ,SUM,W
      CHARACTER*3 STR3
C
C index 3:  ipol= 1,2,3  ==  (+),(-),(z)
C
C
      CHECK = .TRUE.
      CHECK = .FALSE.
      NPOL = 3
      NMVEC = MIN(4,NMVECMAX)
C
      IF ( IECURR.EQ.1 .AND. IEPATH.EQ.1 ) THEN
C
         DO IMV = 1,NMVEC
            DO IPOL = 1,NPOL
               DO IT = 1,NT
                  DO IL = 1,NL
                     MVEVDL0(IL,IT,IPOL,IMV) = C0
                     MVEVIL(IL,IT,IPOL,IMV) = C0
                     BMVEVDL0(IL,IT,IPOL,IMV) = C0
                     BMVEVIL(IL,IT,IPOL,IMV) = C0
                  END DO
               END DO
            END DO
         END DO
C
      END IF
C
      NOSWF = NT*NKM
      NMVEC = 3
C
CF90--------------------------------------------------------------------
      ALLOCATE (JF(NRMAX,2,NKMMAX),JG(NRMAX,2,NKMMAX),STAT=IT)
      IF ( IT.NE.0 ) STOP '      < CALCMVEC > : allocate JF/JG '
      ALLOCATE (ZF(NRMAX,2,NKMMAX),ZG(NRMAX,2,NKMMAX),STAT=IT)
      IF ( IT.NE.0 ) STOP '      < CALCMVEC > : allocate ZF/ZG '
CF90--------------------------------------------------------------------
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
C
         M = NKMMAX
         N = NKMQ(IQAT(1,IT))
         IM = IMT(IT)
         IM = 1
         JTOP = JRWS(IM)
C
         LMAX = NL - 1
C
         CALL CINIT(NKMMAX*NKMMAX*3*NMVECMAX,MEREG)
         CALL CINIT(NKMMAX*NKMMAX*3*NMVECMAX,MEIRR)
C
C=======================================================================
C                      calculate matrix elements
C=======================================================================
C
         DO IKM = 1,N
            READ (NFILCBWF,REC=IKM+(IT-1)*NKM) ITI,LI,MJ,NSOL,STR3,
     &            (KAPCB,IKMCB(K,IKM),(ZG(I,K,IKM),ZF(I,K,IKM),I=1,JTOP)
     &            ,K=1,NSOL)
            IF ( IT.NE.ITI ) STOP ' IT(INI) <> IT  in <MENABIRR>'
            IF ( STR3.NE.'REG' ) STOP 'WFT(INI) <> REG in <MENABIRR>'
C
            READ (NFILCBWF,REC=IKM+(IT-1)*NKM+NOSWF) ITI,LI,MJ,NSOL,
     &            STR3,
     &            (KAPCB,IKMCB(K,IKM),(JG(I,K,IKM),JF(I,K,IKM),I=1,JTOP)
     &            ,K=1,NSOL)
            NSOLCB(IKM) = NSOL
            IF ( IT.NE.ITI ) STOP ' IT(INI) <> IT  in <MENABIRR>'
            IF ( STR3.NE.'IRR' ) STOP 'WFT(INI) <> IRR in <MENABIRR>'
            IF ( IPRINT.GT.3 ) WRITE (1337,*) ITI,LI,MJ,NSOL,STR3,KAPCB
         END DO
C

         DO IKMT2 = 1,N
C
            DO IKMT1 = IKMLLIM1(IKMT2),IKMLLIM2(IKMT2)
C
               DO IMVEC = 1,NMVEC
                  DO K2 = 1,NSOLCB(IKMT2)
                     J2 = IKMCB(K2,IKMT2)
                     DO K1 = 1,NSOLCB(IKMT1)
                        J1 = IKMCB(K1,IKMT1)
                        DO IPOL = 1,NPOL
                           IF ( ABS(AMEMVEC(J1,J2,IPOL,IMVEC)).GT.1E-8 )
     &                          GOTO 10
                        END DO
                     END DO
                  END DO
               END DO
C -------------------------------------- all angular matrix elements = 0
               GOTO 20
C ---------------------------------- non-0 angular matrix elements found
C ------------------------------------- calculate radial matrix elements
C
 10            CONTINUE
               CALL CINTABR(ZG(1,1,IKMT1),ZG(1,1,IKMT2),ZGZG,
     &                      ZF(1,1,IKMT1),ZF(1,1,IKMT2),ZFZF,
     &                      R2DRDI(1,IM),NSOLCB(IKMT1),NSOLCB(IKMT2),
     &                      JTOP,NRMAX)
C
               CALL CINTABR(ZG(1,1,IKMT1),JG(1,1,IKMT2),ZGJG,
     &                      ZF(1,1,IKMT1),JF(1,1,IKMT2),ZFJF,
     &                      R2DRDI(1,IM),NSOLCB(IKMT1),NSOLCB(IKMT2),
     &                      JTOP,NRMAX)
C
C -------------------------------------- calculate total matrix elements
C
               DO K2 = 1,NSOLCB(IKMT2)
                  IKM2 = IKMCB(K2,IKMT2)
                  IMKM2 = IMKMTAB(IKM2)
C
                  DO K1 = 1,NSOLCB(IKMT1)
                     IKM1 = IKMCB(K1,IKMT1)
                     IMKM1 = IMKMTAB(IKM1)
C
                     DO IMV = 1,NMVEC
                        DO IPOL = 1,NPOL
                           MEREG(IKMT1,IKMT2,IPOL,IMV)
     &                        = MEREG(IKMT1,IKMT2,IPOL,IMV)
     &                        + AMEMVEC(IKM1,IKM2,IPOL,IMV)*ZGZG(K1,K2)
                        END DO
                     END DO
C
                  END DO
C
                  DO IMV = 1,NMVEC
                     DO IPOL = 1,NPOL
                        MEREG(IKMT2,IKMT2,IPOL,IMV)
     &                     = MEREG(IKMT2,IKMT2,IPOL,IMV)
     &                     - AMEMVEC(IMKM2,IMKM2,IPOL,IMV)*ZFZF(K2,K2)
                     END DO
                  END DO
C
               END DO
C
               IF ( IKMT1.EQ.IKMT2 ) THEN
                  DO K2 = 1,NSOLCB(IKMT2)
                     IKM2 = IKMCB(K2,IKMT2)
                     IMKM2 = IMKMTAB(IKM2)
                     DO K1 = 1,NSOLCB(IKMT1)
                        IKM1 = IKMCB(K1,IKMT1)
                        IMKM1 = IMKMTAB(IKM1)
C
                        DO IMV = 1,NMVEC
                           DO IPOL = 1,NPOL
                              MEIRR(IKMT1,IKMT2,IPOL,IMV)
     &                           = MEIRR(IKMT1,IKMT2,IPOL,IMV)
     &                           + AMEMVEC(IKM1,IKM2,IPOL,IMV)
     &                           *ZGJG(K1,K2)
                           END DO
                        END DO
C
                     END DO
C
                     DO IMV = 1,NMVEC
                        DO IPOL = 1,NPOL
                           MEIRR(IKMT2,IKMT2,IPOL,IMV)
     &                        = MEIRR(IKMT2,IKMT2,IPOL,IMV)
     &                        - AMEMVEC(IMKM2,IMKM2,IPOL,IMV)
     &                        *ZFJF(K2,K2)
                        END DO
                     END DO
C
                  END DO
               END IF
C
 20         END DO
C
         END DO
C
         IF ( CHECK ) THEN
            DO I = 1,NKM
               DO J = IKMLLIM1(I),IKMLLIM2(I)
                  SUM = ABS(MEZZ(I,J,IT,1)) + ABS(MEZJ(I,J,IT,1))
                  SUM = SUM + ABS(MEREG(I,J,3,1)) + ABS(MEIRR(I,J,3,1))
                  IF ( SUM.GT.1D-8 ) THEN
                     WRITE (1337,*) ' spin '
                     WRITE (1337,'(2i3,2e17.8,2x,2e17.8)') I,J,
     &                      MEZZ(I,J,IT,2),MEZJ(I,J,IT,2)
                     WRITE (1337,'(6x,2e17.8,2x,2e17.8)') MEREG(I,J,3,1)
     &                     ,MEIRR(I,J,3,1),
     &                      (MEZZ(I,J,IT,2)-MEREG(I,J,3,1)),
     &                      (MEZJ(I,J,IT,2)-MEIRR(I,J,3,1))
                     WRITE (1337,*) ' orb '
                     WRITE (1337,'(2i3,2e17.8,2x,2e17.8)') I,J,
     &                      MEZZ(I,J,IT,3),MEZJ(I,J,IT,3)
                     WRITE (1337,'(6x,2e17.8,2x,2e17.8)') MEREG(I,J,3,2)
     &                     ,MEIRR(I,J,3,2),
     &                      (MEZZ(I,J,IT,3)-MEREG(I,J,3,2)),
     &                      (MEZJ(I,J,IT,3)-MEIRR(I,J,3,2))
                  END IF
               END DO
            END DO
         END IF
C
C=======================================================================
         DO IMV = 1,NMVEC
            DO IPOL = 1,NPOL
C
               MVEVD(IT,IPOL,IMV) = 0.0D0
               MVEVI(IT,IPOL,IMV) = 0.0D0
C
               IF ( .NOT.SPLITSS ) THEN
                  BMVEVD(IT,IPOL,IMV) = 0.0D0
                  BMVEVI(IT,IPOL,IMV) = 0.0D0
               END IF
C
               CALL ZGEMM('N','N',N,N,N,CPRE,MEREG(1,1,IPOL,IMV),M,
     &                    TAUT(1,1,IT),M,C0,W1,M)
               CWGT = -1D0
               DO J = 1,N
                  CALL ZAXPY(N,-CPRE,MEIRR(1,J,IPOL,IMV),1,W1(1,J),1)
                  CALL ZCOPY(N,TAUT(1,J,IT),1,W2(1,J),1)
                  CALL ZAXPY(N,CWGT,TSST(1,J,IT),1,W2(1,J),1)
               END DO
               CALL ZGEMM('N','N',N,N,N,CPRE,MEREG(1,1,IPOL,IMV),M,W2,M,
     &                    C0,W3,M)
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
               DO L = 0,LMAX
                  IL = L + 1
C
                  IF ( IREL.LE.1 ) THEN
                     I0 = L*L
                     DO MUE = 1,2*L + 1
                        MVEVDM(IL,MUE,IPOL,IMV) = W1(I0+MUE,I0+MUE)
                        BMVEVDM(IL,MUE,IPOL,IMV) = W3(I0+MUE,I0+MUE)
                     END DO
                  ELSE
                     I0 = 2*L*L + L + L
                     MUETOP = 2*L + 2
                     DO MUE = 1,MUETOP
                        MVEVDM(IL,MUE,IPOL,IMV) = W1(I0+MUE,I0+MUE)
                        BMVEVDM(IL,MUE,IPOL,IMV) = W3(I0+MUE,I0+MUE)
                     END DO
                     I0 = 2*(L-1)*L + L + L
                     DO MUE = 2,MUETOP - 1
                        MVEVDM(IL,MUE,IPOL,IMV)
     &                     = MVEVDM(IL,MUE,IPOL,IMV)
     &                     + W1(I0+MUE-1,I0+MUE-1)
                        BMVEVDM(IL,MUE,IPOL,IMV)
     &                     = BMVEVDM(IL,MUE,IPOL,IMV)
     &                     + W3(I0+MUE-1,I0+MUE-1)
                     END DO
                  END IF
C
                  MVEVDL(IL,IPOL,IMV) = 0.0D0
C
                  IF ( .NOT.SPLITSS ) BMVEVDL(IL,IPOL,IMV) = 0.0D0
C
                  IF ( IREL.GT.1 ) THEN
                     MUETOP = 2*L + 2
                  ELSE
                     MUETOP = 2*L + 1
                  END IF
C
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
                  DO MUE = 1,MUETOP
                     MVEVDL(IL,IPOL,IMV) = MVEVDL(IL,IPOL,IMV)
     &                  + MVEVDM(IL,MUE,IPOL,IMV)
C
                     IF ( .NOT.SPLITSS ) BMVEVDL(IL,IPOL,IMV)
     &                    = BMVEVDL(IL,IPOL,IMV)
     &                    + BMVEVDM(IL,MUE,IPOL,IMV)
                  END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
                  MVEVIL(IL,IT,IPOL,IMV) = MVEVIL(IL,IT,IPOL,IMV)
     &               + WE*MVEVDL(IL,IPOL,IMV)
                  MVEVDL0(IL,IT,IPOL,IMV) = MVEVDL(IL,IPOL,IMV)
C
                  IF ( .NOT.SPLITSS ) THEN
                     BMVEVIL(IL,IT,IPOL,IMV) = BMVEVIL(IL,IT,IPOL,IMV)
     &                  + WE*BMVEVDL(IL,IPOL,IMV)
                     BMVEVDL0(IL,IT,IPOL,IMV) = BMVEVDL(IL,IPOL,IMV)
                  END IF
C
                  IF ( IGRID.NE.4 .OR. IECURR.LE.NETAB ) THEN
C
                     MVEVD(IT,IPOL,IMV) = MVEVD(IT,IPOL,IMV)
     &                  + MVEVDL(IL,IPOL,IMV)
                     MVEVI(IT,IPOL,IMV) = MVEVI(IT,IPOL,IMV)
     &                  + MVEVIL(IL,IT,IPOL,IMV)
C
                     IF ( .NOT.SPLITSS ) THEN
                        BMVEVD(IT,IPOL,IMV) = BMVEVD(IT,IPOL,IMV)
     &                     + BMVEVDL(IL,IPOL,IMV)
                        BMVEVI(IT,IPOL,IMV) = BMVEVI(IT,IPOL,IMV)
     &                     + BMVEVIL(IL,IT,IPOL,IMV)
                     END IF
                  END IF
C
               END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
            END DO
         END DO
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
CF90--------------------------------------------------------------------
      DEALLOCATE (JF,JG,ZF,ZG,STAT=IT)
      IF ( IT.NE.0 ) STOP '      < CALCMVEC > : deallocate JF/JG/ZF/ZG'
CF90--------------------------------------------------------------------
      IF ( SPLITSS .AND. ((IEPATH.EQ.1) .AND. (IECURR.EQ.NETAB)) ) THEN
         DO IMV = 1,NMVEC
            DO IPOL = 1,NPOL
               DO IT = 1,NT
                  DO IL = 1,NL
                     BMVEVDL0(IL,IT,IPOL,IMV) = MVEVDL0(IL,IT,IPOL,IMV)
                  END DO
               END DO
            END DO
         END DO
      END IF
C
C=======================================================================
      IF ( (IGRID.GE.6) .OR. (IECURR.NE.NETAB) .OR. (IEPATH.NE.NEPATH) )
     &     RETURN
C
C     this part of the original Munich subroutine has been moved to
C     < mvecglobal > 
C     main2 --> tbkkr2 --> mvecglobal -- see makefile2    
C
      END
