C*==decipotbas.f    processed by SPAG 6.05Rc at 14:16 on  3 Dec 2004
      SUBROUTINE DECIPOTBAS(IHOST,IQOFF,ITOFF,NQ,NT,RBASIS,QMTET,QMPHI,
     &                      NOQ,KAOEZ,ZAT,IQAT,CONC,IRWS,IPAN,IRCUT,RR,
     &                      DRDI,VISP,NSPIN,KREL,SOLVER,SOCSCL,CSCL,
     &                      VTREL,BTREL,IRMD,IPAND,NEMBD1,NTMAX,NSPIND,
     &                      LMAXD)
C **********************************************************************
C *                                                                    *
C * reads in the potential data for the host atoms from the potential  *
C * file 'decimate.pot'                                                *
C *                                        v.popescu - munich, Dec 04  *
C *                                                                    *
C * Note: so far, only  SPHERICAL case implemented                     *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Arguments
      INTEGER IHOST,IPAND,IQOFF,IRMD,ITOFF,KREL,LMAXD,NEMBD1,NQ,NSPIN,
     &        NSPIND,NT,NTMAX
      CHARACTER*10 SOLVER
      DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL),NTMAX),CONC(NTMAX),
     &                 CSCL(KREL*LMAXD+1,KREL*NTMAX+(1-KREL)),
     &                 DRDI(IRMD,NTMAX),QMPHI(NEMBD1),QMTET(NEMBD1),
     &                 RBASIS(3,NEMBD1),RR(IRMD,NTMAX),
     &                 SOCSCL(KREL*LMAXD+1,KREL*NTMAX+(1-KREL)),
     &                 VISP(IRMD,NTMAX*NSPIND),
     &                 VTREL(IRMD*KREL+(1-KREL),NTMAX),ZAT(NTMAX)
      INTEGER IPAN(NTMAX),IQAT(NEMBD1,NTMAX),IRCUT(0:IPAND,NTMAX),
     &        IRWS(NTMAX),KAOEZ(NEMBD1,NEMBD1),NOQ(NEMBD1)
C     ..
C     .. Locals
      INTEGER I,IH,IHF,IL,IPOT1,IPOT2
      INTEGER NINT
      DOUBLE PRECISION RMT(NTMAX),RWS(NTMAX)
      CHARACTER*3 TXTT(NT)
C ......................................................................
C --> read basis
C
      DO IH = 1,NQ
         IHF = IH + IQOFF
         READ (36+IHOST,99001) IL,(RBASIS(I,IHF),I=1,3)
         IF ( IH.NE.IL ) STOP ' Inconsistent data '
         WRITE (1337,99001) IL,(RBASIS(I,IHF),I=1,3)
      END DO
      READ (36+IHOST,*)
      WRITE (1337,99003)
      DO IH = 1,NQ
         IHF = IH + IQOFF
         READ (36+IHOST,FMT=99002) IL,QMTET(IHF),QMPHI(IHF),NOQ(IHF),
     &                             (KAOEZ(I,IHF),I=1,NOQ(IHF))
         IF ( IH.NE.IL ) STOP ' Inconsistent data '
         WRITE (1337,99004) IH,QMTET(IHF),QMPHI(IHF),NOQ(IHF),
     &                   (KAOEZ(I,IHF),I=1,NOQ(IHF))
      END DO
      WRITE (1337,99005)
      IF ( KREL.EQ.1 ) READ (36+IHOST,'(7X,A10)') SOLVER
      READ (36+IHOST,*)
C
C --> read atoms
C
      DO IH = 1,NT
         IHF = IH + ITOFF
         READ (36+IHOST,99006) IL,ZAT(IHF),IQAT(1,IHF),CONC(IHF),
     &                         IRWS(IHF),IPAN(IHF),
     &                         (IRCUT(I,IHF),I=0,IPAN(IHF))
         IF ( IH.NE.IL ) STOP ' Inconsistent data '
C
         IF ( KREL.EQ.1 ) THEN
            READ (36+IHOST,99007) SOCSCL(1,IHF),CSCL(1,IHF)
            DO IL = 2,LMAXD + 1
               SOCSCL(IL,IHF) = SOCSCL(1,IHF)
               CSCL(IL,IHF) = CSCL(1,IHF)
            END DO
         END IF
C
      END DO
      DO IH = 1,NT
         IHF = IH + ITOFF
         IPOT1 = (IHF-1)*NSPIN + 1
         IPOT2 = IPOT1 + 1
         READ (36+IHOST,*)
         READ (36+IHOST,99008) IL,TXTT(IH),RMT(IHF),RWS(IHF)
         IF ( IH.NE.IL ) STOP ' Inconsistent data '
         WRITE (1337,99009) IH,TXTT(IH),NINT(ZAT(IHF)),CONC(IHF),
     &                   IRWS(IHF),RWS(IHF)
         IF ( KREL.EQ.0 ) THEN
            READ (36+IHOST,*)
            DO I = 1,IRWS(IHF)
               READ (36+IHOST,99010) RR(I,IHF),DRDI(I,IHF),
     &                               (VISP(I,IL),IL=IPOT1,IPOT2)
            END DO
         ELSE
            READ (36+IHOST,'(7X,I3)') IL
            IRWS(IHF) = IRWS(IHF) - IL
            IRCUT(IPAN(IHF),IHF) = IRWS(IHF)
            READ (36+IHOST,*)
            DO I = 1,IRWS(IHF)
               READ (36+IHOST,99010) RR(I,IHF),DRDI(I,IHF),VTREL(I,IHF),
     &                               BTREL(I,IHF)
            END DO
         END IF
      END DO
C
99001 FORMAT (9X,I3,3F12.8)
99002 FORMAT (7X,I3,2(7X,F9.4),7X,I3,7X,8I3)
99003 FORMAT (9X,39('-'),/,9X,'   THETA   ','   PHI   ','OCC',' IT')
99004 FORMAT (9X,I3,2(F9.4),I3,8I3)
99005 FORMAT (9X,39('-'),/,10X,'ATOMS',/,15X,'Z   CONC  IWS    RWS')
99006 FORMAT (7X,I3,7X,F4.0,7X,I3,7X,F7.4,/,17X,I4,7X,I3,7X,6I4)
99007 FORMAT (17X,F10.6,7X,D13.6)
99008 FORMAT (7X,I3,1X,A3,2(/,7X,F12.8))
99009 FORMAT (9X,I3,1X,A3,I3,F7.4,I4,F10.6)
99010 FORMAT (1P,4D20.12)
      END
