c 13.10.95 ***************************************************************
      SUBROUTINE MIXSTR(RMSAVQ,RMSAVM,LPOT,LMPOT,
     +                  I1,NSPIN,ITER,RFPI,FPI,IPF,
     +                  MIXING,FCM,IRC,IRMIN,R,DRDI,VONS,VISP,VINS,
C                       new parameters after inc.p removal
     &                  naez, irmd, irnsd)
c ************************************************************************
      IMPLICIT NONE

      INTEGER naez
      INTEGER irmd
      INTEGER irnsd

C     INTEGER LMPOTD,IRMIND
C     PARAMETER (LMPOTD= (LPOTD+1)**2,
C    +          IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION FCM,FPI,MIXING,RFPI,RMSAVM,RMSAVQ
      INTEGER IPF,ITER,LMPOT,LPOT,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD,NAEZ),R(IRMD,NAEZ),
     +                 VINS((IRMD-IRNSD):IRMD,LMPOT,2),VISP(IRMD,2),
     +                 VONS(IRMD,LMPOT,2)
      INTEGER IRC(NAEZ),IRMIN(NAEZ)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FAC,RMSERM,RMSERQ,VMN,VNM,VNP,VOLDM,VOLDP,VPN
      INTEGER IH,IHP1,IRC1,IRMIN1,J,LM,I1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD,REAL,SQRT
C     ..
c
c---> final construction of the potentials
c     attention : the spherical averaged potential is the lm=1
c                     component of vons times sqrt(4 pi).
c
c     first mixing scheme : straight mixing
c---> determination of the root mean sqare error
c

        IF (NSPIN.EQ.2) THEN
c         IH = 2*I1 - 1
c         IHP1 = IH + 1
          IH = 1
          IHP1 = IH + 1
        ELSE
c         IH = I1
c         IHP1 = IH
          IH = 1
          IHP1 = IH
        END IF

C       ok, that means that for
C       NSPIN=1 (or not 2)  -> IH=1 and IHP1=2
C       NSPIN=2             -> IH=1 and IHP1=1
C       IHP1 means IH plus 1
C       :-P

        IRC1 = IRC(I1)
        RMSERQ = 0.0D0
        RMSERM = 0.0D0
        FAC = 0.5D0/RFPI
c
        DO 10 J = 1,IRC1
          VNP = FAC* (VONS(J,1,IH)+VONS(J,1,IHP1))
          VNM = FAC* (VONS(J,1,IH)-VONS(J,1,IHP1))
          VOLDP = 0.5D0* (VISP(J,IH)+VISP(J,IHP1))
          VOLDM = 0.5D0* (VISP(J,IH)-VISP(J,IHP1))
          RMSERQ = RMSERQ + 2.0D0*REAL(1+MOD(J,2))*R(J,I1)*R(J,I1)*
     +             DRDI(J,I1)* (VNP-VOLDP)* (VNP-VOLDP)
          RMSERM = RMSERM + 2.0D0*REAL(1+MOD(J,2))*R(J,I1)*R(J,I1)*
     +             DRDI(J,I1)* (VNM-VOLDM)* (VNM-VOLDM)
          VPN = VOLDP + MIXING* (VNP-VOLDP)
          VMN = VOLDM + FCM*MIXING* (VNM-VOLDM)
          VONS(J,1,IHP1) = VPN - VMN
          VONS(J,1,IH) = VPN + VMN
   10   CONTINUE
C
        RMSERQ = RMSERQ/ (R(IRC1,I1)**3)
        RMSERM = RMSERM/ (R(IRC1,I1)**3)
        RMSAVQ = RMSAVQ + RMSERQ
        RMSAVM = RMSAVM + RMSERM

        IF (NSPIN.EQ.2) THEN
!          WRITE (IPF,FMT=9000) I1,SQRT(RMSERQ),SQRT(RMSERM)
          WRITE(IPF, REC=(ITER-1)*2*NAEZ + 2*I1-1 ) RMSERQ,RMSERM
        ELSE
 !         WRITE (IPF,FMT=9020) I1,SQRT(RMSERQ)
          WRITE(IPF, REC=(ITER-1)*2*NAEZ + 2*I1-1 ) RMSERQ
        END IF

        IF (LPOT.GT.0) THEN

          RMSERQ = 0.0D0
          RMSERM = 0.0D0
          IRMIN1 = IRMIN(I1)
          DO 30 LM = 2,LMPOT
            DO 20 J = IRMIN1,IRC1
              VNP = 0.5D0* (VONS(J,LM,IH)+VONS(J,LM,IHP1))
              VNM = 0.5D0* (VONS(J,LM,IH)-VONS(J,LM,IHP1))
              VOLDP = 0.5D0* (VINS(J,LM,IH)+VINS(J,LM,IHP1))
              VOLDM = 0.5D0* (VINS(J,LM,IH)-VINS(J,LM,IHP1))
              RMSERQ = RMSERQ + 2.0D0*REAL(1+MOD(J,2))*R(J,I1)*R(J,I1)*
     +                 DRDI(J,I1)* (VNP-VOLDP)* (VNP-VOLDP)
              RMSERM = RMSERM + 2.0D0*REAL(1+MOD(J,2))*R(J,I1)*R(J,I1)*
     +                 DRDI(J,I1)* (VNM-VOLDM)* (VNM-VOLDM)
              VPN = VOLDP + MIXING* (VNP-VOLDP)
              VMN = VOLDM + FCM*MIXING* (VNM-VOLDM)
              VONS(J,LM,IHP1) = VPN - VMN
              VONS(J,LM,IH) = VPN + VMN
   20       CONTINUE
   30     CONTINUE
          RMSERQ = RMSERQ/ (R(IRC1,I1)**3)/FPI
          RMSERM = RMSERM/ (R(IRC1,I1)**3)/FPI
          RMSAVQ = RMSAVQ + RMSERQ
          RMSAVM = RMSAVM + RMSERM

          IF (NSPIN.EQ.2) THEN
!            WRITE (IPF,FMT=9010) I1,SQRT(RMSERQ),SQRT(RMSERM)
            WRITE(IPF, REC=(ITER-1)*2*NAEZ + 2*I1) RMSERQ,RMSERM
          ELSE
!            WRITE (IPF,FMT=9030) I1,SQRT(RMSERQ)
            WRITE(IPF, REC=(ITER-1)*2*NAEZ + 2*I1) RMSERQ
          END IF

        END IF

 9000 FORMAT (5x,' rms-error for atom',i3,1x,':','v+ + v- = ',1p,d11.4,
     +       2x,',',2x,'v+ - v- = ',1p,d11.4)
 9010 FORMAT (5x,' rms-error non spherical contribution for atom ',i3,
     +       1x,':','v+ + v- = ',1p,d11.4,02x,',',2x,'v+ - v- = ',1p,
     +       d11.4)
 9020 FORMAT (5x,' rms-error for atom',i3,1x,':','v+ + v- = ',1p,d11.4)
 9030 FORMAT (5x,' rms-error non spherical contribution for atom ',i3,
     +       1x,':','v+ + v- = ',1p,d11.4)
      END
