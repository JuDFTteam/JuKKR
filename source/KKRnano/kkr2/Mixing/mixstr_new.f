c 13.10.95 ***************************************************************
      SUBROUTINE MIXSTR_NEW(RMSAVQ,RMSAVM,LPOT,LMPOT,
     +                  NSPIN,ITER,RFPI,FPI,
     +                  MIXING,FCM,IRC1,IRMIN1,R,DRDI,VONS,VISP,VINS,
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
      INTEGER ITER,LMPOT,LPOT,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD),R(IRMD),
     +                 VINS((IRMD-IRNSD):IRMD,LMPOT,2),VISP(IRMD,2),
     +                 VONS(IRMD,LMPOT,2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FAC,RMSERM,RMSERQ,VMN,VNM,VNP,VOLDM,VOLDP,VPN
      INTEGER IH,IHP1,IRC1,IRMIN1,J,LM
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
          IH = 1
          IHP1 = IH + 1
        ELSE
          IH = 1
          IHP1 = IH
        END IF

C       ok, that means that for
C       NSPIN=1 (or not 2)  -> IH=1 and IHP1=2
C       NSPIN=2             -> IH=1 and IHP1=1
C       IHP1 means IH plus 1
C       :-P

        RMSERQ = 0.0D0
        RMSERM = 0.0D0
        FAC = 0.5D0/RFPI
c
        DO 10 J = 1,IRC1
          VNP = FAC* (VONS(J,1,IH)+VONS(J,1,IHP1))
          VNM = FAC* (VONS(J,1,IH)-VONS(J,1,IHP1))
          VOLDP = 0.5D0* (VISP(J,IH)+VISP(J,IHP1))
          VOLDM = 0.5D0* (VISP(J,IH)-VISP(J,IHP1))
          RMSERQ = RMSERQ + 2.0D0*REAL(1+MOD(J,2))*R(J)*R(J)*
     +             DRDI(J)* (VNP-VOLDP)* (VNP-VOLDP)
          RMSERM = RMSERM + 2.0D0*REAL(1+MOD(J,2))*R(J)*R(J)*
     +             DRDI(J)* (VNM-VOLDM)* (VNM-VOLDM)
          VPN = VOLDP + MIXING* (VNP-VOLDP)
          VMN = VOLDM + FCM*MIXING* (VNM-VOLDM)
          VONS(J,1,IHP1) = VPN - VMN
          VONS(J,1,IH) = VPN + VMN
   10   CONTINUE
C
        RMSERQ = RMSERQ/ (R(IRC1)**3)
        RMSERM = RMSERM/ (R(IRC1)**3)
        RMSAVQ = RMSAVQ + RMSERQ
        RMSAVM = RMSAVM + RMSERM

        IF (LPOT.GT.0) THEN

          RMSERQ = 0.0D0
          RMSERM = 0.0D0

          DO 30 LM = 2,LMPOT
            DO 20 J = IRMIN1,IRC1
              VNP = 0.5D0* (VONS(J,LM,IH)+VONS(J,LM,IHP1))
              VNM = 0.5D0* (VONS(J,LM,IH)-VONS(J,LM,IHP1))
              VOLDP = 0.5D0* (VINS(J,LM,IH)+VINS(J,LM,IHP1))
              VOLDM = 0.5D0* (VINS(J,LM,IH)-VINS(J,LM,IHP1))
              RMSERQ = RMSERQ + 2.0D0*REAL(1+MOD(J,2))*R(J)*R(J)*
     +                 DRDI(J)* (VNP-VOLDP)* (VNP-VOLDP)
              RMSERM = RMSERM + 2.0D0*REAL(1+MOD(J,2))*R(J)*R(J)*
     +                 DRDI(J)* (VNM-VOLDM)* (VNM-VOLDM)
              VPN = VOLDP + MIXING* (VNP-VOLDP)
              VMN = VOLDM + FCM*MIXING* (VNM-VOLDM)
              VONS(J,LM,IHP1) = VPN - VMN
              VONS(J,LM,IH) = VPN + VMN
   20       CONTINUE
   30     CONTINUE
          RMSERQ = RMSERQ/ (R(IRC1)**3)/FPI
          RMSERM = RMSERM/ (R(IRC1)**3)/FPI
          RMSAVQ = RMSAVQ + RMSERQ
          RMSAVM = RMSAVM + RMSERM

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
