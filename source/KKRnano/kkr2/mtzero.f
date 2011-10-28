c ************************************************************************
      SUBROUTINE MTZERO(LMPOT,IH,NSPIN,VONS,Z,R,DRDI,IMT,IRCUT,
     +                IPAN,ICELL,LMSP,IFUNM,THETAS,IRWS,VAV0,VOL0,
C                     new input parameters after inc.p removal
     &                irmd, irid, nfund, ipand)
      implicit none
c ************************************************************************
c
c     determine muffin tin zero and shift potential to muffin tin zero
c
c     for spin polarized calculations muffin tin zero is related to
c         the average of the 2 spins 
c
c                                            may,2000 (new version)
c
c-----------------------------------------------------------------------
C     INTEGER LMPOTD
C     PARAMETER (LMPOTD= (LPOTD+1)**2)

      INTEGER irmd
      INTEGER irid
      INTEGER nfund
      INTEGER ipand

C     .. Scalar Arguments ..
      INTEGER LMPOT,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION 
     +     DRDI(IRMD,*),
     +     R(IRMD,*),
     +     THETAS(IRID,NFUND,*),
     +     VONS(IRMD,LMPOT,2),
     +     Z(*)
      INTEGER 
     +     IFUNM(*),IMT(*),IPAN(*),IRCUT(0:IPAND,*),IRWS(*),
     +     LMSP(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FPI,RFPI,VAV0,VOL0,ZZOR
      INTEGER ICELL,IFUN,IH,IMT1,IPAN1,IPOT,IR,IRC1,IRH,IS,LM,II
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION V1(IRMD),V2(IRMD),VAV1(2),VOL1(2)
C     ..
C     .. External Subroutines ..
      LOGICAL TEST,OPT
      EXTERNAL SIMP3,SIMPK,TEST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
C     ..
      FPI = 16.0D0*ATAN(1.0D0)
      RFPI = SQRT(FPI)
C
        DO 10 IR = 1,IRMD
          V1(IR) = 0.0D0
          V2(IR) = 0.0D0
   10   CONTINUE
C
C
C
        DO IS = 1,NSPIN
C
C
C
           IPOT = NSPIN* (IH-1) + IS
           IPAN1 = IPAN(IH)
           IMT1 = IMT(IH)
           
           IF (IPAN1.EQ.1) THEN
c     
c---  >     muffin tin or atomic sphere calculation
c     
              IRC1 = IRWS(IH)
              DO 20 IR = IMT1,IRC1
                 V2(IR) = FPI*R(IR,IH)**2
                 ZZOR = 2.0D0*Z(IH)/R(IR,IH)
                 V1(IR) = (VONS(IR,1,IS)/RFPI-ZZOR)*V2(IR)
 20           CONTINUE
              
              CALL SIMP3(V1,VAV1(IS),IMT1,IRC1,DRDI(1,IH))
              CALL SIMP3(V2,VOL1(IS),IMT1,IRC1,DRDI(1,IH))
              
           ELSE                 ! (IPAN1.EQ.1)
c     
c---  >     full potential calculation
c     
              IRC1 = IRCUT(IPAN1,IH)
              IMT1 = IMT(IH) 
              DO 30 IR = IMT1 + 1,IRC1
                 V2(IR) = R(IR,IH)**2*THETAS(IR-IMT1,1,ICELL)*RFPI
                 ZZOR = 2.0D0*Z(IH)/R(IR,IH)
                 V1(IR) = (VONS(IR,1,IS)/RFPI-ZZOR)*V2(IR)
 30           CONTINUE
              DO 50 LM = 2,LMPOT
                 IF (LMSP(LM).GT.0) THEN
                    IFUN = IFUNM(LM)
                    
                    DO 40 IR = IMT1 + 1,IRC1
                       IRH = IR - IMT1
                       V1(IR) = V1(IR) + R(IR,IH)**2*VONS(IR,LM,IS)*
     +                      THETAS(IRH,IFUN,ICELL)
 40                 CONTINUE
                    
                 END IF
                 
 50           CONTINUE          ! LM = 2,LMPOT
              
              CALL SIMPK(V1,VAV1(IS),IPAN1,IRCUT(0,IH),DRDI(1,IH))
              CALL SIMPK(V2,VOL1(IS),IPAN1,IRCUT(0,IH),DRDI(1,IH))
              
           END IF               ! (IPAN1.EQ.1)
C
C
C
        END DO                  ! SPIN LOOP
C
C
C
        IF (NSPIN.EQ.1) THEN
        VAV1(2) = VAV1(1)
        VOL1(2) = VOL1(1)
        END IF
        
        VAV0 = VAV0 + (VAV1(1)+VAV1(2))/2.D0
        VOL0 = VOL0 + (VOL1(1)+VOL1(2))/2.D0
      
      RETURN
      
      END
