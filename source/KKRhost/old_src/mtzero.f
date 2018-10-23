c ************************************************************************
      SUBROUTINE MTZERO(LMPOT,NATYP,CONC,NSPIN,V,VBC,Z,R,DRDI,IMT,IRCUT,
     +                  IPAN,NTCELL,LMSP,IFUNM,THETAS,IRWS,ESHIFT,
     +                  ISHIFT,NSHELL,LSURF)
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
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ESHIFT,VBC(*)
      INTEGER ISHIFT,LMPOT,NATYP,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION 
     +     DRDI(IRMD,*),CONC(NATYPD),
     +     R(IRMD,*),
     +     THETAS(IRID,NFUND,*),
     +     V(IRMD,LMPOTD,*),
     +     Z(*)
      INTEGER 
     +     IFUNM(NATYPD,*),IMT(*),IPAN(*),IRCUT(0:IPAND,*),IRWS(*),
     +     LMSP(NATYPD,*),
     +     NTCELL(*),NSHELL(0:NSHELD)
      LOGICAl LSURF
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FPI,RFPI,VAV0,VOL0,ZZOR
      INTEGER ICELL,IFUN,IH,IMT1,IPAN1,IPOT,IR,IRC1,IRH,IS,LM
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

      VAV0 = 0.0D0
      VOL0 = 0.0D0
      VAV1(1) = 0.D0
      VAV1(2) = 0.D0
      VOL1(1) = 0.D0
      VOL1(2) = 0.D0
      DO 60 IH = 1,NATYP

        DO 10 IR = 1,IRMD
          V1(IR) = 0.0D0
          V2(IR) = 0.0D0
   10   CONTINUE
        DO IS = 1,NSPIN
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
                 V1(IR) = (V(IR,1,IPOT)/RFPI-ZZOR)*V2(IR)
 20           CONTINUE
              
              CALL SIMP3(V1,VAV1(IS),IMT1,IRC1,DRDI(1,IH))
              CALL SIMP3(V2,VOL1(IS),IMT1,IRC1,DRDI(1,IH))
              
           ELSE                 ! (IPAN1.EQ.1)
c     
c---  >     full potential calculation
c     
              IRC1 = IRCUT(IPAN1,IH)
              ICELL = NTCELL(IH)
              IMT1 = IMT(IH) 
              DO 30 IR = IMT1 + 1,IRC1
                 V2(IR) = R(IR,IH)**2*THETAS(IR-IMT1,1,ICELL)*RFPI
                 ZZOR = 2.0D0*Z(IH)/R(IR,IH)
                 V1(IR) = (V(IR,1,IPOT)/RFPI-ZZOR)*V2(IR)
 30           CONTINUE
              DO 50 LM = 2,LMPOT
                 IF (LMSP(ICELL,LM).GT.0) THEN
                    IFUN = IFUNM(ICELL,LM)
                    
                    DO 40 IR = IMT1 + 1,IRC1
                       IRH = IR - IMT1
                       V1(IR) = V1(IR) + R(IR,IH)**2*V(IR,LM,IPOT)*
     +                      THETAS(IRH,IFUN,ICELL)
 40                 CONTINUE
                    
                 END IF
                 
 50           CONTINUE          ! LM = 2,LMPOT
              
              CALL SIMPK(V1,VAV1(IS),IPAN1,IRCUT(0,IH),DRDI(1,IH))
              CALL SIMPK(V2,VOL1(IS),IPAN1,IRCUT(0,IH),DRDI(1,IH))
              
           END IF               ! (IPAN1.EQ.1)
           
        END DO                  ! SPIN LOOP
        IF (NSPIN.EQ.1) THEN
        VAV1(2) = VAV1(1)
        VOL1(2) = VOL1(1)
        END IF

        
c     19.5.99   Nikos
c     This way it is compatible with old kkr and tb-kkr
        if(lsurf.and.(ih.eq.1)) 
     &       write(1337,*) 'Vacancies are ignored for VBC'
        
        IF (LSURF.AND.(Z(IH).LT.1.d0)) GOTO 60
        VAV0 = VAV0 + CONC(IH)*NSHELL(IH)*(VAV1(1)+VAV1(2))/2.D0
        VOL0 = VOL0 + CONC(IH)*NSHELL(IH)*(VOL1(1)+VOL1(2))/2.D0
 60   CONTINUE                  ! IH = 1,NATYP
      IF (.NOT.(OPT('DECIMATE'))) THEN ! added 10.11.99 to fix vbc
         VBC(1) = 0.0D0
         IF (ABS(VAV0).GT.1D-10) VBC(1) = -VAV0/VOL0
         IF (ISHIFT.GT.0) VBC(1) = VBC(1) + ESHIFT
      END IF
      
      WRITE (1337,FMT=9000) VOL0,VAV0,VBC(1)
      VBC(2) = VBC(1) 
c     
c---  > shift potential to muffin tin zero
c
      DO 90 IS = 1,NSPIN
         DO 80 IH = 1,NATYP
            IPOT = NSPIN* (IH-1) + IS
            DO 70 IR = 1,IRCUT(IPAN(IH),IH)
               V(IR,1,IPOT) = V(IR,1,IPOT) + RFPI*VBC(IS)
 70         CONTINUE
            
 80      CONTINUE
 90   CONTINUE
      
      RETURN
      
 9000 FORMAT ('  VOL INT.',F16.9,'  VAV INT.',F16.9,'  VMT ZERO',F16.9)
 9010 FORMAT ('  ATOM ',I4,' VMT ZERO :',F16.9)
      END
