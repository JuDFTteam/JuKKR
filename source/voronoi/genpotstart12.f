      SUBROUTINE GENPOTSTART(NSPIN,IFILE,I13,INS,NBEGIN,NATOMS,ZATOM,
     &                 SITEAT,KSHAPE,IDSHAPE,VOLUMECL,LPOT,
     &                 AOUT_ALL,RWSCL,RMTCL,RMTCORE,MESHN,XRN,DRN,
     &                 IRWS,IRNS,
     &                 ALATNEW,QBOUND,KXC,TXC,LJELL)
c ******************************************************
c * This subroutine reads a general potential format
c * file. and interpolates to the new mesh 
c ******************************************************
      use mod_splint, only: splint_real
      implicit none
c#@# KKRtags: VORONOI potential initialization input-output
      include 'inc.geometry'
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND,INSLPD
      PARAMETER (IRMIND=IRMD-IRNSD,INSLPD= (IRNSD+1)*LMPOTD)
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
      LOGICAL LJELL
C     ..
C     .. Scalar Arguments ..
      REAL*8           ALAT,C,EFERMI,HFIELD,VCONST,ALATNEW,
     &                 QBOUND
      INTEGER IFILE,IINFO,INS,IPE,IPF,IPFE,
     +        KHFELD,KSHAPE,KVREL,KWS,
     +        LMAX,LPOT,
     +        NBEG,NEND,NSPIN,KXC,NBEGIN
C     ..
C     .. Array Arguments ..
      REAL*8           ZATOM(NATYPD),VOLUMECL(*),RWSCL(*),RMTCL(*),
     &                 AOUT_ALL(NATYPD),
     &                 RMTCORE(*),XRN(IRID,NSHAPED),DRN(IRID,NSHAPED)
      INTEGER IRNS(NATYPD),IRWS(NATYPD),ITITLE(20),
     +        LCORE(20),NCORE,IDSHAPE(*),SITEAT(*)
C     ..
C     .. Local Scalars ..
      REAL*8           A1,B1,EA,EFNEW,S1,Z1,DUMMY,RMAX,RMTNW1,RMT1,
     &                 VBC(2),VCON,ECORE1(20),MAXA,AOUT,BOUT,RMTOUT,
     &                 PARSUM,PARSUMDERIV,R0,DIST,DR,RMAXOUT,RMTNEW
      INTEGER I,IA,ICELL,ICORE,IFUN,IH,IMT1,INEW,IO,IPAN1,IR,IRC1,IRI,
     +        IRMINM,IRMINP,IRNS1P,IRT1P,IRWS1,ISAVE,ISPIN,ISUM,
     +        J,NATOMS,IPOT,IRNSTOT,MESHN(NATYPD),MESHN0,ID,
     +        L,LM,LM1,LMPOT,LMPOTP,IRNSOUT,IRMTOUT,IRWSOUT,
     +        N,NCELL,NFUN,NR,IAT,IRNS1,NCORE1,LCORE1(20),IRC,
     &        I1,I2,ISITE
      LOGICAL TEST,POTLM(LMPOTD)
C     ..
C     .. Local Arrays ..
      REAL*8           U(IRMD),DRDI(IRMD),ECORE(20),
     &           RMESH(IRMD),RWS,VINS(IRMIND:IRMD,LMPOTD),
     &                 VM2Z(IRMD),VINSOUT(IRMIND:IRMD,LMPOTD),
     &                 VM2ZOUT(IRMD),VM2ZB(IRMD),ROUT(IRMD),
     &             VINSB(IRMD,LMPOTD),DRDIOUT(IRMD),            
     &             WORK(IRMD,LMPOTD),RA(IRMD) 
      CHARACTER*40 BANER,TEXT,I13
      CHARACTER*4 ELEM_NAME 
      CHARACTER*124 TXC(3)
c     --------------------------------------------------------------
      WRITE(6,*) ' ****  READING  POTENTIAL  **** '
      OPEN(IFILE,STATUS='OLD',FILE=I13,ERR=1000)
      OPEN(19,STATUS='UNKNOWN',FILE='output.pot')

      DO I2=1,LMPOTD
         DO I1=IRMIND,IRMD
            VINS(I1,I2) = 0.D0
         END DO
      END DO
      DO I1=1,IRMD
         VM2Z(I1) = 0.D0
      END DO  

      DO IAT = NBEGIN,NATOMS
         ISITE = SITEAT(IAT)
         ID = IDSHAPE(ISITE)

         WRITE(*,FMT='(A$,I6,A$,I6,A$,I6,A1)') 
     &              'Generating potential for atom',IAT,
     &              ' at site',ISITE,
     &              ' with shape',ID
         WRITE(*,*) ' '

         DO ISPIN = 1,NSPIN

            DO LM=1,LMPOTD
               POTLM(LM) =.FALSE.
            END DO
            IPOT =  NSPIN* (IAT-1) + ISPIN    
            READ (IFILE,FMT=8000) BANER
            WRITE(6,999) BANER(2:23)
 999        FORMAT ('###',A22,'###')
            IF (BANER(2:23).NE.'GENERAL POTENTIAL MESH') THEN
               WRITE(6,*) '  Input potential is not in the '
               WRITE(6,*) '  GENERAL POTENTIAL MESH format'
               WRITE(6,*) '* It cannot be interpolated  *'
               STOP
            END IF
            READ(IFILE,8010) ELEM_NAME,Z1
            IF (NSPIN.EQ.1) THEN
               WRITE (6,8011) ELEM_NAME,Z1
            ELSEIF (ISPIN.EQ.1) THEN
               WRITE (6,8012) ELEM_NAME,Z1
            ELSEIF (ISPIN.EQ.2) THEN
               WRITE (6,8013) ELEM_NAME,Z1
            END IF
            READ(IFILE,8020) TEXT
            READ(IFILE,8030) ALAT,RMAX,RMTNW1,RMT1
              !WRITE(6,*) ALAT,RMAX,RMTNW1,RMT1
            READ(IFILE,8040) NR,IMT1,IRNS1
              !WRITE(6,*) NR,IMT1,IRNS1
            READ(IFILE,8050) A1,B1
              !WRITE(6,*) A1,B1
            READ(IFILE,8060) EFERMI,VBC(ISPIN),VCON
              !WRITE(6,*) EFERMI,VBC(ISPIN),VCON
            READ(IFILE,8070) NCORE1,LMPOT
              !WRITE(6,*) NCORE1,LMPOT
            IF (NCORE1.GE.1) READ (IFILE,FMT=9140) (LCORE1(ICORE),
     +           ECORE1(ICORE),ICORE=1,NCORE1)
           ! WRITE(6,*) (LCORE1(ICORE),
    !+           ECORE1(ICORE),ICORE=1,NCORE1)

            IF (IRNS1.EQ.0) THEN
c     
c---  >       store only the spherically averaged potential
c     (in mt or as - case)
c     this is done always for the host
c     
               READ (IFILE,9051) (VM2Z(IR),IR=1,NR)
            ELSE ! (IRNS1.EQ.0)
c     
c---  >     store the full potential , but the non spherical contribution
c     only from irns1 up to irws1 ;
c     remember that the lm = 1 contribution is multiplied
c     by a factor 1/sqrt(4 pi)
c     
               READ (IFILE,9160)  IRT1P,IRNS1P,LMPOTP,ISAVE
               IRMINP = IRT1P - IRNS1P
               IRMINM = MAX(IRMINP,IRMIND)
               READ (IFILE,FMT=9100) (VM2Z(IR),IR=1,NR)
               IF (LMPOTP.GT.1) THEN
                  LM1 = 2
                  DO 50 LM = 2,LMPOTP
                     IF (LM1.NE.1) THEN
                        IF (ISAVE.EQ.1) THEN
                           READ (IFILE,FMT=9090) LM1
                        ELSE
                           LM1 = LM
                        END IF
                        
                        IF (LM1.GT.1) THEN
                           
                           READ (IFILE,FMT=9100) (U(IR),IR=IRMINP,NR)
                           
                           IF (LM1.LE.LMPOT) THEN
                              POTLM(LM1) = .TRUE.
                              DO 40 IR = IRMINM,NR
                                 VINS(IR,LM1) = U(IR)
 40                           CONTINUE
                           END IF
                           
                        END IF
                        
                     END IF
                     
 50               CONTINUE
               END IF   ! (LMPOTP.GT.1)


            END IF  ! (IRNS1.EQ.0)
c     
c     Now create mesh information
c         
            !write(6,*) ' info info info : potential read in!'
            IRWS1 =  NR
            if (irns1.ne.0) then   ! bug corrected 5.10.2000 
               IRWS1  = IMT1
            END IF 
                             
c     
c---  > generate radial mesh - potential only is stored in potential card
c     
            RMESH(1) = 0.0D0
            DRDI(1) = A1*B1
            DO 70 IR = 2,IRWS1
               EA = EXP(A1*REAL(IR-1))
               RMESH(IR) = B1* (EA-1.0D0)
               DRDI(IR) = A1*B1*EA
 70         CONTINUE          
            IF (IRNS1.NE.0) THEN
               MESHN0 = NR - IMT1
               DIST = RMAX - RMTNW1
               DR = DIST/MESHN0
               DO 80 IRI = 1,MESHN0
                  IR = IRI + IMT1
                  RMESH(IR) = RMTNW1 + DR*FLOAT(IRI)
                  DRDI(IR) = DR
 80            CONTINUE
            END IF
c
c     
c     The input mesh is constructed now Construct the output mesh
c     
            ID = IDSHAPE(IAT)           
            ROUT(1) = 0.D0
            AOUT = AOUT_ALL(ISITE)
            RMAXOUT = RWSCL(ID) 
            RMTOUT  = RMTCL(ID)
            IRWSOUT = IRWS(ISITE)
            IRMTOUT = IRWS(ISITE) - MESHN(ID)
            IRNSOUT = IRNS(ISITE)  ! 22.1.12 Changed from IRNS(ID) to IRNS(IAT)

            IF (KSHAPE.EQ.0) THEN
               BOUT = RMAXOUT / (EXP(AOUT*REAL(IRWSOUT-1))-1.0D0)
               DO IR=2,IRWSOUT
                  EA = EXP(AOUT*REAL(IR-1))
                  ROUT(IR) = BOUT* (EA-1.0D0)
                  DRDIOUT(IR) = AOUT*BOUT*EA
               END DO 
               DO I=1,IRWSOUT
                IF (ROUT(I).LT.RMTOUT) IRMTOUT = I
               END DO
               IF (MOD(IRMTOUT,2).EQ.0) IRMTOUT = IRMTOUT+1
               RMTNEW = ROUT(IRMTOUT)
               RMAXOUT = ROUT(IRWSOUT)
            ELSE
               BOUT = RMTOUT /  (EXP(AOUT*REAL(IRMTOUT-1))-1.0D0)
               DO IR=2,IRMTOUT
                  EA = EXP(AOUT*REAL(IR-1))
                  ROUT(IR) = BOUT* (EA-1.0D0)
                  DRDIOUT(IR) = AOUT*BOUT*EA
               END DO
               DO IRI=1,MESHN(ID)
                  IR = IRI + IRMTOUT
                  ROUT(IR) = ALATNEW*XRN(IRI,ID)   ! scaling is always 1.0d0
                  DRDIOUT(IR) = ALATNEW*DRN(IRI,ID)
               ENd DO
               RMTNEW = ROUT(IRMTOUT)
               RMAXOUT = ROUT(IRWSOUT)
            END IF
c
c
c  Ok now interpolate
c
                   
            MAXA = 1.D35
            CALL SPLINE(IRMD,RMESH,VM2Z,NR,MAXA,MAXA,VM2ZB)  
            IF (INS.GT.0) THEN
               DO LM1=1,LMPOTD
                  IF (POTLM(LM1)) THEN
                     IRNSTOT = NR - IRMINM ! same as IRNS1
c     map it
                     DO I=1,IRNS1
                        WORK(I,LM1) = VINS(NR - IRNS1 + I - 1,LM1)
                        RA(I) = RMESH(NR - IRNS1 + I - 1)
                     END DO
                     CALL SPLINE(IRMD,RA,WORK(1,LM1),IRNSTOT,MAXA,MAXA,
     &                    VINSB(1,LM1))
                  END IF
               ENDDO            ! LM1
            END IF
C
c OK with spline
c
            VM2ZOUT(1) = VM2Z(1)
            DO IR = 2,IRWSOUT
               R0 = ROUT(IR)
               CALL splint_real(RMESH,VM2Z,VM2ZB,NR,R0,PARSUM,
     &                           PARSUMDERIV)
               VM2ZOUT(IR) = PARSUM
            END DO
c
        
c
            IF (INS.GT.0) THEN
               IRC = IRWSOUT - IRNSOUT
               DO LM1=1,LMPOTD
                  IF (POTLM(LM1)) THEN
                     DO IR = IRC+1,IRWSOUT
                        R0 = ROUT(IR)
                        CALL splint_real(RA,WORK(1,LM1),
     &                       VINSB(1,LM1),IRNSTOT,R0,PARSUM,PARSUMDERIV)
                        VINSOUT(IR,LM1) = PARSUM
                     END DO
                  END IF
               END DO
            END IF

c           All interpolation ok now write
            CALL RITESONE(19,ISPIN,Z1,ALATNEW,RMTOUT,RMTNEW,RMAXOUT,
     +                    ITITLE,ROUT,DRDIOUT,VM2ZOUT,IRWSOUT,AOUT,BOUT,
     &                    TXC,KXC,INS,IRNSOUT,
     +                    LPOT,VINSOUT,QBOUND,IRWSOUT,KSHAPE,EFERMI,VBC,
     &                    ECORE1,LCORE1,NCORE1,ELEM_NAME,NSPIN)
c
c
c Next atom or next spin
c
         END DO ! ISPIN = 1,NSPIN
      END DO    ! IAT = NBEGIN,NATOMS

   
      RETURN

 1000 WRITE(6,*) 'Error read file......... ',I13
      WRITE(6,*) 'Error occured on atom... ',iat
      WRITE(6,*) 'Will use jellium starting potentials.'
      LJELL = .TRUE.
      RETURN

 8000 format (A40)
 8010 format (3X,A4,26X,F8.3)
 8011 Format ('#  ',A4,'POTENTIAL             Z = ',F8.3)
 8012 format ('#  ',A4,'POTENTIAL SPIN UP     Z=  ',F8.3)
 8013 format ('#  ',A4,'POTENTIAL SPIN DOWN   Z=  ',F8.3)
 8020  FORMAT(A40)
 8030 format (4f12.8)
 8040 format (1X,3I6)
 8050 format (2D15.8)
 8060 format (3f12.8)
 8070 format (2I5)
 9051 FORMAT (1p,4d20.12)
 9000 FORMAT (7a4,6x,'  exc:',a24,3x,a10)
 9010 FORMAT (3f12.8)
 9020 FORMAT (f10.5,/,f10.5,2f15.10)
 9030 FORMAT (i3,/,2d15.8,/,2i2)
 9140 FORMAT (i5,1p,d20.11)
 9040 FORMAT (f10.5,/,f10.5,2f15.10)
 9050 FORMAT (i3,/,2d15.8,/,2i2)
 9060 FORMAT (1p,2d15.6,1p,d15.8)
 9160 FORMAT (10i5)
 9061 FORMAT (1p,5d15.8)
 9070 FORMAT (i5,1p,d20.11)
c 9080 FORMAT (10x,20a4)
 9080 FORMAT (' < ',20a4)
 9081 FORMAT (' <#',20a4)
 9090 FORMAT (10i5)
 9100 FORMAT (1p,4d20.13)
        END








