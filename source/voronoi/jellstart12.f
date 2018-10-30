      SUBROUTINE JELLSTART(NSPIN,IFILE,I13,INS,NBEGIN,NATOMS,
     &                 ZATOM,SITEAT,KSHAPE,IDSHAPE,
     &                 VOLUMECL,LPOT,
     &                AOUT_ALL,RWSCL,RMTCL,RMTCORE,MESHN,XRN,DRN,THETAS,
     &                 LMIFUN,NFUN,IRWS,IRNS,
     &                 ALATNEW,QBOUND,KXC,TXC)
c ******************************************************
c * This subroutine reads a jellium potential from the database
c * file. and interpolates to the new mesh 
c ******************************************************
      implicit none
c#@# KKRtags: VORONOI initialization potential core-electrons
      include 'inc.geometry'
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IBMAXD
      PARAMETER (IBMAXD=(LMAXD1+1)*(LMAXD1+1))
      INTEGER IRMIND,INSLPD
      PARAMETER (IRMIND=IRMD-IRNSD,INSLPD= (IRNSD+1)*LMPOTD)
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
      INTEGER IRMDJJ
      PARAMETER(IRMDJJ=1501)
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
     &                 RMTCORE(*),XRN(IRID,NSHAPED),DRN(IRID,NSHAPED),
     &                 THETAS(IRID,IBMAXD,NSHAPED)
      INTEGER IRNS(NATYPD),IRWS(NATYPD),ITITLE(20),
     &        LCORE(20),NCORE,IDSHAPE(*),LMIFUN(IBMAXD,NSHAPED),
     &        NFUN(NSHAPED),SITEAT(*)
C     ..
C     .. Local Scalars ..
      REAL*8           A1,B1,EA,EFNEW,S1,Z1,DUMMY,RMAX,RMTNW1,RMT1,
     &                 VBC(2),VCON,ECORE1(20),MAXA,AOUT,BOUT,RMTOUT,
     &                 PARSUM,PARSUMDERIV,R0,DIST,DR,RMAXOUT,RMTNEW
      REAL*8           RWS0,BR,ZA,ZVALI,EINF,AR,AMSH,RFPI,ZZOR
      INTEGER I,IA,ICELL,ICORE,IFUN,IH,IMT1,INEW,IO,IPAN1,IR,IRC1,IRI,
     +        IRMINM,IRMINP,IRNS1P,IRT1P,IRWS1,ISAVE,ISPIN,ISUM,
     +        J,NATOMS,IPOT,IRNSTOT,MESHN(NATYPD),MESHN0,ID,
     +        L,LM,LM1,LMPOT,LMPOTP,IRNSOUT,IRMTOUT,IRWSOUT,
     +        N,NCELL,NR,IAT,IRNS1,NCORE1,LCORE1(20),IRC,NZ,
     &        IRS1,NSEC,NZVALI,NC,II,I1,I2,ISITE
      LOGICAL TEST,POTLM(LMPOTD)
C     ..
C     .. Local Arrays ..
      REAL*8           U(IRMD),DRDI(IRMDJJ),ECORE(20),
     &           RMESH(IRMDJJ),RWS,VINS(IRMIND:IRMD,LMPOTD),
     &                 VM2Z(IRMDJJ),VINSOUT(IRMIND:IRMD,LMPOTD),
     &                 VM2ZOUT(IRMD),VM2ZB(IRMDJJ),ROUT(IRMD),
     &             VINSB(IRMD,LMPOTD),DRDIOUT(IRMD),            
     &             WORK(IRMD,LMPOTD),RA(IRMD) 
      CHARACTER*40 BANER,TEXT,I13
      CHARACTER*4 ELEM_NAME,AAAA,TRAN,exte 
      CHARACTER*124 TXC(3)
      CHARACTER*4 ELEM_FILE(0:113)
      CHARACTER*26 ATOMPOT
      CHARACTER*2 TXTC(20)
      CHARACTER*20 DATA1
       DATA ELEM_FILE/'Vac0',
     & 'H_01','He02','Li03','Be04','B_05','C_06','N_07','O_08',
     & 'F_09','Ne10',
     & 'Na11','Mg12','Al13','Si14','P_15','S_16','Cl17','Ar18',
     & 'K_19','Ca20','Sc21','Ti22',
     & 'V_23','Cr24','Mn25','Fe26','Co27','Ni28',
     & 'Cu29','Zn30',
     & 'Ga31','Ge32','As33','Se34','Br35','Kr36','Rb37','Sr38',
     & 'Y_39','Zr40',
     & 'Nb41','Mo42','Tc43','Ru44','Rh45','Pd46','Ag47','Cd48',
     & 'In49','Sn50',
     & 'Sb51','Te52','I_53','Xe54','Cs55','Ba56','La57','Ce58',
     & 'Pr59','Nd60',
     & 'Pm61','Sm62','Eu63','Gd64','Tb65','Dy66','Ho67','Er68',
     & 'Tm69','Yb70',
     & 'Lu71','Hf72','Ta73','W_74','Re75','Os76','Ir77','Pt78',
     & 'Au79','Hg80',
     & 'Tl81','Pb82','Bi83','Po84','At85','Rn86','Fr87','Ra88',
     & 'Ac89','Th90',
     & 'Pa91','U_92','Np93','Pu94','Am95','Cm96','Bk97','Cf98',
     & 'Es99','Fm__',
     & 'Md__','No__','Lr__','Rf__','Db__','Sg__','Bh__','Hs__',
     & 'Mt__','Uun_',
     & 'Uuu_','Uub_','NoE_'/
c     --------------------------------------------------------------
      

      RFPI = 4.D0*DSQRT(DATAN(1.D0))

      WRITE(6,*) ' ****  JELLSTART POTENTIALS **** '
      WRITE(6,*) 'From atom No.',NBEGIN,'to atom No.',NATOMS

      OPEN(19,STATUS='UNKNOWN',FILE='output.pot')
      DO I2=1,LMPOTD
         DO I1=IRMIND,IRMD
            VINS(I1,I2) = 0.D0
         END DO
      END DO
      DO I1=1,IRMDJJ
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

         DO ISPIN=1,NSPIN
            DO LM=1,LMPOTD
               POTLM(LM) =.FALSE.
            END DO
            IPOT =  NSPIN* (IAT-1) + ISPIN    

c Find out what atom is needed            
c            
            NZ = ZATOM(IAT)
            IF (((NZ.GE.24.AND.NZ.LE.28).OR.(NZ.GE.57.AND.NZ.LE.70))
     &                                          .AND.ISPIN.EQ.2) THEN
            ATOMPOT = 'ElementDataBase/'//ELEM_FILE(NZ)//'.pots2'
            ELSE
            ATOMPOT = 'ElementDataBase/'//ELEM_FILE(NZ)//'.pot  '
            END IF
            WRITE(6,*) 'Using database ....: ',ATOMPOT
            OPEN(21,STATUS='OLD',FILE=ATOMPOT,ERR=1010)
c           IRWS1 =  NR
            
c --------------------------------------------------------------------
          
            EFERMI = .409241D+00
            VBC(1)    = .500D0
            VBC(2)    = .500D0
c           read potential from jellium
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
            READ(21,141) BANER,AAAA,AAAA
            READ(21,142) RWS0,S1,IRS1,BR
            READ(21,142) ZA,ZVALI,NSEC,EINF
c     Calculate number of core states
            NZ = ZA             ! make integer
            NZVALI = ZVALI      ! make integer
            NC = NZ - NZVALI
            NCORE = 0
            IF (NC.EQ.2 ) NCORE = 1 ! 1s
            IF (NC.EQ.4 ) NCORE = 2 ! 1s2s
            IF (NC.EQ.10) NCORE = 3 ! 1s2s2p
            IF (NC.EQ.12) NCORE = 4 ! 1s2s2p3s
            IF (NC.EQ.18) NCORE = 5 ! 1s2s2p3s3p
            IF (NC.EQ.28) NCORE = 6 ! 1s2s2p3s3p3d
            IF (NC.EQ.30) NCORE = 7 ! 1s2s2p3s3p3d4s
            IF (NC.EQ.36) NCORE = 8 ! 1s2s2p3s3p3d4s4p
            IF (NC.EQ.46) NCORE = 9 ! 1s2s2p3s3p3d4s4p4d
            IF (NC.EQ.48) NCORE = 10 ! 1s2s2p3s3p3d4s4p4d4s
            IF (NC.EQ.54) NCORE = 11 ! 1s2s2p3s3p3d4s4p4d4s4p
            IF (NC.EQ.68) NCORE = 12 ! 1s2s2p3s3p3d4s4p4d4s4p4f
            IF (NC.EQ.78) NCORE = 13 ! 1s2s2p3s3p3d4s4p4d4s4p4f5d  
            IF (NC.EQ.80) NCORE = 14 ! 1s2s2p3s3p3d4s4p4d4s4p4f5d6s
            IF (NC.EQ.86) NCORE = 15 ! 1s2s2p3s3p3d4s4p4d4s4p4f5d6s4p 
            WRITE(6,*) '*************************************'
            WRITE(6,*) '   Potential Interpolation Program   '
            WRITE(6,*) '   Using the Jellium Database v1.0   '
            WRITE(6,*) '*************************************'
            WRITE(6,163) EFERMI
            WRITE(6,161) ZA  
            WRITE (6,162) NCORE
            READ(21,133)(LCORE(NC),TXTC(NC),ECORE(NC),NC=1,NCORE)
            IF (TEST('verb0   '))
     &              WRITE(6,*) ' ** Position of the Core States ** '
            DO I=1,NCORE
               IF (TEST('verb0   ')) 
     &              WRITE(6,135) LCORE(I),TXTC(I),ECORE(I)
               IF (TXTC(I).eq.'s ') LCORE(I) = 0
               IF (TXTC(I).eq.'p ') LCORE(I) = 1
               IF (TXTC(I).eq.'d ') LCORE(I) = 2
               IF (TXTC(I).eq.'f ') LCORE(I) = 3
            END DO
            NCORE1 = NCORE
            DO I=1,NCORE1
               LCORE1(I) = LCORE(I)
               ECORE1(I) = ECORE(I)
            END DO              
            WRITE(6,134) EINF
            WRITE(6,*) '**********************************************'  
c 
            READ(21,131)(VM2Z(II),II=1,IRS1)
            READ(21,132)TRAN
            CLOSE (21)
c     make mesh r0
            AR = LOG(S1/BR+1.D0)/FLOAT(IRS1-1)
            EA=EXP(AR)
            AMSH=1.D0
            RMESH(1)=0.D0
            DRDI(1)=BR*AR
            DO 1 I=2,IRS1
               AMSH=AMSH*EA
               RMESH(I)=BR*(AMSH-1.D0)
               DRDI(I)=DRDI(1)*AMSH
 1          CONTINUE
            WRITE(6,*) 'Jellium Potential Read In ',IRS1,' points'
            NR = IRS1
            Z1 = ZA
C     
 131        FORMAT(4D15.8)
 132        FORMAT(A4)
 133        FORMAT(I3,A2,D15.8)
 135        FORMAT(I3,A2,F15.6,' Ry')
 134        FORMAT('All other states are above :',F8.4,' Ry in Energy')
 141        FORMAT(3X,A40,3X,A4,3X,A4)
 142        FORMAT(7X,F8.4,7X,F8.4,7X,I5,7X,F8.4)
 161        FORMAT('Potential Atomic Number :',F7.2)
 162        FORMAT('Number of Core   States :',I4)
 163        FORMAT('Jellium Fermi Energy :',F10.5,' Ry')
c --------------------------------------------------------------------
c     
c     The input mesh has been constructed. Now construct the output mesh.
c     
            ROUT(1) = 0.D0
            AOUT = AOUT_ALL(ISITE)   
            RMAXOUT = RWSCL(ID) 
            RMTOUT  = RMTCL(ID)
            IRWSOUT = IRWS(ISITE)
            IRMTOUT = IRWS(ISITE) - MESHN(ID)
            IRNSOUT = IRNS(ISITE)


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
               END DO
               RMTNEW = ROUT(IRMTOUT)
               RMAXOUT = ROUT(IRWSOUT)
            END IF
c
c  Ok now interpolate
c
                   
            MAXA = 1.D35
            CALL SPLINE(IRMDJJ,RMESH,VM2Z,NR,MAXA,MAXA,VM2ZB)            
c
c OK with spline
c
            VM2ZOUT(1) = VM2Z(1)
            DO IR = 2,IRWSOUT
               R0 = ROUT(IR)
               CALL SPLINT(RMESH,VM2Z,VM2ZB,NR,R0,PARSUM,PARSUMDERIV)
               VM2ZOUT(IR) = PARSUM
            END DO
c

        
c
            IF (INS.GT.0) THEN
               IRC = IRWSOUT - IRNSOUT
               DO LM1=1,LMPOTD
                  DO IR = IRC,IRWSOUT  
                     VINSOUT(IR,LM1) = 0.d0                   
                  END DO
               END DO
            END IF

! Convolute with shapes
            IF (KSHAPE.GT.0) THEN 
               LMPOT = (LPOT + 1)**2
               DO IRI=1,MESHN(ID)
                  IR = IRI + IRMTOUT
                  ZZOR = 2.D0 * ZATOM(IAT) / ROUT(IR) 
                  DO IFUN = 2,NFUN(ID)
                     LM1 = LMIFUN(IFUN,ID)
                     IF (LM1.LE.LMPOT) THEN
                        VINSOUT(IR,LM1) = ! (VM2ZOUT(IR) - ZZOR) * 
     &                                    THETAS(IRI,IFUN,ID) 
                     ENDIF
                  ENDDO
                  VM2ZOUT(IR) = (VM2ZOUT(IR) - ZZOR ) * 
     &                           THETAS(IRI,1,ID) / RFPI + ZZOR ! Because the ratial solver adds -2*Z/R by default without shape convolution
                                                                ! and assumes that the spher. component is without sqrt(4pi) factor.
               END DO
            ENDIF

            CALL RITESONE(19,ISPIN,Z1,ALATNEW,RMTOUT,RMTNEW,RMAXOUT,
     &                    ITITLE,ROUT,DRDIOUT,VM2ZOUT,IRWSOUT,AOUT,BOUT,
     &                    TXC,KXC,INS,IRNSOUT,
     &                    LPOT,VINSOUT,QBOUND,IRWSOUT,KSHAPE,EFERMI,VBC,
     &                    ECORE1,LCORE1,NCORE1,ELEM_FILE(NZ),NSPIN) 
c
c Next atom or next spin
c
         END DO  ! ISPIN=1,NSPIN
      END DO     ! IAT = NBEGIN,NATOMS
   
      RETURN


 1000  WRITE(6,*) 'Error read file......... ',I13
       WRITE(6,*) 'Error occured on atom... ',iat
       STOP
 1010  WRITE(6,*) ' Error in JELLSTART '
       WRITE(6,*) ' Potential.............',ELEM_FILE(NZ)
       WRITE(6,*) ' Does not exist in the database'
       STOP
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








