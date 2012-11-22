      SUBROUTINE RINPUT99(BRAVAIS,ALAT,RBASIS,ABASIS,BBASIS,CBASIS,
     &           CLS,NCLS,E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3,
     +           NSTEPS,IMIX,MIXING,QBOUND,FCM,ITDBRY,
     +           NTCELL,NAEZ,IRM,Z,
     +           NREF,
     +           ICST,IFILE,IPE,IPF,IPFE,
     +           KHFELD,KPRE,KTE,
     +           KVMAD,KVREL,KXC,LMAX,LMPOT,LPOT, 
     +           NSPIN,
     +           REFPOT,
     +           INTERVX,INTERVY,INTERVZ,
     +           HFIELD,
     +           VBC,VCONST,
     +           I13,I19,
     +           RCUTZ,RCUTXY,RCUTJIJ,JIJ,RCUTTRC,
     +           LDAU,
     +           RMTREF,KFORCE,
     +           IGUESS,BCP,QMRBOUND,LCARTESIAN,RMAX,GMAX,
     &           LMAXD, IRNSD, TRC, LPOTD, NSPIND,
     &           IRMD, NAEZD)

      IMPLICIT NONE
c
c

C     .. Local Arrays ..
      CHARACTER*4 TSPIN(3)
      CHARACTER*43 TKCOR(0:3),TVREL(0:2)
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Array Arguments ..
      INTEGER IRNS_dummy,KFGdummy(4),LMXCdummy,
     &        NTCELL(naezd),CLS(*),REFPOT(naezd)

      DOUBLE PRECISION Z(*),MTFACdummy,VBC(*),RBASIS(3,*),RMTREF(*)
      DOUBLE PRECISION BRAVAIS(3,3)

      CHARACTER*24 TXC(4)
      CHARACTER*80 UIO

      CHARACTER*40 I12,I13,I19,I40
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION E1,E2,FCM,HFIELD,MIXING,QBOUND,TK,
     +                 VCONST,ABASIS,BBASIS,CBASIS,
     +                 RCUTZ,RCUTXY,RCUTJIJ,RCUTTRC
      INTEGER ICST,IFILE,IGF,IMIX,
     +        IPE,IPF,IPFE,
     +        IRM,ISHIFT,ITDBRY,KCOR,
     +        KFROZN,KHFELD,KPRE,KTE,KVMAD,
     +        KVREL,KWS,KSHAPE,KXC,LMAX,LMPOT,LPOT,KFORCE,
     +        NPNT1,NPNT2,NPNT3,NPOL,NSPIN,IGUESS,BCP
      INTEGER NSTEPS,NAEZ
      DOUBLE PRECISION ALAT,QMRBOUND
      double precision temp
      INTEGER INTERVX,INTERVY,INTERVZ,NREF,NCLS
      LOGICAL LINIPOL,JIJ,LDAU

      LOGICAL :: LCARTESIAN
      DOUBLE PRECISION :: RMAX, GMAX
      INTEGER, intent(in) :: LMAXD, IRNSD, TRC, LPOTD, NSPIND,
     &                       IRMD, NAEZD
                                 ! atom types located at a given site
C-----------------------------------------------------------------------
      
C     
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BRYMIX,STRMIX,DVEC(10)
      INTEGER I,IL,J,IER
      CHARACTER*43 TSHAPE
      LOGICAL OPT
c
      CHARACTER*8 TESTC(16),OPTC(8),VERSION
      COMMON /TESTC/TESTC
      COMMON /OPTC/OPTC
C     ..
C     .. Data statements ..
      DATA VERSION /'Jun 2010'/
      DATA TSPIN/'non-','    ','    '/
      DATA TSHAPE/' exact cell treatment (shape correction)  '/
      DATA TVREL/
     +     ' non relativistic calculation              ',
     +     ' s.r.a. calculation                        ',
     +     ' fully relativistic calculation            '/
      DATA TKCOR/
     +     ' frozen core approximation                 ',
     +     ' core relaxation s.r.a.                    ',
     +     ' core relaxation nonsra                    ',
     +     ' core relaxation                           '/

c------------ array set up and definition of input parameter -----------
c
      KWS = 2
      KSHAPE = 2
      TXC(1) = ' Morruzi,Janak,Williams '
      TXC(2) = ' von Barth,Hedin        '
      TXC(3) = ' Vosko,Wilk,Nusair      '
      TXC(4) = ' GGA PW91               '
C
      WRITE (6,2004) VERSION
C
C read RUNNING options
C
      CALL IoInput('RUNOPT    ',UIO,1,7,IER)
                   READ (UNIT=UIO,FMT=980)(OPTC(I),I=1,8)
C
C read TEST options
C
      CALL IoInput('TESTOPT   ',UIO,1,7,IER)
                     READ(UNIT=UIO,FMT=980)(TESTC(i),i=1,8)
      CALL IoInput('TESTOPT   ',UIO,2,7,IER)
                     READ(UNIT=UIO,FMT=980)(TESTC(8+i),i=1,8)

      IL=1
      CALL IoInput('NSTEPS    ',UIO,1,7,IER)
                         READ (UNIT=UIO,FMT=*) NSTEPS
      CALL IoInput('NSPIN     ',UIO,1,7,IER)
                         READ (UNIT=UIO,FMT=*) NSPIN
      write(6,2011) nsteps
      write(6,2104)
      write(6,2010) nspin
      write(6,2104)
C
      CALL IoInput('NAEZ      ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) naez
      if (NAEZ.gt.NAEZD) then
               write(6,*) ' set NAEZD to at least ',NAEZ
               stop ' in < RINPUT99 > '
      end if
C
      OPEN(77,FILE='atominfo',FORM='formatted')
      DO I=1,NAEZ

                           READ (UNIT=77,FMT=*)    Z(I),
     +                        LMXCdummy,
     +                       (KFGdummy(J),J=1,4),
     +                        CLS(I),
     +                        REFPOT(I),
     +                        NTCELL(I),
     +                        MTFACdummy,
     +                        IRNS_dummy,
     +                        temp

c     E.R. check for possible out of bounds error
      IF (REFPOT(I) < 1 .or. REFPOT(I) > naezd) then
        WRITE(*,*) "Error in startb1."
        WRITE(*,*) "check atominfo for atom ", I
        STOP
      ENDIF

      RMTREF(REFPOT(I)) = temp

      END DO
      CLOSE (77)
 
      write(6,2028) NAEZ
c
c---> read input
c
      IL=1
      CALL IoInput('LMAX      ',UIO,0,7,IER)
                    READ (UNIT=UIO,FMT=*) LMAX
      CALL IoInput('EMIN      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) E1
      CALL IoInput('EMAX      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) E2

      CALL IoInput('TEMPR     ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) TK

      CALL IoInput('NPOL      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) NPOL
      CALL IoInput('NPT1      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) NPNT1
      CALL IoInput('NPT2      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) NPNT2
      CALL IoInput('NPT3      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) NPNT3

      CALL IoInput('IFILE     ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) ifile
      CALL IoInput('IPE       ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) ipe
      CALL IoInput('ISHIFT    ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) ishift
      IF ( OPT('rigid-ef') .OR. OPT('DECIMATE') ) THEN
        ISHIFT = 2
        WRITE(6,*) ' Rigid Fermi Energy, ISHIFT is set to ',ISHIFT
      END IF
C
      CALL IoInput('IRM       ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) irm
      CALL IoInput('ICST      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) icst

      CALL IoInput('KCOR      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kcor
      CALL IoInput('KVREL     ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kvrel
      CALL IoInput('KHFIELD   ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) khfeld
      CALL IoInput('KEXCOR    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kxc    
c -------------------------------------------------
      CALL IoInput('KTE       ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kte
      CALL IoInput('KPRE      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kpre
      CALL IoInput('KVMAD     ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kvmad
c
      CALL IoInput('IMIX      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) imix
      CALL IoInput('IGREENFUN ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) igf
      CALL IoInput('ITDBRY    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) itdbry
      CALL IoInput('STRMIX    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) strmix
      CALL IoInput('FCM       ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) fcm
      CALL IoInput('QBOUND    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) qbound
      CALL IoInput('QMRBOUND  ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) qmrbound
      CALL IoInput('IGUESS    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) iguess
      CALL IoInput('BCP       ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) bcp
      CALL IoInput('BRYMIX    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) brymix
      CALL IoInput('HFIELD    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) hfield
      CALL IoInput('VCONST    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) vconst
      CALL IoInput('CARTESIAN ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) lcartesian

C
C----------------------------------------------------------------------
C
C --> determination of properties at Fermi level
C
C     has some effect in EMESHT
C     IGF = 1 has no effect
      IF ( OPT('GF-EF   ') ) THEN
         IGF = 1
         IF (NPOL.GT.0) NPOL = 0
         IF (NPOL.LT.0) THEN 
            NPNT1 = 0
            NPNT3 = 0
         END IF
         NPNT2 = 1
      END IF
C ----------------------------------------------------------------------
C     Force calculation 18.5.2000
C
      KFORCE = 0
      IER = 0
      CALL IoInput('KFORCE   ',UIO,1,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) KFORCE

      KFROZN = KCOR
      IF (KCOR.EQ.0) KCOR = 2
c ------------------------------------------------------------------------
      WRITE (6,9210) LMAX
      WRITE (6,9301)
      WRITE (6,9220) E1,E2,TK
      WRITE (6,9302)
      WRITE (6,9230) NPOL,NPNT1,NPNT2,NPNT3
      WRITE (6,9304)
      WRITE (6,9303)
      WRITE (6,9250) IFILE,IPE,ISHIFT
      WRITE (6,9305)
      WRITE (6,9260) KSHAPE,IRM,ICST
      WRITE (6,9309)
      WRITE (6,9270) KCOR,KVREL,KWS,KHFELD,KXC
      WRITE (6,9306)
      WRITE (6,9330) KTE,KPRE,KVMAD
      WRITE (6,9309)
      WRITE (6,9290) IMIX
      WRITE (6,9304)
      WRITE (6,9300) ITDBRY
      WRITE (6,9307)
      WRITE (6,9310) STRMIX,FCM,QBOUND
      WRITE (6,9302)
      WRITE (6,9320) BRYMIX
      WRITE (6,9308)
      WRITE (6,9280) HFIELD,VCONST
c ------------------------------------------------------------------------
c
c
      IPF = 6
      IPFE = IPF + 3
C
      IF (IMIX.GT.2) THEN
        FCM = 1.0D0
        MIXING = BRYMIX
      ELSE
        MIXING = STRMIX
      END IF
c
      IF (IMIX.GE.6) WRITE (6,FMT=9110) (IMIX-5),ITDBRY - 1
c
      WRITE (6,FMT=9090) MIXING,QBOUND
c
      LPOT  = MIN(2*LMAX,LPOTD)
      LMPOT = (LPOT+1)* (LPOT+1)
c
      WRITE (6,FMT=9020) LMAX,LMAXD,NAEZ,NAEZD,IRM,IRMD,NSPIN,NSPIND

c      IF (LMAX.GT.LMAXD .OR. NAEZ.GT.NAEZD .OR. IRM.GT.IRMD .OR.
c     +    NSPIN.GT.NSPIND) CALL RCSTOP('18      ')

      WRITE (6,FMT=9130)
      WRITE (6,FMT=9140)
c     DO I = 1,NAEZ
c       WRITE (6,FMT=9150) I,IRNS(I),IRNSD
c       IF (IRNS(I).GT.IRNSD) CALL RCSTOP('19      ')
c     ENDDO


      IF (LMAX.NE.LMAXD) THEN
        WRITE (6,FMT=9120)

        CALL RCSTOP('20      ')

      END IF

      WRITE (6,FMT=9130)
c
c
c
c
      IF (KHFELD.EQ.1) WRITE (6,FMT=9030) HFIELD
      if (kvrel.le.1 ) then 
          WRITE (6,FMT=9050) TSPIN(NSPIN)
      else
          write (6,fmt=9050) tspin(nspin+1)
      end if
      WRITE (6,FMT=9170) TVREL(KVREL)
      WRITE (6,FMT=9170) TKCOR(KFROZN)
      WRITE (6,FMT=9170) TSHAPE

      WRITE (6,FMT=9100) TXC(KXC+1)
      WRITE (6,FMT=9080)

c
      VBC(1) = VCONST
      VBC(2) = VBC(1)

c
c
c----------------------------------------------------------------------
c
      OPEN(77,FILE='rbasis',FORM='formatted')
      DO I=1,NAEZ
             READ (UNIT=77,FMT=*) (RBASIS(J,I), J=1,3)
      ENDDO                         
      CLOSE (77)
c
c----------------------------------------------------------------------
c
      CALL IoInput('RCLUSTZ   ',UIO,IL,7,IER)
                READ (UNIT=UIO,FMT=*) RCUTZ
      CALL IoInput('RCLUSTXY  ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) RCUTXY
      WRITE(6,*) 'Parameters used for the cluster calculation'
      if (abs(rcutz-rcutxy).lt.1.d-4) then
      write(6,*) 'Clusters inside spheres with radius R = ',rcutz
      else
      write(6,*) 'Clusters inside cylinders with '
      write(6,*) 'Rz = ',rcutz,' Rxy = ',rcutxy
      end if
      write(6,2104)
      write(6,2018)                 ! rbasis
      write(6,2101)
      do i=1,naez
        write(6,2025) i,(rbasis(j,i),j=1,3)
      enddo
c-------------------------------------------------------------
c JIJ calculation switched on by logical switch LJIJ
c-------------------------------------------------------------
      CALL IoInput('LJIJ      ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) JIJ
      IF (JIJ) THEN
      CALL IoInput('RCUTJIJ  ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) RCUTJIJ
      WRITE(6,*) 'Radius for Jij calculation is RCUTJIJ=',RCUTJIJ
      ENDIF
c-------------------------------------------------------------
c in case of active truncation (TRC=1) read radius of truncation
c zone RCUTTRC
c-------------------------------------------------------------
      IF (TRC.EQ.1) THEN
      CALL IoInput('RCUTTRC  ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) RCUTTRC
      WRITE(6,*) 'Radius for Truncation is RCUTTRC=',RCUTTRC
      ENDIF
c-------------------------------------------------------------
c
c read flag in inputcard: LOGICAL LLDAU
c
      CALL IoInput('LLDAU ',UIO,IL,7,IER) 
      READ (UNIT=UIO,FMT=*) LDAU

C ----------------------------------------------------------------------
C Read in the bravais vectors (normalised to alat)
C Notation: BRAVAIS(J,I) J=x,y,z I=1,2,3
C ----------------------------------------------------------------------
C
      DO I=1,3
        DO J=1,3
            BRAVAIS(J,I) = 0D0
         END DO
      END DO

      IER = 0
      DO I = 1,3
         CALL IOINPUT('BRAVAIS   ',UIO,I,7,IER)
         READ (UNIT=UIO,FMT=*) (BRAVAIS(J,I),J=1,3)
      END DO


      CALL IoInput('BASISCALE ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) (DVEC(I),I=1,3)
      ABASIS = DVEC(1)
      BBASIS = DVEC(2) 
      CBASIS = DVEC(3)

      CALL IoInput('ALATBASIS ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) ALAT
      CALL IoInput('BZDIVIDE  ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) INTERVX,INTERVY,INTERVZ
      
      WRITE(6,2019) ABASIS,BBASIS,CBASIS 
      WRITE(6,2107)
      WRITE(6,2014) ALAT
      WRITE(6,2104)
      WRITE(6,2015) INTERVX,INTERVY,INTERVZ 
      WRITE(6,2102)

      IER = 0
      CALL IOINPUT('RMAX      ',UIO,1,7,IER)
      READ (UNIT=UIO,FMT=*) RMAX
      IF ( IER.NE.0 ) STOP ' ERROR: Invalid RMAX setting in the input'
C
      CALL IOINPUT('GMAX      ',UIO,1,7,IER)
      READ (UNIT=UIO,FMT=*) GMAX
      IF ( IER.NE.0 ) STOP ' ERROR: Invalid RMAX setting in the input'
c ------------------------------------------------------------------------

      NCLS = 0
      NREF = 0
c
      DO I=1,NAEZ 
        NCLS = MAX(NCLS,CLS(I)) 
      ENDDO
      DO I=1,NAEZ
        NREF = MAX(NREF,REFPOT(I)) 
      ENDDO
c
      WRITE(6,2016) NCLS,NREF
      WRITE(6,2110)
      WRITE(6,2103)

c ------------------------------------------------------------------------
C
C
      WRITE(6,62) (OPTC(I),I=1,8)                                              
 62   FORMAT(79('-')/' EXECUTION OPTIONS:'/1X,A8,7('//',A8)/79('-'))
      WRITE(6,52) (TESTC(I),I=1,16)                                    
 52   FORMAT(79('-')/' TEST OPTIONS:'/2(1X,A8,7('//',A8)/)/79('-'))
 980  FORMAT(8A8)
C
C ---------------------------------------------------------------------
C 
      IL=1
      CALL IoInput('FILES     ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT='(A40)')  I12
      CALL IoInput('FILES     ',UIO,IL+1,7,IER)
                      READ (UNIT=UIO,FMT='(A40)')  I13
      CALL IoInput('FILES     ',UIO,IL+2,7,IER)
                      READ (UNIT=UIO,FMT='(A40)')  I40
      CALL IoInput('FILES     ',UIO,IL+3,7,IER)
                      READ (UNIT=UIO,FMT='(A40)')  I19

      write(6,*) 'I12="',I12,'"'
      write(6,*) 'I13="',I13,'"'
      write(6,*) 'I40="',I40,'"'
      write(6,*) 'I19="',I19,'"'
      write(6,2100) 
      WRITE(6,2110)
      WRITE(6,*) ' >>>>>>>>> RINPUT99 EXITS NOW <<<<<<<<<< '
      RETURN
C
C *********************************************Input-End ********
C
 2004 FORMAT( /80(1H=)/
     & '|',78X,'|'/
     & '|',1X,'KKRnano ',69X,'|'/
     & '|',1X,'Massively Parallel Screened Korringa-Kohn-Rostoker ',
     &         'Electronic Structure Code',1X,'|'/
     & '|',1X,'for Bulk',69X,'|'/
     & '|',78X,'|'/
     & '|',1X,'established Juelich 2008',26X,
     & 'Version : ',A8,9X,'|'/
     & '|',78X,'|'/80(1H=))
 2010 FORMAT(' NSPIN '/I4)
 2011 FORMAT(' NSTEPS'/I4)
 2014 FORMAT('          ALAT = ',F15.8)
 2015 FORMAT('   INTERVX   INTERVY   INTERVZ'/3I10)
 2016 FORMAT('    NCLS    NREF   '/,2I8)
 2018 FORMAT(' RBASIS'/,
     &     'SITE                BASIS VECTORS                 ')
 2019 FORMAT('         ABASIS         BBASIS         CBASIS'/3F15.8)
 2025 FORMAT((i4,3F15.8))
 2028 FORMAT(' NAEZ ',/,I8)
C ------------------------------------------------------------------------
 2100 FORMAT(79(1H-))
 2101 format(   3(1H-),1H+  , 3(14(1H-),1H+),  30(1H-))
 2102 format( 3(9(1H-),1H+) ,49(1H-))
 2103 FORMAT(10(3(1H-),1H+) ,39(1H-))
 2104 format(   3(1H-),1H+  ,75(1H-))
 2107 format( 3(14(1H-),1H+),34(1H-))
 2110 format( 3(7(1H-),1H+) ,55(1H-))
 9000 FORMAT (I2,3X,4I5)
 9010 FORMAT (1X,I1,1X,4I1)
 9020 FORMAT (/,33x,'check of dimension-data consistency',/,33x,
     +       35 ('-'),/,40x,'lmax   : (',i6,',',i6,')',/,40x,
     +       'NAEZ  : (',i6,',',i6,')',/,40x,'irm    : (',i6,',',i6,
     +       ')',/,40x,'nspin  : (',i6,',',i6,')',/)
 9030 FORMAT (1x,10 ('*'),' external magnetic field applied hfield=',
     +       f8.5)
 9040 FORMAT (3f12.7)
 9050 FORMAT (20x,a4,'spin polarized calculation')
 9060 FORMAT (8i4)
 9070 FORMAT (1x,20x,' calculation with',a8,'-potential')
 9080 FORMAT (1x,79 ('*'))
 9090 FORMAT (' mixing factor used           :',f15.6,/,
     +        ' convergence quality required :',1p,d15.2)
 9100 FORMAT (1x,20x,a24,'exchange-correlation potential')
 9110 FORMAT (/,20x,'broyden"s method # :',i3,
     +       ' is used up to iteration-      ',/,20x,'depth :',i3,
     +       '  then jacobian is fixed and potential      ',/,20x,
     +       'is updated using that jacobian')
 9120 FORMAT (13x,' in case of calculating non - spherical wavefcts ',
     +       'the parameter lmaxd has to be set equal lmax ')
 9130 FORMAT (/)
 9140 FORMAT (20x,'full potential calculation ',
     +       '- cut off of non spherical potential',/,' >',/)
 9150 FORMAT (31x,'representive atom no.',i3,' irns :',i5,' irnsd :',i5)
 9160 FORMAT (21x,a43,/,21x,' using',i3,'-th. born approximation ')
 9170 FORMAT (21x,a43)
 9180 FORMAT (2i5)
 9190 FORMAT (3f12.7,/,4i4)
 9200 FORMAT (3i4,1f12.7)
 9210 FORMAT (' lmax'/,i4)
 9220 FORMAT ('          E1          E2          TK'/,3f12.6)
 9230 FORMAT ('   NPOL  NPNT1  NPNT2  NPNT3'/,4i7)
 9250 FORMAT ('  IFILE    IPE ISHIFT'/,3i7)
 9260 FORMAT (' KSHAPE    IRM    ICST'/,3i7)
 9270 FORMAT ('   KCOR  KVREL    KWS   KHFELD    KXC'/,5i7)
 9280 FORMAT (' external magnetic hfield     :',f15.4/,
     +        ' VCONST                       :',f15.6)
 9290 FORMAT ('   IMIX   '/,i7)
 9300 FORMAT (' ITDBRY'/,i7)
 9310 FORMAT ('      STRMIX        FCM       QBOUND'/,3f12.6)
 9320 FORMAT ('      BRYMIX'/,f12.6)
 9330 FORMAT ('    KTE   KPRE   KVMAD'/,3i7)
 9301 format(   3(1H-),1H+  ,75(1H-))
 9302 format( 3(11(1H-),1H+),43(1H-))
 9303 format(3(6(1H-),1H+) ,58(1H-))
 9304 format(4(6(1H-),1H+) ,51(1H-))
 9305 format(3(6(1H-),1H+),11(1H-),1H+ ,46(1H-))
 9306 format(6(6(1H-),1H+) ,37(1H-))
 9307 format(6(1H-),1H+,72(1H-))
 9308 format(11(1H-),1H+,67(1H-))
 9309 format(5(6(1H-),1H+) ,44(1H-))
C
      END
