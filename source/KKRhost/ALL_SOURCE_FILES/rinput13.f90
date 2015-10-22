      SUBROUTINE RINPUT13(ALAT,RBASIS,ABASIS,BBASIS,CBASIS,CLS,NCLS,& 
     &           EMIN,EMAX,TK,NPOL,NPNT1,NPNT2,NPNT3,&
     &           EBOTSEMI,EMUSEMI,TKSEMI,NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,&
     &           FSEMICORE,ESHIFT,&
     &           NSTEPS,IMIX,MIXING,QBOUND,FCM,ITDBRY,&
     &           IRNS,FPRADIUS,NTCELL,NAEZ,NEMB,KAOEZ,IRM,ZAT,&
     &           NINEQ,NREF,NREFD,&
     &           ICST,IFILE,IGF,INS,INSREF,IPE,IPF,IPFE,&
     &           KCOR,KEFG,KFROZN,KHFIELD,KHYP,KPRE,KSHAPE,KTE,&
     &           KFG,KVMAD,KVREL,KWS,KXC,LAMBDA_XC, &
     &           LMAX,LMMAX,LMPOT,LPOT,NATYP,NSPIN,&
     &           LMXC,TXC,ICC,REFPOT,&
     &           IPRCOR,IRNUMX,ISHIFT,ITCCOR,&
     &           INTERVX,INTERVY,INTERVZ,&
     &           HFIELD,COMPLX,&
     &           KMT,MTFAC,VBC,VCONST,LINIPOL,INIPOL,IXIPOL,LRHOSYM,&
     &           MMIN,MMAX,SINN,SOUT,RIN,ROUT,M2,I12,I13,I19,I25,I40,&
     &           NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,    &
     &           TLEFT,TRIGHT,LINTERFACE,RCUTZ,RCUTXY,RMTREF,RMTREFAT,&
     &           KFORCE,KMROT,QMTET,QMPHI,NCPA,ICPA,ITCPAMAX,CPATOL,&   
     &           NOQ,IQAT,CONC,SOLVER,SOCSCL,CSCL,KREL,SOCSCALE,&
     &           LOPT,UEFF,JEFF,EREFLDAU,KREADLDAU,&
     &           LMAXD,LPOTD,NSPIND,NAEZD,NATYPD,NEMBD,NPRINCD,&
     &           IRMD,IRNSD,NPAN_LOG,NPAN_EQ,NCHEB,R_LOG,IVSHIFT,&
     &           TOLRDIF,LLY,DELTAE,&
     &           LCARTESIAN,BRAVAIS,RMAX,GMAX)
      IMPLICIT NONE
!     ..
!     .. Parameters
      DOUBLE PRECISION CVLIGHT
      PARAMETER (CVLIGHT=274.0720442D0)
!     ..
!     .. Scalar arguments ..
      INTEGER  KREL,LMAXD,LPOTD,NSPIND,NAEZD,NATYPD,NEMBD,NPRINCD,IRMD,IRNSD,NREFD

!     .. Local Arrays ..
      CHARACTER*4 TSPIN(3)
      CHARACTER*8 TKWS(3)
      CHARACTER*43 TINS(0:3),TKCOR(0:3),TVREL(0:2)
      CHARACTER*2  SOCII(-2:-1)
! NOTE of VP : there should be some crosscheck of competing optons
!              e.g., XCPL and CONDUCT cannot be done simultaneously
!              neither SOC1 and SOC2 manipulation etc.
!     ..
!     .. External Subroutines ..
      EXTERNAL RCSTOP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MIN
!     ..
!     .. Array Arguments ..
      INTEGER IRNS(*),KFG(4,*),LMXC(*),NTCELL(*),CLS(*),REFPOT(*)
      INTEGER INIPOL(*),IXIPOL(*)
      DOUBLE PRECISION BRAVAIS(3,3)
      DOUBLE PRECISION ZAT(*),MTFAC(*),VBC(*),RBASIS(3,*),RMTREF(NREFD),RMTREFAT(NAEZD+NEMBD)
      DOUBLE PRECISION FPRADIUS(NATYPD)
      DOUBLE PRECISION TRIGHT(3,NEMBD+1),TLEFT(3,NEMBD+1)
      DOUBLE PRECISION ZPERLEFT(3),ZPERIGHT(3)
!     variables for spin-orbit/speed of light scaling
      DOUBLE PRECISION SOCSCL(KREL*LMAXD+1,KREL*NATYPD+(1-KREL))
      DOUBLE PRECISION SOCSCALE(NATYPD)
      DOUBLE PRECISION CSCL(KREL*LMAXD+1,KREL*NATYPD+(1-KREL))
      CHARACTER*24 TXC(4)
      CHARACTER*256 UIO  ! NCOLIO=256
      CHARACTER*10 SOLVER
      CHARACTER*40 I12,I13,I19,I25,I40
!     ..
!     .. Scalar Arguments ..
      DOUBLE PRECISION EMIN,EMAX,ESHIFT,FCM,HFIELD,MIXING,QBOUND,TK,&
     &       VCONST,ABASIS,BBASIS,CBASIS,RCUTZ,RCUTXY,LAMBDA_XC,TOLRDIF,RMAX,GMAX
      INTEGER ICC,ICST,IFILE,IGF,IMIX,INS,INSREF,&
     &        IPE,IPF,IPFE,IPOTOU,IPRCOR,&
     &        IRM,IRNUMX,ISHIFT,ITCCOR,ITDBRY,KCOR,&
     &        KEFG,KFROZN,KHFIELD,KHYP,KPRE,KSHAPE,KTE,KVMAD,&
     &        KVREL,KWS,KXC,LMAX,LMMAX,LMPOT,LPOT,KFORCE,&
     &        NATYP,NPNT1,NPNT2,NPNT3,NPOL,NSPIN,&
     &        NPAN_LOG,NPAN_EQ,NCHEB
      INTEGER NPOLSEMI,N1SEMI,N2SEMI,N3SEMI
      DOUBLE PRECISION FSEMICORE,EBOTSEMI,EMUSEMI,TKSEMI,R_LOG
      INTEGER NSTEPS,KMT,NAEZ,NEMB
      INTEGER NINEQ
      DOUBLE PRECISION ALAT
      INTEGER MMIN,MMAX,SINN,SOUT,RIN,ROUT
      INTEGER INTERVX,INTERVY,INTERVZ,NREF,NCLS
      INTEGER NLBASIS,NRBASIS,NLEFT,NRIGHT,NDIM      
      INTEGER LLY
      DOUBLE COMPLEX DELTAE  ! LLY Energy difference for numerical derivative
      LOGICAL LINIPOL,LRHOSYM,COMPLX,LINTERFACE,LCARTESIAN,LNEW,LATOMINFO
!----------------------------------------------------------------
!     CPA variables. Routine has been modified to look for
!     the token ATOMINFOC and only afterwards, if not found, for the
!     old token ATOMINFO. The only necessary extra information 
!     required is the site IQAT(IATOM) on which the atom IATOM 
!     is located and the occupancy (concentration) CONC(IATOM). 
!     The rest of CPA variables are deduced from these two. 
!     The tolerance for the CPA-cycle and the number of CPA iterations
!     can be modified adding the token <CPAINFO> in the input file.
!
      INTEGER NCPA,ICPA(NAEZD)   ! ncpa = 0/1 CPA flag
                                 ! icpa = 0/1 site-dependent CPA flag
      INTEGER ITCPAMAX           ! max. number of CPA iterations
      REAL*8  CPATOL             ! convergency tolerance for CPA-cycle
      INTEGER NOQ(NAEZD)         ! number of diff. atom types located
                                 ! on a given site
      INTEGER IQAT(NATYPD) ! the site on which an atom is located
      INTEGER KAOEZ(NATYPD,NAEZD+NEMBD) 
                                 ! atom types located at a given site
      REAL*8 CONC(NATYPD)        ! concentration of a given atom 

      CHARACTER*3 CPAFLAG(0:1)
      REAL*8 SUM
      INTEGER IO,IA,IQ,IPRINT
!-----------------------------------------------------------------------
!     Variables storing the magnetization direction information.
!     QMTET/QMPHI(NAEZD) give the angles to which the magnetic moment
!     on a given site is rotated against the z-axis. Default values
!     0.0 and 0.0, i.e., magnetic moment parallel to the z-axis.
!     The angles are read in after the token RBASISANG is found 
!     (sought in input file prior to RBASIS token)
!
!   *  KMROT                                                           *
!   *  0: no rotation of the magnetisation                             *
!   *  1: individual rotation of the magnetisation for every site      *
!   ( see also the routine < FINDGROUP > and ff)
!
      INTEGER KMROT
      REAL*8 QMTET(NAEZD),QMPHI(NAEZD)
! ----------------------------------------------------------------------
! LDA+U
      INTEGER LOPT(NATYPD),KREADLDAU
      DOUBLE PRECISION EREFLDAU(NATYPD),UEFF(NATYPD),JEFF(NATYPD)
! LDA+U
! ----------------------------------------------------------------------
! IVSHIFT test option
      INTEGER IVSHIFT
      LOGICAL TEST,OPT
      EXTERNAL TEST,OPT
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION BRYMIX,STRMIX,TX,TY,TZ,DVEC(10)
      INTEGER I,IL,J,IER,IER2,I1,II,IR,M2,IDOSEMICORE
      CHARACTER*43 TSHAPE
      DOUBLE PRECISION SOSCALE,CTLSCALE
      INTEGER IMANSOC(NATYPD),NASOC,ISP(NATYPD)
      LOGICAL MANSOC,MANCTL
!
      CHARACTER*8 TESTC(32),OPTC(32),VERSION
      COMMON /TESTC/TESTC
      COMMON /OPTC/OPTC
!     ..
!     .. Data statements ..
      DATA VERSION /'May 2015'/
!       DATA VERSION /'Feb 2005'/
      DATA TSPIN/'non-','    ','    '/
      DATA TSHAPE/' exact cell treatment (shape correction)  '/
      DATA TVREL/&
     &     ' non relativistic calculation              ',&
     &     ' s.r.a. calculation                        ',&
     &     ' fully relativistic calculation            '/
      DATA TKCOR/&
     &     ' frozen core approximation                 ',&
     &     ' core relaxation s.r.a.                    ',&
     &     ' core relaxation nonsra                    ',&
     &     ' core relaxation                           '/
      DATA TINS/' spherical averaged input potential        ',&
     &     ' non spherical input potential for cluster ',&
     &     ' non spherical input potential for cluster ',&
     &     ' non spherical input potential             '/
      DATA TKWS/' full mt','   ws   ',' full ws'/
!
      DATA CPAFLAG/' NO','YES'/
      DATA SOCII/'xy','zz'/ 
!     ..
!
!------------ array set up and definition of input parameter -----------
!
      TXC(1) = ' Morruzi,Janak,Williams '
      TXC(2) = ' von Barth,Hedin        '
      TXC(3) = ' Vosko,Wilk,Nusair      '
      TXC(4) = ' GGA PW91               '

      IPRINT = 0
      WRITE (1337,2004) VERSION

      OPEN(111,FILE='inputcard_generated.txt') ! Write out found or assumed values

      RMTREFAT(:) = -1.D0 ! Signals the need for later calculation
      RMTREF(:) = -1.D0
      NEMB = 0


!
! read RUNNING options
!
      CALL IoInput('RUNOPT          ',UIO,1,7,IER)
      IF (IER.NE.0) THEN 
         WRITE(111,*) 'RUNOPT not found'
      ELSE
         READ (UNIT=UIO,FMT=980)(OPTC(I),I=1,8)
         WRITE(111,FMT='(A6)') 'RUNOPT'
         WRITE(111,FMT=980)  (OPTC(I),I=1,8)
      ENDIF
!
! read TEST options
!
      CALL IoInput('TESTOPT         ',UIO,1,7,IER)
      IF (IER.NE.0) THEN 
         WRITE(111,*) 'TESTOPT not found'
      ELSE
         READ(UNIT=UIO,FMT=980)(TESTC(i),i=1,8)
         CALL IoInput('TESTOPT         ',UIO,2,7,IER)
         READ(UNIT=UIO,FMT=980)(TESTC(8+i),i=1,8)
         WRITE(111,FMT='(A7)') 'TESTOPT'
         WRITE(111,FMT=980)(TESTC(i),i=1,8)
         WRITE(111,FMT=980)(TESTC(8+i),i=1,8)
      ENDIF


!==========================================================
! Begin lattice structure definition


      CALL IoInput('ALATBASIS       ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) ALAT
         WRITE(111,*) 'ALATBASIS=',ALAT
      ELSE
         WRITE(111,*) 'ALATBASIS not found in inputcard'
         WRITE(*,*) 'rinput13: ALATBASIS not found in inputcard'
         STOP 'rinput13: ALATBASIS not found in inputcard'
      ENDIF
      ! move this writeout back to line 330 where A,B,CBASIS is set
      !WRITE(6,2019) ABASIS,BBASIS,CBASIS 
      !WRITE(6,2107)
      !WRITE(6,2014) ALAT


      ! Set 2-d or 3-d geometry
      LINTERFACE = .FALSE.
      CALL IoInput('INTERFACE       ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN 
         READ (UNIT=UIO,FMT=*) LINTERFACE
         WRITE(111,*) 'INTERFACE=',LINTERFACE
      ELSE
         WRITE(111,*) 'Default INTERFACE= ',LINTERFACE
      ENDIF

      NDIM = 3
      IF (LINTERFACE) NDIM = 2
      IF (.NOT.LINTERFACE.AND..NOT.OPT('SUPRCELL')) THEN 
         WRITE(*,*) '3D-calculation, adding run-option "full inv" for full inversion.'
         CALL ADDOPT('full inv')
      ENDIF


      WRITE(111,*) 'Bravais vectors in units of ALAT'
      BRAVAIS(1:3,1:3) = 0D0
      DO I = 1,NDIM
         CALL IOINPUT('BRAVAIS         ',UIO,I,7,IER)
         IF (IER.NE.0) STOP 'RINPUT: BRAVAIS NOT FOUND'
         READ (UNIT=UIO,FMT=*) (BRAVAIS(J,I),J=1,NDIM)
      END DO
      WRITE(111,FMT='(A7)') 'BRAVAIS'
      DO I = 1,NDIM
         WRITE(111,*) (BRAVAIS(J,I),J=1,NDIM)
      ENDDO


      CALL IoInput('NAEZ            ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) NAEZ
         WRITE(111,*) 'NAEZ=',NAEZ
      ELSE
         WRITE(111,*) 'NAEZ not found'
         STOP 'NAEZ not found in <RINPUT13>'
      ENDIF
      IF (NAEZ.GT.NAEZD) THEN
         WRITE(6,*) ' set NAEZD to at least ',NAEZ
         STOP ' in < RINPUT13 > '
      END IF


      LCARTESIAN = .FALSE.
      IER = 0
      CALL IOINPUT('CARTESIAN       ',UIO,1,7,IER)
      IF ( IER.EQ.0 ) THEN 
         READ (UNIT=UIO,FMT=*) LCARTESIAN
         WRITE(111,*) 'CARTESIAN= ',LCARTESIAN
      ELSE
         WRITE(111,*) 'Default CARTESIAN= ',LCARTESIAN
      ENDIF


      ! Basis atoms
      WRITE(111,FMT='(A16)') '<RBASIS>        '
      DO I=1,NAEZ
         CALL IoInput('<RBASIS>        ',UIO,I,7,IER)
         IF (IER.EQ.0) THEN
            READ (UNIT=UIO,FMT=*) (RBASIS(J,I), J=1,3)
            WRITE(111,FMT='(3E24.12)') (RBASIS(J,I), J=1,3)
         ELSE
            IER=0
            CALL IoInput('RBASIS          ',UIO,I,7,IER)
            IF (IER.EQ.0) THEN
               READ (UNIT=UIO,FMT=*) (RBASIS(J,I), J=1,3)
               WRITE(111,FMT='(3E24.12)') (RBASIS(J,I), J=1,3)
            ELSE
               WRITE(*,*) 'RINPUT13: Keyword <RBASIS> or RBASIS not found. Stopping.'
               STOP 'RINPUT13: RBASIS'
            ENDIF
         ENDIF
      ENDDO                         ! I=1,NAEZ
      CALL IDREALS(RBASIS(1,1),3*NAEZ,IPRINT)



      DVEC(1:3) = 1.D0
      CALL IoInput('BASISCALE       ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) (DVEC(I),I=1,3)
         WRITE(111,FMT='(A10,3E12.4)') 'BASISCALE=',DVEC(1:3)
      ELSE
         WRITE(111,FMT='(A18,3E12.4)') 'Default BASISCALE=',DVEC(1:3)
      ENDIF

      CALL IDREALS(DVEC(1),3,IPRINT)
      ABASIS = DVEC(1)
      BBASIS = DVEC(2) 
      CBASIS = DVEC(3)

      WRITE(1337,2019) ABASIS,BBASIS,CBASIS 
      WRITE(1337,2107)
      WRITE(1337,2014) ALAT

      !----------------------------------------------------------------------
      ! Begin read left- and right-host information in 2d-case.
      ! Set up the embeding positions

      IF (LINTERFACE) THEN

         WRITE(1337,9410)

         NRIGHT = 10
         CALL IoInput('NRIGHTHO        ',UIO,1,7,IER)
         IF (IER.EQ.0) THEN
            READ (UNIT=UIO,FMT=*) NRIGHT
            WRITE(111,*) 'NRIGHTHO=',NRIGHT
         ELSE
            WRITE(111,*) 'Default NRIGHTHO=',NRIGHT
         ENDIF

         NLEFT = 10
         CALL IoInput('NLEFTHOS        ',UIO,1,7,IER)
         IF (IER.EQ.0) THEN
            READ (UNIT=UIO,FMT=*) NLEFT
            WRITE(111,*) 'NLEFTHOS=',NLEFT 
         ELSE
            WRITE(111,*) 'Default NLEFTHOS=',NLEFT
         ENDIF

         CALL IoInput('<NLBASIS>       ',UIO,1,7,IER)
         IF (IER.NE.0) THEN
            WRITE(1337,*) 'rinput13: <NLBASIS> not found in inputcard'
            IER = 0
            CALL IoInput('NLBASIS         ',UIO,1,7,IER)
            IF (IER.NE.0) THEN
               WRITE(*,*) 'rinput13: NLBASIS also not found in inputcard'
               STOP 'rinput13: NLBASIS not found in inputcard'
            ENDIF
         ENDIF
         IF (IER.EQ.0) THEN
            READ (UNIT=UIO,FMT=*) NLBASIS
            WRITE(111,*) '<NLBASIS>=',NLBASIS
         ENDIF

         CALL IoInput('<NRBASIS>       ',UIO,1,7,IER)
         IF (IER.NE.0) THEN
            WRITE(1337,*) 'rinput13: <NRBASIS> not found in inputcard'
            IER = 0
            CALL IoInput('NRBASIS         ',UIO,1,7,IER)
            IF (IER.NE.0) THEN
               WRITE(*,*) 'rinput13: NRBASIS also not found in inputcard'
               STOP 'rinput13: NRBASIS not found in inputcard'
            ENDIF
         ENDIF
         IF (IER.EQ.0) THEN
            READ (UNIT=UIO,FMT=*) NRBASIS
            WRITE(111,*) '<NRBASIS>=',NRBASIS
         ENDIF


         NEMB = NLBASIS + NRBASIS
         WRITE(1337,*) 'Number of embedded atoms NEMB=NLBASIS + NRBASIS=',NEMB
         IF(NEMB.GT.NEMBD) THEN
           write(6,*) 'Please, increase the parameter nembd (',nembd,') in inc.p to',nemb
            STOP 'ERROR in NEMBD.'
         ENDIF


         IER = 0
         ! Check if the keywords exist for old/new treatment of left and right host
         CALL IoInput('LEFTBASIS       ',UIO,1,7,IER)
         IF (IER.EQ.0) THEN
            LNEW = .FALSE.
         ELSE
            LNEW = .TRUE.
            IER = 0
            CALL IoInput('<RBLEFT>        ',UIO,1,7,IER)
         ENDIF
         IF (IER.NE.0) THEN
            WRITE(*,*) 'rinput13: LEFTBASIS or <RBLEFT> not found in inputcard'
            STOP 'rinput13: LEFTBASIS or <RBLEFT> not found in inputcard'
         ENDIF
         IER = 0
         CALL IoInput('RIGHBASIS       ',UIO,1,7,IER)
         IF (IER.EQ.0) THEN
            LNEW = .FALSE.
         ELSE
            LNEW = .TRUE.
            IER = 0
            CALL IoInput('<RBRIGHT>       ',UIO,1,7,IER)
         ENDIF
         IF (IER.NE.0) THEN
            WRITE(*,*) 'rinput13: RIGHBASIS or <RBRIGHT> not found in inputcard'
            STOP 'rinput13: RIGHBASIS or <RBRIGHT> not found in inputcard'
         ENDIF


         ! In leftbasis and rightbasis, kaoez is used only in decimation case.
         ! Then it indicates the correspondence of the atom-coordinate given
         ! by leftbasis and rightbasis to the left- and right-host t-matrix read in 
         ! by decimaread. For the slab case, kaoez is not used in the embedded positions.
         IF (LNEW) THEN

            WRITE(111,FMT='(A82)') '<RBLEFT>                                                      <RMTREFL>   <KAOEZL>'
            DO I=1,NLBASIS
               CALL IoInput('<RBLEFT>        ',UIO,I,7,IER)
               READ (UNIT=UIO,FMT=*) (TLEFT(I1,I),I1=1,3)
               KAOEZ(1,NAEZ+I) = I            ! Default
               CALL IoInput('<KAOEZL>        ',UIO,I,7,IER)
               IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) KAOEZ(1,NAEZ+I)
               CALL IoInput('<RMTREFL>       ',UIO,I,7,IER)
               IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) RMTREFAT(NAEZ+I)
               WRITE (111,FMT='(3E20.12,3X,F9.6,3X,I5)') (TLEFT(I1,I),I1=1,3),RMTREFAT(NAEZ+I),KAOEZ(1,NAEZ+I)
            ENDDO
            WRITE(111,FMT='(A82)') '<RBRIGHT>                                                     <RMTREFR>   <KAOEZL>'
            DO I=1,NRBASIS
               CALL IoInput('<RBRIGHT>       ',UIO,I,7,IER)
               READ (UNIT=UIO,FMT=*) (TRIGHT(I1,I),I1=1,3)
               KAOEZ(1,NAEZ+NLBASIS+I) = I     ! Default
               CALL IoInput('<KAOEZR>        ',UIO,I,7,IER)
               IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) KAOEZ(1,NAEZ+NLBASIS+I)
               CALL IoInput('<RMTREFR>       ',UIO,I,7,IER)
               IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) RMTREFAT(NAEZ+NLBASIS+I)
               WRITE (111,FMT='(3E20.12,3X,F9.6,3X,I5)') (TRIGHT(I1,I),I1=1,3),RMTREFAT(NAEZ+NLBASIS+I),KAOEZ(1,NAEZ+NLBASIS+I)
            ENDDO

         ELSE ! (LNEW) now old-style input

            DO I=1,NLBASIS
               CALL IoInput('LEFTBASIS       ',UIO,I,7,IER)
               READ (UNIT=UIO,FMT=*) (TLEFT(I1,I),I1=1,3),II,IR
               KAOEZ(1,NAEZ+I) = II    ! changed 1.11.99
               REFPOT(NAEZ+I) = IR    
            END DO
            DO I=1,NRBASIS
               CALL IoInput('RIGHBASIS       ',UIO,I,7,IER)
               READ (UNIT=UIO,FMT=*) (TRIGHT(I1,I),I1=1,3),II,IR
               KAOEZ(1,NAEZ+NLBASIS+I) = II  ! changed 1.11.99
               REFPOT(NAEZ+NLBASIS+I) = IR  
            END DO
         ENDIF

         CALL IDREALS(TLEFT,3*(NEMBD+1),IPRINT)
         CALL IDREALS(TRIGHT,3*(NEMBD+1),IPRINT)


         ! Put The additional atoms in the "embeding" positions

         DO I=1,NLBASIS
            RBASIS(1:3,NAEZ+I) = TLEFT(1:3,I)
         END DO
         DO I=1,NRBASIS
            RBASIS(1:3,NAEZ+NLBASIS+I) = TRIGHT(1:3,I)
         END DO 
         !cccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! In RBASIS we have first the basis atoms or the interface
         ! atoms then the left host then the right host the host
         ! goes in the NEMB positions 
         !
         ! IN CASE OF CPA the host is treated as an effective 
         ! CPA medium, that is, there is only one kind of atom
         ! occupying a crystallographic site.
         !
         !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
         CALL IoInput('ZPERIODL        ',UIO,1,7,IER)
         IF (IER.NE.0) THEN
            WRITE(*,*) 'rimput13: ZPERIODL not found in inputcard'
            STOP 'rimput13: ZPERIODL not found in inputcard'
         ELSE
            READ (UNIT=UIO,FMT=*) (ZPERLEFT(I1),I1=1,3)
            WRITE(111,FMT='(A9,3E20.12)') 'ZPERIODL=',(ZPERLEFT(I1),I1=1,3)
         ENDIF
         CALL IDREALS(ZPERLEFT(1),3,IPRINT)

         CALL IoInput('ZPERIODR        ',UIO,1,7,IER)
         IF (IER.NE.0) THEN
            WRITE(*,*) 'rimput13: ZPERIODR not found in inputcard'
            STOP 'rimput13: ZPERIODR not found in inputcard'
         ELSE
            READ (UNIT=UIO,FMT=*) (ZPERIGHT(I1),I1=1,3)
            WRITE(111,FMT='(A9,3E20.12)') 'ZPERIODR=',(ZPERIGHT(I1),I1=1,3)
         ENDIF
         CALL IDREALS(ZPERIGHT(1),3,IPRINT)

         WRITE(1337,9430) NLEFT,NLBASIS
         WRITE(1337,9440) NRIGHT,NRBASIS
         WRITE(1337,9450) (ZPERLEFT(i1),I1=1,3)
         WRITE(1337,9460) (ZPERIGHT(i1),I1=1,3)
         WRITE(1337,9465)
         WRITE(1337,9470)
         DO I=NLEFT,1,-1
            DO I1=NLBASIS,1,-1
            tx = TLEFT(1,i1) + (I-1)*ZPERLEFT(1)
            ty = TLEFT(2,i1) + (I-1)*ZPERLEFT(2)
            tz = TLEFT(3,i1) + (I-1)*ZPERLEFT(3)
            WRITE(1337,9420) (I-1)*NLBASIS+i1, tx,ty,tz,KAOEZ(1,I1)
            END DO 
         END DO
          WRITE(1337,9475)
         DO I=1,NAEZ
            WRITE(1337,9420) I, (RBASIS(I1,I),I1=1,3)
         END DO
          WRITE(1337,9480)
          DO I=1,NRIGHT
            DO I1=1,NRBASIS
            tx = TRIGHT(1,i1) + (I-1)*ZPERIGHT(1)
            ty = TRIGHT(2,i1) + (I-1)*ZPERIGHT(2)
            tz = TRIGHT(3,i1) + (I-1)*ZPERIGHT(3) 
            WRITE(1337,9420) (I-1)*NRBASIS+i1,tx,ty,tz,KAOEZ(1,I1)
            END DO 
         END DO  

      END IF                  ! LINTERFACE      
      ! End read left- and right-host information in 2d-case.
      !----------------------------------------------------------------------



! End lattice structure definition
!==========================================================
!==========================================================
! Begin atom type information


      ! Different atoms
      NATYP = NAEZ
      CALL IoInput('NATYP           ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) NATYP     
         WRITE(111,*) 'NATYP=',NATYP
      ELSE
         WRITE(111,*) 'Default NATYP= ',NAEZ
      ENDIF
      IF (NATYP.GT.NATYPD) THEN
         WRITE(6,*) 'RINPUT13: NATYP > NATYPD',NATYP,NATYPD
         STOP ' IN < RINPUT13 > '
      END IF
      IF (NATYP.LT.NAEZ) THEN
         WRITE(6,*) 'RINPUT13: NATYP < NAEZ ',NATYP,NAEZ
         STOP ' IN < RINPUT13 > '
      END IF


      ! although NSPIND is fixed to 1 in REL mode,
      ! NSPIN should be used as 1 or 2 at this stage 
      ! to indicate a non- or spin-polarised potential 
      ! that has to be read in. NSPIN is set to 1 before
      ! being passed to the subsequent programs.
      ! < TESTDIM > has been accordingly modified
      CALL IoInput('NSPIN           ',UIO,1,7,IER)
      IF (IER.NE.0) THEN
         WRITE(111,*) 'NSPIN not found'
         STOP 'NSPIN not found'
      ELSE
         READ (UNIT=UIO,FMT=*) NSPIN
         WRITE(111,*) 'NSPIN=',NSPIN
      ENDIF          

      WRITE(1337,2010) NSPIN
      WRITE(1337,2104)

      ! Atomic number
      ZAT(1:NATYP) = -1.D0  ! Negative value signals read-in from pot-file
      CALL IoInput('<ZATOM>         ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         WRITE(111,'(A10)') '<ZATOM>   '
         DO I = 1,NATYP
            CALL IoInput('<ZATOM>         ',UIO,I,7,IER)
            IF (IER.EQ.0) THEN
               READ (UNIT=UIO,FMT=*) ZAT(I)
               WRITE(111,FMT='(F6.3)') ZAT(I)
            ENDIF
         ENDDO
      ELSE
         WRITE(111,*) 'zatom will be read in from pot-file'
      ENDIF


      DO I=1,NAEZ
         ICPA(I) = 0
         NOQ(I) = 1
      END DO
      NCPA = 0


      DO I = 1,NAEZ
         KAOEZ(1,I) = I       ! default
         IQAT(I) = I          ! Basis-Site of atom I
      ENDDO
      IF (NATYP.EQ.NAEZ) CONC(1:NATYP) = 1.D0

      ! CPA calculation, read concentrations
      IF (NATYP.GT.NAEZ) THEN 

         NCPA = 1
         NOQ(1:NAEZ) = 0 ! re-initialize

         IER = 0
         IER2 = 0
         CALL IoInput('<SITE>          ',UIO,1,7,IER)
         CALL IoInput('<CPA-CONC>      ',UIO,1,7,IER2)
         IF (IER.NE.0.OR.IER2.NE.0) THEN
            WRITE(1337,*) '<SITE> or <CPA-CONC> not found, will search for ATOMINFOC'
         ELSE


            WRITE(111,FMT='(A18)') '<SITE>  <CPA-CONC>'
            DO I = 1,NATYP
               CALL IoInput('<SITE>          ',UIO,I,7,IER)
               READ (UNIT=UIO,FMT=*) IQAT(I)
               CALL IoInput('<CPA-CONC>      ',UIO,I,7,IER)
               READ (UNIT=UIO,FMT=*) CONC(I)
               WRITE(111,FMT='(I5,4X,E16.8)') IQAT(I),CONC(I)
            ENDDO
            
            DO I = 1,NATYP
               IQ = IQAT(I)
               NOQ(IQ) = NOQ(IQ) + 1
               IF ( NOQ(IQ) .GT. 1 ) ICPA(IQ) = 1
               KAOEZ(NOQ(IQ),IQ) = I
            ENDDO

            DO IQ=1,NAEZ
               SUM = 0D0
               IF (NOQ(IQ).LT.1) THEN
                  WRITE(6,*) 'RINPUT13: CPA: SITE',IQ,'HAS NO ASSIGNED ATOM'
                  STOP 'RINPUT13: CPA'
               ENDIF
               DO IO=1,NOQ(IQ)
                  SUM = SUM + CONC(KAOEZ(IO,IQ))
               END DO
               IF ( ABS(SUM-1.D0).GT.1D-6) THEN
                  WRITE(6,*) ' SITE ', IQ, ' CONCENTRATION <> 1.0 !'
                  WRITE(6,*) ' CHECK YOUR <ATOMINFO-CPA> INPUT '
                  STOP       ' IN <RINPUT99>'
               END IF
            END DO

         ENDIF


      ENDIF ! (NATYP.GT.NAEZ)



! End atom type information
!==========================================================
!==========================================================
! Begin relativistic treatment information

      KCOR = 2
!      CALL IoInput('KCOR      ',UIO,1,7,IER)
!                      READ (UNIT=UIO,FMT=*) kcor

      KVREL = 1   ! 0=Schroedinger / 1=SRA / 2=Dirac
      CALL IoInput('KVREL           ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) kvrel
         WRITE(111,*) 'KVREL= ',KVREL
      ELSE
         WRITE(111,*) 'Default KVREL= ',KVREL
      ENDIF


      IF (OPT('NEWSOSOL')) THEN ! Spin-orbit
         IF ( OPT('NEWSOSOL') .AND. (NSPIN.NE.2) ) STOP ' set NSPIN = 2 for SOC solver in inputcard'
         NPAN_LOG = 30
         NPAN_EQ = 30
         NCHEB = 10
         R_LOG = 0.1D0
         CALL IoInput('NPAN_LOG        ',UIO,1,7,IER)
         IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) NPAN_LOG
         CALL IoInput('NPAN_EQ         ',UIO,1,7,IER)
         IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) NPAN_EQ
         CALL IoInput('NCHEB           ',UIO,1,7,IER)
         IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) NCHEB
         CALL IoInput('R_LOG           ',UIO,1,7,IER)
         IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) R_LOG
         WRITE(111,*) 'NPAN_LOG= ',NPAN_LOG
         WRITE(111,*) 'NPAN_EQ= ',NPAN_EQ
         WRITE(111,*) 'NCHEB= ',NCHEB
         WRITE(111,*) 'R_LOG= ',R_LOG
      ENDIF

      SOCSCALE(1:NATYPD) = 1.D0   ! Spin-orbit scaling
      CALL IoInput('<SOCSCL>        ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) (SOCSCALE(I1),I1=1,NATYP)
         WRITE(111,FMT='(A10,50E10.2)') '<SOCSCL>= ',(SOCSCALE(I1),I1=1,NATYP)
      ELSE
         WRITE(111,FMT='(A18,50E10.2)') 'Default <SOCSCL>= ',(SOCSCALE(I1),I1=1,NATYP)
      ENDIF



! End relativistic treatment information
!==========================================================
!==========================================================
! Begin cell control




      IRNS(1:NATYP) = -1        ! Negative value signals to use FPRADIUS
      FPRADIUS(1:NATYP) = -1.D0 ! Negative value signals to use IRNS from pot-file (sub. startb1)
      CALL IoInput('<FPRADIUS>      ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         WRITE(111,'(A10)') '<FPRADIUS>'
         DO I = 1,NATYP
            CALL IoInput('<FPRADIUS>      ',UIO,I,7,IER)
            IF (IER.EQ.0) THEN
               READ (UNIT=UIO,FMT=*) FPRADIUS(I)
            ENDIF
            WRITE(111,FMT='(F6.3)') FPRADIUS(I)
         ENDDO
      ELSE
         WRITE(111,*) 'fpradius will be read in from pot-file'
      ENDIF

!

      INS = 1
      CALL IoInput('INS             ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) ins
         WRITE(111,*) 'INS= ',INS
      ELSE
         WRITE(111,*) 'Default INS= ',INS
      ENDIF


      KSHAPE = 2
      IF (INS.EQ.0) KSHAPE = 0
      CALL IoInput('KSHAPE          ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) KSHAPE
         WRITE(111,*) 'KSHAPE= ',KSHAPE
      ELSE
         WRITE(111,*) 'Default KSHAPE= ',KSHAPE
      ENDIF

      IF ( (KREL.EQ.1).AND.(KSHAPE.NE.0) ) THEN
         WRITE(1337,*) ' WARNING : KSHAPE set to ZERO for REL case'
         WRITE(111,*) ' WARNING : kshape set to ZERO for REL case'
         KSHAPE = 0
      END IF


      
      ! Read cell information
      WRITE(1337,*) 'Cell information <SHAPE>:'
      WRITE(111,FMT='(A16)') '<SHAPE>         '
      DO I = 1,NATYP
         NTCELL(I) = IQAT(I) ! Default: Different shape function per atom
         CALL IoInput('<SHAPE>         ',UIO,I,7,IER)
         IF (IER.EQ.0) THEN
            READ (UNIT=UIO,FMT=*) NTCELL(I)
            WRITE(111,FMT='(I6)') NTCELL(I)
         ENDIF
      ENDDO
      
! End cell control
!==========================================================
!==========================================================
! Begin exchange correlation treatment information

      KXC = 2 ! 0=vBH 1=MJW 2=VWN 3=PW91
      CALL IoInput('KEXCOR          ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) kxc    
         WRITE(111,*) 'KEXCOR= ',KXC
      ELSE
         WRITE(111,*) 'Default KEXCOR= ',KXC
      ENDIF

      ! Scale magnetic moment (0 < Lambda_XC < 1,  0=zero moment, 1= full moment)
      LAMBDA_XC = 1.D0 
      CALL IoInput('LAMBDA_XC       ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) LAMBDA_XC
         WRITE(111,*) 'LAMBDA_XC= ',LAMBDA_XC
      ELSE
         WRITE(111,*) 'Default LAMBDA_XC= ',LAMBDA_XC
      ENDIF

      ! LDA+U treatment
      ! -> Initialise UEFF,JEFF,LOPT,EREFLDAU for all atoms
      LOPT(1:NATYP) = -1        !  not perform lda+u (default)
      UEFF(1:NATYP) = 0.D0
      JEFF(1:NATYP) = 0.D0
      EREFLDAU(1:NATYP) = 0.5D0
      IF (OPT('LDA+U   ')) THEN


         !Check for LDA+U consistency -- if INS=0 suppress it
         IF ((INS.EQ.0) ) THEN
            WRITE (1337,*)
            WRITE (1337,*)&
                 &        ' WARNING: LDA+U should be used only in NON-SPHERICAL',&
                 &        ' case (INS=1) '
            WRITE (1337,*) ' Running option LDA+U will be ignored'
            WRITE (1337,*)
            DO I=1,32
               IF (OPTC(I)(1:8).EQ.'LDA+U   ') OPTC(I)='        '
            END DO
         END IF

         ! -> get number of atoms for lda+u:

         IER = 0
         CALL IoInput('NAT_LDAU        ',UIO,1,7,IER)
         IF ( IER.NE.0 ) THEN 
            NASOC = NATYP
         ELSE 
            READ (UNIT=UIO,FMT=*) NASOC
            IF ( NASOC.GT.NATYP ) STOP ' main0: NAT_LDAU > NATYP'
         END IF

         ! -> read in UEFF,JEFF,LOPT,EREFLDAU for the desired atoms

         IL = 0
         DO I=1,NASOC
            IER = 0
            CALL IoInput('LDAU_PARA       ',UIO,I,7,IER)
            IF ( IER.EQ.0 ) THEN
               READ (UNIT=UIO,FMT=*) I1,LOPT(I1),UEFF(I1),JEFF(I1),EREFLDAU(I1)
               IL = IL + 1
            END IF
         ENDDO
         IF ( IL.NE.NASOC ) THEN
            WRITE(6,*) ' ERROR: LDA+U invoked for ',NASOC,' atoms'
            WRITE(6,*) '        Some (all) parameters are missing in the input-file'
            STOP 
         END IF
         KREADLDAU = 0
         IER = 0
         CALL IoInput('KREADLDAU       ',UIO,1,7,IER)
         IF ( IER.EQ.0 ) READ (UNIT=UIO,FMT=*) KREADLDAU


      END IF

! End exchange correlation treatment information
!==========================================================



!==========================================================
! Begin external field control


      KHFIELD = 0
      HFIELD = 0.D0
      CALL IoInput('HFIELD          ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) HFIELD
         IF (HFIELD.NE.0.D0) KHFIELD = 1
         WRITE(111,*) 'HFIELD= ',HFIELD
      ELSE
         WRITE(111,*) 'Default HFIELD= ',HFIELD
      ENDIF

      VCONST = 0.D0
      CALL IoInput('VCONST          ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) vconst
         WRITE(111,*) 'VCONST= ',VCONST
      ELSE
         WRITE(111,*) 'Default VCONST= ',VCONST
      ENDIF


      IF (TEST('atptshft')) THEN
        write(1337,*) 'READ IN IVSHIFT'
        CALL IoInput('IVSHIFT         ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) ivshift
      ENDIF


      ! Initial polarization
      LINIPOL = .FALSE.
      CALL IoInput('LINIPOL         ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) LINIPOL
         WRITE(111,*) 'LINIPOL= ',LINIPOL
      ELSE
         WRITE(111,*) 'Default: LINIPOL= ',LINIPOL
      ENDIF

      INIPOL(1:NATYPD) = 0
      IF (LINIPOL) THEN
      INIPOL(1:NATYPD) = 1
        CALL IoInput('XINIPOL         ',UIO,1,7,IER)
        IF (IER.EQ.0) THEN
           READ (UNIT=UIO,FMT=*) (inipol(I),I=1,natyp) 
           WRITE(111,FMT='(A10,80I2)') 'XINIPOL=  ',(INIPOL(I),I=1,NATYP)
        ELSE
           WRITE(111,FMT='(A18,80I2)') 'Default XINIPOL=  ',(INIPOL(I),I=1,NATYP)
        ENDIF
      ENDIF


      WRITE (1337,2021) (INIPOL(I),I=1,NATYP) 
      WRITE(1337,2103)




! End external field control
!==========================================================
!==========================================================
! Begin Green function calculation control (diag./non-diag)


      IGF = 0
      CALL IoInput('IGREENFUN       ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) igf
         WRITE(111,*) 'IGREENFUN= ',IGF
      ELSE
         WRITE(111,*) 'Default IGREENFUN= ',IGF
      ENDIF
      IF (OPT('KKRFLEX ')) IGF = 1

      ICC = 0
      CALL IoInput('ICC             ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) ICC
         WRITE(111,*) 'ICC= ',ICC
      ELSE
         WRITE(111,*) 'Default ICC= ',ICC
      ENDIF
      IF (OPT('KKRFLEX ')) ICC = 1
      IF ( ( OPT('XCPL    ') ).OR.( OPT('CONDUCT ') ) ) ICC = -1

      IF ( ICC.NE.0 .AND. IGF.EQ.0 ) IGF = 1
      IF ( ICC.EQ.0 .AND. IGF.NE.0 ) ICC = -1


! End Green function calculation control (diag./non-diag)
!==========================================================
! =============================================================
! Begin accuracy parameters

      ! Angular momentum cutoff
      CALL IoInput('LMAX            ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) LMAX
         WRITE(111,*) 'LMAX=',LMAX
      ELSE
         STOP 'LMAX not found'
      ENDIF


      ! Brilloun zone mesh
      INTERVX = 10
      INTERVY = 10
      INTERVZ = 10
      CALL IoInput('BZDIVIDE        ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) INTERVX,INTERVY,INTERVZ
         WRITE(111,FMT='(A9,3I5)') 'BZDIVIDE=',INTERVX,INTERVY,INTERVZ
      ELSE
          WRITE(111,FMT='(A17,3I5)') 'Default BZDIVIDE=',INTERVX,INTERVY,INTERVZ
      ENDIF
      WRITE(1337,2104)
      WRITE(1337,2015) INTERVX,INTERVY,INTERVZ 
      WRITE(1337,2102)


      ! Energy contour
      NPOL = 7
!      IF (OPT('dos     ').OR.OPT('DOS     ')) NPOL = 0
      CALL IoInput('NPOL            ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) NPOL
         WRITE(111,*) 'NPOL=',NPOL
      ELSE
         WRITE(111,*) 'Default NPOL=',NPOL
      ENDIF

      CALL IoInput('EMIN            ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) EMIN
         WRITE(111,*) 'EMIN= ',EMIN
      ELSE IF (NPOL.EQ.0) THEN
         EMIN = -1.D0
         WRITE(111,*) 'Default for DOS: EMIN= ',EMIN
      ELSE
         WRITE(1337,*) 'Error in rinput13: EMIN not found'
         WRITE(111,*) 'Error in rinput13: EMIN not found'
         STOP 'Error in rinput13: EMIN not found'
      ENDIF


      EMAX = 1.D0
      CALL IoInput('EMAX            ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) EMAX
         WRITE(111,*) ' EMAX=',EMAX
      ELSE
         WRITE(111,*) 'Default  EMAX=',EMAX
      ENDIF

      TK = 800.D0
      CALL IoInput('TEMPR           ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) TK
         WRITE(111,*) 'TEMPR=',TK
      ELSE
         WRITE(111,*) 'Default TEMPR=',TK
      ENDIF

      NPNT1 = 3
      IF (NPOL.EQ.0) NPNT1 = 0
      CALL IoInput('NPT1            ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) NPNT1
         WRITE(111,*) ' NPT1=',NPNT1
      ELSE
         WRITE(111,*) 'Default  NPT1=',NPNT1
      ENDIF


      NPNT2 = NINT((EMAX-EMIN)*20.D0) ! 20 pts/Ryd
      IF (NPOL.EQ.0) NPNT2 = NINT((EMAX-EMIN)*100.D0) ! For dos, 100 pts/Ryd
      CALL IoInput('NPT2            ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) NPNT2
         WRITE(111,*) ' NPT2=',NPNT2
      ELSE
         WRITE(111,*) 'Default  NPT2=',NPNT2
      ENDIF

      NPNT3 = 3
      IF (NPOL.EQ.0) NPNT3 = 0
      CALL IoInput('NPT3            ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) NPNT3
         WRITE(111,*) ' NPT3=',NPNT3
      ELSE
         WRITE(111,*) 'Default  NPT3=',NPNT3
      ENDIF



      ! -> semicore 
      ! initialise variables
      IDOSEMICORE = 0
      EBOTSEMI = EMIN
      EMUSEMI = EBOTSEMI
      NPOLSEMI = 0
      N1SEMI = 0
      N2SEMI = 0
      N3SEMI = 0
      FSEMICORE = 1.D0

      IER = 0
      IF ( OPT('SEMICORE') ) THEN
         CALL IoInput('EBOTSEMI        ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                       READ (UNIT=UIO,FMT=*) EBOTSEMI
         CALL IoInput('EMUSEMI         ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                      READ (UNIT=UIO,FMT=*) EMUSEMI

         ! -> EMUSEMI < EBOT
         IF ( EMUSEMI.GE.EMIN ) GOTO 99800
         CALL IoInput('TKSEMI          ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                      READ (UNIT=UIO,FMT=*) TKSEMI

         CALL IoInput('NPOLSEMI        ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                      READ (UNIT=UIO,FMT=*) NPOLSEMI
         CALL IoInput('N1SEMI          ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                      READ (UNIT=UIO,FMT=*) N1SEMI
         CALL IoInput('N2SEMI          ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                      READ (UNIT=UIO,FMT=*) N2SEMI
         CALL IoInput('N3SEMI          ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                      READ (UNIT=UIO,FMT=*) N3SEMI
         CALL IoInput('FSEMICORE       ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                      READ (UNIT=UIO,FMT=*) FSEMICORE
         IDOSEMICORE = 1
99800    CONTINUE
         IF ( IDOSEMICORE.EQ.0 ) THEN 
            WRITE (1337,*)
            WRITE (1337,*) ' WARNING: SEMICORE used',&
     &           ' with incomplete/incorrect contour description'
            WRITE (1337,*) ' Running option SEMICORE will be ignored'
            WRITE (111,*)
            WRITE (111,*) ' WARNING: SEMICORE used',&
     &           ' with incomplete/incorrect contour description'
            WRITE (111,*) ' Running option SEMICORE will be ignored'
            DO I=1,32
               IF (OPTC(I)(1:8).EQ.'SEMICORE') OPTC(I)='        '
            END DO
         END IF
      END IF


      ! CPA convergence parameters
      CPATOL = 1D-4
      ITCPAMAX = 20
      CALL IOINPUT('CPAINFO         ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) CPATOL,ITCPAMAX
      ELSE
         WRITE(111,*) 'Default cpainfo:'
      ENDIF
      WRITE(111,FMT='(A7)') 'CPAINFO'
      WRITE(111,FMT='(E12.4,I5)') CPATOL,ITCPAMAX


      !==========================================================
      ! Begin screening cluster information

      RCUTZ = 11.D0/ALAT  ! Default 11 Bohr radii
      CALL IoInput('RCLUSTZ         ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) RCUTZ
         WRITE(111,*) 'RCLUSTZ=',RCUTZ
      ELSE
         WRITE(111,*) 'Default RCLUSTZ=',RCUTZ
      ENDIF


      RCUTXY = RCUTZ
      CALL IoInput('RCLUSTXY        ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) RCUTXY
         WRITE(111,*) 'RCLUSTXY=',RCUTXY
      ELSE
         WRITE(111,*) 'Default RCLUSTXY=',RCUTXY
      ENDIF

      WRITE(1337,*) 'Parameters used for the cluster calculation'
      if (abs(rcutz-rcutxy).lt.1.d-4) then
      write(1337,*) 'Clusters inside spheres with radius R = ',rcutz
      else
      write(1337,*) 'Clusters inside cylinders with '
      write(1337,*) 'Rz = ',rcutz,' Rxy = ',rcutxy
      end if
      write(1337,2104)
      write(1337,2018)                 ! rbasis
      write(1337,2101)
      do i=1,naez
        write(1337,2025) i,(rbasis(j,i),j=1,3),&
     &       QMTET(I),QMPHI(I),ICPA(I),NOQ(I),(KAOEZ(J,I),J=1,NOQ(I))
      enddo         





      DO I=1,NAEZ
         CALL IoInput('<RMTREF>        ',UIO,I,7,IER)
         IF (IER.EQ.0) THEN
            READ (UNIT=UIO,FMT=*) RMTREFAT(I)
         ENDIF
      ENDDO
      IF (IER.EQ.0) THEN
         WRITE(111,FMT='(A18)') '        <RMTREF>  '
      ELSE
         WRITE(111,FMT='(A18)') 'Default <RMTREF>  '
      ENDIF
      DO I=1,NAEZ
         WRITE(111,FMT='(9X,F9.6)') RMTREFAT(I)
      ENDDO


      ! End screening cluster information
      !==========================================================


      ! Number of Born iterations
      ICST = 2
      CALL IoInput('ICST            ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) ICST
         WRITE(111,*) 'ICST=',ICST
      ELSE
         WRITE(111,*) 'Default ICST=',ICST
      ENDIF


      ! Usage of Lloyd's formula
      LLY = 0 ! LLY Default=0 : do not apply Lloyds formula
      IF (OPT('LLOYD   ').OR.OPT('Lloyd   ').OR.OPT('lloyd   ')) LLY = 1
      CALL IoInput('<LLOYD>         ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) LLY
         WRITE(111,*) '<LLOYD>=',LLY
      ELSE
         WRITE(111,*) 'Default <LLOYD>=',LLY  
      ENDIF
      IF (LLY.NE.0) WRITE(1337,*) 'Applying Lloyds formula, LLY=',LLY

      DELTAE = (1.D-5,0.D0) ! Difference for numer. derivative in Lloyds formula
      CALL IoInput('<DELTAE>        ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) DELTAE 
         WRITE(111,*) '<DELTAE>=',DELTAE 
      ELSE
         WRITE(111,*) 'Default <DELTAE>=',DELTAE   
      ENDIF



! End accuracy parameters
!=============================================================
!==========================================================
! Begin old-type of ATOMINFO
      LATOMINFO = .FALSE.
      ! Initialize all clusters to 1
      CLS(1:NAEZD+NEMBD) = 1
      WRITE(1337,*) 'ATOMINFOC or ATOMINFO:'
      DO I=1,NATYP
         CALL IoInput('ATOMINFOC       ',UIO,I+1,7,IER)
         IA = 1
         IF ( IER.EQ.0 ) THEN
            LATOMINFO = .TRUE.
                           READ (UNIT=UIO,FMT=*)    ZAT(I),&
     &                        LMXC(I),&
     &                        (KFG(J,I),J=1,4),&
     &                        J,&
     &                        IER,&
     &                        NTCELL(I),&
     &                        MTFAC(I),&
     &                        IRNS(I),&
     &                        RMTREF(IER),IQAT(I),CONC(I)
            IQ = IQAT(I)
            REFPOT(IQ) = IER
            RMTREFAT(I) = RMTREF(IER)
            CLS(IQ) = J
            NOQ(IQ) = NOQ(IQ) + 1
            IF ( NOQ(IQ) .GT. 1 ) THEN
                ICPA(IQ) = 1
                NCPA = 1
            END IF
            KAOEZ(NOQ(IQ),IQ) = I
         ELSE
            IER = 0
            CALL IoInput('ATOMINFO        ',UIO,I+1,7,IER)
            IF (IER.EQ.0) THEN
               LATOMINFO = .TRUE.
                           READ (UNIT=UIO,FMT=*)    ZAT(I),&
     &                        LMXC(I),&
     &                        (KFG(J,I),J=1,4),&
     &                        J,&
     &                        REFPOT(I),&
     &                        NTCELL(I),&
     &                        MTFAC(I),&
     &                        IRNS(I),&
     &                        RMTREF(REFPOT(I))
                IQAT(I) = I
                RMTREFAT(I) = RMTREF(REFPOT(I))
                CLS(I) = J
                CONC(I) = 1.D0
                NOQ(I) = 1
                KAOEZ(1,I) = I
            ENDIF
         END IF
      END DO

      ! If old-style ATOMINFO is present, and if a 2-dim calculation is performed,
      ! and also if the RMTREF of the "outside region" is not read in explicitly
      ! (LNEW is false) then assign the RMTREF of the outside region according to
      ! the already-read-in REFPOT under LEFTBASIS  and RIGHBASIS.
      IF (LATOMINFO.AND.LINTERFACE.AND..NOT.LNEW) THEN
         DO I = NAEZ + 1,NAEZ + NEMB
            RMTREFAT(I) = RMTREF(REFPOT(I))
         ENDDO
      ENDIF

      NCLS = 0
      NREF = 0

      ! Determine total number of clusters
      DO I=1,NATYP 
        NCLS = MAX(NCLS,CLS(IQAT(I))) 
      ENDDO

      ! Determine total number of different reference potentials
      DO I=1,NAEZ + NEMB
        NREF = MAX(NREF,REFPOT(I)) 
      ENDDO

      !in line 1792  this is done: NINEQ = NAEZ, so here NINEQ is still undefinded
      !so we move this writeout back
      ! 
      !WRITE(6,2016) NCLS,NREF,NINEQ
      !WRITE(6,2110)
      !WRITE(6,2103)

      DO IQ=1,NAEZ
         SUM = 0D0
         IF (NOQ(IQ).LT.1) THEN
            WRITE(6,*) 'RINPUT13: CPA: SITE',IQ,'HAS NO ASSIGNED ATOM'
            STOP 'RINPUT13: CPA'
         ENDIF
         DO IO=1,NOQ(IQ)
            SUM = SUM + CONC(KAOEZ(IO,IQ))
         END DO
         IF ( ABS(SUM-1.D0).GT.1D-6) THEN
            WRITE(6,*) ' SITE ', IQ, ' CONCENTRATION <> 1.0 !'
            WRITE(6,*) ' CHECK YOUR <ATOMINFO-CPA> INPUT '
            STOP       ' IN <RINPUT99>'
         END IF
      END DO

! End old-type of ATOMINFO
!==========================================================



      ! Write out atominfo
      WRITE(1337,2028) NATYP
      WRITE(1337,2104)
      WRITE(1337,1029) ( &
     &     ZAT(I),&
     &     LMXC(I),&
     &     (KFG(J,I),J=1,4),&
     &     CLS(IQAT(I)),&
     &     REFPOT(IQAT(I)),&
     &     NTCELL(I),&
     &     MTFAC(I),&
     &     IRNS(I),&
     &     IQAT(I),CONC(I),I=1,NATYP)
      WRITE(1337,2108)
      WRITE(1337,2104)


! =============================================================
! Begin SCF convergence control

      NSTEPS = 1
      CALL IoInput('NSTEPS          ',UIO,1,7,IER)
      IF (IER.NE.0) THEN 
         WRITE(111,*) 'Default NSTEPS=',NSTEPS
      ELSE
         READ (UNIT=UIO,FMT=*) NSTEPS
      ENDIF
      IF (NPOL.EQ.0) THEN
         NSTEPS = 1
         WRITE(1337,*) 'NPOL=0, setting NSTEPS to 1'
      ENDIF
      IF (IGF.NE.0) THEN
         NSTEPS = 1
         WRITE(1337,*) 'IGF.NE.0, setting NSTEPS to 1'
      ENDIF
      IF (ICC.NE.0) THEN
         NSTEPS = 1
         WRITE(1337,*) 'ICC.NE.0, setting NSTEPS to 1'
      ENDIF
      IF (OPT('XCPL    ')) THEN
         NSTEPS = 1
         WRITE(1337,*) 'RUNOPT XCPL used, setting NSTEPS to 1'
      ENDIF
      IF (OPT('KKRFLEX ')) THEN
         NSTEPS = 1
         WRITE(1337,*) 'RUNOPT KKRFLEX used, setting NSTEPS to 1'
      ENDIF



      WRITE(1337,2011) NSTEPS
      WRITE(1337,2104)

      IMIX = 0
      CALL IoInput('IMIX            ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) imix
         WRITE(111,*) 'IMIX= ',IMIX
      ELSE
         WRITE(111,*) 'Default IMIX= ',IMIX
      ENDIF
      IF (NPOL.EQ.0) THEN
         WRITE(1337,*) 'NPOL=0, setting IMIX= 0'
         IMIX = 0
      ENDIF

      STRMIX = 0.01D0
      CALL IoInput('STRMIX          ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) STRMIX
         WRITE(111,*) 'STRMIX= ',STRMIX
      ELSE
         WRITE(111,*) 'Default STRMIX= ',STRMIX
      ENDIF
      IF (NPOL.EQ.0) THEN
         WRITE(1337,*) 'NPOL=0, setting STRMIX= 0.'
         STRMIX = 0
      ENDIF


      ITDBRY = 40
      CALL IoInput('ITDBRY          ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) itdbry
         WRITE(111,*) 'ITDBRY= ',ITDBRY
      ELSE
         WRITE(111,*) 'Default ITDBRY= ',ITDBRY
      ENDIF

      FCM = 20.D0
      CALL IoInput('FCM             ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) fcm
         WRITE(111,*) 'FCM= ',FCM
      ELSE
         WRITE(111,*) 'Default FCM= ',FCM
      ENDIF

      QBOUND = 1.D-7
      CALL IoInput('QBOUND          ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) qbound
         WRITE(111,*) 'QBOUND= ',QBOUND
      ELSE
         WRITE(111,*) 'Default QBOUND= ',QBOUND
      ENDIF


      BRYMIX = 0.01D0
      CALL IoInput('BRYMIX          ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) brymix
         WRITE(111,*) 'BRYMIX= ',BRYMIX
      ELSE
         WRITE(111,*) 'Default BRYMIX= ',BRYMIX
      ENDIF


      CALL IOINPUT('RMAX            ',UIO,1,7,IER)
      IF ( IER.NE.0 ) STOP 'rinput13: RMAX not in the inputcard'
      READ (UNIT=UIO,FMT=*) RMAX
      WRITE(111,*) 'RMAX= ',RMAX

      CALL IOINPUT('GMAX            ',UIO,1,7,IER)
      IF ( IER.NE.0 ) STOP 'rinput13: GMAX not in the inputcard'
      READ (UNIT=UIO,FMT=*) GMAX
      WRITE(111,*) 'GMAX= ',GMAX



! End SCF convergence control
! =============================================================
! ======================================================================
! Begin file name definitions

      
      IL=1
      CALL IoInput('FILES           ',UIO,IL,7,IER)
      IF (IER.EQ.0) THEN
         CALL IoInput('FILES           ',UIO,IL,7,IER)
         READ (UNIT=UIO,FMT='(A40)')  I12
         CALL IoInput('FILES           ',UIO,IL+1,7,IER)
         READ (UNIT=UIO,FMT='(A40)')  I13
         CALL IoInput('FILES           ',UIO,IL+2,7,IER)
         READ (UNIT=UIO,FMT='(A40)')  I40
         CALL IoInput('FILES           ',UIO,IL+3,7,IER)
         READ (UNIT=UIO,FMT='(A40)')  I19
         CALL IoInput('FILES           ',UIO,IL+4,7,IER)
         READ (UNIT=UIO,FMT='(A40)')  I25
      ELSE
         I13 = 'potential                               ' ! 40 chars
         I19 = 'shapefun                                ' ! 40 chars
         I25 = 'scoef                                   ' ! 40 chars
         I12 = '                                        ' ! 40 chars (not used)
         I40 = '                                        ' ! 40 chars (not used)
      ENDIF



      WRITE(1337,*) 'I12="',I12,'"'
      WRITE(1337,*) 'I13="',I13,'"'
      WRITE(1337,*) 'I40="',I40,'"'
      WRITE(1337,*) 'I19="',I19,'"'
      WRITE(1337,*) 'I25="',I25,'"'

! End file name definitions
! ======================================================================


      IFILE = 13
      CALL IoInput('<IFILE>         ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) ifile
         WRITE(111,*) '<IFILE>= ',IFILE
      ELSE
         WRITE(111,*) 'Default <IFILE>= ',IFILE
      ENDIF


      IPE = 1    ! Used to print out in calrmt      


      ISHIFT = 0
      CALL IoInput('ISHIFT          ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) ishift
         WRITE(111,*) 'ISHIFT= ',ISHIFT
      ELSE
         WRITE(111,*) 'Default ISHIFT= ',ISHIFT
      ENDIF
      IF (  OPT('rigid-ef').OR. OPT('DECIMATE') ) THEN
        ISHIFT = 2
        WRITE(1337,*) ' Rigid Fermi Energy, ISHIFT is set to ',ISHIFT
        WRITE(111,*) ' Rigid Fermi Energy, ishift is set to ',ISHIFT
      END IF
      IF ( TEST('no-neutr').OR.OPT('no-neutr') )   THEN
         ISHIFT = 1
         WRITE(1337,*) 'No charge neutrality required, ISHIFT is set to',ISHIFT
         WRITE(111,*) 'No charge neutrality required, ISHIFT is set to',ISHIFT
      ENDIF


      ESHIFT = 0.D0
      IRM = IRMD                ! never used
      INSREF = 0
      KWS = 2
      KHYP = 0
!      CALL IoInput('KHYPERF   ',UIO,1,7,IER)
!                      READ (UNIT=UIO,FMT=*) khyp  



      TOLRDIF = 0.5D0 ! Set free GF to zero for r<tolrdif (a.u.)(vir. atoms)
      CALL IoInput('<TOLRDIF>       ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) TOLRDIF
         WRITE(111,*) '<TOLRDIF>=',TOLRDIF
      ELSE
         WRITE(111,*) 'Default <TOLRDIF>=',TOLRDIF  
      ENDIF



! -------------------------------------------------
      KTE = 1
!     CALL IoInput('KTE       ',UIO,1,7,IER)
!                      READ (UNIT=UIO,FMT=*) kte

      KPRE = 1
!     CALL IoInput('KPRE      ',UIO,1,7,IER)
!                      READ (UNIT=UIO,FMT=*) kpre

      KEFG = 0
!     CALL IoInput('KEFG      ',UIO,1,7,IER)
!                      READ (UNIT=UIO,FMT=*) kefg

      KVMAD = 0
      CALL IoInput('KVMAD           ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) kvmad
         WRITE(111,*) 'KVMAD= ',KVMAD
      ELSE
         WRITE(111,*) 'Default KVMAD= ',KVMAD
      ENDIF

!



!----------------------------------------------------------------------
!
! --> determination of properties at Fermi level
!
      IF ( OPT('GF-EF   ') ) THEN
         IGF = 1
         IF (NPOL.GT.0) NPOL = 0
         IF (NPOL.LT.0) THEN 
            NPNT1 = 0
            NPNT3 = 0
         END IF
         NPNT2 = 1
      END IF

      IF (OPT('DOS-EF  ')) THEN
        NPOL = 0
        NPNT2 = 1
      END IF
! ----------------------------------------------------------------------
! ---------------------------------------------------------------------
!     

      KFORCE = 0
      IF (INS.GT.0) THEN      
         CALL IoInput('KFORCE          ',UIO,1,7,IER)
         IF (IER.EQ.0) THEN
            READ (UNIT=UIO,FMT=*) KFORCE
            WRITE(111,*) 'KFORCE= ',KFORCE 
         ELSE
            WRITE(111,*) 'Default KFORCE= ',KFORCE 
         ENDIF            
      END IF

      KFROZN = KCOR
      IF (KCOR.EQ.0) KCOR = 2




! ------------------------------------------------------------------------
      WRITE (1337,9210) LMAX
      WRITE (1337,9301)
      WRITE (1337,9220) EMIN,EMAX,TK
      WRITE (1337,9302)
      WRITE (1337,9230) NPOL,NPNT1,NPNT2,NPNT3
      WRITE (1337,9304)
      WRITE (1337,9303)
      WRITE (1337,9250) IFILE,IPE,ISHIFT,ESHIFT
      WRITE (1337,9305)
      WRITE (1337,9260) KSHAPE,IRM,INS,ICST,INSREF
      WRITE (1337,9309)
      WRITE (1337,9270) KCOR,KVREL,KWS,KHYP,KHFIELD,KXC
      WRITE (1337,9306)
      WRITE (1337,9330) KTE,KPRE,KEFG,KVMAD
      WRITE (1337,9309)
      WRITE (1337,9290) IMIX,IGF,ICC
      WRITE (1337,9304)
      WRITE (1337,9300) ITDBRY
      WRITE (1337,9307)
      WRITE (1337,9310) STRMIX,FCM,QBOUND
      WRITE (1337,9302)
      WRITE (1337,9320) BRYMIX
      WRITE (1337,9308)
      WRITE (1337,9280) HFIELD,VCONST
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------


      IPF = 1337
      IPFE = IPF + 3

      IF (OPT('SEARCHEF')) THEN
         IMIX=0
         MIXING=0.0d0
         STRMIX=MIXING
         ITDBRY=1
         QBOUND=1.0d-10
         WRITE(1337,'(1X,A)') 'Option SEARCHEF used overriding INPUT for'
         WRITE(1337,'(1X,A)') 'IMIX,MIX,QBOUND,ITDBRY: 0, 0.0, 1E-10, 1'
         WRITE(1337,*)
      ENDIF

      IF (IMIX.GT.2) THEN
        FCM = 1.0D0
        MIXING = BRYMIX
      ELSE
        MIXING = STRMIX
      END IF

      IF (IMIX.GE.6) WRITE (1337,FMT=9110) (IMIX-5),ITDBRY - 1

      WRITE (1337,FMT=9090) MIXING,QBOUND
!--------------------------------------------------------
      WRITE (1337,FMT=9091) CPAFLAG(NCPA)
      IF (NCPA.NE.0) WRITE(1337,9092) ITCPAMAX,CPATOL
!--------------------------------------------------------

      LMMAX = (LMAX+1)**2
      LPOT  = MIN(2*LMAX,LPOTD)
      LMPOT = (LPOT+1)* (LPOT+1)

      WRITE (1337,FMT=9020) LMAX,LMAXD,NATYP,NATYPD,IRM,IRMD,NSPIN,NSPIND


      IF (INS.GT.0) THEN
        WRITE (1337,FMT=9130)
        WRITE (1337,FMT=9140)
        DO 20 I = 1,NATYP
          WRITE (1337,FMT=9150) I,IRNS(I),IRNSD

          IF (IRNS(I).GT.IRNSD) CALL RCSTOP('19      ')

   20   CONTINUE

        IF (LMAX.NE.LMAXD) THEN
          WRITE (1337,FMT=9120)

          CALL RCSTOP('20      ')

        END IF

      END IF


      WRITE (1337,FMT=9130)




      IF (KHFIELD.EQ.1) WRITE (1337,FMT=9030) HFIELD
      IF (KVREL.LE.1 ) THEN 
          WRITE (1337,FMT=9050) TSPIN(NSPIN)
      ELSE
          WRITE (1337,FMT=9050) TSPIN(NSPIN+1)
      END IF
      WRITE (1337,FMT=9170) TVREL(KVREL)
      WRITE (1337,FMT=9170) TKCOR(KFROZN)
      IF (KSHAPE.EQ.0) THEN
        WRITE (1337,FMT=9070) TKWS(KWS+1)

      ELSE
        WRITE (1337,FMT=9170) TSHAPE
      END IF

      WRITE (1337,FMT=9100) TXC(KXC+1)
      IF (INS.GT.0) WRITE (1337,FMT=9160) TINS(INS),ICST
      WRITE (1337,FMT=9080)


      VBC(1) = VCONST
      VBC(2) = VBC(1)

      LRHOSYM = .FALSE.
      CALL IoInput('LRHOSYM         ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) lrhosym
         WRITE(111,*) 'LRHOSYM= ',LRHOSYM
      ELSE
         WRITE(111,*) 'Default LRHOSYM= ',LRHOSYM
      ENDIF


      IF ( (NCPA.NE.0).AND.LRHOSYM ) THEN
         WRITE(1337,*) ' WARNING : CHARGE SYMMETRISATION NOT ALLOWED FOR CPA '
         WRITE(1337,*) '        YOUR SETTING IN INPUT FILE IS OVERRIDDEN'
         WRITE(111,*) ' WARNING : CHARGE SYMMETRISATION NOT ALLOWED FOR CPA '
         WRITE(111,*) '    YOUR SETTING IN INPUT FILE IS OVERRIDDEN'
         LRHOSYM = .FALSE.
      END IF
           

      IF (LRHOSYM) THEN

        CALL IoInput('IXIPOL          ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) (ixipol(I),I=1,natyp) 
        write (1337,2022) (ixipol(i),i=1,natyp) 
        write (1337,2103)
        DO I=1,NATYP
          IF ( IXIPOL(I).NE.0 .AND. ABS(IXIPOL(ABS(IXIPOL(I)))).NE.I) THEN
            write(6,*) 'Error in IXIPOL at atom ',I,'.'
            stop 'IXIPOL'
          END IF
        END DO
      ELSE
        DO I=1,NATYP
          IXIPOL(I) = 0
        END DO
        write (1337,2022) (ixipol(i),i=1,natyp) 
        write (1337,2103)
      END IF
      write(1337,2023) NAEZ,NEMB
      write(1337,2110)


      NINEQ = NAEZ
      WRITE(1337,2016) NCLS,NREF,NINEQ
      WRITE(1337,2110)
      WRITE(1337,2103)


!----------------------------------------------------------------------

      KMROT = 0

      DO I=1,NAEZ

         ! --->  atoms equivalent by inversional symmetry

         QMTET(I)=0D0
         QMPHI(I)=0D0
         IER = 0 
         CALL IoInput('RBASISANG       ',UIO,I,7,IER)

         IF( IER.EQ.0 ) THEN
            READ (UNIT=UIO,FMT=*) (RBASIS(J,I), J=1,3),QMTET(I),QMPHI(I)
            IF( ABS(QMTET(I)) .GT. 1D-6 ) KMROT = 1
            IF( ABS(QMPHI(I)) .GT. 1D-6 ) KMROT = 1
         ENDIF
      ENDDO                         ! I=1,NAEZ
      CALL IDREALS(RBASIS(1,1),3*NAEZ,IPRINT)




!-------------------------------------------------------------

      if (nemb.gt.0) write(1337,*) 
      write(1337,2031) ((rbasis(j,i),j=1,3),i,refpot(i),i=naez+1,naez+nemb)



! ------------------------------------------------------------------------
      IF ( .not. OPT('VIRATOMS') ) THEN
        DO I=1,NAEZ
          DO IO=1,NOQ(I)
            IF (KAOEZ(IO,I).LT.1) STOP 'Error in KAOEZ'
          END DO
        ENDDO
      END IF


! ------------------------------------------------------------------------
      WRITE(1337,2111)

!Check for DECIMATE consistency

      IF (OPT('DECIMATE')) THEN
         IF ( MOD(NPRINCD,NLBASIS).NE.0 ) THEN
            WRITE(6,*) ' Decimation cannot continue '
            WRITE(6,*) 'NPRINCD=',NPRINCD,' NLBASIS=',NLBASIS
            STOP
         END IF
         IF ( MOD(NPRINCD,NRBASIS).NE.0 )  THEN
            WRITE(6,*) ' Decimation cannot continue '
            WRITE(6,*) 'NPRINCD=',NPRINCD,' NRBASIS=',NRBASIS
            STOP
         END IF
      END IF

!Check for ITERMDIR consistency -- if KMROT=0 suppress it

      IF ( (OPT('ITERMDIR')).AND.(KMROT.EQ.0) ) THEN
         WRITE (1337,*)
         WRITE (1337,*)&
     &        ' WARNING: ITERMDIR running option used with collinear/',&
     &        'parallel Oz starting'
         WRITE (1337,*)&
     &        '          system (KMROT = 0 ). Please check token',&
     &        ' RBASISANG in your input'
         WRITE (1337,*) ' Running option ITERMDIR will be ignored'
         WRITE (1337,*)
         DO I=1,32
            IF (OPTC(I)(1:8).EQ.'ITERMDIR') OPTC(I)='        '
         END DO
      END IF

!Check for XCPL consistency 

      MANCTL = ( KMROT.EQ.0 ).AND.( KREL.EQ.0 ).AND.( NPOL.NE.0 ).AND.( NSPIN.GT.1 )
      IF ( (OPT('XCPL    ') ).AND.( .NOT.MANCTL ) ) THEN
         WRITE (1337,*)
         WRITE (1337,*)&
     &        ' WARNING: XCPL running option requires collinear ',&
     &        'magnetic systems, complex'
         WRITE (1337,*)&
     &        '          energy contour (NPOL<>0) in a NON/SCALAR',&
     &        ' relativistic mode (KREL=0)'
         WRITE (1337,*) ' Running option XCPL will be ignored'
         WRITE (1337,*)
         DO I=1,32
            IF (OPTC(I)(1:8).EQ.'XCPL    ') OPTC(I)='        '
         END DO
      END IF


      WRITE(1337,62) (OPTC(I),I=1,8)                                              
 62   FORMAT(79('-')/' EXECUTION OPTIONS:'/1X,A8,7('//',A8)/79('-'))
      WRITE(1337,52) (TESTC(I),I=1,16)                                    
 52   FORMAT(79('-')/' TEST OPTIONS:'/2(1X,A8,7('//',A8)/)/79('-'))
 980  FORMAT(8A8)

! ---------------------------------------------------------------------
! Initialise SOLVER, SOC and CTL parameters in REL case


      IF (KREL.EQ.1) THEN
         SOLVER='BS        '

         CALL IoInput('SOLVER          ',UIO,0,7,IER)
         IF (IER.EQ.0) THEN
              READ (UNIT=UIO,FMT=*) SOLVER
              IF ( SOLVER(1:2) .EQ. 'BS' ) THEN
                 SOLVER = 'BS        '
              ELSE
                 IF( SOLVER .NE. 'ABM-OP    ' ) SOLVER='ABM-OP    '
              END IF
         END IF

         SOCSCL(1:LMAXD+1,1:NATYP) = 1.D0
         CSCL(1:LMAXD+1,1:NATYP) = CVLIGHT
         MANSOC=.FALSE.
         MANCTL=.FALSE.

! ============================================================= SOC-MAN
! For Dirac-ASA

         IF (OPT('SOC     ')) THEN
            CALL IOInput('SOSCALE         ',UIO,0,7,IER)
            IF (IER.EQ.0) THEN
               READ (UNIT=UIO,FMT=*) SOSCALE
               IF (SOSCALE.GT.-2.5D0) THEN 
                  IF (SOSCALE.GE.0.0D0) THEN           ! SOC-I
                     SOLVER='ABM-SOC   '
                     MANSOC=.TRUE.
                  ELSE                                  ! SOC-II
                     SOLVER       = 'ABM-SOC-II'
                     MANSOC=.TRUE.
                     DO I=1,NATYP
                        SOCSCL(1:LMAXD+1,I) = SOSCALE
                     END DO
                     WRITE(1337,99010) SOCII(NINT(SOSCALE))
                  END IF
               ELSE
                  WRITE(1337,99001) '< SOC >'
                  WRITE(1337,99003)
               END IF
            ELSE
               WRITE(1337,99002) '< SOC >'
               WRITE(1337,99003)
            END IF

            IF ( MANSOC .AND. (SOSCALE.GE.0D0) ) THEN
               IMANSOC(1:NATYP) = 1

! ---> now look for a possible include/exclude list (SOCLIST= +/- NASOC)
! ---> if SOCLIST is not found, ALL the atoms will have SOC modified with
! ---> SOSCALE (+NASOC=only NASOC atoms, -NASOC=all but these NASOC atoms)
!      Note that this is allowed only for SOC-I manipulation
!              
               CALL IOInput('SOCLIST         ',UIO,0,7,IER)
               IF (IER.EQ.0) THEN
                  READ(UNIT=UIO,FMT=*) NASOC,(ISP(I),I=1,ABS(NASOC))
                  
                  IF (NASOC.NE.0) THEN
                     IF (NASOC.LT.0) THEN ! exclude this atoms
                        DO I=1,-NASOC
                           IMANSOC(ISP(I)) = 0
                        END DO
                     ELSE
                        IMANSOC(1:NATYP) = 0
                        DO I=1,NASOC
                           IMANSOC(ISP(I)) = 1
                        END DO
                     END IF
                  END IF
               END IF
     
               WRITE(1337,2100)
               DO I=1,NATYP
                  IF (IMANSOC(I).EQ.1) THEN
                     SOCSCL(1:LMAXD+1,I)=SOSCALE
                  END IF
               END DO
               WRITE(1337,99004)
               IF (NASOC.EQ.0) WRITE(1337,99005)
               IF (NASOC.GT.0) THEN
                  WRITE(1337,99006)
                  WRITE(1337,99008) (ISP(I),I=1,NASOC)
               END IF
               IF (NASOC.LT.0) THEN
                  WRITE(1337,99007)
                  WRITE(1337,99008) (ISP(I),I=1,ABS(NASOC))
               END IF
               WRITE(1337,99009) SOSCALE
               WRITE(1337,2100)
            END IF
         END IF

! ============================================================= SOC-MAN

         WRITE(1337,'('' SOLVER used for the DIRAC equation : '',2X,A)') SOLVER
         WRITE(1337,2100)

! ============================================================= CTL-MAN

         IF (OPT('CSCALE  ')) THEN
            CALL IOInput('CTLSCALE        ',UIO,0,7,IER)
            IF (IER.EQ.0) THEN
               READ (UNIT=UIO,FMT=*) CTLSCALE
               IF (CTLSCALE.GE.1D-12) THEN 
                  MANCTL=.TRUE.
               ELSE
                  WRITE(1337,99001) '< CSCALE >'
                  WRITE(1337,99011)
               END IF
            ELSE
               WRITE(1337,99002) '< CSCALE >'
               WRITE(1337,99011)
            END IF

            IF (MANCTL) THEN
               CSCL(1:LMAXD+1,1:NATYP) = CSCL(1:LMAXD+1,1:NATYP)/DSQRT(CTLSCALE)
               WRITE(1337,99012)
               WRITE(1337,99005)
               WRITE(1337,99009) 1.D0/DSQRT(CTLSCALE)
            END IF
            WRITE(1337,2100)
         END IF

! ============================================================= CTL-MAN
      END IF
! ================================================================ LDA+U



      WRITE(1337,2100) 
      WRITE(1337,2040) KMROT
      WRITE(1337,2110)
      WRITE(1337,*) ' >>>>>>>>> RINPUT13 EXITS NOW <<<<<<<<<< '


      CLOSE(111) ! Close file inputcard_generated.txt

      RETURN
! *********************************************Input-End ********
 1029 FORMAT((F4.0,I4,4x,4I1,3I4,F8.4,I4,I5,1x,f8.5))
! ------------------------------------------------------------------------
 2004 FORMAT( /79(1H*)/&
     &     '*',77X,'*'/&
     &     '*',10X,'Screened Korringa-Kohn-Rostoker ',&
     &             'Electronic Structure Code',10X,'*'/&
     &     '*',27X,'for Bulk and Interfaces',27X,'*'/&
     &     '*',77X,'*'/&
     &     '*',4X,'Juelich-Munich 2001 - 2015',25X,&
     &     'Version : ',A8,4X,'*'/&
     &     '*',77X,'*'/79(1H*))
 2010 FORMAT(' NSPIN '/I4)
 2011 FORMAT(' NSTEPS'/I4)
 2013 FORMAT('      M2    MMIN    MMAX    SINN',&
     &       '    SOUT     RIN    ROUT'/7I8)
 2014 FORMAT('          ALAT = ',F15.8)
 2015 FORMAT('   INTERVX   INTERVY   INTERVZ'/3I10)
 2016 FORMAT('    NCLS    NREF   NINEQ'/,3I8)
 2018 FORMAT(' RBASIS'/,&
     &     'SITE                BASIS VECTORS                 ',&
     &     'THETA   PHI CPA OCC KAOEZ')
 2019 FORMAT('         ABASIS         BBASIS         CBASIS'/3F15.8)
 2021 FORMAT(' INIPOL'/,(10I4))
 2022 FORMAT(' IXIPOL'/,(10I4))
 2023 FORMAT('    NAEZ    NEMB  '/,2I8)
 2025 FORMAT((i4,3F15.8,2F6.1,2(1x,I3),4I3))
 2028 FORMAT(' NATYP '/,I4/,&
     &     '   Z lmx     KFG cls pot ntc  MTFAC irns SITE  CONC')
 2029 FORMAT(' KMT   '/,I4)
 2031 FORMAT((3F15.8,2I6))
 2032 FORMAT(' NTCELLR'/,(10I4))
 2040 FORMAT(' KMROT'/,4I8)
! ------------------------------------------------------------------------
 2100 FORMAT(79(1H-))
 2101 format(   3(1H-),1H+  , 3(14(1H-),1H+),  30(1H-))
 2102 format( 3(9(1H-),1H+) ,49(1H-))
 2103 FORMAT(10(3(1H-),1H+) ,39(1H-))
 2104 format(   3(1H-),1H+  ,75(1H-))
 2107 format( 3(14(1H-),1H+),34(1H-))
 2108 format( 2(3(1H-),1H+),  7(1H-),1H+,      3(3(1H-),1H+),&
     &          7(1H-),1H+,   3(1H-),1H+,      39(1H-))
 2110 format( 3(7(1H-),1H+) ,55(1H-))
 2111 format( 7(7(1H-),1H+) ,23(1H-))
 9020 FORMAT (/,33x,'check of dimension-data consistency',/,33x,&
     &       35 ('-'),/,40x,'lmax   : (',i6,',',i6,')',/,40x,&
     &       'natyp  : (',i6,',',i6,')',/,40x,'irm    : (',i6,',',i6,&
     &       ')',/,40x,'nspin  : (',i6,',',i6,')',/)
 9030 FORMAT (1x,10 ('*'),' external magnetic field applied hfield=',&
     &       f8.5)
 9050 FORMAT (20x,a4,'spin polarized calculation')
 9070 FORMAT (1x,20x,' calculation with',a8,'-potential')
 9080 FORMAT (1x,79 ('*'))
 9090 FORMAT (' mixing factor used           :',f15.6,/,&
     &        ' convergence quality required :',1p,d15.2)
 9091 FORMAT (' make use of CPA algorithm    :',1x,a14)
 9092 FORMAT ('         max. iterations      :',i15,/,&
     &        '         req. CPA convergency :',1p,d15.2)
 9100 FORMAT (1x,20x,a24,'exchange-correlation potential')
 9110 FORMAT (/,20x,'broyden"s method # :',i3,&
     &       ' is used up to iteration-      ',/,20x,'depth :',i3,&
     &       '  then jacobian is fixed and potential      ',/,20x,&
     &       'is updated using that jacobian')
 9120 FORMAT (13x,' in case of calculating non - spherical wavefcts ',&
     &       'the parameter lmaxd has to be set equal lmax ')
 9130 FORMAT (/)
 9140 FORMAT (20x,'full potential calculation ',&
     &       '- cut off of non spherical potential',/,' >',/)
 9150 FORMAT (31x,'representive atom no.',i3,' irns :',i5,' irnsd :',i5)
 9160 FORMAT (21x,a43,/,21x,' using',i3,'-th. born approximation ')
 9170 FORMAT (21x,a43)
 9210 FORMAT (' lmax'/,i4)
 9220 FORMAT ('          EMIN        EMAX        TK'/,3f12.6)
 9230 FORMAT ('   NPOL  NPNT1  NPNT2  NPNT3'/,4i7)
 9240 FORMAT (' IRNUMX ITCCOR IPRCOR'/,3i7)
 9250 FORMAT ('  IFILE    IPE ISHIFT ESHIFT'/,3i7,f12.6)
 9260 FORMAT (' KSHAPE    IRM    INS   ICST INSREF'/,5i7)
 9270 FORMAT ('   KCOR  KVREL    KWS   KHYP KHFIELD   KXC'/,6i7)
 9280 FORMAT (' external magnetic hfield     :',f15.4/,&
     &        ' VCONST                       :',f15.6)
 9290 FORMAT ('   IMIX    IGF    ICC'/,3i7)
 9300 FORMAT (' ITDBRY'/,i7)
 9310 FORMAT ('      STRMIX        FCM       QBOUND'/,3f12.6)
 9320 FORMAT ('      BRYMIX'/,f12.6)
 9330 FORMAT ('    KTE   KPRE   KEFG  KVMAD '/,5i7)
 9301 format(   3(1H-),1H+  ,75(1H-))
 9302 format( 3(11(1H-),1H+),43(1H-))
 9303 format(3(6(1H-),1H+) ,58(1H-))
 9304 format(4(6(1H-),1H+) ,51(1H-))
 9305 format(3(6(1H-),1H+),11(1H-),1H+ ,46(1H-))
 9306 format(6(6(1H-),1H+) ,37(1H-))
 9307 format(6(1H-),1H+,72(1H-))
 9308 format(11(1H-),1H+,67(1H-))
 9309 format(5(6(1H-),1H+) ,44(1H-))
 9410 format('*** SLAB - INTERFACE CALCULATION ***'/)
 9420 format(I5,3F14.8,I5)
 9430 format('Number of LEFT  Host Layers : ',I5,' with ',I5,' basis')
 9440 format('Number of RIGHT Host Layers : ',I5,' with ',I5,' basis')
 9450 format('Left  side periodicity : ',3F10.5)
 9460 format('Right side periodicity : ',3F10.5)
 9465 format('    Geommetry used : '/,&
     &       ' ATOM       TX          TY          TZ ')
 9470 format('--------------- Left  Host -------------- ')
 9475 format('---------------   S L A B  -------------- ')
 9480 format('--------------- Right Host -------------- ')
99001 FORMAT(/,1X,&
     &     "WARNING: Option ",A," used with an INVALID ",&
     &     "scaling parameter.")
99002 FORMAT(/,1X,&
     &     "WARNING: Option ",A," found but NO value given for the",&
     &     " scaling parameter.")
99003 FORMAT(15X,'++++++++++   SOC option will be IGNORED   ++++++++++',&
     &     /,1X,'Please use SOCSCALE= XXX (real>-2.5) in the inputcard',&
     &     ' to make your option valid ',/)
99004 FORMAT(1X,'The SOC will be SCALED',$)
99005 FORMAT(' for ALL the atoms in the unit cell.')
99006 FORMAT(' for the FOLLOWING atoms in the unit cell :')
99007 FORMAT(' for all the atoms in the unit cell EXCLUDING :')
99008 FORMAT(1X,6(2X,I3))
99009 FORMAT(1X,'Scaling factor = ',1P,D9.2)
99010 FORMAT(1X,'The SOC is manipulated',' -- part of the SOC kept: ',A)
99011 FORMAT(15X,'+++++++++  CSCALE option will be IGNORED  ++++++++++',&
     &     /,1X,'Please use CTLSCALE= X (real>=1D-12) in the inputcard',&
     &     ' to make your option valid ',/)
99012 FORMAT(1X,'The CLIGHT will be SCALED',$)

      END SUBROUTINE RINPUT13
!---------------------------------------------------------------------
!---------------------------------------------------------------------

SUBROUTINE ADDOPT(STRING)
IMPLICIT NONE
INTEGER NOPTD
PARAMETER (NOPTD=32)
CHARACTER*8 STRING
CHARACTER*8 OPTC(NOPTD)
INTEGER II
LOGICAL OPT,LADDED
EXTERNAL OPT
COMMON /OPTC/OPTC

IF (.NOT.OPT('        ')) THEN
   WRITE(*,*) 'Error in ADDOPT for ',STRING,' : No free slots in array OPTC.'
   STOP 'Error in ADDOPT: No free slots in array OPTC.'
ENDIF

IF (.NOT.OPT(STRING)) THEN
   II = 1
   DO WHILE (II.LE.NOPTD)
      IF (OPTC(II).EQ.'        ') THEN
         OPTC(II) = STRING
         II = NOPTD + 1
      ENDIF
      II = II + 1
   ENDDO
ENDIF

END SUBROUTINE ADDOPT
