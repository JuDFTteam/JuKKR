      SUBROUTINE READINPUT(BRAVAIS,LCARTESIAN,RBASIS,ABASIS,BBASIS,CBASIS, &
     &           DX,DY,DZ, &
     &           ALATC,BLATC,CLATC, &
     &           IRNS,NAEZ,NEMB,KAOEZ,IRM,ZAT,SITEAT, &
     &           INS,KSHAPE, &
     &           LMAX,LMMAX,LPOT,  &
     &           NATYP,NSPIN, &
     &           NMIN,NRAD,NSMALL, &
     &           TOLHS,TOLVDIST,TOLAREA, &
     &           KXC,TXC, &
     &           I13, &
     &           NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,   &  
     &           TLEFT,TRIGHT,LINTERFACE,RCUTZ,RCUTXY,RMTCORE, &
     &           LMTREF,RMTREF,SIZEFAC,NFACELIM, EFSET)   
      use mod_version_info, only: serialnr
!#@# KKRtags: VORONOI input-output
      implicit none
      include 'inc.geometry'
!     .. Local Arrays ..
      CHARACTER*4 TSPIN(2)
      CHARACTER*8 TKWS(3)
      CHARACTER*43 TINS(0:3),TKCOR(0:3),TVREL(0:2)
!     ..
!     .. External Subroutines ..
      EXTERNAL RCSTOP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MIN
!     ..
!     .. Array Arguments ..
      INTEGER IRNS(*),KFG(4,NATYPD),LMXC(NATYPD)
      INTEGER NTCELL(NATYPD),CLS(NATYPD)
      INTEGER REFPOT(NATYPD+NEMBD)
      INTEGER KAOEZ(*)
      INTEGER NBSHIFT ! How many basis atoms have been shifted from polyhedron center
      INTEGER SITEAT(NATYPD),NOQ(NAEZD)
      REAL*8  ZAT(*),MTFAC(NATYPD),BRAVAIS(3,3),RBASIS(3,*),RMTCORE(*),RMTREF(*)
      REAL*8  DX(*),DY(*),DZ(*),RBUNSHIFT(3,NATYPD),RAUX(3)  ! Unshifted positions
      REAL*8  TRIGHT(3,*),TLEFT(3,*),ZPERLEFT(3),ZPERIGHT(3) 
      REAL*8  MTWGHT(0:NTOTD)
      REAL*8  SIZEFAC(-NLEMBD*NEMBD:NTOTD)
      CHARACTER*124 TXC(3)
      CHARACTER*256 UIO
      CHARACTER*40 I12,I13,I19,I25,I40
      character*80 text
!     ..
!     .. Scalar Arguments ..
      REAL*8        ALAT,E1,E2,ESHIFT,FCM,HFIELD,MIXING,QBOUND,TK, &
     &       VCONST,ABASIS,BBASIS,CBASIS,RCUTZ,RCUTXY,RMTREFDEF
      REAL*8 TOLHS,TOLVDIST,TOLAREA
      REAL*8 EFSET ! set Fermi level to this value
      INTEGER ICC,ICST,IFILE,IGF,IMIX,INS, &
     &        IPE,IPF,IPFE,IPOTOU,IPRCOR,ISITE, &
     &        IRM,IRNUMX,ISHIFT, &
     &        KPRE,KSCOEF,KSHAPE, &
     &        KVREL,KWS,KXC,LMAX,LMMAX,LMPOT,LPOT,MD, &
     &        NATYP,NPNT1,NPNT2,NPNT3,NPOL,NSPIN,INDX,IAT
      INTEGER NMIN,NSMALL,NRAD,NFACELIM,NBR
      INTEGER NSTEPS,KMT,NAEZ,NEMB
      INTEGER NINEQ,NEMBZ,NZ,CENTEROFINV(3)
      REAL*8        ALATC,BLATC,CLATC
      INTEGER MMIN,MMAX,SINN,SOUT,RIN,ROUT
      INTEGER INTERVX,INTERVY,INTERVZ,NREF,NCLS
      INTEGER NLBASIS,NRBASIS,NLEFT,NRIGHT       
      LOGICAL COMPLX,LINTERFACE,LMTREF,LNEW,LCARTESIAN 
!     ..
!     .. Local Scalars ..
      REAL*8        BRYMIX,STRMIX,E3,TX,TY,TZ
      INTEGER I,IL,IP,J,IER,I1,IC,II,M2,IQ
      CHARACTER*43 TSHAPE
!
      CHARACTER*8 TESTC(16),OPTC(8)
      COMMON /TESTC/TESTC
      COMMON /OPTC/OPTC
!     ..
!     .. Data statements ..
      DATA RMTREFDEF/2.0/ ! default rmt-ref
      DATA TSPIN/'non-','    '/
      DATA TSHAPE/' exact cell treatment (shape correction)  '/
      DATA TINS/' spherical averaged input potential        ', &
     &     ' non spherical input potential for cluster ', &
     &     ' non spherical input potential for cluster ', &
     &     ' non spherical input potential             '/
      DATA TKWS/' full mt','   ws   ',' full ws'/ 
!     ..
!
!------------ array set up and definition of input parameter -----------
!
      TXC(1) = ' Morruzi,Janak,Williams  #serial: ' // serialnr 
      TXC(2) = ' von Barth,Hedin         #serial: ' // serialnr
      TXC(3) = ' Vosko,Wilk,Nusair       #serial: ' // serialnr

      OPEN(111,FILE='inputcard_generated.txt') ! Write out found or assumed values



      IL=1
      I=1

      KXC = 1

! =============================================================================
! Begin Structure
      WRITE(*,*) 'Begin Structure'

! set interface/bulk (2D/3D) mode
      LINTERFACE = .FALSE.
      CALL IoInput('INTERFACE       ',UIO,IL,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) LINTERFACE
         WRITE(111,*) 'INTERFACE= ',LINTERFACE
      ENDIF

! Read in the bravais vectors (normalized to alatc)
! Notation: BRAVAIS(J,I) J=x,y,z I=1,2,3
! If the third bravais vector is zero, then surface (2-dimentional) geometry
! is implied.
!

      
      CALL IoInput('BRAVAIS         ',UIO,I,7,IER)
      IF (IER.NE.0) THEN
         WRITE(*,*) 'readinput: error: BRAVAIS not found, stopping.'
         STOP 'readinput: error: BRAVAIS not found'
      ENDIF

      NBR = 3
      IF (LINTERFACE) NBR = 2

      BRAVAIS(:,:) = 0.D0
      DO I=1,NBR
         CALL IoInput('BRAVAIS         ',UIO,I,7,IER)
              READ (UNIT=UIO,FMT=*) (BRAVAIS(J,I), J=1,3)
      ENDDO

      IF (BRAVAIS(1,3).EQ.0.D0.AND.BRAVAIS(2,3).EQ.0.D0.AND.BRAVAIS(3,3).EQ.0.D0) THEN
         LINTERFACE=.TRUE.
         WRITE(*,*) 'Surface geometry'
      ENDIF

      WRITE(111,FMT='(A10)') 'BRAVAIS   '
      DO I=1,NBR
         WRITE(111,FMT='(3E16.8)') BRAVAIS(1:3,I)
      ENDDO


      LCARTESIAN= .FALSE.  ! defalt is false, then bravais lattice coordinates are used
      CALL IoInput('CARTESIAN       ',UIO,1,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) LCARTESIAN
      WRITE(*,*) 'readinput: CARTESIAN=',LCARTESIAN
      WRITE(111,*) 'CARTESIAN= ',LCARTESIAN


      CALL IoInput('NAEZ            ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) NAEZ
         WRITE(*,*) 'NAEZ=', NAEZ
      ELSE
         WRITE(*,*) 'readinput12: NAEZ not found, stopping.'
         STOP 'readinput12: NAEZ not found, stopping.'
      ENDIF

      IF(NAEZ.GT.NAEZD) THEN
        WRITE(6,*) 'Please, increase the parameter naezd (',naezd,') in inc.p to',naez
        STOP 'readinput12: ERROR in NAEZD.'
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

      NFACELIM = NFACED   ! Upper limit of allowed number of faces can be smaller than dimension for speedup
      CALL IoInput('<NFACELIM>      ',UIO,1,7,IER)
      IF (IER.EQ.0) THEN
         READ (UNIT=UIO,FMT=*) NFACELIM
         WRITE(111,*) '<NFACELIM>=',NFACELIM
      ENDIF

      IF (NFACELIM.GT.NFACED) THEN
         WRITE(*,*) 'readinput: NFACELIM.GT.NFACED',NFACELIM,NFACED
         STOP 'readinput: NFACELIM.GT.NFACED'
      ENDIF



      RBUNSHIFT(1:3,1:NAEZ) = RBASIS(1:3,1:NAEZ)
      DX(1:NAEZ) = 0.D0
      DY(1:NAEZ) = 0.D0
      DZ(1:NAEZ) = 0.D0
      CALL IoInput('NBSHIFT         ',UIO,0,7,IER)
      IF (IER.EQ.0)  THEN 
         WRITE(*,*) 'Reading unshifted positions...'
         READ (UNIT=UIO,FMT=*) NBSHIFT
         IF (NBSHIFT.GT.NAEZ) STOP 'NBSHIFT > NAEZ'
         WRITE(*,*) '... for',NBSHIFT,' basis atoms.'
         DO I = 1,NBSHIFT
            CALL IoInput('RBUNSHIFT       ',UIO,I,7,IER)
            READ (UNIT=UIO,FMT=*) INDX,(RBUNSHIFT(J,INDX), J=1,3)
         ENDDO
         !     Now put the unshifted positions in array rbasis,
         !     so that the voronoi cells are centered in these.
         !     The extra shift is saved in DX,DY,DZ.
         DO I = 1,NAEZ
            DX(I) = RBASIS(1,I) - RBUNSHIFT(1,I)
            DY(I) = RBASIS(2,I) - RBUNSHIFT(2,I)
            DZ(I) = RBASIS(3,I) - RBUNSHIFT(3,I)
            RBASIS(1:3,I) = RBUNSHIFT(1:3,I)
         ENDDO
      ENDIF

      DO I=1,NAEZ
         KAOEZ(I) = I
      END DO

      ABASIS = 1.D0
      BBASIS = 1.D0
      CBASIS = 1.D0
      CALL IoInput('BASISCALE       ',UIO,IL,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) ABASIS,BBASIS,CBASIS 

      CALL IoInput('ALATBASIS       ',UIO,IL,7,IER)
      IF (IER.EQ.0) THEN 
         READ (UNIT=UIO,FMT=*) ALATC
      ELSE
         WRITE(*,*) 'readinput12: ALATBASIS not found, stopping.'
         STOP 'readinput12: ALATBASIS not found, stopping.'
      ENDIF
      BLATC = 1.D0
      CLATC = 1.D0





      ! =============================================================================
      ! Start Interface mode
      !
      ! --- > Read Left and Right host and set up the embending positions
      !
      !

      NLEFT = 0 ; NRIGHT = 0 ; NEMB = 0 ! in 3d geometry

      IF (LINTERFACE) THEN

         WRITE(6,9410)


         IER = 0
         ! Check if the keywords exist for old/new treatment of left and right host
         CALL IoInput('LEFTBASIS       ',UIO,1,7,IER)
         IF (IER.EQ.0) THEN
            LNEW = .FALSE.
         ELSE
            LNEW = .TRUE. ! New type input, KAOEZ not needed, keywords <RBLEFT>,<RBRIGHT> used.
            IER = 0
            CALL IoInput('<RBLEFT>        ',UIO,1,7,IER)
         ENDIF
         IF (IER.NE.0) THEN
            WRITE(*,*) 'readinput12: LEFTBASIS or <RBLEFT> not found in inputcard'
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
            WRITE(*,*) 'readinput12: RIGHBASIS or <RBRIGHT> not found in inputcard'
            STOP 'readinput12: RIGHBASIS or <RBRIGHT> not found in inputcard'
         ENDIF



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

         IF (NLEFT+NRIGHT.GT.NLEMBD) THEN
            WRITE(*,*) 'readinput12: NLEFTHOS+NRIGHTHO=',NLEFT+NRIGHT,' exceeds dimension NLEMBD=',NLEMBD
            STOP 'readinput12: NLEFTHOS+NRIGHTHO exceeds dimension NLEMBD'
         ENDIF

         CALL IoInput('<NLBASIS>       ',UIO,IL,7,IER)
         IF (IER.EQ.0) THEN
            READ (UNIT=UIO,FMT=*) NLBASIS
         ELSE
            IER=0
            CALL IoInput('NLBASIS         ',UIO,IL,7,IER)
            IF (IER.EQ.0) THEN
               READ (UNIT=UIO,FMT=*) NLBASIS
            ELSE
               WRITE(*,*) 'readinput12:<NLBASIS> or NLBASIS not found in inputcard'
               STOP 'readinput12: NLBASIS not found in inputcard'
            ENDIF
         ENDIF

         CALL IoInput('<NRBASIS>       ',UIO,IL,7,IER)
         IF (IER.EQ.0) THEN
            READ (UNIT=UIO,FMT=*) NRBASIS
         ELSE
            IER=0
            CALL IoInput('NRBASIS         ',UIO,IL,7,IER)
            IF (IER.EQ.0) THEN
               READ (UNIT=UIO,FMT=*) NRBASIS
            ELSE
               WRITE(*,*) 'readinput12:<NRBASIS> or NRBASIS not found in inputcard'
               STOP 'readinput12: NLBASIS not found in inputcard'
            ENDIF
         ENDIF


         ! Information is enought to define NEMB
         NEMB = NLBASIS + NRBASIS

         IF(NEMB.GT.NEMBD) THEN
            write(6,*) 'Please, increase the parameter nembd (',nembd,') in inc.p to',nemb
            STOP 'ERROR in NEMBD.'
         ENDIF


         IF (LNEW) THEN

            DO I=1,NLBASIS
               CALL IoInput('<RBLEFT>        ',UIO,I,7,IER)
               READ (UNIT=UIO,FMT=*) (TLEFT(I1,I),I1=1,3)
               KAOEZ(NAEZ+I) = I            ! Default
               CALL IoInput('<KAOEZL>        ',UIO,I,7,IER)
               IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) KAOEZ(NAEZ+I)
               CALL IoInput('<RMTREFL>       ',UIO,I,7,IER)
               IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) RMTREF(NAEZ+I)
               WRITE (111,FMT='(3E20.12,3X,F9.6,3X,I5)') (TLEFT(I1,I),I1=1,3),RMTREF(NAEZ+I),KAOEZ(NAEZ+I)
            ENDDO
            WRITE(111,FMT='(A82)') '<RBRIGHT>                                                     <RMTREFL>   <KAOEZL>'
            DO I=1,NRBASIS
               CALL IoInput('<RBRIGHT>       ',UIO,I,7,IER)
               READ (UNIT=UIO,FMT=*) (TRIGHT(I1,I),I1=1,3)
               KAOEZ(NAEZ+NLBASIS+I) = I     ! Default
               CALL IoInput('<KAOEZR>        ',UIO,I,7,IER)
               IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) KAOEZ(NAEZ+NLBASIS+I)
               CALL IoInput('<RMTREFR>       ',UIO,I,7,IER)
               IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) RMTREF(NAEZ+NLBASIS+I)
               WRITE (111,FMT='(3E20.12,3X,F9.6,3X,I5)') (TRIGHT(I1,I),I1=1,3),RMTREF(NAEZ+NLBASIS+I),KAOEZ(NAEZ+NLBASIS+I)
            ENDDO


         ELSE   

            DO I=1,NLBASIS
               TLEFT(:,I) = 0.D0
               CALL IoInput('LEFTBASIS       ',UIO,I,7,IER)
               READ (UNIT=UIO,FMT=*) (TLEFT(I1,I),I1=1,3),II
               KAOEZ(NAEZ+I) = II    ! changed 1.11.99
               write(6,*) 'this must be 1',KAOEZ(NAEZ+I),NAEZ+I

               RMTREF(NAEZ+I) = RMTREFDEF
               CALL IoInput('<RMTREFL>       ',UIO,I,7,IER) ! referen-pot rmt
               IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) RMTREF(NAEZ+I)

            END DO
            DO I=1,NRBASIS
               TRIGHT(:,I) = 0.D0
               CALL IoInput('RIGHBASIS       ',UIO,I,7,IER)
               READ (UNIT=UIO,FMT=*) (TRIGHT(I1,I),I1=1,3),II
               KAOEZ(NAEZ+NLBASIS+I) = II  ! changed 1.11.99

               RMTREF(NAEZ+NLBASIS+I) = RMTREFDEF
               CALL IoInput('<RMTREFR>       ',UIO,I,7,IER) ! referen-pot rmt
               IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) RMTREF(NAEZ+NLBASIS+I)
               
            END DO
         ENDIF


         ! Put The additional atoms in the "embending" positions

         DO I=1,NLBASIS
           DO I1=1,3
             RBASIS(I1,NAEZ+I) = TLEFT(I1,I)
           END DO
         END DO
         DO I=1,NRBASIS
           DO I1=1,3
             RBASIS(I1,NAEZ+NLBASIS+I) = TRIGHT(I1,I)
           END DO
         END DO 
         !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! In RBASIS we have first the basis atoms or the interface
         ! atoms then the left host then the right host the host
         ! goes in the NEMB positions 
         !
         !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
         CALL IoInput('ZPERIODL        ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) (ZPERLEFT(i1),I1=1,3)
         CALL IoInput('ZPERIODR        ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) (ZPERIGHT(i1),I1=1,3) 
         WRITE(6,9430) NLEFT,NLBASIS
         WRITE(6,9440) NRIGHT,NRBASIS
         WRITE(6,9450) (ZPERLEFT(i1),I1=1,3)
         WRITE(6,9460) (ZPERIGHT(i1),I1=1,3)
         WRITE(6,9465)
         WRITE(6,9470)
         DO I=NLEFT,1,-1
            DO I1=NLBASIS,1,-1
            tx = TLEFT(1,i1) + (I-1)*ZPERLEFT(1)
            ty = TLEFT(2,i1) + (I-1)*ZPERLEFT(2)
            tz = TLEFT(3,i1) + (I-1)*ZPERLEFT(3)
            WRITE(6,9420) (I-1)*NLBASIS+i1, tx,ty,tz
            END DO 
         END DO
          WRITE(6,9475)
         DO I=1,NAEZ
            WRITE(6,9420) I, (RBASIS(I1,I),I1=1,3)
         END DO
          WRITE(6,9480)
          DO I=1,NRIGHT
            DO I1=1,NRBASIS
            tx = TRIGHT(1,i1) + (I-1)*ZPERIGHT(1)
            ty = TRIGHT(2,i1) + (I-1)*ZPERIGHT(2)
            tz = TRIGHT(3,i1) + (I-1)*ZPERIGHT(3) 
            WRITE(6,9420) (I-1)*NRBASIS+i1,tx,ty,tz
            END DO 
         END DO  

      END IF                    ! LINTERFACE  

      ! End Interface mode
      ! =============================================================================


      ! =============================================================================
      ! Begin define volume weights
      WRITE(*,*) 'Begin define volume weights'
      ! returned to main program in array sizefac that is also read in from atominfo below

      ! sizefac runs over all atoms, incl. left/right embedded atoms in all layers 
      ! up to nleft/nright, and incl. all imp. atoms (set in subroutine that sets the
      ! impurity details).
      ! However, mtwght goes not over all left/right embedded atoms in all layers 
      ! but only over all left/right embedded basis atoms.
      SIZEFAC(:) = 1.D0      ! Default
 
      ! Wished muffin-tin (units of alat)
      DO ISITE=1,NAEZ
         MTWGHT(ISITE) = DSQRT(SIZEFAC(ISITE)) ! Initialize
         CALL IoInput('<MTWAL>         ',UIO,ISITE,7,IER)
         IF (IER.EQ.0)  READ (UNIT=UIO,FMT=*) MTWGHT(ISITE)
         IF (IER.EQ.0)  WRITE(*,*) 'Read in wished MT radius in ALAT:',MTWGHT(ISITE)
      ENDDO
      ! Wished muffin-tin (atomic units, overrides units of alat)
      DO ISITE=1,NAEZ
         CALL IoInput('<MTWAU>         ',UIO,ISITE,7,IER)
         IF (IER.EQ.0)  THEN 
            READ (UNIT=UIO,FMT=*) MTWGHT(ISITE)
            MTWGHT(ISITE) = MTWGHT(ISITE)/ALATC
            WRITE(*,*) 'Read in wished MT radius in AU:',MTWGHT(ISITE)
         ENDIF
      ENDDO


      ! Same for left/right embedded atoms
      IF (LINTERFACE) THEN
         IL = 0
         DO ISITE = NAEZ + 1,NAEZ + NLBASIS
            MTWGHT(ISITE) = 1.D0 ! Initialize
            IL = IL + 1
            CALL IoInput('<LFMTWAL>       ',UIO,IL,7,IER)
            IF (IER.EQ.0)  READ (UNIT=UIO,FMT=*) MTWGHT(ISITE)
            IF (IER.EQ.0)  WRITE(*,*) 'Read in wished MT radius in ALAT:',MTWGHT(ISITE)
         ENDDO
         IL = 0
         DO ISITE = NAEZ + NLBASIS + 1,NAEZ + NLBASIS + NRBASIS
            MTWGHT(ISITE) = 1.D0 ! Initialize
            IL = IL + 1
            CALL IoInput('<RTMTWAL>       ',UIO,IL,7,IER)
            IF (IER.EQ.0)  READ (UNIT=UIO,FMT=*) MTWGHT(ISITE)
            IF (IER.EQ.0)  WRITE(*,*) 'Read in wished MT radius in ALAT:',MTWGHT(ISITE)
         ENDDO

         ! Wished muffin-tin (atomic units, overrides units of alat)
         IL = 0
         DO ISITE = NAEZ + 1,NAEZ + NLBASIS
            IL = IL + 1
            CALL IoInput('<LFMTWAU>       ',UIO,IL,7,IER)
            IF (IER.EQ.0)  THEN 
               READ (UNIT=UIO,FMT=*) MTWGHT(ISITE)
               MTWGHT(ISITE) = MTWGHT(ISITE)/ALATC
               WRITE(*,*) 'Read in wished MT radius in AU:',MTWGHT(ISITE)
            ENDIF
         ENDDO
         IL = 0
         DO ISITE = NAEZ + NLBASIS + 1,NAEZ + NLBASIS + NRBASIS
            IL = IL + 1
            CALL IoInput('<RTMTWAU>       ',UIO,IL,7,IER)
            IF (IER.EQ.0)  THEN 
               READ (UNIT=UIO,FMT=*) MTWGHT(ISITE)
               MTWGHT(ISITE) = MTWGHT(ISITE)/ALATC
               WRITE(*,*) 'Read in wished MT radius in AU:',MTWGHT(ISITE)
            ENDIF
         ENDDO
         IL = 1
      ENDIF


      ! Re-define size-factor as muffin-tin squared (units of alat)
      DO ISITE=1,NAEZ
         SIZEFAC(ISITE) = MTWGHT(ISITE)**2
      ENDDO

      ! Define weights for left/right embedded atoms
      IF (LINTERFACE) THEN
         ISITE = 0
         DO II = 1,NLEFT
            DO I1 = NAEZ + 1,NAEZ + NLBASIS
               ISITE = ISITE + 1
               SIZEFAC(-ISITE) = MTWGHT(I1)**2
            ENDDO
         ENDDO
         DO II = 1,NRIGHT
            DO I1 = NAEZ + NLBASIS + 1,NAEZ + NLBASIS + NRBASIS
               ISITE = ISITE + 1
               SIZEFAC(-ISITE) = MTWGHT(I1)**2
            ENDDO
         ENDDO
      ENDIF
      WRITE(*,*) 'End Define volume weights'
      ! End Define volume weights
      ! =============================================================================

      
      WRITE(*,*) 'End Structure'
! End Structure    
! =============================================================================

! =============================================================================
! Start clusters setup
      CALL IoInput('RCLUSTZ         ',UIO,IL,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) RCUTZ
      IF (IER.NE.0) THEN 
         WRITE(*,*) 'readinput12: RCLUSTZ not found, stopping.'
         STOP 'readinput12: RCLUSTZ not found, stopping.'
      END IF

      RCUTXY = RCUTZ
      CALL IoInput('RCLUSTXY        ',UIO,IL,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) RCUTXY

      ! RMT of ref-pot in TB-KKR to determine TB-KKR clusters
      LMTREF = .FALSE.
      DO ISITE=1,NAEZ
         RMTREF(ISITE) = RMTREFDEF ! Initialize
         CALL IoInput('<RMTREF>        ',UIO,ISITE,7,IER)
         IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) RMTREF(ISITE)
         IF (IER.EQ.0) LMTREF = .TRUE.
         RMTREFDEF = RMTREF(1)  ! Re-define default in case it is only included for 1st atom.
      ENDDO


      WRITE(6,*) 'Parameters Used for the cluster calculation'
      if (abs(rcutz-rcutxy).lt.1.d-4) then
         write(6,*) 'Clusters inside spheres with radius R = ',rcutz
      else
         write(6,*) 'Clusters inside cylindels with '
         write(6,*) 'Rz = ',rcutz,' Rxy = ',rcutxy
      end if
      write(6,2018)                 ! rbasis
      write(6,2025) ((rbasis(j,i),j=1,3),i,i=1,naez)
      if (nemb.gt.0) write(6,*) 
      write(6,2031) ((rbasis(j,i),j=1,3),i,refpot(kaoez(i)),i=naez+1,naez+nemb)

! End clusters setup
! =============================================================================

! =============================================================================
! Start chemistry

      IL=1
      NSPIN = 1
      CALL IoInput('NSPIN           ',UIO,1,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) NSPIN
      IF (IER.NE.0) WRITE(*,*) 'NSPIN not found, setting NSPIN=1.'
      IF (NSPIN.NE.1.AND.NSPIN.NE.2) THEN
         WRITE(*,*) 'readinput: NSPIN not 1 or 2',NSPIN
         STOP 'readinput: incomprehensible value for NSPIN'
      ENDIF



      NATYP = NAEZ
      CALL IoInput('NATYP           ',UIO,1,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) NATYP
      IF (IER.NE.0) WRITE(*,*) 'NATYP not found, setting NATYP=NAEZ.'
      IF (NATYP.LT.NAEZ) THEN 
         WRITE(*,*) 'readinput: NATYP.LT.NAEZ',NATYP,NAEZ
         STOP 'readinput: NATYP.LT.NAEZ'
      ENDIF

      WRITE(6,8000) NAEZ,NATYP,NSPIN
 8000 FORMAT('NAEZ=',I5,' ; NATYP=',I5' ; NSPIN=',I2)


      RMTCORE(1:NAEZ) = 1.D10     ! A large number to be reduced to touching MT later
      ZAT(1:NAEZ) = 29.D0         ! Default is copper.


      CALL IoInput('ATOMINFO        ',UIO,I+1,7,IER)
      IF (IER.EQ.0) THEN
         DO I=1,NAEZ
            CALL IoInput('ATOMINFO        ',UIO,I+1,7,IER)
                IF (IER.EQ.0) READ (UNIT=UIO,FMT=*)    ZAT(I), &
     &                           LMXC(I), &
     &                          (KFG(J,I),J=1,4), &
     &                           CLS(I), &
     &                           REFPOT(I), &
     &                           NTCELL(I), &
     &                           MTFAC(I), &
     &                           IRNS(I),RMTCORE(I),SIZEFAC(I)
         END DO
       ENDIF

      ! Read ZAT (overrides ATOMINFO read-in)
      DO I = 1,NATYP
        CALL IoInput('<ZATOM>         ',UIO,I,7,IER)
         IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) ZAT(I)
      ENDDO
      ! Read RMTCORE (overrides ATOMINFO read-in)
      DO I = 1,NAEZ
        CALL IoInput('<RMTCORE>       ',UIO,I,7,IER)
         IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) RMTCORE(I)
      ENDDO

      DO IAT = 1,NATYPD
         SITEAT(IAT) = IAT  ! For mapping in jellstart if there is no cpa mode
      ENDDO

! CPA mode:
      IF (NATYP.GT.NAEZ) THEN ! CPA mode, read site for each atom

         SITEAT(1:NATYPD) = 0  ! reset the site index per atom, it will be read in in case of cpa.

         IER = 0
         CALL IoInput('<SITE>          ',UIO,1,7,IER)
         IF (IER.NE.0) THEN
            WRITE(*,*) 'readinput: <SITE> in CPA mode not found, stopping.'
            STOP 'readinput: <SITE> in CPA mode not found.'
         ENDIF

         WRITE(111,FMT='(A18)') '<SITE>  '
         DO IAT = 1,NATYP
            CALL IoInput('<SITE>          ',UIO,IAT,7,IER)
            READ (UNIT=UIO,FMT=*) SITEAT(IAT)
            WRITE(111,FMT='(I5)') SITEAT(IAT)
            IF (SITEAT(IAT).GT.NAEZ.OR.SITEAT(IAT).LT.1) THEN
               WRITE(*,*) 'readinput: CPA SITE ill-defined',IAT,SITEAT(IAT)
               STOP 'readinput: CPA SITE ill-defined'
            ENDIF
         ENDDO

         NOQ(1:NAEZD) = 0 
         DO IAT = 1,NATYP
            ISITE = SITEAT(IAT)
            NOQ(ISITE) = NOQ(ISITE) + 1
         ENDDO

         DO ISITE=1,NAEZ
            IF (NOQ(ISITE).LT.1) THEN
               WRITE(6,*) 'readinput: CPA: SITE',ISITE,'HAS NO ASSIGNED ATOM'
               STOP 'readinput: CPA, atoms not properly assigned.'
            ENDIF
         END DO

      ENDIF
! End CPA mode        


      ! read in whished value of Fermi level (core state energies of starting potential is shifted accordingly)
      EFSET = -1.0d0 ! default value -1 signals no shift (EF=0.4...)
      CALL IoInput('EFSET           ',UIO,1,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) EFSET
      WRITE(*,*) 'readinput: EFSET=',EFSET
      WRITE(111,*) 'EFSET= ',EFSET


      WRITE(6,2028) NATYP
      WRITE(6,2104)
      WRITE(6,1029) ( &
     &     ZAT(I), &
     &     LMXC(I), &
     &     (KFG(J,I),J=1,4), &
     &     CLS(I), &
     &     REFPOT(I), &
     &     NTCELL(I), &
     &     MTFAC(I), &
     &     IRNS(I),RMTCORE(I),SIZEFAC(I),I=1,NATYP)
      WRITE(6,2108)
      WRITE(6,2104)

! End  chemistry
! =============================================================================


! =============================================================================
! Start control

      LMAX = 3
      CALL IoInput('LMAX            ',UIO,0,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) LMAX
      IF (IER.NE.0) WRITE(*,*) 'LMAX not found, setting LMAX=3.'

      KSHAPE = 2 ! Default
      CALL IoInput('KSHAPE          ',UIO,1,7,IER)
      IF (IER.EQ.0)   READ (UNIT=UIO,FMT=*) kshape

      IRM = 484 ! Default
      CALL IoInput('IRM             ',UIO,1,7,IER)
      IF (IER.EQ.0)   READ (UNIT=UIO,FMT=*) irm

      INS = 1 ! Default
      CALL IoInput('INS             ',UIO,1,7,IER)
      IF (IER.EQ.0)   READ (UNIT=UIO,FMT=*) ins

      CALL IoInput('NMIN            ',UIO,1,7,IER)
      IF (IER.EQ.0)   READ (UNIT=UIO,FMT=*) NMIN
      CALL IoInput('NRAD            ',UIO,1,7,IER)
      IF (IER.EQ.0)   READ (UNIT=UIO,FMT=*) NRAD
      NRAD = MAX(NMIN,NRAD)
      CALL IoInput('NSMALL          ',UIO,1,7,IER)
      IF (IER.EQ.0)   READ (UNIT=UIO,FMT=*) NSMALL
      NSMALL = MAX(NMIN,NSMALL)
      
      ! Tolerance for voronoi construction, defaults in maindriver data
      CALL IoInput('<TOLHS>         ',UIO,1,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) TOLHS

      CALL IoInput('<TOLVD>         ',UIO,1,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) TOLVDIST

      CALL IoInput('<TOLAREA>       ',UIO,1,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) TOLAREA


      IF (INS.GT.0) THEN
         WRITE (6,FMT=9130)
         WRITE (6,FMT=9140)
         DO I = 1,NATYP
            WRITE (6,FMT=9150) I,IRNS(I),IRNSD

            IF (IRNS(I).GT.IRNSD) CALL RCSTOP('19      ')

         ENDDO

         IF (LMAX.GT.LMAXD) THEN
            WRITE (6,FMT=9120)

            CALL RCSTOP('20      ')

         END IF

      END IF

      LMMAX = (LMAX+1)**2
      LPOT  = MIN(2*LMAX,LPOTD)
      LMPOT = (LPOT+1)* (LPOT+1)

! End control
! =============================================================================


      write(6,2101)
      write(6,2012) (kaoez(i),i=1,naez)
!
! ------------------------------------------------------------------------
      do i=1,naez
        if (kaoez(i).lt.1) STOP 'Error in KAOEZ'
      enddo
      NCLS = 0
      NREF = 0
!
      DO I=1,NATYP 
        NCLS = MAX(NCLS,CLS(I)) 
      ENDDO
      DO I=1,NATYP 
        NREF = MAX(NREF,REFPOT(I)) 
      ENDDO
 
      WRITE(6,2016) NCLS,NREF
! ------------------------------------------------------------------------

 
! =============================================================================
! Start read run and test options and files
 
           
      CALL IoInput('TESTOPT         ',UIO,1,7,IER)
                     READ(UNIT=UIO,FMT=980)(TESTC(i),i=1,8)
      CALL IoInput('TESTOPT         ',UIO,2,7,IER)
                     READ(UNIT=UIO,FMT=980)(TESTC(8+i),i=1,8)
      WRITE(6,52) (TESTC(I),I=1,16)                                                 
 52   FORMAT(79('-')/' TEST OPTIONS:'/2(1X,A8,7('//',A8)/)/79('-'))
 
 
      CALL IoInput('RUNOPT          ',UIO,1,7,IER)
                   READ (UNIT=UIO,FMT=980)(OPTC(i),i=1,8)
      WRITE(6,62) (OPTC(i),i=1,8)                                              
 62   FORMAT(79('-')/' EXECUTION OPTIONS:'/1X,A8,7('//',A8)/79('-'))
 980  FORMAT(8A8)
 
      IL=1
      I12='                                        '
      CALL IoInput('FILES           ',UIO,IL,7,IER)
                     IF(IER.EQ.0) READ (UNIT=UIO,FMT='(A40)')  I12
      I13='                                        '
      CALL IoInput('FILES           ',UIO,IL+1,7,IER)
                     IF(IER.EQ.0) READ (UNIT=UIO,FMT='(A40)')  I13
      I40='                                        '
      CALL IoInput('FILES           ',UIO,IL+2,7,IER)
                     IF(IER.EQ.0) READ (UNIT=UIO,FMT='(A40)')  I40
      I19='                                        '
      CALL IoInput('FILES           ',UIO,IL+3,7,IER)
                     IF(IER.EQ.0) READ (UNIT=UIO,FMT='(A40)')  I19

      write(6,*) 'I12="',I12,'"'
      write(6,*) 'I13="',I13,'"'
      write(6,*) 'I19="',I19,'"'

! End read run and test options and files
! =============================================================================

      CLOSE(111)

      WRITE(6,2110)
      WRITE(6,*) ' >>>>>>>>> RINPUT99 EXITS NOW <<<<<<<<<< '

      RETURN
! *********************************************Input-End ********
 1000 FORMAT(10I4)
 1001 FORMAT(3F12.9)
 1002 FORMAT(80A1)
 1003 FORMAT(10L4)
 1004 FORMAT(I4)
 1005 FORMAT(4I4)
 1006 FORMAT(3F12.9,I4)
 1029 FORMAT((F4.0,I4,4x,4I1,3I4,F8.4,I4,E12.4,F8.4))
! ------------------------------------------------------------------------
 2001 FORMAT(/,' ATOM',I3/)
 2002 FORMAT(/,1X,80(A1)///' INPUT DATA :  (NBG=1: LDA; NBG=2: LSDA)')
 2004 FORMAT( 79(1H=)/ &
     &     'I',77X,'I'/ &
     &     'I',4X,' Screened Korringa-Kohn-Rostoker ', &
     &             'Electronic Structure Code',10X,'I'/ &
     &     'I',4X,'for Bulk and Interfaces',20X,'I'/ &
     &     'I',77X,'I',/,'I',4X,'        ',A8,57X,'I'/ &
     &     'I',77X,'I',/,'I',4X,'Version ',A8,57X,'I'/ &
     &     'I',77X,'I'/ &
     &     'I',77X,'I',/,79(1H=),/,/,1X,80A1 )
 2010 FORMAT(' NSPIN '/I4)
 2011 FORMAT(' NSTEPS'/I4)
 2012 FORMAT(' KAOEZ '/,(10I4))
 2013 FORMAT('      M2    MMIN    MMAX    SINN', &
     &       '    SOUT     RIN    ROUT'/7I8)
 2014 FORMAT('          ALATC          BLATC          CLATC'/3F15.8)
 2015 FORMAT('   INTERVX   INTERVY   INTERVZ'/3I10)
 2016 FORMAT('    NCLS    NREF   '/,3I8)
 2017 FORMAT(' COMPLX'/(3X,L1))
 2018 FORMAT(' RBASIS')
 2019 FORMAT('         ABASIS         BBASIS         CBASIS'/3F15.8)
 2020 FORMAT(' NPRINCD  NLAYER'/,2I8)
 2021 FORMAT(' INIPOL'/,(10I4))
 2022 FORMAT(' IXIPOL'/,(10I4))
 2023 FORMAT('    NAEZ    NEMB   NEMBZ'/,3I8)
 2024 FORMAT(' NZ    '/,I4)
 2025 FORMAT((3F15.8,I6))
 2026 FORMAT(' COFINV'/,3F15.8)
 2027 FORMAT(' LATT  '/,I4)
 2028 FORMAT &
     &(' NATYP '/,I4/,'   Z lmx     KFG cls pot ntc  MTFAC irns   MT Size')
 2030 FORMAT(' EQINV '/,(10I4))
 2031 FORMAT((3F15.8,2I6))
 2040 FORMAT('  IFULLD  ISPARD  ISLABD'/,4I8)
 2050 FORMAT(' Dimension and Input Data CHECK')
! ------------------------------------------------------------------------
 2100 FORMAT(79(1H-))
 2101 format(   3(1H-),1H+  , 3(14(1H-),1H+),  30(1H-))
 2102 format( 3(9(1H-),1H+) ,49(1H-))
 2103 FORMAT(10(3(1H-),1H+) ,39(1H-))
 2104 format(   3(1H-),1H+  ,75(1H-))
 2105 format( 3(3(1H-),1H+) ,67(1H-))
 2106 format( 5(3(1H-),1H+) ,59(1H-))
 2107 format( 3(14(1H-),1H+),34(1H-))
 2108 format( 2(3(1H-),1H+),  7(1H-),1H+,      3(3(1H-),1H+), &
     &          7(1H-),1H+,   3(1H-),1H+,      39(1H-))
 2109 format( 5(7(1H-),1H+) ,39(1H-))
 2110 format( 3(7(1H-),1H+) ,55(1H-))
 2111 format( 7(7(1H-),1H+) ,23(1H-))
 2112 format( 2(7(1H-),1H+) ,63(1H-))
 9000 FORMAT (I2,3X,4I5)
 9010 FORMAT (1X,I1,1X,4I1)
 9020 FORMAT (/,33x,'check of dimension-data consistency',/,33x, &
     &       35 ('-'),/,40x,'lmax   : (',i6,',',i6,')',/,40x, &
     &       'natyp  : (',i6,',',i6,')',/,40x,'irm    : (',i6,',',i6, &
     &       ')',/,40x,'nspin  : (',i6,',',i6,')',/)
 9030 FORMAT (1x,10 ('*'),' external magnetic field applied hfield=', f8.5)
 9040 FORMAT (3f12.7)
 9050 FORMAT (20x,a4,'spin polarized calculation')
 9060 FORMAT (8i4)
 9070 FORMAT (1x,20x,' calculation with',a8,'-potential')
 9080 FORMAT (1x,79 ('*'))
 9090 FORMAT (' mixing factor used           :',f15.6,/, &
     &        ' convergence quality required :',1p,d15.2)
 9100 FORMAT (1x,20x,a24,'exchange-correlation potential')
 9110 FORMAT (/,20x,'broyden"s method # :',i3, &
     &       ' is used up to iteration-      ',/,20x,'depth :',i3, &
     &       '  then jacobian is fixed and potential      ',/,20x, &
     &       'is updated using that jacobian')
 9120 FORMAT (13x,' in case of calculating non - spherical wavefcts ', &
     &       'the parameter lmaxd has to be less-equal lmax !')
 9130 FORMAT (/)
 9140 FORMAT (20x,'full potential calculation ', &
     &       '- cut off of non spherical potential',/,' >',/)
 9150 FORMAT (31x,'representive atom no.',i3,' irns :',i5,' irnsd :',i5)
 9160 FORMAT (21x,a43)
 9170 FORMAT (21x,a43)
 9180 FORMAT (2i5)
 9190 FORMAT (3f12.7,/,4i4)
 9200 FORMAT (3i4,1f12.7)
 9210 FORMAT (' lmax'/,i4)
 9220 FORMAT ('          E1          E2          TK'/,3f12.6)
 9230 FORMAT ('   NPOL  NPNT1  NPNT2  NPNT3'/,4i7)
 9240 FORMAT (' IRNUMX ITCCOR IPRCOR'/,3i7)
 9250 FORMAT ('  IFILE    IPE ISHIFT ESHIFT'/,3i7,f12.6)
 9260 FORMAT (' KSHAPE    IRM    INS   '/,3i7)
 9280 FORMAT (' external magnetic hfield     :',f15.4/, &
     &        ' VCONST                       :',f15.6)
 9290 FORMAT ('   IMIX IPOTOU    IGF    ICC'/,4i7)
 9300 FORMAT (' ITDBRY'/,i7)
 9310 FORMAT ('      STRMIX        FCM       QBOUND'/,3f12.6)
 9320 FORMAT ('      BRYMIX'/,f12.6)
 9330 FORMAT ('    KTE   KPRE   KEFG  KVMAD KSCOEF'/,5i7)
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
 9420 format(I5,3F12.6)
 9430 format('Number of LEFT  Host Layers : ',I5,' with ',I5,' basis')
 9440 format('Number of RIGHT Host Layers : ',I5,' with ',I5,' basis')
 9450 format('Left  side periodicity : ',3F10.5)
 9460 format('Right side periodicity : ',3F10.5)
 9465 format('    Geommetry used : '/, &
     &       ' ATOM       TX          TY          TZ ')
 9470 format('--------------- Left  Host -------------- ')
 9475 format('---------------   S L A B  -------------- ')
 9480 format('--------------- Right Host -------------- ')
   END SUBROUTINE READINPUT





