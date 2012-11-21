c ************************************************************************
      SUBROUTINE STARTB1(IFILE,IPF,IPFE,IPE,KHFELD,
     &                   NBEG,NEND,
     &                   RMTNEW,RMT,ITITLE,HFIELD,IMT,IRC,VCONST,
     &                   IRNS,LPOT,NSPIN,IRMIN,NTCELL,IRCUT,IPAN,
     &                   THETAS,IFUNM,NFU,LLMSP,LMSP,EFERMI,
     &                   VBC,RWS,LCORE,NCORE,DRDI,
     &                   R,ZAT,A,B,IRWS,INIPOL,IINFO,
     &                   IPAND, IRID, NFUND, IRMD, NCELLD,
     &                   NAEZD, IRNSD)
c ************************************************************************
c   reads the input potentials
c
c    units :       ry - units for energy
c                  the lattice constant and all other lengths
c                                           given in bohr units
c                  the planck constant h/2pi=1
c                  the electron charge e=sqrt(2)
c                  the electron mass m=1/2
c                  the speed of light c = 2/alpha = 274.0720442
c                      with the fein structure constant alpha
c
c    remember that the input potentials do not include the electro-
c             static contribution of the nucleus of the cell itself
c             this has to be added explicitly !
c
c   as input is used: lmax=maximum angular momentum
c                    nbeg .. nend=number of different atoms
c
c
c     in case of shape corrections this routine  reads from unit 19
c     a suitable radial  mesh 'xrn',its derivate 'drn' and the shape
c     functions 'thetas' .          thus, the region from the muffin
c     tin to the circumscribed  sphere radii is divided  into  'npan'
c     pannels, each one containing 'nm(ipan)' points in order to take
c     care of the  discontinuities of the shape-function  derivative.
c     at the output one obtains :
c            llmsp (icell,ifun)       = integer array giving the com-
c                                       posite  index  lm=l*(l+1)+m+1
c                                       of the ifun-th shape function
c            lmsp  (icell,lm)         = (0,1)  if the lm-th component
c                                       is vanishing or not
c            nfu   (icell)            = number  of   shape   function
c                                       components in cell 'icell'
c
c
c     modified for bandstructure code
c
c                                 b.drittler nov. 1989
c-----------------------------------------------------------------------
C     .. Parameters ..
      IMPLICIT NONE

C      include 'inc.p'
C      from inc.p
      INTEGER, INTENT(IN) :: IPAND
      INTEGER, INTENT(IN) :: IRID
      INTEGER, INTENT(IN) :: NFUND
      INTEGER, INTENT(IN) :: IRMD
      INTEGER, INTENT(IN) :: NCELLD
      INTEGER, INTENT(IN) :: NAEZD
      INTEGER, INTENT(IN) :: IRNSD

C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT,CVLIGHT,EFERMI,HFIELD,VBC(*),VCONST
      INTEGER IFILE,IINFO,IPE,IPF,IPFE,
     &        KHFELD,
     &        LPOT,
     &        NBEG,NEND,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),B(*),DRDI(IRMD,*),ECORE(20,2),
     &                 R(IRMD,*),RMT(*),RMTNEW(*),
     &                 RWS(*),THETAS(IRID,NFUND,NCELLD),
C                      VINS(IRMIND,LMPOTD,2)
     &                 VINS(IRMD-IRNSD:IRMD,(LPOT+1)**2,2),
     &                 VISP(IRMD,2),
     &                 ZAT(*),ZATINFO(NAEZD)

C             IFUNM(LMXSPD,NAEZD)
      INTEGER IFUNM((2*LPOT+1)**2,NAEZD),
     &        IMT(*),INIPOL(*),IPAN(*),
     &        IRC(*),IRCUT(0:IPAND,*),
     &        IRMIN(*),IRNS(*),IRWS(*),ITITLE(20,*),
     &        LCORE(20,*),LLMSP(NFUND,NAEZD),
     &        LMSP((2*LPOT+1)**2,NAEZD),
     &        NCORE(*),NFU(NCELLD),NTCELL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,B1,EA,EFNEW
      INTEGER I,IA,ICELL,ICORE,IFUN,IH,IMT1,INEW,IO,IPAN1,IR,IRI,
     &        IRMINM,IRMINP,IRNS1P,IRT1P,IRWS1,ISAVE,ISPIN,ISUM,
     &        J,
     &        LM,LM1,LMPOT,LMPOTP,
     &        N,NCELL,NFUN,NR
      LOGICAL TEST
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DRN(IRID,NCELLD),SCALE(NCELLD),U(IRMD),
     &                 XRN(IRID,NCELLD)

      INTEGER MESHN(NCELLD),NM(IPAND,NCELLD),NPAN(NCELLD)
C     ..
C     .. External Subroutines ..
      EXTERNAL CALRMT,RINIT,TEST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ANINT,EXP,LOG,MAX,MOD,REAL,SQRT

      INTEGER ISHAPE


      INTEGER LMPOTD
      INTEGER IRMIND
      INTEGER LMXSPD
      INTEGER LRECPOT

C     ..
c-----------------------------------------------------------------------
c
      ISHAPE = 0

      IRMIND= IRMD-IRNSD
      LMPOTD= (LPOT+1)**2
      LMXSPD= (2*LPOT+1)**2

C     I/O Record length for potential file
      LRECPOT=8*(LMPOTD*(IRNSD+1)+IRMD+20)

c ---> output of radial mesh information
c
      IO = 0
      IF (IINFO.NE.0 .AND. TEST('RMESH   ')) IO = 1
c
c---> set speed of light
c
      CVLIGHT = 274.0720442D0

C
C intitialize potential arrays
C
      VINS = 0.0D0
      VISP = 0.0D0
      ECORE = 0.0D0

C     more initialisations
      THETAS = 0.0d0

c-----------------------------------------------------------------------
c
c---> read radial mesh information of the shape functions and
c     shape functions THETAS in the first iteration - if needed
c
      IF (ISHAPE.EQ.0) THEN
        ISHAPE = 1
        READ (19,FMT=9000) NCELL
        WRITE (6,FMT=*) '  ncell : ',NCELL,NCELLD
c
        IF(NCELL.NE.NCELLD) THEN
          WRITE(6,*) 'Please, change the parameter ncelld (',NCELLD,
     +         ') in inc.p to',NCELL
          STOP 'STARTB - NCELLD'
        ENDIF
c
        READ (19,FMT=9010) (SCALE(ICELL),ICELL=1,NCELL)
        DO 30 ICELL = 1,NCELL
          READ (19,FMT=9000) NPAN(ICELL),MESHN(ICELL)
c
          IF(NPAN(ICELL)+1.GT.IPAND) THEN
            WRITE(6,*) 'Please, change the parameter ipand (',IPAND,
     +           ') in inc.p to',NPAN(ICELL)+1
            STOP 'STARTB - IPAND'
          ENDIF
c
          IF(MESHN(ICELL).GT.IRID) THEN
            WRITE(6,*) 'Please, change the parameter irid (',IRID,
     +           ') in inc.p to',MESHN(ICELL)
            STOP 'STARTB - IRID'
          ENDIF
c
          READ (19,FMT=9000) (NM(IPAN1,ICELL),IPAN1=2,NPAN(ICELL)+1)
          READ (19,FMT=9010) (XRN(IR,ICELL),DRN(IR,ICELL),IR=1,
     +      MESHN(ICELL))

          READ (19,FMT=9000) NFU(ICELL)
          NFUN = NFU(ICELL)
          WRITE (6,FMT=*) '  nfun  : ',NFUN,NFUND
c
          IF(NFUN.GT.NFUND) THEN
            WRITE(6,*) 'Please, change the parameter nfund (',NFUND,
     +           ') in inc.p to',NFUN
            STOP 'STARTB - NFUND'
          ENDIF
c
          DO 10 LM = 1,LMXSPD
            LMSP(LM,ICELL) = 0
   10     CONTINUE

          DO 20 IFUN = 1,NFUN
            READ (19,FMT=9000) LM
            LLMSP(IFUN,ICELL) = LM
            LMSP(LM,ICELL) = 1
            IFUNM(LM,ICELL) = IFUN
            READ (19,FMT=9010) (THETAS(N,IFUN,ICELL),N=1,MESHN(ICELL))
   20     CONTINUE

   30   CONTINUE
      END IF 
c-----------------------------------------------------------------------
c
      LMPOT = (LPOT+1)* (LPOT+1)
C
      OPEN (66,ACCESS='direct',RECL=LRECPOT*2,FILE='vpotnew',
     +     FORM='unformatted')
C
      DO IH = NBEG,NEND
        ZATINFO(IH) = ZAT(IH)
      ENDDO

      DO 150 IH = NBEG,NEND
        DO 140 ISPIN = 1,NSPIN
          I = NSPIN* (IH-1) + ISPIN
C
          IF (IFILE.NE.0) THEN
            IRCUT(0,IH) = 0
            ICELL = NTCELL(IH)
            IPAN(IH) = 1 + NPAN(ICELL)
c
c---> read title of potential card
c
            READ (IFILE,FMT=9020) (ITITLE(IA,I),IA=1,20)
            IF (IINFO.NE.0.AND.I.LE.8) THEN 
                WRITE (6,FMT=9080) (ITITLE(IA,I),IA=1,20)
            END IF
c
c---  >read muffin-tin radius , lattice constant and new muffin radius
c      (new mt radius is adapted to the given radial mesh)
c
            READ (IFILE,FMT=9030) RMT(IH),ALAT,RMTNEW(IH)
c
c---> read nuclear charge , lmax of the core states ,
c     wigner seitz radius , fermi energy and energy difference
c     between electrostatic zero and muffin tin zero
c
            READ (IFILE,FMT=9040) ZAT(IH),RWS(IH),EFNEW,VBC(ISPIN)
c
c---> if efermi .eq. 0 use value from in5 removed E.R.
c
c            IF (EFNEW.NE.0.0D0 .AND. I.EQ.1) EFERMI = EFNEW
            IF (I.EQ.1) EFERMI = EFNEW
c
c---> read : number of radial mesh points
c     (in case of ws input-potential: last mesh point corresponds
c     to ws-radius, in case of shape-corrected input-potential
c     last mesh point of the exponential mesh corresponds to
c     mt-radius/nevertheless this point is always in the array
c     irws(ih)),number of points for the radial non-muffin-tin
c     mesh  needed for shape functions, the constants a and b
c     for the radial exponential mesh : r(i) = b*(exp(a*(i-1))-1)
c     the no. of different core states and some other stuff
c
            READ (IFILE,FMT=9050) IRWS(IH),A(IH),B(IH),NCORE(I),INEW
c
            NR = IRWS(IH)

            IF (NR.GT.IRMD) THEN
              write(6,*) 'Increase parameter IRMD in ''inc.p''',
     +             ' to a value .ge. ',NR,' (= IRWS(',IH,')).'
              STOP 'STARTB1 - IRWS'
            END IF
c
c---> read the different core states : l and energy

C           check: e.r.
            if (ncore(I).GT.20) THEN
              write(*,*) "Error: More than 20 core states."
              STOP
            endif

c
            IF (NCORE(I).GE.1) THEN
                DO ICORE=1,NCORE(I)
                    READ (IFILE,FMT=9070) LCORE(ICORE,I),
     &                   ECORE(ICORE,ISPIN)
                END DO
            END IF
c     
c--->  read full potential - the non spherical contribution from irmin
c      to irt - remember that the lm = 1 contribution is multiplied by
c      1/sqrt(4 pi)
c
            READ (IFILE,FMT=9090) IRT1P,IRNS1P,LMPOTP,ISAVE
            IRMINP = IRT1P - IRNS1P

c     do some checks to get the numerous problems with inconsistencies
c     under control  e.r.

c           crash if IRMINP < IRMIND
c           crash if IRT1P > IRMD
c           crash if LMPOTP inconsistent

            if (irminp < irmind) then
              write (*,*) "startb1 error: IRMINP < IRMIND"
              write (*,*) "potential entry: ", IH
              stop
            endif

            if (irt1p > irmd) then
              write (*,*) "startb1 error: IRTP1 > IRMS"
              write (*,*) "potential entry: ", IH
              stop
            endif

            if (lmpotp /= (LPOT+1)**2) then
              write(*,*) "startb1 error: lmpotp /= (LPOT+1)**2"
              write (*,*) "potential entry: ", IH
              stop
            endif

            if (irnsd < irns1p) then
              write(*,*) "startb1 error: irnsd < irnsp1"
              write (*,*) "potential entry: ", IH
              stop
            endif

c           assign irns here
            IRNS(IH) = IRNS1P

            IRMINM = MAX(IRMINP,IRMIND)
            READ (IFILE,FMT=9100) (VISP(IR,ISPIN),IR=1,NR)
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
                    DO 40 IR = IRMINM,NR
                        VINS(IR,LM1,ISPIN) = U(IR)
   40                 CONTINUE
                    END IF

                  END IF

                END IF

   50         CONTINUE

            END IF

            IRWS1 = IRWS(IH)
c
c---> redefine new mt-radius in case of shape corrections
c
            RMTNEW(IH) = SCALE(ICELL)*ALAT*XRN(1,ICELL)
            IMT1 = ANINT(LOG(RMTNEW(IH)/B(IH)+1.0D0)/A(IH)) + 1

C     E.R. try to set correct muffin-tin index
C     there are IRWS (not IRMD, IRWS <= IRMD)
C     radial mesh points and MESHN(ICELL) points belong
C     to the outer muffin-tin region where shape functions apply
C     Therefore IMT1 must be:
C           IMT1 = IRMD - MESHN(ICELL)
            IMT1 = IRWS1 - MESHN(ICELL)

C     TODO FIXME CHECK: therefore IRWS1 must be =IRMD ??? No!

c
c---> for proper core treatment imt must be odd
c     shift potential by one mesh point if imt is even
c
            IF (MOD(IMT1,2).EQ.0) THEN
              IMT1 = IMT1 + 1
              DO 60 IR = IMT1,2,-1
                VISP(IR,ISPIN) = VISP(IR-1,ISPIN)
   60         CONTINUE
            END IF
c
            IMT(IH) = IMT1
            B(IH) = RMTNEW(IH)/ (EXP(A(IH)*REAL(IMT1-1))-1.0D0)
c
c---> generate radial mesh - potential only is stored in potential card
c     INEW = 1
c     p. zahn, jan. 99
c
            A1 = A(IH)
            B1 = B(IH)
            R(1,IH) = 0.0D0
            DRDI(1,IH) = A1*B1

C all the mesh points are first filled with an exponential mesh
C then the outer meshpoints are overwritten by a equidist. mesh
C in next code block E.R.

            DO 70 IR = 2,IRWS1
              EA = EXP(A1*REAL(IR-1))
              R(IR,IH) = B1* (EA-1.0D0)
              DRDI(IR,IH) = A1*B1*EA
   70       CONTINUE

c
c---> fill cell-type depending mesh points in the non-muffin-tin-region
c
            DO 80 IRI = 1,MESHN(ICELL)
              IR = IRI + IMT1
              R(IR,IH) = SCALE(ICELL)*ALAT*XRN(IRI,ICELL)
              DRDI(IR,IH) = SCALE(ICELL)*ALAT*DRN(IRI,ICELL)
   80       CONTINUE

            RWS(IH) = R(IRWS1,IH)
c
            CALL CALRMT(IPF,IPFE,IPE,IMT(IH),ZAT(IH),RMT(IH),RWS(IH),
     +                  RMTNEW(IH),ALAT,DRDI(1,IH),A(IH),B(IH),IRWS1,
     +                  R(1,IH),IO)
c
            IRCUT(1,IH) = IMT(IH)

CDEBUG E.R.
            if (IRMD < IRWS1) then
              write(*,*) "ERROR in STARTB1"
              write(*,*) "IRMD>=IRWS1 failed."
              write(*,*) "Atom, IRMD, IRWS1"
              write(*,*) IH, IRMD, IRWS1
              stop
            end if
CDEBUG


CDEBUG E.R.
C           if ((IRMD - IMT1) /= MESHN(ICELL)) then
C             write(*,*) "ERROR in STARTB1"
C             write(*,*) "Assertion IRMD - IMT1 == MESHN(ICELL) failed."
C             write(*,*) "Atom, IRMD, IMT1, MESHN(ICELL)"
C             write(*,*) IH, IRMD, IMT1, MESHN(ICELL)
C             stop
C           end if
CDEBUG
            ISUM = IMT(IH)
            DO 90 IPAN1 = 2,IPAN(IH)
              ISUM = ISUM + NM(IPAN1,ICELL)
              IRCUT(IPAN1,IH) = ISUM
   90       CONTINUE
            NR = ISUM
c
            IRC(IH) = IRCUT(IPAN(IH),IH)

CDEBUG E.R.
            if (IRC(IH) > IRMD) then
              write(*,*) "STARTB1: Assertion IRC(IH) <= IRMD failed."
              stop
            end if
CDEBUG
c
c---> fill array irmin in case of full potential
c
            IRMIN(IH) = NR - IRNS(IH)

CDEBUG E.R.
      if (IRMIN(IH) < (IRMD - IRNSD)) then
        write(*,*) "STARTB1: Assertion IRMIN(IH)>=(IRMD - IRNSD) fail."
        stop
      end if
CDEBUG

c
c--->  first iteration : shift all potentials (only for test purpose)
            DO 120 J = 1,NR
              VISP(J,ISPIN) = VISP(J,ISPIN) + VCONST
  120       CONTINUE

          END IF                    ! (IFILE.NE.0)

c
          IF (KHFELD.EQ.1 ) THEN
c          IF (KHFELD.EQ.1 .AND. NSPIN.EQ.2) THEN
c
c--->       maybe apply a magnetic field
c
            write(6,*) 'atom',ih,'spin',ispin,'shifted by',
     +           -REAL(2*ISPIN-3)*HFIELD*INIPOL(IH)
            DO 130 J = 1,IRCUT(IPAN(IH),IH)
              VISP(J,ISPIN) = VISP(J,ISPIN) - 
     +                        REAL(2*ISPIN-3)*HFIELD*INIPOL(IH)
  130       CONTINUE
          END IF

  140   CONTINUE                    ! ISPIN = 1,NSPIN
C
C write to unformatted file 'vpotnew'
C
        WRITE(66,REC=IH) VINS,VISP,ECORE
C
  150 CONTINUE                      ! IH = NBEG,NEND
C
      DO IH = NBEG,NEND
        IF (ZATINFO(IH).NE.ZAT(IH)) THEN
          WRITE(6,'(79(1H=),/,15X,A,I5,A)') 'trouble on atom ',IH,' ...'
          STOP ' ERROR: inconsistency of Z in atominfo & potential'
        ENDIF
      ENDDO
C
      CLOSE(66)

      RETURN


 9000 FORMAT (16i5)
 9010 FORMAT (4d20.12)
 9020 FORMAT (20a4)
 9030 FORMAT (3f12.8)
 9040 FORMAT (f10.5,/,f10.5,2f15.10)
 9050 FORMAT (i3,/,2d15.8,/,2i2)
 9070 FORMAT (i5,1p,d20.11)
 9080 FORMAT (' <#',20a4)
 9090 FORMAT (10i5)
 9100 FORMAT (1p,4d20.13)
      END                           ! STARTB1
