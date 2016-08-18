C **********************************************************************
      SUBROUTINE STARTB1(IFILE,IPF,IPFE,IPE,KVREL,KWS,KHFELD,LMAX,
     +                  NBEG,NEND,ALAT,RMTNEW,RMT,ITITLE,HFIELD,IMT,IRC,
     +                  VCONST,INS,IRNS,FPRADIUS,LPOT,NSPIN,VINS,IRMIN,
     +                  KSHAPE,NTCELL,IRCUT,IPAN,THETAS,IFUNM,NFU,LLMSP,
     +                  LMSP,EFERMI,VBC,CVLIGHT,DROR,RS,S,VM2Z,RWS,
     +                  ECORE,LCORE,NCORE,DRDI,R,Z,A,B,IRWS,INIPOL,
     +                  IINFO,LMPOTD,IRMIND,IRMD,LMXSPD,IPAND,IRID,
     +                  IRNSD,LMAXD,NATYPD,NCELLD,NFUND,NSPOTD,IVSHIFT)
C **********************************************************************
c   reads the input potentials
c
c    units :       ry - units for energy
c                  the lattice constant and all other lengths
c                                           given in bohr units
c                  the planck constant h/2pi=1
c                  the electron charge e=sqrt(2)
c                  the electron mass m=1/2
c                  the speed of light c = 2/alpha = 274.0720442
c                      with the fine structure constant alpha
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
C **********************************************************************
C     .. Parameters ..
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER LMPOTD,IRMIND,IRMD,LMXSPD
      INTEGER IPAND,IRID,IRNSD,LMAXD,NATYPD,NCELLD,NFUND,NSPOTD
      INTEGER IFILE,IINFO,INS,IPE,IPF,IPFE,KHFELD,KSHAPE,KVREL,KWS,
     +        LMAX,LPOT,NBEG,NEND,NSPIN,IVSHIFT
      DOUBLE PRECISION ALAT,CVLIGHT,EFERMI,HFIELD,VCONST
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),B(*),DRDI(IRMD,*),DROR(IRMD,*),ECORE(20,*),
     +                 R(IRMD,*),RMT(*),RMTNEW(*),RS(IRMD,0:LMAXD,*),
     +                 RWS(*),S(0:LMAXD,*),THETAS(IRID,NFUND,*),VBC(*),
     +                 VINS(IRMIND:IRMD,LMPOTD,*),VM2Z(IRMD,*),Z(*),
     &                 FPRADIUS(NATYPD)
      INTEGER IFUNM(NATYPD,*),IMT(*),INIPOL(*),IPAN(*),
     +        IRC(*),IRCUT(0:IPAND,*),
     +        IRMIN(*),IRNS(*),IRWS(*),ITITLE(20,*),
     +        LCORE(20,*),LLMSP(NATYPD,*),LMSP(NATYPD,*),
     +        NCORE(*),NFU(*),NTCELL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,B1,EA,EFNEW,S1,Z1
      INTEGER I,IA,ICELL,ICORE,IFUN,IH,IMT1,INEW,IO,IPAN1,IR,IRC1,IRI,
     +        IRMINM,IRMINP,IRNS1P,IRT1P,IRWS1,ISAVE,ISPIN,ISUM,
     +        J,L,LM,LM1,LMPOT,LMPOTP,N,NCELL,NFUN,NR
      INTEGER INSLPD
      LOGICAL TEST
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DRN(IRID,NCELLD),SCALE(NCELLD),U(IRMD),
     +                 XRN(IRID,NCELLD)
      DOUBLE PRECISION VSPSME(IRMD) ! dummy for potcut IMPURITY-compatible
      INTEGER MESHN(NCELLD),NM(IPAND,NCELLD),NPAN(NCELLD)
C     ..
C     .. External Subroutines ..
      EXTERNAL CALRMT,POTCUT,RINIT,TEST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC NINT,EXP,LOG,MAX,MOD,DBLE,SQRT
C     ..
C     .. Data statement ..
      INTEGER ISHAPE
      DATA ISHAPE / 0 /
C     ..
c-----------------------------------------------------------------------
c
c ---> output of radial mesh information
c
      IO = 0
      IF (IINFO.NE.0 .AND. TEST('RMESH   ')) IO = 1
c
c---> set speed of light
c
      CVLIGHT = 274.0720442D0
      INSLPD= (IRNSD+1) * LMPOTD * NSPOTD
      CALL RINIT(INSLPD,VINS(IRMIND,1,1))
c-----------------------------------------------------------------------
c
c---> read radial mesh information of the shape functions and
c     shape functions THETAS in the first iteration - if needed
c
      IF ((KSHAPE.NE.0) .AND. (ISHAPE.EQ.0)) THEN
        ISHAPE = 1
        READ (19,FMT=9000) NCELL
        WRITE (1337,FMT=*) '  ncell : ',NCELL,NCELLD
c
        IF(NCELL.GT.NCELLD) THEN
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
          WRITE (1337,FMT=*) '  nfun  : ',NFUN,NFUND
c
          IF(NFUN.GT.NFUND) THEN
            WRITE(6,*) 'Please, change the parameter nfund (',NFUND,
     +           ') in inc.p to',NFUN
            STOP 'STARTB - NFUND'
          ENDIF
c
          DO 10 LM = 1,LMXSPD
            LMSP(ICELL,LM) = 0
   10     CONTINUE

          DO 20 IFUN = 1,NFUN
            READ (19,FMT=9000) LM
            LLMSP(ICELL,IFUN) = LM
            LMSP(ICELL,LM) = 1
            IFUNM(ICELL,LM) = IFUN
            READ (19,FMT=9010) (THETAS(N,IFUN,ICELL),N=1,MESHN(ICELL))
   20     CONTINUE

   30   CONTINUE
      END IF                        ! ((KSHAPE.NE.0) .AND. (IFILE.NE.0))
c-----------------------------------------------------------------------
c
      LMPOT = (LPOT+1)* (LPOT+1)
      DO 150 IH = NBEG,NEND
        DO 140 ISPIN = 1,NSPIN
          I = NSPIN* (IH-1) + ISPIN

          IF (IFILE.NE.0) THEN
            IRCUT(0,IH) = 0
            IF (INS.NE.0) THEN
c p.z.            IF (KSHAPE.NE.0) THEN
              ICELL = NTCELL(IH)
              IPAN(IH) = 1 + NPAN(ICELL)

            ELSE

              IPAN(IH) = 1
            END IF
c
c---> read title of potential card
c
            READ (IFILE,FMT=9020) (ITITLE(IA,I),IA=1,20)
            IF (IINFO.NE.0) THEN 
              IF (INS.EQ.0) THEN
                WRITE (1337,FMT=9080) (ITITLE(IA,I),IA=1,20)
              ELSE
                WRITE (1337,FMT=9081) (ITITLE(IA,I),IA=1,20)
              END IF
            END IF
c
c---  >read muffin-tin radius , lattice constant and new muffin radius
c      (new mt radius is adapted to the given radial mesh)
c
            READ (IFILE,FMT=*) RMT(IH),ALAT,RMTNEW(IH)
            !READ (IFILE,FMT=9030) RMT(IH),ALAT,RMTNEW(IH)
c
c---> read nuclear charge , lmax of the core states ,
c     wigner seitz radius , fermi energy and energy difference
c     between electrostatic zero and muffin tin zero
c
c            READ (IFILE,FMT=9040) Z(IH),RWS(IH),EFNEW,VBC(ISPIN)
            READ (IFILE,*) Z1,RWS(IH),EFNEW,VBC(ISPIN)
            IF (Z(IH).LT.0.D0) Z(IH) = Z1
            IF (Z1.NE.Z(IH).AND.Z(IH).GE.0.D0) THEN
               WRITE(*,*) 'Warning: For atom ',IH,
     &         ': ZATOM different in inputcard and in potential.',
     &          Z(IH),Z1
            ENDIF
c
c---> if efermi .eq. 0 use value from in5
c
            IF (EFNEW.NE.0.0D0 .AND. I.EQ.1) EFERMI = EFNEW
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
            !READ (IFILE,FMT=9050) IRWS(IH),A(IH),B(IH),NCORE(I),INEW
            READ (IFILE,FMT=*) IRWS(IH)
            READ (IFILE,FMT=*) A(IH),B(IH)
            READ (IFILE,FMT=*) NCORE(I),INEW
c
            NR = IRWS(IH)

            IF (NR.GT.IRMD) THEN
              write(6,*) 'Increase parameter IRMD in ''inc.p''',
     +             ' to a value .ge. ',NR,' (= IRWS(',IH,')).'
              STOP 'STARTB1 - IRWS'
            END IF
c
c---> read the different core states : l and energy
c
            IF (NCORE(I).GE.1) THEN
                DO ICORE=1,NCORE(I)
                    READ (IFILE,FMT=9070) LCORE(ICORE,I),
     &                   ECORE(ICORE,I)
                END DO
            END IF
c     
            IF (INS.LT.1) THEN
c
c--->  read radial mesh points, its derivative, the spherically averaged
c      charge density and the input potential without the nuclear pot.
c
              IF (INEW.EQ.0) THEN
                READ (IFILE,FMT=9060) (R(IR,IH),DRDI(IR,IH),VM2Z(IR,I),
     +               IR=1,NR)
              ELSE
                READ (IFILE,FMT=*) (VM2Z(IR,I),IR=1,NR)
              END IF
            ELSE                    ! (INS.LT.1)
c
c--->  read full potential - the non spherical contribution from irmin
c      to irt - remember that the lm = 1 contribution is multiplied by
c      1/sqrt(4 pi)
c
              READ (IFILE,FMT=9090) IRT1P,IRNS1P,LMPOTP,ISAVE
              IRMINP = IRT1P - IRNS1P
              IRMINM = MAX(IRMINP,IRMIND)
              READ (IFILE,FMT=9100) (VM2Z(IR,I),IR=1,NR)
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
                          VINS(IR,LM1,I) = U(IR)
   40                   CONTINUE
                      END IF

                    END IF

                  END IF

   50           CONTINUE

              END IF

            END IF                  ! (INS.LT.1)
c

            IRWS1 = IRWS(IH)
c
c---> redefine new mt-radius in case of shape corrections
c
            IF (INS.NE.0) THEN
c p.z.      IF (KSHAPE.NE.0) THEN
              RMTNEW(IH) = SCALE(ICELL)*ALAT*XRN(1,ICELL)
              IMT1 = NINT(LOG(RMTNEW(IH)/B(IH)+1.0D0)/A(IH)) + 1
c
c---> for proper core treatment imt must be odd
c     shift potential by one mesh point if imt is even
c
              IF (MOD(IMT1,2).EQ.0) THEN
                IMT1 = IMT1 + 1
                DO 60 IR = IMT1,2,-1
                  VM2Z(IR,I) = VM2Z(IR-1,I)
   60           CONTINUE
              END IF
c
              IMT(IH) = IMT1
              B(IH) = RMTNEW(IH)/ (EXP(A(IH)*DBLE(IMT1-1))-1.0D0)
            END IF                  ! (KSHAPE.NE.0)
c
c---> generate radial mesh - potential only is stored in potential card
c     INEW = 1
c     p. zahn, jan. 99
c
            A1 = A(IH)
            B1 = B(IH)
            R(1,IH) = 0.0D0
            DRDI(1,IH) = A1*B1
            DO 70 IR = 2,IRWS1
              EA = EXP(A1*DBLE(IR-1))
              R(IR,IH) = B1* (EA-1.0D0)
              DRDI(IR,IH) = A1*B1*EA
              DROR(IR,IH) = A1/ (1.0D0-1.0D0/EA)
   70       CONTINUE

c
c---> fill cell-type depending mesh points in the non-muffin-tin-region
c
            IF (INS.NE.0) THEN
c p.z.      IF (KSHAPE.NE.0) THEN
              DO 80 IRI = 1,MESHN(ICELL)
                IR = IRI + IMT1
                R(IR,IH) = SCALE(ICELL)*ALAT*XRN(IRI,ICELL)
                DRDI(IR,IH) = SCALE(ICELL)*ALAT*DRN(IRI,ICELL)
                DROR(IR,IH) = DRDI(IR,IH)/R(IR,IH)
   80         CONTINUE
            END IF

            RWS(IH) = R(IRWS1,IH)
c
c---> kshape.eq.0 : calculate new rmt adapted to exp. mesh
c
            CALL CALRMT(IPF,IPFE,IPE,IMT(IH),Z(IH),RMT(IH),RWS(IH),
     +                  RMTNEW(IH),ALAT,DRDI(1,IH),A(IH),B(IH),IRWS1,
     +                  R(1,IH),IO,INS)
c p.z. +                  R(1,IH),IO,KSHAPE)
c
            IF (INS.GT.0) THEN
c p.z.            IF (KSHAPE.GT.0) THEN
              IRCUT(1,IH) = IMT(IH)
              ISUM = IMT(IH)
              DO 90 IPAN1 = 2,IPAN(IH)
                ISUM = ISUM + NM(IPAN1,ICELL)
                IRCUT(IPAN1,IH) = ISUM
   90         CONTINUE
              NR = ISUM
              IF (IRT1P.NE.NR) THEN
                 WRITE(*,*) 'STARTB1: Error: IRT1P.NE.NR',IRT1P,NR,
     &                      ' for atom',IH
                 STOP 'STARTB1: IRT1P.NE.NR'
              ENDIF
            ELSE                    ! (KSHAPE.GT.0)

              NR = IRWS(IH)
              IF (KWS.GE.1) THEN
                IRCUT(1,IH) = IRWS1

              ELSE
                IRCUT(1,IH) = IMT(IH)
              END IF

            END IF                  ! (KSHAPE.GT.0)
c
            IRC(IH) = IRCUT(IPAN(IH),IH)
c
c---> fill array irmin in case of full potential
c
            IF (INS.NE.0)  THEN
               IF (FPRADIUS(IH).GE.0.D0) THEN
                  IRMIN(IH) =  MIN( 
     &                 FLOOR(LOG(FPRADIUS(IH)/B(IH)+1.0D0)/A(IH)) + 1 , 
     &                 IMT(IH)          )
                  IRNS(IH) = NR - IRMIN(IH)
               ELSE IF (IRNS(IH).GE.MESHN(ICELL)) THEN
                  IRMIN(IH) = NR - IRNS(IH)
               ELSE
                  IRNS(IH) = IRNS1P
                  IRMIN(IH) = NR - IRNS(IH)
               ENDIF   
            ENDIF
c
c---> generate arrays for the calculation of the wave functions
c
            Z1 = Z(IH)
            DO 110 L = 0,LMAX
              IF (KVREL.GE.1) THEN
                 S1 = SQRT(DBLE(L*L+L+1)-4.0D0*Z1*Z1/ (CVLIGHT*CVLIGHT))
                 IF (Z1.EQ.0.0D0) S1 = DBLE(L)
              ELSE
                S1 = DBLE(L)
              END IF

              S(L,IH) = S1
              RS(1,L,IH) = 0.0D0
              DO 100 IR = 2,NR
                RS(IR,L,IH) = R(IR,IH)**S1
  100         CONTINUE
  110       CONTINUE                ! L = 0,LMAX
c
c---> cut input potential at rmt if given only at exponential mesh
c
            IF (KSHAPE.EQ.1) THEN
              IMT1 = IMT(IH)
              IRC1 = IRCUT(IPAN(IH),IH)
              CALL POTCUT(IMT1,IRC1,INS,LMPOT,R(1,IH),VM2Z(1,I),VSPSME,
     +                    VINS(IRMIND,1,I),Z(IH),IRMD,IRMIND)
            END IF
c
c--->  first iteration : shift all potentials (only for test purpose)
c      in case of test option 'atptshft' shift only potential of atom at position ivshift
c

            IF (TEST('atptshft').AND.(IH.EQ.IVSHIFT)) THEN
              write(1337,*) 'atptshft',IH,IVSHIFT,VCONST,NR,IRMIN(IH)
              DO 120 J = 1,IRMIN(IH)
                VM2Z(J,I) = VM2Z(J,I) + VCONST
  120         CONTINUE
            ELSEIF(.NOT.TEST('atptshft')) THEN
              DO 121 J = 1,NR
                VM2Z(J,I) = VM2Z(J,I) + VCONST
  121         CONTINUE
            ENDIF

          END IF                    ! (IFILE.NE.0)

c
          IF (KSHAPE.EQ.0 .AND. KWS.EQ.0) THEN
c
c---> in case of a mt calculation cut potential at mt radius
c
            IMT1 = IMT(IH)
            IRWS1 = IRWS(IH)
            CALL POTCUT(IMT1,IRWS1,INS,LMPOT,R(1,IH),VM2Z(1,I),VSPSME,
     +                  VINS(IRMIND,1,I),Z(IH),IRMD,IRMIND)

          END IF                    ! KSHAPE.EQ.0 .AND. KWS.EQ.0
c
!           move this to end of main0 so that non-spherical part is also shifted
!           IF (KHFELD.EQ.1 ) THEN
! c          IF (KHFELD.EQ.1 .AND. NSPIN.EQ.2) THEN
! c
! c--->       maybe apply a magnetic field
! c
!             write(6,*) 'atom',ih,'spin',ispin,'shifted by',
!      +           -DBLE(2*ISPIN-3)*HFIELD*INIPOL(IH)
!             DO 130 J = 1,IRCUT(IPAN(IH),IH)
!               VM2Z(J,I) = VM2Z(J,I) - DBLE(2*ISPIN-3)*HFIELD*INIPOL(IH)
!   130       CONTINUE
!           END IF

  140   CONTINUE                    ! ISPIN = 1,NSPIN

  150 CONTINUE                      ! IH = NBEG,NEND

      IF (INS.NE.0) THEN
         I = 0
         DO IH = NBEG,NEND
            IF (IRMIN(IH).LT.IRMIND) THEN
               WRITE(*,*) 'IRMIN < IRMIND for atom',IH
               WRITE(*,*) IRMIN(IH),IRMIND
               WRITE(*,*) 'Increase dimension IRNSD'
               I = 1
            ENDIF
         ENDDO
         IF (I.NE.0) STOP 'stop startb1 IRNS IRNSD'
      ENDIF


      RETURN


 9000 FORMAT (16i5)
 9010 FORMAT (4d20.12)
 9020 FORMAT (20a4)
 9030 FORMAT (3f12.8)
 9040 FORMAT (f10.5,/,f10.5,2f15.10)
 9050 FORMAT (i3,/,2d15.8,/,2i2)
 9060 FORMAT (1p,2d15.6,1p,d15.8)
 9070 FORMAT (i5,1p,d20.11)
c 9080 FORMAT (10x,20a4)
 9080 FORMAT (' < ',20a4)
 9081 FORMAT (' <#',20a4)
 9090 FORMAT (10i5)
 9100 FORMAT (1p,4d20.13)
      END                           ! STARTB1
