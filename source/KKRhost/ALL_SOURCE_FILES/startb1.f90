!-------------------------------------------------------------------------------
! SUBROUTINE: STARTB1
!> @brief Reads the input potentials
!> @details   units :
!> - Rydbergs - units for energy
!> - The lattice constant and all other lengths given in bohr units
!> - The planck constant \f$\frac{h}{2\pi}=1\f$
!> - The electron charge \f$e=\sqrt{2}\f$
!> - The electron mass \f$m=\frac{1}{2}\f$
!> - The speed of light \f$c = \frac{2}{\alpha} = 274.0720442\f$ with the fine structure constant \f$\alpha\f$
!>
!> In case of shape corrections this routine reads from unit 19 a suitable radial mesh 'xrn',its derivate 'drn' and the shape
!> functions 'thetas'. Thus, the region from the muffin-tin to the circumscribed sphere radii is divided  into 'npan'
!> pannels, each one containing 'nm(ipan)' points in order to take care of the  discontinuities of the shape-function  derivative.
!
!> @note remember that the input potentials do not include the electro-static contribution of the nucleus of the cell itself
!> this has to be added explicitly!
!> @note Modified for bandstructure code
!> - Jonathan Chico: Removed inc.p dependencies and rewrote to Fortran90
!> @author B. Drittler
!> @date Nov. 1989
!-------------------------------------------------------------------------------
subroutine STARTB1(IFILE,IPF,IPFE,IPE,KREL,KWS,LMAX,  &
   NBEG,NEND,ALAT,RMTNEW,RMT,ITITLE,IMT,IRC,          &
   VCONST,INS,IRNS,FPRADIUS,LPOT,NSPIN,VINS,IRMIN,    &
   KSHAPE,NTCELL,IRCUT,IPAN,THETAS,IFUNM,NFU,LLMSP,   &
   LMSP,EFERMI,VBC,DROR,RS,S,VM2Z,RWS,                &
   ECORE,LCORE,NCORE,DRDI,R,ZAT,A,B,IRWS,             &
   IINFO,LMPOT,IRMIND,IRM,LMXSPD,IPAND,IRID,          &
   IRNSD,NATYP,NCELLD,NFUND,NSPOTD,IVSHIFT,NPOTD)

   use Constants

   implicit none
   ! ..
   ! .. Input variables
   integer, intent(in) :: IRM    !< Maximum number of radial points
   integer, intent(in) :: KWS    !< 0 (MT), 1(ASA)
   integer, intent(in) :: INS    !< 0 (MT), 1(ASA), 2(Full Potential)
   integer, intent(in) :: IPE    !< Not real used, IPFE should be 0
   integer, intent(in) :: IPF    !< Not real used, IPFE should be 0
   integer, intent(in) :: IPFE   !< Not real used, IPFE should be 0
   integer, intent(in) :: IRID   !< Shape functions parameters in non-spherical part
   integer, intent(in) :: LMAX   !< Maximum l component in wave function expansion
   integer, intent(in) :: LPOT   !< Maximum l component in potential expansion
   integer, intent(in) :: NBEG   !< Starting number for reading the potential
   integer, intent(in) :: NEND   !< Final number for reading the potential
   integer, intent(in) :: KREL   !< Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
   integer, intent(in) :: NPOTD     !< (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
   integer, intent(in) :: NSPIN  !< Counter for spin directions
   integer, intent(in) :: IPAND  !< Number of panels in non-spherical part
   integer, intent(in) :: IRNSD
   integer, intent(in) :: NFUND  !< Shape functions parameters in non-spherical part
   integer, intent(in) :: IFILE  !< Unit specifier for potential card
   integer, intent(in) :: IINFO
   integer, intent(in) :: NATYP  !< Number of kinds of atoms in unit cell
   integer, intent(in) :: LMPOT  !< (LPOT+1)**2
   integer, intent(in) :: IRMIND !< IRM-IRNSD
   integer, intent(in) :: LMXSPD !< (2*LPOT+1)**2
   integer, intent(in) :: NCELLD !< Number of cells (shapes) in non-spherical part
   integer, intent(in) :: NSPOTD !< Number of potentials for storing non-sph. potentials
   integer, intent(in) :: KSHAPE !< Exact treatment of WS cell
   integer, intent(in) :: IVSHIFT
   double precision, intent(in) :: VCONST    !< Potential shift
   integer, dimension(NATYP), intent(in) :: NTCELL !< Index for WS cell
   double precision, dimension(NATYP), intent(in) :: FPRADIUS  !< R point at which full-potential treatment starts
   ! .. In/Out variables
   double precision, intent(inout) :: ALAT      !< Lattice constant in a.u.
   double precision, intent(inout) :: EFERMI    !< Fermi energy
   integer, dimension(NATYP), intent(inout) :: NFU !< number of shape function components in cell 'icell'
   integer, dimension(NATYP), intent(inout) :: IMT !< R point at MT radius
   integer, dimension(NATYP), intent(inout) :: IRC !< R point for potential cutting
   integer, dimension(NATYP), intent(inout) :: IPAN !< Number of panels in non-MT-region
   integer, dimension(NATYP), intent(inout) :: IRNS !< Position of atoms in the unit cell in units of bravais vectors
   integer, dimension(NATYP), intent(inout) :: IRWS   !< R point at WS radius
   integer, dimension(NATYP), intent(inout) :: IRMIN  !< Max R for spherical treatment
   integer, dimension(NPOTD), intent(inout) :: NCORE  !< Number of core states
   integer, dimension(NATYP,LMXSPD), intent(inout)    :: LMSP   !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
   integer, dimension(20,NPOTD), intent(inout)        :: LCORE  !< Angular momentum of core states
   integer, dimension(NATYP,NFUND), intent(inout)     :: LLMSP  !< lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
   integer, dimension(0:IPAND,NATYP), intent(inout)   :: IRCUT !< R points of panel borders
   integer, dimension(NATYP,LMXSPD), intent(inout)    :: IFUNM
   integer, dimension(20,NPOTD), intent(inout)        :: ITITLE !< Titles of the potential card
   double precision, dimension(NATYP), intent(inout)  :: A      !< Constants for exponential R mesh
   double precision, dimension(NATYP), intent(inout)  :: B      !< Constants for exponential R mesh
   double precision, dimension(NATYP), intent(inout)  :: ZAT    !< Nuclear charge
   double precision, dimension(2), intent(inout)      :: VBC    !< Potential constants
   double precision, dimension(NATYP), intent(inout)  :: RMT    !< Muffin-tin radius of true system
   double precision, dimension(NATYP), intent(inout)  :: RWS    !< Wigner Seitz radius
   double precision, dimension(NATYP), intent(inout)  :: RMTNEW !< Adapted muffin-tin radius
   double precision, dimension(0:LMAX,NATYP), intent(inout) :: S
   double precision, dimension(IRM,NATYP), intent(inout)    :: R     !< Radial mesh ( in units a Bohr)
   double precision, dimension(IRM,NATYP), intent(inout)    :: DRDI  !< Derivative dr/di
   double precision, dimension(IRM,NATYP), intent(inout)    :: DROR
   double precision, dimension(IRM,NPOTD), intent(inout)    :: VM2Z
   double precision, dimension(20,NPOTD), intent(inout)     :: ECORE  !< Core energies
   double precision, dimension(IRM,0:LMAX,NATYP), intent(inout)         :: RS
   double precision, dimension(IRMIND:IRM,LMPOT,NSPOTD), intent(inout)  :: VINS   !< Non-spherical part of the potential
   double precision, dimension(IRID,NFUND,NCELLD), intent(inout)        :: THETAS !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
   ! .. Local Scalars
   integer :: INSLPD,LMSHAPEMAX
   integer :: J,L,LM,LM1,LMPOTP,N,NCELL,NFUN,NR
   integer :: IRMINM,IRMINP,IRNS1P,IRT1P,IRWS1,ISAVE,ISPIN,ISUM
   integer :: I,IA,ICELL,ICORE,IFUN,IH,IMT1,INEW,IO,IPAN1,IR,IRC1,IRI
   double precision :: A1,B1,EA,EFNEW,S1,Z1
   logical :: TEST
   ! ..
   ! .. Local Arrays
   integer, dimension(NCELLD) :: NPAN
   integer, dimension(NCELLD) :: MESHN
   integer, dimension(IPAND,NCELLD) :: NM
   double precision, dimension(IRM)    :: U
   double precision, dimension(NCELLD) :: SCALE
   double precision, dimension(IRID)   :: RDUMMY
   double precision, dimension(IRM)    :: VSPSME ! dummy for potcut IMPURITY-compatible
   double precision, dimension(IRID,NCELLD) :: DRN
   double precision, dimension(IRID,NCELLD) :: XRN
   ! ..
   ! .. External Subroutines ..
   external :: CALRMT,POTCUT,RINIT,TEST
   ! ..
   ! .. Intrinsic Functions ..
   intrinsic :: NINT,EXP,LOG,MAX,MOD,DBLE,SQRT
   ! ..
   ! .. Data statement ..
   integer :: ISHAPE
   data ISHAPE / 0 /
   !----------------------------------------------------------------------------
   ! Output of radial mesh information
   !----------------------------------------------------------------------------
   IO = 0
   if (IINFO.NE.0 .AND. TEST('RMESH   ')) IO = 1
   !----------------------------------------------------------------------------
   ! Set speed of light
   !----------------------------------------------------------------------------
   INSLPD= (IRNSD+1) * LMPOT * NSPOTD
   LMSHAPEMAX = (4*LMAX+1)**2
   call RINIT(INSLPD,VINS(IRMIND,1,1))
   !----------------------------------------------------------------------------
   ! Read radial mesh information of the shape functions and
   ! shape functions THETAS in the first iteration - if needed
   !----------------------------------------------------------------------------
   if ((KSHAPE.NE.0) .AND. (ISHAPE.EQ.0)) then
      ISHAPE = 1
      read (19,FMT=9000) NCELL
      write (1337,FMT=*) '  ncell : ',NCELL,NCELLD
      !       check consistency with shape numbers from inputcard
      if(maxval(NTCELL(1:NATYP))>NCELL) then
         write(*,*) 'Found ',NCELL,'shapes in shapefun file but need',  &
         maxval(NTCELL(1:NATYP)),'according to inputcard/default values'
         write(*,*) 'Did you set <SHAPE> correctly in inputcard?'
         stop 'Error consistency shapes from input/shapefun file'
      endif
      !
      if(NCELL.GT.NCELLD) then
         write(6,*) 'Please, change the parameter ncelld (',NCELLD,') in the inputcard to',NCELL
         stop 'STARTB - NCELLD'
      endif
      !
      read (19,FMT=9010) (SCALE(ICELL),ICELL=1,NCELL)
      do 30 ICELL = 1,NCELL
         read (19,FMT=9000) NPAN(ICELL),MESHN(ICELL)
         !
         if(NPAN(ICELL)+1.gt.IPAND) then
            write(6,*) 'Please, change the parameter ipand (',IPAND,') in the inputcard to',NPAN(ICELL)+1
            stop 'STARTB - IPAND'
         endif
         !
         if(MESHN(ICELL).GT.IRID) then
            write(6,*) 'Please, change the parameter irid (',IRID,') in the inputcard to',MESHN(ICELL)
            stop 'STARTB - IRID'
         endif
         !
         read (19,FMT=9000) (NM(IPAN1,ICELL),IPAN1=2,NPAN(ICELL)+1)
         read (19,FMT=9010) (XRN(IR,ICELL),DRN(IR,ICELL),IR=1,MESHN(ICELL))

         read (19,FMT=9000) NFU(ICELL)
         NFUN = NFU(ICELL)
         write (1337,FMT=*) '  nfun  : ',NFUN,NFUND
         !
         if(NFUN.gt.NFUND) then
            write(6,*) 'Please, change the parameter nfund (',NFUND,') in the inputcard to',NFUN
            stop 'STARTB - NFUND'
         endif
         !
         do 10 LM = 1,LMXSPD
            LMSP(ICELL,LM) = 0
         10 continue

         do 20 IFUN = 1,NFUN
            read (19,FMT=9000) LM
            if(LM<=LMSHAPEMAX)then
               LLMSP(ICELL,IFUN) = LM
               LMSP(ICELL,LM) = 1
               IFUNM(ICELL,LM) = IFUN
               read (19,FMT=9010) (THETAS(N,IFUN,ICELL),N=1,MESHN(ICELL))
            else
               read (19,FMT=9010) (RDUMMY(N),N=1,MESHN(ICELL))
            endif
         20 continue

      30 continue
   end if                        ! ((KSHAPE.NE.0) .AND. (IFILE.NE.0))
   !----------------------------------------------------------------------------
   !LMPOT = (LPOT+1)* (LPOT+1)
   do 150 IH = NBEG,NEND
      do 140 ISPIN = 1,NSPIN
         I = NSPIN* (IH-1) + ISPIN

         if (IFILE.NE.0) then
            IRCUT(0,IH) = 0
            if (INS.NE.0) then
               ! p.z.            IF (KSHAPE.NE.0) THEN
               ICELL = NTCELL(IH)
               IPAN(IH) = 1 + NPAN(ICELL)
            else
               IPAN(IH) = 1
            end if
            !-------------------------------------------------------------------
            ! Read title of potential card
            !-------------------------------------------------------------------
            read (IFILE,FMT=9020) (ITITLE(IA,I),IA=1,20)
            if (IINFO.NE.0) then
               if (INS.EQ.0) then
                  write (1337,FMT=9080) (ITITLE(IA,I),IA=1,20)
               else
                  write (1337,FMT=9081) (ITITLE(IA,I),IA=1,20)
               end if
            end if
            !
            !-------------------------------------------------------------------
            ! Read muffin-tin radius , lattice constant and new muffin radius
            ! (new mt radius is adapted to the given radial mesh)
            !-------------------------------------------------------------------
            read (IFILE,FMT=*) RMT(IH),ALAT,RMTNEW(IH)
            !READ (IFILE,FMT=9030) RMT(IH),ALAT,RMTNEW(IH)
            !
            !-------------------------------------------------------------------
            ! Read nuclear charge , lmax of the core states ,
            ! wigner seitz radius , fermi energy and energy difference
            ! between electrostatic zero and muffin tin zero
            !-------------------------------------------------------------------
            !            READ (IFILE,FMT=9040) ZAT(IH),RWS(IH),EFNEW,VBC(ISPIN)
            read (IFILE,*) Z1
            read (IFILE,*) RWS(IH),EFNEW,VBC(ISPIN)

            !             READ (IFILE,*) Z1,RWS(IH),EFNEW,VBC(ISPIN)
            if (ZAT(IH).LT.0.D0) ZAT(IH) = Z1
            if (Z1.NE.ZAT(IH).AND.ZAT(IH).GE.0.D0) then
               write(*,*) 'Warning: For atom ',IH,                   &
               ': ZATOM different in inputcard and in potential.',   &
               ZAT(IH),Z1
            endif
            !-------------------------------------------------------------------
            ! If efermi .eq. 0 use value from in5
            !-------------------------------------------------------------------
            if (EFNEW.NE.0.0D0 .AND. I.EQ.1) EFERMI = EFNEW
            !-------------------------------------------------------------------
            ! Read : number of radial mesh points
            ! (in case of ws input-potential: last mesh point corresponds
            ! to ws-radius, in case of shape-corrected input-potential
            ! last mesh point of the exponential mesh corresponds to
            ! mt-radius/nevertheless this point is always in the array
            ! irws(ih)),number of points for the radial non-muffin-tin
            ! mesh  needed for shape functions, the constants a and b
            ! for the radial exponential mesh : r(i) = b*(exp(a*(i-1))-1)
            ! the no. of different core states and some other stuff
            !-------------------------------------------------------------------
            !READ (IFILE,FMT=9050) IRWS(IH),A(IH),B(IH),NCORE(I),INEW
            read (IFILE,FMT=*) IRWS(IH)
            read (IFILE,FMT=*) A(IH),B(IH)
            read (IFILE,FMT=*) NCORE(I),INEW
            !
            NR = IRWS(IH)

            if (NR.GT.IRM) then
               write(6,*) 'Increase parameter IRM in the inputcard ', &
               ' to a value .ge. ',NR,' (= IRWS(',IH,')).'
               stop 'STARTB1 - IRWS'
            end if
            !-------------------------------------------------------------------
            ! Read the different core states : l and energy
            !-------------------------------------------------------------------
            if (NCORE(I).GE.1) then
               do ICORE=1,NCORE(I)
                  read (IFILE,FMT=9070) LCORE(ICORE,I),ECORE(ICORE,I)
               end do
            end if
            !
            if (INS.lt.1) then
               !----------------------------------------------------------------
               ! Read radial mesh points, its derivative, the spherically averaged
               ! charge density and the input potential without the nuclear pot.
               !----------------------------------------------------------------
               if (INEW.eq.0) then
                  read (IFILE,FMT=9060) (R(IR,IH),DRDI(IR,IH),VM2Z(IR,I),IR=1,NR)
               else
                  read (IFILE,FMT=*) (VM2Z(IR,I),IR=1,NR)
               end if
            else                    ! (INS.LT.1)
               !-------------------------------------------------------------------
               ! Read full potential - the non spherical contribution from irmin
               ! to irt - remember that the lm = 1 contribution is multiplied by
               ! 1/sqrt(4 pi)
               !-------------------------------------------------------------------
               read (IFILE,FMT=9090) IRT1P,IRNS1P,LMPOTP,ISAVE
               IRMINP = IRT1P - IRNS1P
               IRMINM = MAX(IRMINP,IRMIND)
               read (IFILE,FMT=9100) (VM2Z(IR,I),IR=1,NR)
               if (LMPOTP.GT.1) then
                  LM1 = 2
                  do 50 LM = 2,LMPOTP
                     if (LM1.NE.1) then
                        if (ISAVE.EQ.1) then
                           read (IFILE,FMT=9090) LM1
                        else
                           LM1 = LM
                        end if

                        if (LM1.GT.1) then
                           read (IFILE,FMT=9100) (U(IR),IR=IRMINP,NR)
                           if (LM1.LE.LMPOT) then
                              do 40 IR = IRMINM,NR
                                 VINS(IR,LM1,I) = U(IR)
                              40 continue
                           end if
                        end if
                     end if
                  50 continue
               end if
            end if                  ! (INS.LT.1)
            !
            IRWS1 = IRWS(IH)
            !----------------------------------------------------------------------
            ! Redefine new mt-radius in case of shape corrections
            !----------------------------------------------------------------------
            if (INS.NE.0) then
               ! p.z.      IF (KSHAPE.NE.0) THEN
               RMTNEW(IH) = SCALE(ICELL)*ALAT*XRN(1,ICELL)
               IMT1 = NINT(LOG(RMTNEW(IH)/B(IH)+1.0D0)/A(IH)) + 1
               !-------------------------------------------------------------------
               ! For proper core treatment imt must be odd
               ! shift potential by one mesh point if imt is even
               !-------------------------------------------------------------------
               if (MOD(IMT1,2).EQ.0) then
                  IMT1 = IMT1 + 1
                  do 60 IR = IMT1,2,-1
                     VM2Z(IR,I) = VM2Z(IR-1,I)
                  60 continue
               end if
               !
               IMT(IH) = IMT1
               B(IH) = RMTNEW(IH)/ (EXP(A(IH)*DBLE(IMT1-1))-1.0D0)
            end if                  ! (KSHAPE.NE.0)
            !----------------------------------------------------------------------
            ! Generate radial mesh - potential only is stored in potential card
            ! INEW = 1
            ! p. zahn, jan. 99
            !----------------------------------------------------------------------
            A1 = A(IH)
            B1 = B(IH)
            R(1,IH) = 0.0D0
            DRDI(1,IH) = A1*B1
            do 70 IR = 2,IRWS1
               EA = EXP(A1*DBLE(IR-1))
               R(IR,IH) = B1* (EA-1.0D0)
               DRDI(IR,IH) = A1*B1*EA
               DROR(IR,IH) = A1/ (1.0D0-1.0D0/EA)
            70 continue
            !----------------------------------------------------------------------
            ! Fill cell-type depending mesh points in the non-muffin-tin-region
            !----------------------------------------------------------------------
            if (INS.NE.0) then
               ! p.z.      IF (KSHAPE.NE.0) THEN
               do 80 IRI = 1,MESHN(ICELL)
                  IR = IRI + IMT1
                  R(IR,IH) = SCALE(ICELL)*ALAT*XRN(IRI,ICELL)
                  DRDI(IR,IH) = SCALE(ICELL)*ALAT*DRN(IRI,ICELL)
                  DROR(IR,IH) = DRDI(IR,IH)/R(IR,IH)
               80 continue
            end if

            RWS(IH) = R(IRWS1,IH)
            !----------------------------------------------------------------------
            ! Kshape.eq.0 : calculate new rmt adapted to exp. mesh
            !----------------------------------------------------------------------
            call CALRMT(IPF,IPFE,IPE,IMT(IH),ZAT(IH),RMT(IH),RWS(IH),   &
               RMTNEW(IH),ALAT,DRDI(1,IH),A(IH),B(IH),IRWS1,            &
               R(1,IH),IO,INS)
            ! p.z. +                  R(1,IH),IO,KSHAPE)
            !
            if (INS.GT.0) then
               ! p.z.            IF (KSHAPE.GT.0) THEN
               IRCUT(1,IH) = IMT(IH)
               ISUM = IMT(IH)
               do 90 IPAN1 = 2,IPAN(IH)
                  ISUM = ISUM + NM(IPAN1,ICELL)
                  IRCUT(IPAN1,IH) = ISUM
               90 continue
               NR = ISUM
               if (IRT1P.NE.NR) then
                  WRITE(*,*) 'STARTB1: Error: IRT1P.NE.NR',IRT1P,NR,' for atom',IH
                  stop 'STARTB1: IRT1P.NE.NR'
               endif
            else                    ! (KSHAPE.GT.0)
               NR = IRWS(IH)
               if (KWS.GE.1) then
                  IRCUT(1,IH) = IRWS1
               else
                  IRCUT(1,IH) = IMT(IH)
               end if
            end if                  ! (KSHAPE.GT.0)
            !
            IRC(IH) = IRCUT(IPAN(IH),IH)
            !----------------------------------------------------------------------
            ! Fill array irmin in case of full potential
            !----------------------------------------------------------------------
            if (INS.NE.0)  then
               if (FPRADIUS(IH).GE.0.D0) then
                  IRMIN(IH) = MIN(FLOOR(LOG(FPRADIUS(IH)/B(IH)+1.0D0)/A(IH))+1,IMT(IH))
                  IRNS(IH) = NR - IRMIN(IH)
               else if (IRNS(IH).GE.MESHN(ICELL)) then
                  IRMIN(IH) = NR - IRNS(IH)
               else
                  IRNS(IH) = IRNS1P
                  IRMIN(IH) = NR - IRNS(IH)
               endif
            endif
            !----------------------------------------------------------------------
            ! Generate arrays for the calculation of the wave functions
            !----------------------------------------------------------------------
            Z1 = ZAT(IH)
            do 110 L = 0,LMAX
               if (KREL.GE.1) then
                  S1 = SQRT(DBLE(L*L+L+1)-4.0D0*Z1*Z1/ (CVLIGHT*CVLIGHT))
                  if (Z1.EQ.0.0D0) S1 = DBLE(L)
               else
                  S1 = DBLE(L)
               end if

               S(L,IH) = S1
               RS(1,L,IH) = 0.0D0
               do 100 IR = 2,NR
                  RS(IR,L,IH) = R(IR,IH)**S1
               100 continue
            110 continue                ! L = 0,LMAX
            !-------------------------------------------------------------------
            ! Cut input potential at rmt if given only at exponential mesh
            !-------------------------------------------------------------------
            if (KSHAPE.EQ.1) then
               IMT1 = IMT(IH)
               IRC1 = IRCUT(IPAN(IH),IH)
               call POTCUT(IMT1,IRC1,INS,LMPOT,R(1,IH),VM2Z(1,I),VSPSME,   &
                  VINS(IRMIND,1,I),ZAT(IH),IRM,IRMIND)
            end if
            !-------------------------------------------------------------------
            ! First iteration : shift all potentials (only for test purpose)
            ! in case of test option 'atptshft' shift only potential of atom at position ivshift
            !-------------------------------------------------------------------
            if (TEST('atptshft').AND.(IH.EQ.IVSHIFT)) then
               write(1337,*) 'atptshft',IH,IVSHIFT,VCONST,NR,IRMIN(IH)
               do 120 J = 1,IRMIN(IH)
                  VM2Z(J,I) = VM2Z(J,I) + VCONST
               120 continue
            elseif(.NOT.TEST('atptshft')) then
               do 121 J = 1,NR
                  VM2Z(J,I) = VM2Z(J,I) + VCONST
               121 continue
            endif
         end if                    ! (ifile.ne.0)
         !
         if (KSHAPE.EQ.0 .AND. KWS.EQ.0) then
            !-------------------------------------------------------------------
            ! In case of a mt calculation cut potential at mt radius
            !-------------------------------------------------------------------
            IMT1 = IMT(IH)
            IRWS1 = IRWS(IH)
            call POTCUT(IMT1,IRWS1,INS,LMPOT,R(1,IH),VM2Z(1,I),VSPSME,  &
               VINS(IRMIND,1,I),ZAT(IH),IRM,IRMIND)
         end if                    ! KSHAPE.EQ.0 .AND. KWS.EQ.0
      140 continue                    ! ISPIN = 1,NSPIN
   150 continue                      ! IH = NBEG,NEND

   if (INS.NE.0) then
      I = 0
      do IH = NBEG,NEND
        if (IRMIN(IH).LT.IRMIND) then
            write(*,*) 'IRMIN < IRMIND for atom',IH
            write(*,*) IRMIN(IH),IRMIND
            write(*,*) 'Increase dimension IRNSD'
            I = 1
         endif
      enddo
      if (I.NE.0) stop 'stop startb1 IRNS IRNSD'
   endif

   return

   9000 format (16i5)
   9010 format (4d20.12)
   9020 format (20a4)
   9030 format (3f12.8)
   9040 format (f10.5,/,f10.5,2f15.10)
   9050 format (i3,/,2d15.8,/,2i2)
   9060 format (1p,2d15.6,1p,d15.8)
   9070 format (i5,1p,d20.11)
   ! 9080 format (10x,20a4)
   9080 format (' < ',20a4)
   9081 format (' <#',20a4)
   9090 format (10i5)
   9100 format (1p,4d20.13)
end subroutine STARTB1
