!-------------------------------------------------------------------------------
! SUBROUTINE: RITES
!> @brief this subroutine stores in 'ifile' the necessary results
!> (potentials etc.) to start self-consistency iterations
!
!> @ details Modified for the full potential case - if ins .gt. 0 there
!> is written a different potential card
!> if the sum of absolute values of an lm component of vins (non
!> spher. potential) is less than the given rms error qbound this
!> component will not be stored .
!
!> see to subroutine start , where most of the arrays are described)
!
!> @note modified by B. Drittler  aug. 1988
!> @note Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine RITES(IFILE,NATPS,NATYP,NSPIN,Z,ALAT,RMT,RMTNEW,RWS,&
   ITITLE,R,DRDI,VM2Z,IRWS,A,B,TXC,KXC,INS,IRNS,               &
   LPOT,VINS,QBOUND,IRC,KSHAPE,EFERMI,VBC,ECORE,               &
   LCORE,NCORE,ECOREREL,NKCORE,KAPCORE)

   ! .. Scalar Arguments
   integer, intent(in) :: INS    !< 0 (MT), 1(ASA), 2(Full Potential)
   integer, intent(in) :: KXC    !< Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
   integer, intent(in) :: LPOT   !< Maximum l component in potential expansion
   integer, intent(in) :: IFILE  !< Unit specifier for potential card
   integer, intent(in) :: NATPS
   integer, intent(in) :: NATYP  !< Number of kinds of atoms in unit cell
   integer, intent(in) :: NSPIN  !< Counter for spin directions
   integer, intent(in) :: KSHAPE !< Exact treatment of WS cell
   double precision, intent(in) :: ALAT   !< Lattice constant in a.u.
   double precision, intent(in) :: QBOUND !< Convergence parameter for the potential
   double precision, intent(in) :: EFERMI !< Fermi energy
   ! .. Array Arguments
   double precision, dimension(*), intent(in) :: A       !< Constants for exponential R mesh
   double precision, dimension(*), intent(in) :: B       !< Constants for exponential R mesh
   double precision, dimension(*), intent(in) :: Z
   double precision, dimension(*), intent(in) :: RWS     !< Wigner Seitz radius
   double precision, dimension(2), intent(in) :: VBC     !< Potential constants
   double precision, dimension(*), intent(in) :: RMT     !< Muffin-tin radius of true system
   double precision, dimension(*), intent(in) :: RMTNEW  !< Adapted muffin-tin radius
   double precision, dimension(IRMD,*), intent(in) :: R  !< Radial mesh ( in units a Bohr)
   double precision, dimension(IRMD,*), intent(in) :: VM2Z
   double precision, dimension(IRMD,*), intent(in) :: DRDI     !< Derivative dr/di
   double precision, dimension(20,*), intent(in)   :: ECORE    !< Core energies !(2), 22.5,2000
   double precision, dimension(IRMIND:IRMD,LMPOTD,*), intent(in) :: VINS   !< Non-spherical part of the potential
   !----------------------------------------------------------------------------
   integer, dimension(20,NATYP), intent(in) :: NKCORE
   integer, dimension(20,2*NATYP), intent(in) :: KAPCORE
   double precision, dimension(KREL*20+(1-KREL),2*NATYP), intent(in) :: ECOREREL  ! relativistic core energies
   !----------------------------------------------------------------------------
   integer, dimension(*), intent(in) :: IRC     !< R point for potential cutting
   integer, dimension(*), intent(in) :: IRNS    !< Position of atoms in the unit cell in units of bravais vectors
   integer, dimension(*), intent(in) :: IRWS    !< R point at WS radius
   integer, dimension(*), intent(in) :: NCORE   !< Number of core states
   integer, dimension(20,*), intent(in) :: ITITLE
   integer, dimension(20,*), intent(in) :: LCORE   !< Angular momentum of core states
   character(len=124), dimension(*), intent(in) :: TXC
   ! .. Local Scalars
   integer :: I,ICORE,IH,INEW,IP,IR,IRMIN,IRNS1,IS,ISAVE,J,LM,LMNR,LMPOT,NCORE1,NR
   double precision :: A1,B1,RMAX,RMT1,RMTNW1,RV,SIGN,SUM,Z1
   ! .. Local Arrays
   integer, dimension(20) :: LCORE1
   double precision, dimension(20) :: ECORE1
   double precision, dimension(IRMD) :: DRADI
   double precision, dimension(IRMD) :: RA
   double precision, dimension(IRMD) :: VM2ZA
   double precision, dimension(20,2) :: ECORE2
   character(len=3), dimension(4) :: TXTK
   character(len=1), dimension(0:3) :: TXTL
   !     ..
   DATA TXTL/'s','p','d','f'/
   DATA TXTK/'1/2','3/2','5/2','7/2'/
   !     ..
   ! -------------------------------------------------------------------
   ISAVE = 1
   INEW  = 1
   !
   !
   LMPOT = (LPOT+1)* (LPOT+1)
   do IH = 1,NATYP
      do IS = 1,NSPIN
         if (IS.EQ.NSPIN) then
            SIGN = 1.0D0
         else
            SIGN = -1.0D0
         end if
         IP = NSPIN* (IH-1) + IS

         RMT1 = RMT(IH)
         RMTNW1 = RMTNEW(IH)
         Z1 = Z(IH)
         RMAX = RWS(IH)
         if (KSHAPE.EQ.0) then
            NR = IRWS(IH)
         else
            NR = IRC(IH)
         end if

         IRNS1 = IRNS(IH)
         IRMIN = NR - IRNS1
         A1 = A(IH)
         B1 = B(IH)
         NCORE1 = NCORE(IP)
         !
         do J = 1,NR
            RA(J) = R(J,IH)
            DRADI(J) = DRDI(J,IH)
            !-------------------------------------------------------------------
            ! Store only lm=1 component of the potential
            !-------------------------------------------------------------------
            VM2ZA(J) = VM2Z(J,IP)
         enddo ! J
         !
         open(IFILE, FILE='out_potential', FORM='formatted')
         write (IFILE,FMT=9000) (ITITLE(I,IP),I=1,7),TXC(KXC+1)
         write (IFILE,FMT=9010) RMT1,ALAT,RMTNW1
         write (IFILE,FMT=9020) Z1,RMAX,EFERMI,VBC(IS)
         write (IFILE,FMT=9030) NR,A1,B1,NCORE1,INEW
         !
         if (NCORE1.GE.1) then
            !
            if (KREL.EQ.0) then
               do J = 1,NCORE1
                  LCORE1(J) = LCORE(J,IP)
                  ECORE1(J) = ECORE(J,IP)
               end do
               write (IFILE,FMT=9040) (LCORE1(ICORE),ECORE1(ICORE),ICORE=1,NCORE1)
            else
               do J = 1,NCORE1
                  LCORE1(J) = LCORE(J,IP)
                  ECORE2(J,1) = ECOREREL(J,2*IH-1)
                  ECORE2(J,2) = ECOREREL(J,2*IH)
               end do
               !----------------------------------------------------------------
               ! independent of spin, the \mu-averaged relativistic core energies
               ! are written out for \kappa = -l-1,l
               ! format compatible with the non-(scalar) relativistic mode
               ! however, the next read in has no meaning for the REL core-solver
               ! a detailed output of the core energies is supplied by < CORE >
               !----------------------------------------------------------------
               do ICORE=1,NCORE1
                  write (IFILE,FMT=9041) LCORE1(ICORE),     &
                     (ECORE2(ICORE,I+1),TXTL(LCORE1(ICORE)),&
                     TXTK(IABS(KAPCORE(ICORE,2*IH-1+I))),   &
                     I=0,NKCORE(ICORE,IH)-1)
               end do
            end if
            !
         end if
         !
         if (INS.EQ.0 .OR. (IH.LT.NATPS.AND.INS.LE.2)) then
            !-------------------------------------------------------------------
            ! store only the spherically averaged potential
            ! (in mt or as - case)
            ! this is done always for the host
            !-------------------------------------------------------------------
            if (INEW.EQ.0) then
               write (IFILE,FMT=9050)(RA(IR),DRADI(IR),VM2ZA(IR),IR=1,NR)
            else
               write (IFILE,FMT=9051) (VM2ZA(IR),IR=1,NR)
            end if
         else
            !-------------------------------------------------------------------
            ! store the full potential , but the non spherical contribution
            ! only from irns1 up to irws1 ;
            ! remember that the lm = 1 contribution is multiplied
            ! by a factor 1/sqrt(4 pi)
            !-------------------------------------------------------------------
            write (IFILE,FMT=9060) NR,IRNS1,LMPOT,ISAVE
            write (IFILE,FMT=9070) (VM2ZA(IR),IR=1,NR)
            if (LPOT.GT.0) then
               LMNR = 1
               do LM = 2,LMPOT
                  SUM = 0.0D0
                  do IR = IRMIN,NR
                     RV = VINS(IR,LM,IP)*RA(IR)
                     SUM = SUM + RV*RV*DRADI(IR)
                  enddo ! IR
                  if (SQRT(SUM).GT.QBOUND) then
                     LMNR = LMNR + 1
                     write (IFILE,FMT=9060) LM
                     write (IFILE,FMT=9070) (VINS(IR,LM,IP),IR=IRMIN,NR)
                  end if
               enddo !LM
               !----------------------------------------------------------------
               ! Write a one to mark the end
               !----------------------------------------------------------------
               if (LMNR.LT.LMPOT) write (IFILE,FMT=9060) ISAVE
            end if
         end if
      enddo !IS
   enddo !IH

   close(IFILE)

   9000 format (7a4,6x,'  exc:',a124,3x,a10)
   9010 format (3f12.8)
   9020 format (f10.5,/,f10.5,2f20.15)
   9030 format (i3,/,2d15.8,/,2i2)
   9040 format (i5,1p,d20.11)
   9041 format (i5,2(1p,d20.11,2x,a1,a3))
   9050 format (1p,2d15.6,1p,d15.8)
   9051 format (1p,4d20.12)
   9060 format (10i5)
   9070 format (1p,4d20.13)
end subroutine RITES
