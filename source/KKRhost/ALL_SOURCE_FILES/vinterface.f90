!*==vinterface.f    processed by SPAG 6.05Rc at 16:57 on  1 Feb 2002
!-------------------------------------------------------------------------------
! SUBROUTINE: VINTERFACE
!> @brief This is calculating the intra-atomic contibution of the potential in
!>  the case of an interface taking into account the bulk potential on
!>  the two sides.
!
!> @details It uses the structure dependent matrices AVMAD which are calculated
!>  once in the subroutine MADELUNG2D() and saved in the DA-file
!>  avmad.unformatted ( May 2004)
!>
!>  For each site in a layer the summation in all other layers is split
!>  into three parts: within the slab, over the NLEFT*NLBASIS left host
!>  sites and over the NRIGHT*NRBASIS right host sites, the last two
!>  steps only in case of decimation run
!
!-------------------------------------------------------------------------------
!> @note
!> - Adapted for the case of more atoms on the same site, summation is
!>  done over the occupants of that site, the charge is weighted with
!>  the appropriate concentration of the occupant  V. Popescu feb. 2002
!-------------------------------------------------------------------------------
!>
!> - Impurity-program adopted feb. 2004 (according to N. Papanikalou)
!>
!> - Jonathan Chico Feb. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine VINTERFACE(CMOM,CMINST,LPOT,NSPIN,NLAYERS,NATYP,V,ZAT,&
      R,IRWS,IRCUT,IPAN,KSHAPE,NOQ,KAOEZ,IQAT,&
      CONC,CATOM,ICC,HOSTIMP,&
      NLBASIS,NLEFT,NRBASIS,NRIGHT,&
      CMOMHOST,CHRGNT,VINTERS)

   use Constants

   implicit none

   ! .. Input variables ..
   integer, intent(in) :: ICC    !< Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
   integer, intent(in) :: IRM    !< Maximum number of radial points
   integer, intent(in) :: LPOT   !< Maximum l component in potential expansion
   integer, intent(in) :: NAEZ   !< Number of atoms in unit cell
   integer, intent(in) :: LMPOT  !< (LPOT+1)**2
   integer, intent(in) :: NSPIN  !< Counter for spin directions
   integer, intent(in) :: NATYP  !< Number of kinds of atoms in unit cell
   integer, intent(in) :: NPOTD  !< (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
   integer, intent(in) :: NLEFT  !< Number of repeated basis for left host to get converged electrostatic potentials
   integer, intent(in) :: NEMBD1 !< NEMB+1
   integer, intent(in) :: NRIGHT !< Number of repeated basis for right host to get converged electrostatic potentials
   integer, intent(in) :: KSHAPE !< Exact treatment of WS cell
   integer, intent(in) :: NLAYERS
   integer, intent(in) :: NLBASIS   !< Number of basis layers of left host (repeated units)
   integer, intent(in) :: NRBASIS   !< Number of basis layers of right host (repeated units)
   double precision, intent(in) :: CHRGNT
   integer, dimension(NAEZ), intent(in)      :: NOQ    !< Number of diff. atom types located
   integer, dimension(NATYP), intent(in)     :: IRWS   !< R point at WS radius
   integer, dimension(NATYP), intent(in)     :: IPAN   !< Number of panels in non-MT-region
   integer, dimension(NATYP), intent(in)     :: IQAT   !< The site on which an atom is located on a given site
   integer, dimension(0:NATYP), intent(in)   :: HOSTIMP
   integer, dimension(0:IPAND,NATYP), intent(in)         :: IRCUT    !< R points of panel borders
   integer, dimension(NATYP,NAEZ+NEMBD1-1), intent(in)   :: KAOEZ    !< Kind of atom at site in elem. cell
   double precision, dimension(NATYP), intent(in) :: ZAT    !< Nuclear charge
   double precision, dimension(NATYP), intent(in) :: CONC   !< Concentration of a given atom
   double precision, dimension(NATYP), intent(in) :: CATOM
   double precision, dimension(IRM,NATYP), intent(in)    :: R        !< Radial mesh ( in units a Bohr)
   double precision, dimension(LMPOT,NATYP), intent(in)  :: CMOM     !< LM moment of total charge
   double precision, dimension(LMPOT,NATYP), intent(in)  :: CMINST   !< charge moment of interstitial
   double precision, dimension(LMPOT,NEMBD1), intent(in) :: CMOMHOST !< Charge moments of each atom of the (left/right) host
   ! .. In/out variables
   double precision, dimension(IRM,LMPOT,NPOTD), intent(inout) :: V

   ! .. Local variables
   integer :: ILEFT,IRIGHT
   integer :: I,IATOM,IB,IH1,ILAY1,ILAY2,IO2
   integer :: IPOT,IRS1,ISPIN,IT1,IT2,L,LM,LM2,LMPOT,M
   integer :: LRECAMAD,IREC,NLEFTOFF,NRIGHTOFF,NLEFTALL,NRIGHTALL
   double precision :: CM1,FPI
   logical :: OPT,TEST,LREAD
   double precision, dimension(LMPOT)  :: AC
   double precision, dimension(LMPOT)  :: CM
   double precision, dimension(2)      :: CHARGE
   double precision, dimension(NAEZ)   :: MONOPOL
   double precision, dimension(LMPOT,LMPOT)  :: AVMAD
   double precision, dimension(LMPOT,NAEZ)   :: VINTERS
   ! .. Intrinsic Functions ..
   intrinsic ATAN,SQRT
   ! .. External Functions/Subroutines
   external OPT,TEST

   if(TEST('flow    ')) write (1337,*) '>>>>>> Vinterface'

   inquire(FILE='avmad.unformatted',EXIST=LREAD) ! ewald2d

   if (LREAD) then
      LRECAMAD = WLENGTH*2*LMPOTD*LMPOTD
      open (69,ACCESS='direct',RECL=LRECAMAD,FILE='avmad.unformatted',FORM='unformatted')
   else
      LRECAMAD = WLENGTH*2*LMPOTD*LMPOTD + WLENGTH*2*LMPOTD
      open (69,ACCESS='direct',RECL=LRECAMAD,FILE='abvmad.unformatted',FORM='unformatted')
   endif

   write(1337,FMT=99001)
   write(1337,FMT=99002)

   FPI = 4.D0*PI
   LMPOT = (LPOT+1)**2

   if ( OPT('DECIMATE') ) then
      !-------------------------------------------------------------------------
      ! Setup the charges to put in the ghost layers in the case of
      ! decimation technique to achieve charge neutrality
      !-------------------------------------------------------------------------
      CHARGE(1) = -CHRGNT/(2.D0*SQRT(FPI))
      CHARGE(2) = -CHRGNT/(2.D0*SQRT(FPI))
      !
      NLEFTOFF = NLAYERS * NLAYERS                     ! record offsets
      NRIGHTOFF = NLEFTOFF + NLAYERS * NLEFT * NLBASIS ! left and right
      NLEFTALL = NLEFT * NLBASIS
      NRIGHTALL = NRIGHT * NRBASIS
   end if
   !----------------------------------------------------------------------------
   !                   START CALCULATION IN THE LAYERS
   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------
   ! Loop over atoms in slab
   !----------------------------------------------------------------------------
   !
   do IT1 = 1,NATYP
      !-------------------------------------------------------------------------
      ! Take a site occupied by IT1
      !-------------------------------------------------------------------------
      ILAY1 = IQAT(IT1)
      !
      if ( KSHAPE.NE.0 ) then
         IRS1 = IRCUT(IPAN(IT1),IT1)
      else
         IRS1 = IRWS(IT1)
      end if
      !
      do LM = 1,LMPOT
         AC(LM) = 0.D0
      end do
      !-------------------------------------------------------------------------
      ! 1.  Summation in all layers in the slab
      !-------------------------------------------------------------------------
      do ILAY2 = 1,NLAYERS
         IREC = ILAY2 + NLAYERS*(ILAY1-1)
         read(69,REC=IREC) AVMAD

         !----------------------------------------------------------------------
         ! Keep the monopole term -- Hoshino is doing (SMONOPOL(I) -SMONOPOL(0))
         !----------------------------------------------------------------------
         if ( ILAY1.EQ.ILAY2 ) MONOPOL(ILAY1) = AVMAD(1,1)
         !----------------------------------------------------------------------
         ! Loop over all occupants of site ILAY2
         !----------------------------------------------------------------------
         do IO2 = 1,NOQ(ILAY2)
            IT2 = KAOEZ(IO2,ILAY2)
            !
            do LM = 1,LMPOT
               CM(LM) = CMOM(LM,IT2)
               !----------------------------------------------------------------
               ! Add contribution of interstial in case of shapes
               !----------------------------------------------------------------
               if ( KSHAPE.NE.0 ) CM(LM) = CM(LM) + CMINST(LM,IT2)
            end do
            CM(1) = CM(1) - ZAT(IT2)/SQRT(FPI)
            !
            do LM = 1,LMPOT
               do LM2 = 1,LMPOT
                  AC(LM) = AC(LM) + AVMAD(LM,LM2)*CM(LM2)*CONC(IT2)
               end do
            end do
         end do
         !----------------------------------------------------------------------
         ! Loop over all occupants of site ILAY2
         !----------------------------------------------------------------------
      end do                 ! ILAY2 loop in all interface planes
      !-------------------------------------------------------------------------
      do ILAY2 = 1,NLAYERS
         !----------------------------------------------------------------------
         ! Loop over all occupants of site ILAY2
         !----------------------------------------------------------------------
         do IO2 = 1,NOQ(ILAY2)
            IT2 = KAOEZ(IO2,ILAY2)
            !
            CM1 = CMOM(1,IT2)
            if ( KSHAPE.NE.0 ) CM1 = CM1 + CMINST(1,IT2)
            !
            CM1 = CM1 - ZAT(IT2)/SQRT(FPI)
            AC(1) = AC(1) - MONOPOL(ILAY1)*CM1*CONC(IT2)
            !
         end do
         !----------------------------------------------------------------------
         ! Loop over all occupants of site ILAY2
         !----------------------------------------------------------------------
      end do
      !-------------------------------------------------------------------------
      ! Correction: charge neutrality is imposed (see P. Lang)
      !-------------------------------------------------------------------------
      if ( OPT('DECIMATE') ) then
         !----------------------------------------------------------------------
         ! 2.  Summation in the LEFT bulk side
         !----------------------------------------------------------------------
         !----------------------------------------------------------------------
         ! Loop over all occupants of LEFT host
         !----------------------------------------------------------------------
         ILEFT = 0
         do IH1 = 1,NLEFT
            do IB = 1,NLBASIS
               ILEFT = ILEFT + 1
               IREC = ILEFT + NLEFTALL*(ILAY1-1) + NLEFTOFF
               read(69,REC=IREC) AVMAD
               !
               IATOM = IB
               do LM = 1,LMPOT
                  do LM2 = 1,LMPOT
                     AC(LM) = AC(LM)+ AVMAD(LM,LM2)*CMOMHOST(LM2,IATOM)
                  end do
               end do
               !
               if ( (IH1.EQ.1) .AND. (IB.EQ.1) ) then
                  AC(1) = AC(1) +(AVMAD(1,1)-MONOPOL(ILAY1))*CHARGE(1)
               endif
            end do
         end do
         !----------------------------------------------------------------------
         if ( ILEFT.NE.NLEFTALL) then
            write(6,*) ' < VINTERFACE > : index error ','ILEFT <> NLEFT*NLBASIS'
            stop
         end if
         !----------------------------------------------------------------------
         ! 3.  Summation in the RIGHT bulk side
         !----------------------------------------------------------------------
         !----------------------------------------------------------------------
         ! Loop over all occupants of RIGHT host
         !----------------------------------------------------------------------
         IRIGHT = 0
         do IH1 = 1,NRIGHT
            do IB = 1,NRBASIS
               IRIGHT = IRIGHT + 1
               IREC = IRIGHT + NRIGHTALL*(ILAY1-1) + NRIGHTOFF
               read(69,REC=IREC) AVMAD
               !
               IATOM = NLBASIS + IB
               do LM = 1,LMPOT
                  do LM2 = 1,LMPOT
                     AC(LM) = AC(LM)+ AVMAD(LM,LM2)*CMOMHOST(LM2,IATOM)
                  end do
               end do
               !
               if ( (IH1.EQ.1) .AND. (IB.EQ.1) ) then
                  AC(1) = AC(1) +(AVMAD(1,1)-MONOPOL(ILAY1))*CHARGE(2)
               endif
            end do
         end do
         !----------------------------------------------------------------------
         if ( IRIGHT.NE.NRIGHTALL ) then
            write(6,*) ' < VINTERFACE > : index error ','IRIGHT <> NRIGHT*NRBASIS'
            stop
         end if
      end if                 ! (OPT(DECIMATE)
      !-------------------------------------------------------------------------
      write (1337,FMT=99003) IT1,(CATOM(IT1)-ZAT(IT1)),(AC(1)/SQRT(4.D0*PI)),(AC(3)/SQRT(4.D0*PI))
      !-------------------------------------------------------------------------
      ! Loop over spins of atom IT1
      !-------------------------------------------------------------------------
      do ISPIN = 1,NSPIN
         !----------------------------------------------------------------------
         ! Determine the right potential number
         !----------------------------------------------------------------------
         IPOT = NSPIN*(IT1-1) + ISPIN
         !----------------------------------------------------------------------
         ! In the case of l=0 : r(1)**l is not defined
         !----------------------------------------------------------------------
         V(1,1,IPOT) = V(1,1,IPOT) + AC(1)
         !
         do L = 0,LPOT
            do M = -L,L
               LM = L*L + L + M + 1
               do I = 2,IRS1
                  V(I,LM,IPOT) = V(I,LM,IPOT) + (-R(I,IT1))**L*AC(LM)
               end do
            end do
         end do
      end do
      !-------------------------------------------------------------------------
      ! This part (ICC.GT.0) should be presumably reconsidered for impurity
      ! calculation in host-CPA case
      !-------------------------------------------------------------------------
      if ( ICC.GT.0  .or. OPT('KKRFLEX ')) then
         do L = 0,LPOT
            so M = -L,L
               LM = L*L + L + M + 1
               VINTERS(LM,ILAY1) = AC(LM)
            end do
         end do
      end if
      !-------------------------------------------------------------------------
   end do
   !----------------------------------------------------------------------------
   close (69)
   write(1337,'(15X,45(1H-),/)')
   write(1337,'(79(1H=))')
   if ( (ICC==0) .and. (.not.OPT('KKRFLEX ')) ) return
   !----------------------------------------------------------------------------
   ! Now Prepare output for Impurity calculation
   !----------------------------------------------------------------------------
   open (91,FILE='intercell_ref',STATUS='unknown',FORM='formatted')
   write(1337,*)
   write(1337,*) '                     ','Writing intercell potential for impurity'
   write(1337,'(/,20X,55(1H-))')
   write(1337,99004) HOSTIMP(0),LMPOT
   write(1337,'(20X,55(1H-),/,35X,"  i host lm  Vint")')
   do I=1,HOSTIMP(0)
      write(1337,*)
      LM = 1
      write(1337,'(35X,I4,I4,I3,1X,F10.6)') I, HOSTIMP(I),LM,VINTERS(LM,HOSTIMP(I))
      do LM=2,9
         write (1337,'(43X,I3,1X,F10.6)') LM,VINTERS(LM,HOSTIMP(I))
      end do
      write(1337,'(20X,55(1H-))')
   end do
   write(1337,'(79(1H=),/)')

   write(91,99005) HOSTIMP(0),LMPOT
   do I=1,HOSTIMP(0)
      write(91,99006) (VINTERS(LM,HOSTIMP(I)),LM=1,LMPOT)
   end do
   close(91)
   !
   return
   !
   99001 format (79(1H=),/,25X,' INTERFACE MADELUNG POTENTIALS ')
   99002 format (/,15X,' ATOM ','  Delta_Q  ','   MONOPOLE       DIPOLE',/,15X,45(1H-))
   99003 format (15X,i4,2X,F10.6,1X,1P,D13.6,1X,1P,D13.6)
   99004 format (22X,I4,' host atoms, LMPOT = ',I2,' output up to LM = 9')
   99005 format (3I6)
   99006 format (4D20.10)
end subroutine VINTERFACE
