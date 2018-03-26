!-------------------------------------------------------------------------------
! SUBROUTINE: TESTDIM
!> @brief Testing the dimension of several arrays
!> @note Jonathan Chico: Some of these tests seem unnecessary with the changes done to the
!> inc.p
!-------------------------------------------------------------------------------
subroutine TESTDIM(NSPIN,NAEZ,NEMB,NATYP,LMAX,IRM,INS,INSREF,NREF,IRNS,    &
      NCLS,NLAYER,KREL,LMAXD,NSPIND,NCLSD,NPRINCD,KNOSPH,IRM,IRNSD,KORBIT)

   implicit none

   integer, intent(in) :: INS       !< 0 (MT), 1(ASA), 2(Full Potential)
   integer, intent(in) :: IRM       !< Maximum number of radial points
   integer, intent(in) :: NEMB      !< Number of 'embedding' positions
   integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
   integer, intent(in) :: NAEZ      !< Number of atoms in unit cell
   integer, intent(in) :: NREF      !< Number of diff. ref. potentials
   integer, intent(in) :: KREL      !< Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
   integer, intent(in) :: NCLS      !< Number of reference clusters
   integer, intent(in) :: NSPIN     !< Counter for spin directions
   integer, intent(in) :: NATYP     !< Number of kinds of atoms in unit cell
   integer, intent(in) :: NCLSD     !< Maximum number of different TB-clusters
   integer, intent(in) :: IRNSD
   integer, intent(in) :: NLAYER    !< Number of principal layer
   integer, intent(in) :: INSREF    !< INS for reference pot. (usual 0)
   integer, intent(in) :: KNOSPH    !< switch for spherical/non-spherical (0/1) program.
   integer, intent(in) :: KORBIT    !< Spin-orbit/non-spin-orbit (1/0) added to the Schroedinger or SRA equations. Works with FP. KREL and KORBIT cannot be both non-zero.
   integer, intent(in) :: NSPIND    !< KREL+(1-KREL)*(NSPIN+1)
   integer, intent(in) :: NPRINCD   !< Number of principle layers, set to a number >= NRPINC in output of main0
   integer, dimension(NATYP), intent(in) :: IRNS   !< Position of atoms in the unit cell in units of bravais vectors
   !
   integer :: STOP_MARK
   integer :: I,J
   logical :: TEST,OPT
   external :: TEST,OPT
   !
   ! ---> dimension tests
   !
   write(1337,2050)
   !
   stop_mark=0
   if((NSPIN.gt.NSPIND).and.(krel.eq.0)) then
      write(6,*) 'There is an inconsistenciy between spin polarised calculation and relativistic options'
      stop_mark=1
   endif
   if(max(ins,insref).gt.knosph) then
      write(6,*) 'Please, change the parameter insd in',' the inputcard to',max(ins,insref)
      stop_mark=1
   endif
   J=1
   do I=1,NATYP
      J=MAX(J,IRNS(I))
   enddo
   if (INS.eq.0 .and. J.gt.1) then
      write(6,*) 'IRNS(*) is set to 1 in case of ','spherical potential treatment.'
      do I=1,NATYP
         IRNS(I) = 1
      enddo
      J = 1
   end if
   !
   if (J.gt.IRSND) then
      write(6,*) 'Please, change the parameter irnsd in',' the inputcard to',j
      stop_mark=1
   endif
   !
   if ( .not. OPT('VIRATOMS') ) then
      if(nref.gt.natyp) then
         write(6,*) 'There are some inconsistencies in the input file./', &
            ' nref(=',nref,') is greater than natyp (=',natyp,').'
         stop_mark=1
      endif
   end if

   if((KREL.eq.1).and.(KORBIT.eq.1)) then
      write(6,*) 'Full relativistic for ASA and new SO solver',&
         'KREL',KREL,'KORBIT',KORBIT
      stop_mark=1
   endif

   if(.not.OPT('NEWSOSOL').and.KORBIT.EQ.1) then
      write(6,*) 'Option NEWSOSOL not found, change KORBIT in the inputcard from',  &
         KORBIT,'to 0'
      stop_mark=1
   endif

   if(OPT('NEWSOSOL').and.KORBIT.eq.0) then
      write(6,*) 'Using option NEWSOSOL, change KORBIT in the inputcard from',   &
         KORBIT,'to 1'
      stop_mark=1
   endif
   !----------------------------------------------------------------------------
   !  OPT 'WIRE' is only useful with OPT 'full inv' or
   !  OPT 'SPARSE  ' because of sparsity of
   !  the KKR matrix ( not tridiagonal like for 2D and 3D systems)
   !----------------------------------------------------------------------------
   if (OPT('WIRE    ') .and..not. (OPT('full inv') .or. OPT('SPARSE  ') )) then
      write(6,*) 'Use option ''full inv'' or ''SPARSE  '' ','for WIRE calculation.'
      stop_mark=1
   end if
   !
   if (OPT('COMPLEX ') .and..not.( OPT('EigenV  ') .or.OPT('wfct    ') .or.   &
      OPT('iso surf')     )   ) then
      write(6,*) 'Use option ''COMPLEX '' only for eigenvalue determination.'
      stop_mark=1
   end if
   !
   if (TEST('CONT    ')) then
      NEMB = 0
      write(6,*) 'No usage of embedding points. NEMB is set to ',NEMB,'.'
   end if
   !
   if (.not.OPT('full inv').and. .not.OPT('SPARSE  ')) then
      !-------------------------------------------------------------------------
      ! Constants for O(N) algorithm for matrix inversion
      !-------------------------------------------------------------------------
      NLAYER=NAEZ/NPRINCD
      WRITE(1337,2020) NPRINCD,NLAYER
      WRITE(1337,2112)
      ! ignore this test if full inversion is done
      IF ( .not.OPT('full inv') ) THEN
         IF (NLAYER*NPRINCD.NE.NAEZ) THEN
            write(6,*) 'NLAYER*NPRINCD ( = ',NLAYER*NPRINCD,').NE.NAEZ ( = ',NAEZ,')'
            stop_mark=1
         END IF
      END IF

   END IF
   !----------------------------------------------------------------------------
   ! STOP IF A DIMENSION ERROR OCCURED
   !----------------------------------------------------------------------------
   if (stop_mark.gt.0) stop 'STOP : Dimension Error.'
   2020 format(' NPRINCD  NLAYER'/,2I8)
   2112 format( 2(7(1H-),1H+) ,63(1H-))
   2050 format(' Dimension and Input Data CHECK')
   return
end subroutine TESTDIM
