!-------------------------------------------------------------------------------
! SUBROUTINE: RHOSYMM
!> @brief Symmetrize the charge densities and magnetic moments of
!> atoms which are magnetic 'antisymmetric'
!> (dependencies in IXIPOL(*))
!
!> @author P. Zahn
!> @date Aug. 1996
!> @note -Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine RHOSYMM(LMPOT,NSPIN,NSTART,NEND,RHO2NS,IXIPOL,IRWS,IRCUT,IPAN,KSHAPE,&
   NATYP,IRM)

   use global_variables

   implicit none
   ! .. Input variables

   integer, intent(in) :: IRM    !< Maximum number of radial points
   integer, intent(in) :: NEND
   integer, intent(in) :: NATYP  !< Number of kinds of atoms in unit cell
   integer, intent(in) :: LMPOT  !< (LPOT+1)**2
   integer, intent(in) :: NSPIN  !< Counter for spin directions
   integer, intent(in) :: KSHAPE !< Exact treatment of WS cell
   integer, intent(in) :: NSTART
   integer, dimension(*), intent(in) :: IPAN !< Number of panels in non-MT-region
   integer, dimension(*), intent(in) :: IRWS !< R point at WS radius
   integer, dimension(*), intent(in) :: IXIPOL  !< Constraint of spin pol.
   integer, dimension(0:IPAND,*), intent(in) :: IRCUT !< R points of panel borders
   ! .. In/Out variables
   double precision, dimension(IRM,LMPOT,NATYP,*), intent(inout) :: RHO2NS !< radial density
   ! .. Local variables
   integer :: I,IATYP,IATYP1,IRC,IRC1,LM
   double precision :: FAC
   ! .. Intrinsic Functions
   intrinsic :: ABS
   !----------------------------------------------------------------------------
   !
   do IATYP = NSTART,NEND
      !
      IATYP1 = ABS(IXIPOL(IATYP))
      !
      FAC = 1.D0
      if (IXIPOL(IATYP).LT.0) FAC = -1.d0
      !
      if (IATYP1.GE.IATYP) then
         write(1337,*) 'Symmetrize atom ',IATYP,' with ',IATYP1,'.'
         if (KSHAPE.NE.0) then
            IRC  = IRCUT(IPAN(IATYP),IATYP)
            IRC1 = IRCUT(IPAN(IATYP1),IATYP1)
         else
            IRC  = IRWS(IATYP)
            IRC1 = IRWS(IATYP1)
         end if
         !
         if (IRC.NE.IRC1) then
            write(6,*) 'Error in RHOSYMM : ***********************'
            write(6,*) 'Radial mesh of atoms ',iatyp,' and ',iatyp1,' are not equal.'
         end if
         !
         do LM = 1,LMPOT
            do I = 1,IRC1
               RHO2NS(I,LM,IATYP,1) = (RHO2NS(I,LM,IATYP,1)+RHO2NS(I,LM,IATYP1,1))/2.d0
               RHO2NS(I,LM,IATYP1,1) = RHO2NS(I,LM,IATYP,1)
               if (NSPIN.GT.1) then
                  RHO2NS(I,LM,IATYP,2) = ( RHO2NS(I,LM,IATYP,2) + &
                     FAC*RHO2NS(I,LM,IATYP1,2) )/2.d0
                  RHO2NS(I,LM,IATYP1,2) = FAC*RHO2NS(I,LM,IATYP,2)
               end if
            enddo ! I =1,IRC1
         enddo ! LM =1,LMPOT
      end if                      ! (IATYP1.GT.IATYP)
   enddo ! IATYP=NSTART,NEND

   return

end subroutine RHOSYMM
