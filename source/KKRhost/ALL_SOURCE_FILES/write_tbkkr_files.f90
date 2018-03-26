!-------------------------------------------------------------------------------
! SUBROUTINE: WRITE_TBKKR_FILES
!> @brief Printing to file the TBKKR files, containing key information of the parameters
!> and structure of the system.
!> @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine WRITE_TBKKR_FILES(LMAX,NEMB,NCLS,NATYP,NAEZ,IELAST,INS,   &
      ALAT,BRAVAIS,RECBV,RBASIS,CLS,NACLS,RCLS,EZOA,ATOM,RR,NSPIN)

   use mod_version_info, only: serialnr

   implicit none
   !interface
   integer, intent(in) :: INS       !< 0 (MT), 1(ASA), 2(Full Potential)
   integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
   integer, intent(in) :: NEMB      !< Number of 'embedding' positions
   integer, intent(in) :: NCLS      !< Number of reference clusters
   integer, intent(in) :: NAEZ      !< Number of atoms in unit cell
   integer, intent(in) :: NATYP     !< Number of kinds of atoms in unit cell
   integer, intent(in) :: NSPIN     !< Counter for spin directions
   integer, intent(in) :: IELAST
   double precision, intent(in) :: ALAT      !< Lattice constant in a.u.
   double precision, dimension(3,3), intent(in) :: RECBV   !< Reciprocal basis vectors
   double precision, dimension(3,3), intent(in) :: BRAVAIS !< Bravais lattice vectors
   double precision, dimension(3,NAEZ+NEMB), intent(in) :: RBASIS   !< Position of atoms in the unit cell in units of bravais vectors
   double precision, dimension(3,NACLSD,NCLSD), intent(in) :: RCLS   !< Real space position of atom in cluster
   double precision, dimension(3,0:NR), intent(in) :: RR    !< Set of real space vectors (in a.u.)
   integer, dimension(NAEZ), intent(in) :: CLS  !< Cluster around atomic sites
   integer, dimension(NCLSD), intent(in) :: NACLS  !< Number of atoms in cluster
   integer, dimension(NACLSD,NAEZ), intent(in) :: EZOA   !< EZ of atom at site in cluster
   integer, dimension(NACLSD,NAEZ), intent(in) :: ATOM   !< Atom at site in cluster
   !.. Local variables
   integer :: I1,I2,J, NACLSMAX
   !     .. External Functions ..
   logical :: OPT
   external :: OPT

   NACLSMAX = 1
   do I1 = 1,NCLS
      if (NACLS(I1).GT.NACLSMAX) NACLSMAX = NACLS(I1)
   enddo

   open(934,FILE='TBkkr_params.txt',FORM='formatted')
   write(934,'(A,A)') '#FILEVERSION= 2'//'   # serial: ',serialnr
   write(934,'(I8,4X,A)') LMAXD, 'lmaxd',   LMAX, 'lmax',         &
                          KORBIT, 'korbit',                       &
                          NSPIN, 'nspin, used as nspind',         &  ! write nspin instead of npsind for program to work with nspin==1 case
                          NR, 'nrd',
                          NEMB, 'nembd',   NEMB, 'nemb',          &
                          NCLS, 'ncls once again, used as nclsd', &
                          NCLS, 'ncls',                           &  ! ncls instead of nclsd for smaller files (see kloopz writeout)
                          NATYP, 'natypd', NATYP, 'natyp',        &
                          NAEZ, 'naezd',   NAEZ, 'naez',          &
                          NACLSMAX, 'naclsmax, used as naclsd',   &  ! naclsmax instead of naclsd for smaller files
                          IELAST, 'ielast',                       &
                          INS, 'ins'
   close(934)

   open(935,FILE='TBkkr_container.txt',FORM='formatted')
   write(935,'(A,A)') '#FILEVERSION= 2'//'   # serial: ',serialnr
   !write out lattice information
   write(935,'(A)') 'alat:'
   write(935,'(ES25.16)') ALAT
   write(935,'(A)') 'bravais:'
   write(935,'(3ES25.16)') ((BRAVAIS(I1,I2),I1=1,3),I2=1,3)
   write(935,'(A)') 'recbv:'
   write(935,'(3ES25.16)') ((RECBV(I1,I2),I1=1,3),I2=1,3)
   write(935,'(A)') 'RBASIS:'
   write(935,'(3ES25.16)') ((RBASIS(J,I1), J=1,3),I1=1,NAEZD+NEMBD)

   !write out cluster information
   write(935,'(A)') 'CLS:'
   write(935,'(1I8)') (CLS(I1),I1=1,NATYPD)
   write(935,'(A)') 'NACLS:'
   write(935,'(1I8)') (NACLS(I1), I1=1,NCLS)
   write(935,'(A)') 'RCLS:'
   do I2=1,NCLS
      do I1=1,NACLSMAX
         write(935,'(3ES25.16)') RCLS(:,I1,I2)
      end do
   end do
   write(935,'(A)') 'EZOA:'
   write(935,'(1I8)') ((EZOA(I1,I2),I1=1,NACLSMAX),I2=1,NAEZD)
   write(935,'(A)') 'ATOM:'
   write(935,'(1I8)') ((ATOM(I1,I2),I1=1,NACLSMAX),I2=1,NAEZD)
   write(935,'(A)') 'RR:'
   do I1=0,NRD
      write(935,'(3ES25.16)') RR(:,I1)
   end do

   close(935)

end subroutine WRITE_TBKKR_FILES
