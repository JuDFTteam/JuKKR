!-------------------------------------------------------------------------------
! SUBROUTINE: writekkrflex
!> @brief Subroutine dealing with the printing of the needed kkrflex files for the
!> realization of an impurity calculation with the KKRFLEX software package.
!> @details It specifically prints the following files:
!> - kkrflex_tmat
!> - kkrflex_intercell_ref
!> - kkrflex_intercell_cmoms
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine writekkrflex(NATOMIMP,NSPIN,IELAST,LMPOT,LMMAXD,ALAT,NATYP,&
      KSHAPE,VBC,ATOMIMP,HOSTIMP,NOQ,ZAT,KAOEZ,CONC,CMOM,CMINST,VINTERS)

   use mod_types, only: t_tgmat
   use mod_wunfiles, only: t_params, read_angles
   use mod_version_info
   use mod_md5sums

   implicit none

   ! .. Input variables
   integer, intent(in) :: NEMB      !< Number of 'embedding' positions
   integer, intent(in) :: NAEZ      !< Number of atoms in unit cell
   integer, intent(in) :: LMPOT     !< (LPOT+1)**2
   integer, intent(in) :: NSPIN     !< Counter for spin directions
   integer, intent(in) :: NATYP     !< Number of kinds of atoms in unit cell
   integer, intent(in) :: LMPOT     !< (LPOT+1)**2
   integer, intent(in) :: KSHAPE    !< Exact treatment of WS cell
   integer, intent(in) :: IELAST
   integer, intent(in) :: LMMAXD    !< (KREL+KORBIT+1)(LMAX+1)^2
   integer, intent(in) :: LMGF0D    !< (LMAX+1)**2
   integer, intent(in) :: NATOMIMP  !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
   double precision, intent(in) :: ALAT   !< Lattice constant in a.u.
   integer, dimension(NAEZ), intent(in)            :: NOQ  !< Number of diff. atom types located
   integer, dimension(NATOMIMP), intent(in)        :: ATOMIMP
   integer, dimension(0:NATYP), intent(in)         :: HOSTIMP
   integer, dimension(NATYP,NAEZ+NEMB), intent(in) :: KAOEZ !< Kind of atom at site in elem. cell
   double precision, dimension(2), intent(in)      :: VBC  !< Potential constants
   double precision, dimension(NATYP), intent(in)  :: ZAT    !< Nuclear charge
   double precision, dimension(NATYP), intent(in)  :: CONC   !< Concentration of a given atom
   double precision, dimension(LMPOT,NATYP), intent(in)  :: CMOM
   double precision, dimension(LMPOT,NATYP), intent(in)  :: CMINST
   double precision, dimension(LMPOT,NAEZ), intent(in)   :: VINTERS
   ! .. Local variables
   integer :: ispin,ie,i1,iatom,irec,i,lm
   double precision, dimension(NATYP) :: THETA, PHI
   double complex, dimension(LMMAXD,LMMAXD) :: TMAT0
   ! .. External Functions
   logical :: OPT
   external OPT

   write(1337,*) 'KKRFLEX WRITEOUT'
   write(1337,*) OPT('KKRFLEX ')

   if ( OPT('KKRFLEX ') ) then
      open (6699,FILE='kkrflex_tmat',STATUS='unknown')
      call version_print_header(6699,'; '//md5sum_potential//'; '//md5sum_shapefun)
      write(6699,*) '#',NATOMIMP,NSPIN,IELAST,LMMAXD,KORBIT
      if (t_tgmat%tmat_to_file) then
         open (69,ACCESS='direct',RECL=WLENGTH*4*LMMAXD*LMMAXD,FILE='tmat',FORM='unformatted')
      end if

      !read in non-collinear angles
      call read_angles(t_params,NATYP,THETA,PHI)

      do IATOM = 1,NATOMIMP
         I1=ATOMIMP(IATOM)
         if (KORBIT.EQ.0) then
            do ISPIN=1,NSPIN
               do IE=1,IELAST
                  if (I1<=NATYP) then
                     IREC = IE+IELAST*(ISPIN-1)+IELAST*NSPIN*(I1-1)
                     if (t_tgmat%tmat_to_file) then
                        read (69,REC=IREC) TMAT0
                     else
                        stop 'WRONG tmat_to_file for KKRFLEX writeout!!'
                        !not correctly read in with this option, to fix this communication is needed
                        !                      tmat0(:,:) = t_tgmat%tmat(:,:,irec)
                     end if
                  else
                     TMAT0=(0.0D0,0.0D0)
                  end if
                  write(6699,'(4I5,50000E14.7)') IATOM,ISPIN,IE,0,TMAT0
               end do !ie=1,ielast
            end do !ispin=1,nspin
         elseif (KORBIT.EQ.1) then
            ISPIN=1
            do IE=1,IELAST
               if (I1<=NATYP) then
                  IREC=IE+IELAST*(I1-1)
                  read(69,REC=IREC) TMAT0
                  !perform here a transformation from the local to the global
                  !spin-frame of reference
                  call ROTATEMATRIX(TMAT0,THETA(I1),PHI(I1),LMGF0D,0)
               else
                  TMAT0=(0d0,0d0)
               endif
               write(6699,'(4I5,50000E14.7)') IATOM,ISPIN,IE,0,TMAT0
            enddo
         endif
      end do
      close(69)
      close(6699)

      open (91,FILE='kkrflex_intercell_ref',STATUS='unknown')
      call version_print_header(91,'; '//md5sum_potential//'; '//md5sum_shapefun)
      write(91,*) '# Intercell potential of each atom'
      write(91,*) '# '
      write(91,*) '# NATOMIMP',NATOMIMP
      write(91,*) '# lmpot',lmpot
      write(91,*) '# KSHAPE',KSHAPE
      write(91,*) '# NATOMIMP, lmpot, ALAT VBC(1), VBC(2)'
      write(91,'(2I5,10F14.7)') NATOMIMP,lmpot,ALAT,VBC(1),VBC(2)
      do IATOM = 1,NATOMIMP      ! Bauer 2011-10-11
         I=ATOMIMP(IATOM)         !
         write(1337,*) 'ac2',I,HOSTIMP(I),lmpot,(VINTERS(LM,I),LM=1,lmpot)
         write(91,'(5000F14.7)') (VINTERS(LM,I),LM=1,lmpot)
      end do
      close(91)

      open (91,FILE='kkrflex_intercell_cmoms',STATUS='unknown')
      call version_print_header(91,'; '//md5sum_potential//'; '//md5sum_shapefun)
      write(91,*) '# Charge moments of each atom in the unit cell'
      write(91,*) '# Values given are CMOM + CMINST'
      write(91,*) '# First colums is the core charge other'
      write(91,*) '# other colums the charge moments'
      write(91,*) '# NATOMIMP',NATOMIMP
      write(91,*) '# lmpot',lmpot
      write(91,*) '# KSHAPE',KSHAPE

      do IATOM = 1,NATOMIMP
         I=ATOMIMP(IATOM)
         I1=KAOEZ(1,I)
         write(1337,*) 'NOQ',I,NOQ(I)
         if (NOQ(I)/=1 .and. NOQ(I)/=0)stop '[vmadelblk] VIRATOMS: NOQ/=1'
         if (NOQ(I)==0) then
           write(91,'(5000F14.7)') 0.0D0,(0.0D0,LM=1,lmpot)
         else
            if ( KSHAPE.NE.0 ) then
               write(91,'(5000F14.7)') ZAT(I1),  &
                  ((CMOM(LM,I1)+ CMINST(LM,I1))*CONC(I1),LM=1,lmpot)
            else
               write(91,'(5000F14.7)') ZAT(I1),(CMOM(LM,I1)*CONC(I1),LM=1,lmpot)
            end if
         end if
      end do
      close(91)
   end if

end subroutine writekkrflex
