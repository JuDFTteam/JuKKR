#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

!-------------------------------------------------------------------------------
!> A wrapper for the subroutine STARTB1
!> Returns: EFERMI (from first entry in potential file)

! determine the record length
subroutine STARTB1_wrapper_new(alat, LPOT,NSPIN, &
                           NTCELL, &
                           EFERMI,ZAT, radius_muffin_tin, &
                           IPAND,IRID,NFUND,IRMD,NCELLD,NAEZD,IRNSD)
  implicit none

  ! Arguments
  double precision, intent(in) :: alat
  INTEGER, INTENT(IN) :: IPAND
  INTEGER, INTENT(IN) :: IRID
  INTEGER, INTENT(IN) :: NFUND
  INTEGER, INTENT(IN) :: IRMD
  INTEGER, INTENT(IN) :: NCELLD
  INTEGER, INTENT(IN) :: NAEZD
  INTEGER, INTENT(IN) :: IRNSD
  DOUBLE PRECISION :: EFERMI
  INTEGER :: LPOT
  INTEGER :: NSPIN
  DOUBLE PRECISION, dimension(*) :: ZAT
  double precision, dimension(naezd), intent(in) :: radius_muffin_tin
  INTEGER, dimension(*) :: NTCELL

  integer :: max_reclen

  ! write atoms file and get maximum record length for potential file
  call write_atoms_file(alat, NSPIN, &
                            NTCELL, &
                            ZAT, radius_muffin_tin, &
                            NAEZD, max_reclen)

end subroutine

!------------------------------------------------------------------------------
!> Write the 'atoms' file (= binary analogue to atominfo - used for parallel reading in kkr2)
!> and determine the record length ('max_reclen') for the binary direct access potential file.
subroutine write_atoms_file(alat, NSPIN, &
                            NTCELL, &
                            ZAT, radius_muffin_tin, &
                            NAEZD, max_reclen)
  use read_formatted_mod
  use BasisAtom_mod
  implicit none

  ! Arguments
  double precision, intent(in) :: alat
  INTEGER :: NSPIN
  DOUBLE PRECISION, dimension(*) :: ZAT
  double precision, dimension(naezd), intent(in) :: radius_muffin_tin
  integer, intent(in) :: NAEZD
  INTEGER, dimension(*) :: NTCELL
  integer, intent(out) :: max_reclen

  integer :: iatom
  integer :: ispin
  integer :: lpot_atom
  integer :: reclen

  type (PotentialEntry) :: entry(2)
  type (BasisAtom) :: atom
  integer, parameter :: UNIT = 13

  max_reclen = 0
  reclen = 0

  call openBasisAtomDAFile(atom, 37, 'atoms')

  open (UNIT, file='potential', status='old', form='formatted')
    do iatom = 1, naezd

      do ispin = 1, nspin
        call create_read_PotentialEntry(entry(ispin), UNIT)

        ! do some consistency checks
        if (abs(ZAT(iatom) - entry(ispin)%header%Z_nuclear) > 1.d-8) then
          write(*,*) "ERROR: Mismatch of nuclear charge between atominfo and potential file for entry: ", iatom
          write(*,*) ZAT(iatom), entry(ispin)%header%Z_nuclear
          STOP
        endif

        if (abs(alat - entry(ispin)%header%alat) > 1.d-8) then
          write(*,*) "WARNING: alat in potential file is not the same as in input.conf."
        endif
      enddo

      lpot_atom = int( sqrt(dble(entry(1)%sblock%lmpot)) - 1.0d0 + 0.01d0)
      CHECKASSERT( (lpot_atom + 1)**2 == entry(1)%sblock%lmpot)

      ! number of radial points must be the same for every spin direction
      CHECKASSERT( entry(1)%sblock%irt1p == entry(nspin)%sblock%irt1p )
      CHECKASSERT( entry(1)%sblock%irns == entry(nspin)%sblock%irns )
      CHECKASSERT( entry(1)%sblock%lmpot == entry(nspin)%sblock%lmpot )

      call createBasisAtom(atom, 1, lpot_atom, nspin, entry(1)%sblock%IRT1P - entry(1)%sblock%IRNS, entry(1)%sblock%IRT1P)  ! create dummy basis atom

      inquire (iolength = reclen) atom%potential%VINS, &
                                  atom%potential%VISP, &
                                  atom%core%ECORE

      max_reclen = max(max_reclen, reclen)

      ! write file 'atoms'
      atom%atom_index = iatom
      atom%cell_index = NTCELL(iatom)
      atom%Z_nuclear = ZAT(iatom)

      atom%radius_muffin_tin = radius_muffin_tin(iatom)

      atom%core%NCORE = 0
      atom%core%LCORE = 0
      atom%core%ITITLE = 0

      do ispin = 1, nspin
        atom%core%NCORE(ispin) = entry(ispin)%csblock%NCORE
        atom%core%LCORE(:, ispin) = entry(ispin)%csblock%LCORE
        atom%core%ITITLE(:, ispin) = entry(ispin)%header%ITITLE
      enddo

      call writeBasisAtomDA(atom, 37, iatom)

      call destroyBasisAtom(atom)

      do ispin = 1, nspin
        call destroy_PotentialEntry(entry(ispin))
      enddo

    enddo ! end loop over atoms
  close(13)

  call closeBasisAtomDAFile(37)

end subroutine


!------------------------------------------------------------------------------
!> Writes the file 'atoms' with some basic information about each atom
!> like nuclear charge, muffin-tin-radius, nspin, core config...
!subroutine writeAtomData(naezd, lpot, nspin, irmd, irnsd, ntcell, &
!                         radius_muffin_tin, ZAT, ITITLE, NCORE, LCORE)
!  use BasisAtom_mod
!  use read_formatted_mod
!  implicit none
!
!  integer, intent(in) :: naezd
!  integer, intent(in) :: lpot
!  integer, intent(in) :: nspin
!  integer, intent(in) :: irmd
!  integer, intent(in) :: irnsd
!  integer, intent(in) :: ntcell(naezd)
!  double precision, intent(in) :: radius_muffin_tin(naezd)
!  double precision, intent(in) :: ZAT(naezd)
!  integer, intent(in) :: ITITLE(20, naezd*nspin)
!  integer, intent(in) :: NCORE(naezd*nspin)
!  integer, intent(in) :: LCORE(20, naezd*nspin)
!
!  type (BasisAtom) :: atom
!  integer :: ii, ispin, ipot
!  integer :: max_reclen
!
!  call createBasisAtom(atom, 1, lpot, nspin, (irmd-irnsd), irmd)  ! create dummy basis atom
!
!  call openBasisAtomDAFile(atom, 37, 'atoms')
!#ifndef TASKLOCAL_FILES
!  call openBasisAtomPotentialIndexDAFile(atom, 38, 'vpotnew.0.idx')
!#endif
!
!  inquire (iolength = max_reclen) atom%potential%VINS, &
!                                  atom%potential%VISP, &
!                                  atom%core%ECORE
!
!  do ii = 1, naezd
!    atom%atom_index = ii
!    atom%cell_index = NTCELL(ii)
!    atom%Z_nuclear = ZAT(ii)
!    atom%radius_muffin_tin = radius_muffin_tin(ii)
!
!    atom%core%NCORE = 0
!    atom%core%LCORE = 0
!    atom%core%ITITLE = 0
!
!    do ispin = 1, nspin
!      ipot = NSPIN * (ii-1) + ispin
!      atom%core%NCORE(ispin) = NCORE(ipot)
!      atom%core%LCORE(:, ispin) = LCORE(:, ipot)
!      atom%core%ITITLE(:, ispin) = ITITLE(:, ipot)
!    enddo
!
!    call writeBasisAtomDA(atom, 37, ii)
!#ifndef TASKLOCAL_FILES
!    call writeBasisAtomPotentialIndexDA(atom, 38, ii, max_reclen)
!#endif
!
!  enddo
!#ifndef TASKLOCAL_FILES
!  call closeBasisAtomPotentialIndexDAFile(38)
!#endif
!  call closeBasisAtomDAFile(37)
!
!  call destroyBasisAtom(atom)
!
!end subroutine

