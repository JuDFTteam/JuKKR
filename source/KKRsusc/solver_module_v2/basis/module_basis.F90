  module basis
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!         Radial basis functions and spherical harmonics
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Global parameters and calculation options
  use global
! Save reference radial wfn from KKR program
  use ref_wfn_mod
! Save radial basis functions
  use save_wfns_mod
! Read radial basis functions
  use read_wfns_mod
! Make basis for each l and spin
  use new_basis1_mod
! Make basis for each l and both spins
  use new_basis2_mod
! Make basis for all l and both spins
  use new_basis3_mod
! Construct basis by Gram-Schmidt process
  use find_basis_mod
! Construct basis by diagonalizing overlap matrix
  use find_basis2_mod
! Output overlaps of basis functions
  use out_overlap_mod
! Set up the basis for the projected GFs
  use overlaps_gf_mod
! Set up the basis for the density
  use overlaps_susc

  implicit none


  end module basis
