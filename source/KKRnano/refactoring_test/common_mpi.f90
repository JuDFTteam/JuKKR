! This is a replacement for the COMMON block
! named 'MPI' in the old kkrnano, which conflicted with
! the module mpi

module common_mpi
  implicit none
  save

  integer :: MYRANK
  integer :: NROFNODES
end module common_mpi
