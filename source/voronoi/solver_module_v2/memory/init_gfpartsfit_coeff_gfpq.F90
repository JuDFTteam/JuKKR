  subroutine init_gfpartsfit_coeff_gfpq(my_rank,mpi_size,mpi_iebounds)
! Allocates arrays for different parts of projected GF
  use global

  implicit none

  integer(kind=i4b) :: is, il, im, ib, i, lm, ia, ilms, jlms, my_rank, mpi_size
  real(kind=r8b)    :: ram
! --> mpi_bounds
  integer(kind=i4b), intent(in) :: mpi_iebounds(2,0:mpi_size-1) 

! -----------------------------------------------------------------------
!   Storage for structural GF and coefficients in full projection basis
! -----------------------------------------------------------------------
  allocate(pzl(nlmsb,nlms,nasusc,mpi_iebounds(1,my_rank):mpi_iebounds(2,my_rank)))             ! storage for lhs pz in lms with basis labels
  allocate(pzr(nlmsb,nlms,nasusc,mpi_iebounds(1,my_rank):mpi_iebounds(2,my_rank)))             ! storage for rhs pz in lms with basis labels
  allocate(gfpq(nlmsb,nlmsb,nasusc,mpi_iebounds(1,my_rank):mpi_iebounds(2,my_rank)))           ! storage for onsite GF pq in lms with basis labels
  pzl = 0.d0; pzr = 0.d0;  gfpq = 0.d0
  if (isra == 1) then
    allocate(fzl(nlmsb,nlms,nasusc,mpi_iebounds(1,my_rank):mpi_iebounds(2,my_rank)))           ! storage for lhs fz
    allocate(fzr(nlmsb,nlms,nasusc,mpi_iebounds(1,my_rank):mpi_iebounds(2,my_rank)))           ! storage for rhs fz
    allocate(gfps(nlmsb,nlmsb,nasusc,mpi_iebounds(1,my_rank):mpi_iebounds(2,my_rank)))         ! storage for ps
    allocate(gffq(nlmsb,nlmsb,nasusc,mpi_iebounds(1,my_rank):mpi_iebounds(2,my_rank)))         ! storage for fq
    allocate(gffs(nlmsb,nlmsb,nasusc,mpi_iebounds(1,my_rank):mpi_iebounds(2,my_rank)))         ! storage for fs
    fzl = 0.d0; fzr = 0.d0; gfps = 0.d0; gffq = 0.d0; gffs = 0.d0
  end if
! -----------------------------------------------------------------------
! All done!
  end subroutine init_gfpartsfit_coeff_gfpq
