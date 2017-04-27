  subroutine init_basis(my_rank)
! Allocate the storage for wfns and projection coefficients
  use global

  implicit none

  integer(kind=i4b) :: is, il, m, ib, i, lm, ia, ilms, jlms, my_rank
  real(kind=r8b)    :: ram

! -----------------------------------------------------------------------
!                  Basis and projection coefficients
! -----------------------------------------------------------------------
! How much memory is being allocated for basis and projection
! The basis can change from iteration to iteration
  ram = 8.d0*nrmax*nbmax*(nlmax+1)*nsmax*nasusc       ! phiref
  ram = ram + 16.d0*(nbmax+1)*nbmax*(nlmax+1)*nsmax*nasusc  ! pzc, pqc
  if (isra == 1) ram = ram + 16.d0*(3*nbmax+1)*nbmax*(nlmax+1)*nsmax*nasusc  ! fzc, psc, fqc, fsc
  ram = ram/(1024.d0**2)
  if (my_rank == 0) then
    write(*,'(/" init_basis: projection RAM=",f16.3," MB")') ram
  end if ! my_rank
! The basis is currently constructed per l-channel
  allocate(phiref(nrmax,nbmax,0:nlmax,nsmax,nasusc))  ! reference wfns
  phiref = 0.d0; nowfns  = .true.; nobasis = .true.
  allocate(pzc(nbmax,0:nlmax,nsmax,nasusc))           ! coefficients of pz in radial basis
  allocate(pqc(nbmax,nbmax,0:nlmax,nsmax,nasusc))     ! coefficients of pzqz in radial basis
  pzc = 0.d0; pqc = 0.d0
  if (isra == 1) then
    allocate(fzc(nbmax,0:nlmax,nsmax,nasusc))         ! coefficients of fz in radial basis
    allocate(psc(nbmax,nbmax,0:nlmax,nsmax,nasusc))   ! coefficients of pzsz in radial basis
    allocate(fqc(nbmax,nbmax,0:nlmax,nsmax,nasusc))   ! coefficients of fzqz in radial basis
    allocate(fsc(nbmax,nbmax,0:nlmax,nsmax,nasusc))   ! coefficients of fzsz in radial basis
    fzc = 0.d0; psc = 0.d0; fqc = 0.d0; fsc = 0.d0
  end if
  allocate(noregcoeffs(nsmax,nesusc))                             ! were reg coeffs computed?
  allocate(noirrcoeffs(nsmax,nesusc))                             ! were irr coeffs computed?
  allocate(noregcoeffs_soc(nsmax,nsmax,nesusc))                   ! were reg coeffs computed (SOC)?
  allocate(noirrcoeffs_soc(nsmax,nsmax,nesusc))                   ! were irr coeffs computed (SOC)?
  if (lsusc) allocate(mtotsusc(nbmax,nbmax,lmmax,lmmax,nasusc2))  ! total magnetization
  if (lsusc) allocate(mxcsusc(nbmax,nbmax,lmmax,lmmax,nasusc2))   ! magnetization from xc potential
  if (lsusc) allocate(msocsusc(nbmax,nbmax,lmmax,lmmax,nasusc2))   ! magnetization from SOC+Bext potential
  noregcoeffs  = .true.
  noirrcoeffs  = .true.
  noregcoeffs_soc = .true.
  noirrcoeffs_soc = .true.

! All done!
  end subroutine init_basis
