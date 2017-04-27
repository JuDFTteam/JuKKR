  subroutine new_basis1(nb,ia,il,is,my_rank)
! The list of trial wavefunctions is reduced to a minimal basis
! New basis functions with norm below tol are discarded
  use global

  implicit none

! --> reduced number of basis functions
  integer(kind=i4b), intent(out) :: nb
! --> susc atom
  integer(kind=i4b), intent(in)  :: ia
! --> angular momentum
  integer(kind=i4b), intent(in)  :: il
! --> spin channel
  integer(kind=i4b), intent(in)  :: is
! --> mpi 
  integer(kind=i4b), intent(in)  :: my_rank
! -----------------------------------------------------------------
  integer(kind=i4b) :: nr, nb1, ib
  real(kind=r8b)    :: basis(nrmax,nbmax)

! Good morning
  if (noparameters) stop 'new_basis1: run init_param!'
  nr = nrpts(ia)
  nb = iwsusc(il,is,ia)
  if (my_rank == 0) then
    if (lhdio) write(iodb,'("new_basis1:",5i4)') ia, il, is, nr, nb
  end if ! my_rank 
! no wfns for this l and spin -- intentional
  if (nb == 0) then
    nowfns(il,is,ia) = .false.
    nobasis(il,is,ia) = .false.
    return
  end if
! no wfns for this l and spin -- unintentional
  if (nowfns(il,is,ia)) stop 'new_basis1: save wavefunctions first!'
! Store trial basis functions and do GS
  nb1 = nb
  if (nb1 > nbmax) stop 'new_basis1: increase nbmax!'
  basis(1:nr,1:nb1) = phiref(1:nr,1:nb1,il,is,ia)
! Decide how to orthogonalize
  if (ibasismethod == 0) then
!   do nothing
  else if (ibasismethod == 1) then
    call find_basis1(nb1,nbmax,basis,nr,nrmax,drproj(:,ia),basistol,npanat(ia),ircutat(:,ia),.true.,my_rank)
  else if (ibasismethod == 2) then
    call find_basis2(nb1,nbmax,basis,nr,nrmax,drproj(:,ia),basistol,npanat(ia),ircutat(:,ia),.true.,my_rank)
  else
    stop 'new_basis1: unknown basis construction method'
  end if
  phiref(1:nr,1:nb,il,is,ia) = basis(1:nr,1:nb)
! The size of the basis is updated
  nb = nb1
  iwsusc(il,is,ia) = nb
! Basis constructed
  nobasis(il,is,ia) = .false.
! All done!
  end subroutine new_basis1
