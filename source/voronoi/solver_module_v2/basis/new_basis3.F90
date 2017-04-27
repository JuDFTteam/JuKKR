  subroutine new_basis3(nb,ia,my_rank)
! The list of trial wavefunctions is reduced to a minimal basis
! New basis functions with norm below tol are discarded
  use global

  implicit none

! --> reduced number of basis functions
  integer(kind=i4b), intent(out) :: nb
! --> susc atom
  integer(kind=i4b), intent(in)  :: ia
! --> Mpi 
  integer(kind=i4b), intent(in)  :: my_rank
! -----------------------------------------------------------------
  integer(kind=i4b) :: nr, nb1, nb2, ib, ib0, ib1, il, is, ns
  real(kind=r8b)    :: basis(nrmax,nbmax)

! Good morning
  if (noparameters) stop 'new_basis3: run init_param!'
  ns = issusc(ia)
  nr = nrpts(ia)
! Store trial basis functions and do GS
  ib1 = 0
  do il=0,nlmax
    do is=1,ns
      nb1 = iwsusc(il,is,ia)
!     no wfns for this l and spin -- intentional
      if (nb1 == 0) then
        nowfns(il,is,ia) = .false.
        cycle
      end if
!     no wfns for this l and spin -- unintentional
      if (nowfns(il,is,ia)) stop 'new_basis3: save wfns'
      ib0 = ib1 + 1
      ib1 = ib1 + nb1
      if (ib1 > nbmax) stop 'new_basis3: increase nbmax!'
      basis(1:nr,ib0:ib1) = phiref(1:nr,1:nb1,il,is,ia)
      if (my_rank) then 
        if (lhdio) write(iodb,'(3i4)') nb1, ib0, ib1
      end if ! my_rank
    end do
  end do
  nb1 = ib1
  nb  = nb1
  if (my_rank) then
    if (lhdio) write(iodb,'("new_basis3:",5i4)') ia, nlmax, ns, nr, nb
  end if ! my_rank
! No wfns for this atom -- hopefully intentional
  if (nb == 0) then
    nobasis(0:nlmax,1:ns,ia) = .false.
    return
  end if
! Decide how to orthogonalize
  if (ibasismethod == 1) then
    call find_basis1(nb1,nbmax,basis,nr,nrmax,drproj(:,ia),basistol,npanat(ia),ircutat(:,ia),.true.,my_rank)
  else if (ibasismethod == 2) then
    call find_basis2(nb1,nbmax,basis,nr,nrmax,drproj(:,ia),basistol,npanat(ia),ircutat(:,ia),.true.,my_rank)
  else
    stop 'new_basis3: unknown basis construction method'
  end if
! A bit redundant, but probably easier to use
  do il=0,nlmax
    do is=1,ns
      phiref(1:nr,1:nb,il,is,ia) = basis(1:nr,1:nb)
    end do
  end do
! The size of the basis is updated
  nb = nb1
  iwsusc(0:nlmax,1:ns,ia) = nb
! Basis constructed
  nobasis(0:nlmax,1:ns,ia) = .false.
! All done!
  end subroutine new_basis3
