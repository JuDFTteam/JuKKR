  subroutine find_basis2(nb,nbmax,phi,nr,nrmax,wr,tol,numpan,numrcut,printout,my_rank)
! Diagonalization of the overlap matrix
  use global, only: i4b, r8b, c8b, iodb

  implicit none

! --> input number of basis function / output reduced number of bfs
  integer(kind=i4b), intent(inout) :: nb
! --> dimensions of array phi
  integer(kind=i4b), intent(in)    :: nbmax, nrmax
! --> input trial basis functions / output independent bfs only
  real(kind=r8b),    intent(inout) :: phi(nrmax,nbmax)
! --> actual number of points in radial mesh
  integer(kind=i4b), intent(in)    :: nr
! --> radial mesh weights
  real(kind=r8b),    intent(in)    :: wr(nrmax)
! --> tolerance for discarding a basis function
  real(kind=r8b),    intent(in)    :: tol
! --> whether to print something
  logical,           intent(in)    :: printout
! --> Number of panels > 1
  integer(kind=i4b), intent(in)    :: numpan, numrcut(numpan+1)
! --> Mpi
  integer(kind=i4b), intent(in)    :: my_rank
! -----------------------------------------------------------------
  integer(kind=i4b) :: ib, jb, nb1, ikeep(nbmax), ibasis(nbmax)
  real(kind=r8b)    :: rnorm(nbmax), overlaps(nb,nb)
  complex(kind=c8b) :: work(nrmax), norm
  real(kind=r8b)    :: evals(nb), rwork(2*nb*(1 + 3*nb) + 1)
  integer(kind=i4b) :: iwork(5*nb+3), info, lrwork, liwork
  real(kind=r8b)    :: invsqrt(nb,nb), newphi(nrmax,nbmax)
  complex(kind=c8b), external :: radint

  if (nb > nbmax) stop 'find_basis2: nb > nbmax'
  if (nr > nrmax) stop 'find_basis2: nr > nrmax'
  if (my_rank == 0) then
    if (printout) write(iodb,'("find_basis2:",2i4,es12.1)') nr, nb, tol
  end if ! my_rank 
! begin by computing the overlaps
  do jb=1,nb
    do ib=1,nb
      work = phi(:,ib)*phi(:,jb)
      norm = radint(nr,work(1:nr),wr(1:nr),numpan,numrcut)
      overlaps(ib,jb) = real(norm)     
    end do
    if (my_rank == 0) then
      if (printout) write(iodb,'("find_basis2:",i4,1000es12.1)') jb, overlaps(1:nb,jb)
    end if ! my_rank
  end do
! diagonalize; eigenvalues in ascending order
  lrwork  = 2*nb*(1 + 3*nb) + 1
  liwork = 5*nb + 3
  call dsyevd('V','U',nb,overlaps,nb,evals,rwork,lrwork,iwork,liwork,info)
  if (my_rank == 0) then
    if (printout) write(iodb,'("find_basis2: eigen")')
  end if ! my_rank
  do jb=1,nb
    if (sum(overlaps(1:nb,jb)) < 0.d0) overlaps(1:nb,jb) = -overlaps(1:nb,jb)
    if (my_rank == 0) then
      if (printout) write(iodb,'(es12.3," | ",1000f12.6)') evals(jb), overlaps(1:nb,jb)
    end if ! my_rank
  end do
! new basis
  nb1 = 0; newphi = 0.d0
  do jb=1,nb
    if (abs(evals(nb-jb+1)) > tol*evals(nb)) then
      nb1 = nb1 + 1
      do ib=1,nb
        newphi(:,jb) = newphi(:,jb) + phi(:,ib)*overlaps(ib,nb-jb+1)
      end do
      work = newphi(:,jb)*newphi(:,jb)
      norm = radint(nr,work(1:nr),wr(1:nr),numpan,numrcut)
!     this should return the eigenvalues
      if (my_rank == 0) then
        if (printout) write(iodb,'("find_basis2: newphi",i4,2es12.3)') jb, norm
      end if ! my_rank
      newphi(:,jb) = newphi(:,jb)/sqrt(norm)
    end if
  end do
! new overlaps
  do jb=1,nb1
    do ib=1,nb1
      work = newphi(:,ib)*newphi(:,jb)
      norm = radint(nr,work(1:nr),wr(1:nr),numpan,numrcut)
      overlaps(ib,jb) = real(norm)
    end do 
    if (my_rank == 0) then 
      if (printout) write(iodb,'("find_basis2:",i4,1000es12.1)') jb, overlaps(1:nb1,jb)
    end if ! my_rank
  end do
! new non-degenerate basis functions
  nb  = nb1
  phi = newphi
! All done!
  end subroutine find_basis2
