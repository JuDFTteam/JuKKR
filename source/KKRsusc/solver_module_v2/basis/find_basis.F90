  subroutine find_basis(nb,phi,nr,wr,tol,numpan,numrcut)
! Gram-Schmidt with basis reduction

  implicit none

! --> input number of basis function / output reduced number of bfs
  integer(kind=i4b), intent(inout) :: nb
! --> input trial basis functions / output independent bfs only
  real(kind=r8b),    intent(inout) :: phi(nrmax,nbmax)
! --> actual number of points in radial mesh
  integer(kind=i4b), intent(in)    :: nr
! --> radial mesh weights
  real(kind=r8b),    intent(in)    :: wr(nrmax)
! --> tolerance for discarding a basis function
  real(kind=r8b),    intent(in)    :: tol
! --> Number of panels > 1
  integer(kind=i4b), intent(in) :: numpan, numrcut(numpan+1)
! -----------------------------------------------------------------
  integer(kind=i4b) :: ib, jb, nb1, ikeep(nbmax), ibasis(nbmax)
  real(kind=r8b)    :: rnorm(nbmax)
  complex(kind=c8b) :: work(nrmax), norm

  if (nb > nbmax) stop 'find_basis: nb > nbmax'
  if (nr > nrmax) stop 'find_basis: nr > nrmax'
  write(iodb,'("find_basis:",2i4,es12.1)') nr, nb, tol
! begin by normalizing the initial basis function
  do ib=1,nb
    work = phi(:,ib)*phi(:,ib)
    norm = radint(nr,work(1:nr),wr(1:nr),numpan,numrcut)
    phi(:,ib) = phi(:,ib)/sqrt(real(norm))
    write(iodb,'("find_basis:",i4,2es12.1)') ib, norm
  end do
! whether to keep a basis function
  ikeep(1:nb) = 1
! number of non-degenerate basis functions
  nb1 = 0
  do ib=1,nb
!   Gram-Schmidt
    do jb=1,ib-1
      if (ikeep(jb) == 1) then
!     subtract contributions from previous basis function
        work = phi(:,jb)*phi(:,ib)
        norm = radint(nr,work(1:nr),wr(1:nr),numpan,numrcut)
        phi(:,ib) = phi(:,ib) - real(norm)*phi(:,jb)
        write(iodb,'("GS   ",i4,2es12.3)') jb, norm
      end if
    end do
!   norm of current basis function
    work = phi(:,ib)*phi(:,ib)
    norm = radint(nr,work(1:nr),wr(1:nr),numpan,numrcut)
    rnorm(ib) = real(norm)
    write(iodb,'("Norm ",i4,es12.3)') ib, sqrt(rnorm(ib))
    if (rnorm(ib) < max(tol**2,gstol**2)) then
!   degenerate, discard it
      ikeep(ib) = 0
    else
!   one more
      nb1 = nb1 + 1
!   which is number...
      ibasis(nb1) = ib
!   normalize to the WS cell
      phi(:,ib) = phi(:,ib)/sqrt(rnorm(ib))
    end if
  end do
  write(iodb,'("ibasis=",100i4)') ibasis(1:nb1)
! how many non-degenerate basis functions
  nb = nb1
! fetch them to the front of the array
  do ib=1,nb
    phi(:,ib) = phi(:,ibasis(ib))
  end do
! the rest is history
  phi(:,nb+1:nbmax) = 0.d0
! All done!
  end subroutine find_basis

