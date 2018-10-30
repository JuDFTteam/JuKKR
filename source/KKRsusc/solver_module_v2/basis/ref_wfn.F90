  subroutine ref_wfn(ia,il,is,ib,phi,numpan,numrcut)
! Save wavefunction for projection
! phi is computed at real energy without normalization at WS radius
! Should be inside a loop over the energies used for projection
! If using irregular solutions put them at the end of the list
  use global

  implicit none

! --> susc atom
  integer(kind=i4b), intent(in) :: ia
! --> angular momentum
  integer(kind=i4b), intent(in) :: il
! --> spin channel
  integer(kind=i4b), intent(in) :: is
! --> basis function
  integer(kind=i4b), intent(in) :: ib
! --> regular scattering solution
  complex(kind=c8b), intent(in) :: phi(nrmax)
! --> Number panels > 1
  integer(kind=i4b), intent(in) :: numpan, numrcut(numpan+1)
! -----------------------------------------------------------------
  integer(kind=i4b) :: nr0, nr1, nr
  complex(kind=c8b) :: proj(nrmax), work(nrmax), norm
  complex(kind=c8b), external :: radint

! Good morning
  if (noparameters) stop 'ref_wfn: run init_param!'
  if (normesh(ia))  stop 'ref_wfn: save rmesh first!'
  if (ib > nbmax)   stop 'ref_wfn: ib > nbmax!'
! Store trial basis functions
  nr0 = nrpts0(ia); nr1 = nrpts1(ia); nr = nr1 - nr0 + 1
  phiref(1:nr,ib,il,is,ia) = real(phi(nr0:nr1))
! Make overall sign positive
  work = phiref(1:nr,ib,il,is,ia)
  norm = radint(nr,work(1:nr),drproj(1:nr,ia),numpan,numrcut(1:numpan+1))
  if (real(norm) < 0.d0) phiref(1:nr,ib,il,is,ia) = -phiref(1:nr,ib,il,is,ia)
! Normalize to 1
  work = phiref(1:nr,ib,il,is,ia)*phiref(1:nr,ib,il,is,ia)
  norm = radint(nr,work(1:nr),drproj(1:nr,ia),numpan,numrcut(1:numpan+1))
  phiref(1:nr,ib,il,is,ia) = phiref(1:nr,ib,il,is,ia)/sqrt(real(norm))
! Check if the basis is completely stored
  if (ib == iwsusc(il,is,ia)) nowfns(il,is,ia) = .false.
! All done!
  end subroutine ref_wfn
