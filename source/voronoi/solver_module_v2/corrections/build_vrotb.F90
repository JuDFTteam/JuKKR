  subroutine build_vrotb(ia,br,u0,u1,vrotb)
! assembles the SOC potential in the projection basis
  use global, only: i4b, r8b, c8b, iodb, nsmax, nrmax, nrpts, drmesh, phiref, nlmsba, nlmsb, nlms, i2lmsb, i2lms, i2lm, atol, pauli, npanat, ircutat 

  implicit none

! Which atom
  integer(kind=i4b), intent(in)    :: ia
! Initial and final orientations of the potential
  real(kind=r8b),    intent(in)    :: u0(3), u1(3)
! Spherical spin-averaged radial exchange potential
  real(kind=c8b),    intent(in)    :: br(nrmax)
! Rotated potential in the basis (added to)
  complex(kind=c8b), intent(inout) :: vrotb(nlmsb,nlmsb)
! -----------------------------------------------------------------
  integer(kind=i4b) :: is, js, ib, jb, nbi, nbj, i3(3), i2(2)
  integer(kind=i4b) :: jm, jl, im, il, jlm, ilm, i, j
  integer(kind=i4b) :: nr
  real(kind=r8b)    :: dr(nrmax)
  complex(kind=c8b) :: work(nrmax), norm, spin(nsmax,nsmax)
  integer(kind=i4b) :: imax, jmax
  complex(kind=c8b), external :: radint

! use the fact that the basis was constructed per l channel
! potential is spherical -> diagonal in lm
!  maxnorm = 0.d0
  nr = nrpts(ia)
  dr(1:nr) = drmesh(1:nr,ia)
! build the spin matrix
  spin = 0.d0
  do i=1,3
    spin = spin + pauli(:,:,i)*(u1(i) - u0(i))
  end do
  write(iodb,'("spin matrix ia=",i4)') ia
  write(iodb,'(4f10.5)') ((spin(i,j),j=1,nsmax),i=1,nsmax)
! now matrix elements in the basis
  do j=1,nlmsba(ia)
    i3 = i2lmsb(:,j,ia)
    jb = i3(1); jlm = i3(2); js = i3(3)
    i2 = i2lm(:,jlm)
    jm = i2(1); jl = i2(2)
    do i=1,nlmsba(ia)
      i3 = i2lmsb(:,i,ia)
      ib = i3(1); ilm = i3(2); is = i3(3)
      i2 = i2lm(:,ilm)
      im = i2(1); il = i2(2)
!   ----------------------------------------------------------------
!   selection rules
      if (ilm == jlm) then
!     revise the basis
        work(1:nr) = phiref(1:nr,ib,il,is,ia)*br(1:nr)*phiref(1:nr,jb,il,js,ia)
        norm = radint(nr,work,dr,npanat(ia),ircutat(:,ia))
!        if (abs(norm) > maxnorm) then
!          maxnorm = abs(norm)
!          imax = i; jmax = j
!        end if
        vrotb(i,j) = vrotb(i,j) + norm*spin(is,js)!*16.d0*atan(1.d0)*2.d0
      end if
!   ----------------------------------------------------------------
    end do
  end do
!  write(iodb,'(" SOC kernel, ia=",i4)') ia
!  write(iodb,'("vrotb norm:",6i4,2es16.8)') i2lmsb(:,imax,ia), i2lmsb(:,jmax,ia), maxnorm
!  write(iodb,'("  ib ilm  is, jb jlm  js, matrix element")')
!  do j=1,nlmsba(ia)
!    do i=1,nlmsba(ia)
!      if (abs(vrotb(i,j)) > 1.d-8) then
!        write(iodb,'(6i4,2es10.2)') i2lmsb(:,i,ia), i2lmsb(:,j,ia), vrotb(i,j)
!      end if
!    end do
!  end do
! All done
  end subroutine build_vrotb

