  subroutine torque_soc2(ie,ia,vsoc,u0,u1)
! assembles the SOC potential in the projection basis
  use global

  implicit none

! Which energy
  integer(kind=i4b), intent(in)    :: ie
! Which atom
  integer(kind=i4b), intent(in)    :: ia
! Spherical spin-averaged radial SOC potential
  complex(kind=c8b), intent(in)    :: vsoc(nrmax)
! initial and final directions of magnetization
  real(kind=r8b),    intent(in)  :: u0(3), u1(3)
! -----------------------------------------------------------------
  real(kind=r8b),    parameter :: tol = 1.d-6
  integer(kind=i4b) :: is, js, ib, jb, nbi, nbj, i3(3), i2(2)
  integer(kind=i4b) :: jm, jl, im, il, jlm, ilm, i, j
  integer(kind=i4b) :: nr
  real(kind=r8b)    :: dr(nrmax)
  complex(kind=c8b) :: work(nrmax), norm
  integer(kind=i4b) :: imax, jmax
  complex(kind=c8b), external :: radint
  real(kind=r8b)    :: axis(3), axislen

! use the structure of L.S
! also use the fact that the basis was constructed per l channel
!  maxnorm = 0.d0
  nr = nrpts(ia)
  dr(1:nr) = drmesh(1:nr,ia)
! ----------------------------------------------------------------------
! Rotation axis for dR/dtheta
  axis(1) = u0(2)*u1(3) - u0(3)*u1(2)
  axis(2) = u0(3)*u1(1) - u0(1)*u1(3)
  axis(3) = u0(1)*u1(2) - u0(2)*u1(1)
  axislen = sqrt(dot_product(axis,axis))
  if (axislen < tol) then
    axis = (/0.d0,1.d0,0.d0/)
  else
    axis = axis/axislen
  end if
!  if (ie == nescf) write(*,'("axis=",3f6.3)') axis
! ----------------------------------------------------------------------
  torque(:,:,:,ia,ie) = 0.d0
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
      if (il == jl) then
!       revise the basis
        work(1:nr) = phiref(1:nr,ib,il,is,ia)*vsoc(1:nr)*phiref(1:nr,jb,il,js,ia)
        norm = radint(nr,work,dr,npanat(ia),ircutat(:,ia))
!        if (abs(norm) > maxnorm) then
!          maxnorm = abs(norm)
!          imax = i; jmax = j
!        end if
!       theta:
        torque(i,j,1,ia,ie) = torque(i,j,1,ia,ie) + norm*axis(1)*(pauli(is,js,2)*lorb(ilm,jlm,3) - pauli(is,js,3)*lorb(ilm,jlm,2))  &
                                                  + norm*axis(2)*(pauli(is,js,3)*lorb(ilm,jlm,1) - pauli(is,js,1)*lorb(ilm,jlm,3))  &
                                                  + norm*axis(3)*(pauli(is,js,1)*lorb(ilm,jlm,2) - pauli(is,js,2)*lorb(ilm,jlm,1))
!       phi:
        torque(i,j,2,ia,ie) = torque(i,j,2,ia,ie) + norm*(pauli(is,js,1)*lorb(ilm,jlm,2) - pauli(is,js,2)*lorb(ilm,jlm,1))
!       nothing
        torque(i,j,3,ia,ie) = 0.d0
      end if
!   ----------------------------------------------------------------
    end do
  end do
!  write(iodb,'(" SOC kernel, ia=",i4)') ia
!  write(iodb,'("vsocb norm:",6i4,2es16.8)') i2lmsb(:,imax,ia), i2lmsb(:,jmax,ia), maxnorm
!  write(iodb,'("  ib ilm  is, jb jlm  js, matrix element")')
!  do j=1,nlmsba(ia)
!    do i=1,nlmsba(ia)
!      if (abs(vsocb(i,j)) > 1.d-8) then
!        write(iodb,'(6i4,2es10.2)') i2lmsb(:,i,ia), i2lmsb(:,j,ia), vsocb(i,j)
!      end if
!    end do
!  end do
! All done
  end subroutine torque_soc2

