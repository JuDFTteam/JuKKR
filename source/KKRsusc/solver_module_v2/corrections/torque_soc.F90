  subroutine torque_soc(ie,ia,vsocb,u0,u1)
! Matrix elements for torque version of anisotropy
  use global

  implicit none

! energy and atom
  integer(kind=i4b), intent(in)  :: ie, ia
! initial and final directions of magnetization
  real(kind=r8b),    intent(in)  :: u0(3), u1(3)
! SOC potential (choose energy)
  complex(kind=c8b), intent(in)  :: vsocb(nlmsb,nlmsb)
! ----------------------------------------------------------------------
  complex(kind=c8b), parameter :: iu = (0.d0,1.d0)
  real(kind=r8b),    parameter :: tol = 1.d-6
  integer(kind=i4b) :: i3(3), i2(2), ilms, jlms
  integer(kind=i4b) :: i, j, is, js, ib, jb, ilm, jlm
  integer(kind=i4b) :: q, qs, qb, qlm, k
  real(kind=r8b)    :: axis(3), axislen

! Rotation axis
  axis(1) = u0(2)*u1(3) - u0(3)*u1(2)
  axis(2) = u0(3)*u1(1) - u0(1)*u1(3)
  axis(3) = u0(1)*u1(2) - u0(2)*u1(1)
  axislen = sqrt(dot_product(axis,axis))
  if (axislen < tol) then
    torque(:,:,:,ia,ie) = 0.d0
    return
  end if
  axis = axis/axislen
! loop over cartesian Pauli matrices
  do k=1,3
    torque(:,:,k,ia,ie) = 0.d0
    do j=1,nlmsba(ia)
      i3 = i2lmsb(:,j,ia)
      jb = i3(1); jlm = i3(2); js = i3(3)
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        do q=1,nlmsba(ia)
          i3 = i2lmsb(:,q,ia)
          qb = i3(1); qlm = i3(2); qs = i3(3)
!         selection rules
          if (ib == qb .and. ilm == qlm) torque(i,j,k,ia,ie) = torque(i,j,k,ia,ie) + 0.5d0*iu*axis(k)*pauli(is,qs,k)*vsocb(q,j)
          if (jb == qb .and. jlm == qlm) torque(i,j,k,ia,ie) = torque(i,j,k,ia,ie) - 0.5d0*iu*axis(k)*vsocb(i,q)*pauli(qs,js,k)
        end do
      end do
    end do
  end do
! All done!
  end subroutine torque_soc
