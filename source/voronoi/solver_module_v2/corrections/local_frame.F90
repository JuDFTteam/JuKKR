  subroutine local_frame(ia,ja,magdir1,gfij,gf)
! Transforms the ij block from global to local frame
  use global

  implicit none

! Which GF block is it
  integer(kind=i4b), intent(in)    :: ia, ja
! Spin quantization axes
  real(kind=r8b),    intent(in)    :: magdir1(3,nasusc)
! GF block to be transformed to local frame
  complex(kind=c8b), intent(inout) :: gfij(nlmsb,nlmsb)
! Workspace
  complex(kind=c8b), intent(inout) :: gf(nlmsb,nlmsb)
! ----------------------------------------------------------------------
  real(kind=r8b), parameter :: uz(3) = (/0.d0,0.d0,1.d0/)
  complex(kind=c8b) :: spinroti(nsmax,nsmax), spinrotj(nsmax,nsmax)
  integer(kind=i4b) :: i3(3)
  integer(kind=i4b) :: i, j, is, js, ib, jb, ilm, jlm
  integer(kind=i4b) :: p, q, ps, qs, pb, qb, plm, qlm


! Spin rotation matrices from global to local frame
  call spin_rotation(magdir1(:,ia),uz,pauli,spinroti)
  call spin_rotation(magdir1(:,ja),uz,pauli,spinrotj)
! Rotate block to local frame
  gf(1:nlmsba(ia),1:nlmsba(ja)) = gfij(1:nlmsba(ia),1:nlmsba(ja))
  gfij(:,:) = 0.d0
  do j=1,nlmsba(ja)
    i3 = i2lmsb(:,j,ja)
    jb = i3(1); jlm = i3(2); js = i3(3)
    do i=1,nlmsba(ia)
      i3 = i2lmsb(:,i,ia)
      ib = i3(1); ilm = i3(2); is = i3(3)
      do q=1,nlmsba(ja)
        i3 = i2lmsb(:,q,ja)
        qb = i3(1); qlm = i3(2); qs = i3(3)
        do p=1,nlmsba(ia)
          i3 = i2lmsb(:,p,ia)
          pb = i3(1); plm = i3(2); ps = i3(3)
!         selection rules
          if (ib == pb .and. ilm == plm .and. jb == qb .and. jlm == qlm) then
            gfij(i,j) = gfij(i,j) + spinroti(is,ps)*gf(p,q)*conjg(spinrotj(js,qs))
          end if
        end do
      end do
    end do
  end do
! All done!
  end subroutine local_frame
