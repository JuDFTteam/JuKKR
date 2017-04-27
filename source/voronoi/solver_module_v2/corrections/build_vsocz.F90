  subroutine build_vsocz(vsocb,magdir1,vsocz)
! computes the SOC z contribution to the spin splitting
  use global

  implicit none

! longitudinal part of SOC potential
  complex(kind=c8b), intent(in)  :: vsocb(nlmsb,nlmsb,nasusc)
! spin quantization axes
  real(kind=r8b),    intent(in)  :: magdir1(3,nasusc)
! contribution from SOC z to spin splitting
  complex(kind=c8b), intent(out) :: vsocz(nbmax,nbmax,lmmax,lmmax,nasusc2)
! ----------------------------------------------------------------------
  integer(kind=i4b) :: i3(3), ia, ia2, i1, ib1, ilm1, is1, i2, ib2, ilm2, is2

  vsocz = 0.d0
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do i2=1,nlmsba(ia)
    do i1=1,nlmsba(ia)
      i3 = i2lmsb(:,i2,ia)
      ib2 = i3(1); ilm2 = i3(2); is2 = i3(3)
      i3 = i2lmsb(:,i1,ia)
      ib1 = i3(1); ilm1 = i3(2); is1 = i3(3)
      vsocz(ib1,ib2,ilm1,ilm2,ia2) = vsocz(ib1,ib2,ilm1,ilm2,ia2) + sum(magdir1(:,ia)*pauli(is2,is1,1:3))*vsocb(i1,i2,ia)
    end do
    end do
  end do
! All done!
  end subroutine build_vsocz
