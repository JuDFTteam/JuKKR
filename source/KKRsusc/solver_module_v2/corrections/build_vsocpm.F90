  subroutine build_vsocpm(ie,vsocb,magdir1,vsocpm)
! computes the non-local SOC contribution to the spin splitting
! ( 0  V+ ) ( Gup  0  ) ( 0  V+ ) = ( V+ Gdn V-      0     )
! ( V-  0 ) (  0  Gdn ) ( V-  0 )   (     0      V- Gup V+ )
  use global

  implicit none

! which energy
  integer(kind=i4b), intent(in)  :: ie
! transverse part of SOC potential
  complex(kind=c8b), intent(in)  :: vsocb(nlmsb,nlmsb,nasusc)
! spin quantization axes
  real(kind=r8b),    intent(in)  :: magdir1(3,nasusc)
! non-local contribution from SOC to spin splitting
  complex(kind=c8b), intent(out) :: vsocpm(nbmax,nbmax,lmmax,lmmax,nasusc2,nasusc2)
! ----------------------------------------------------------------------
  integer(kind=i4b) :: i3(3), ia, ia2, ja, ja2
  integer(kind=i4b) :: i1, ib1, ilm1, is1, i2, ib2, ilm2, is2, j1, jb1, jlm1, js1, j2, jb2, jlm2, js2
  complex(kind=c8b) :: gfij(nlmsb,nlmsb)


! here the GFs don't include the effects of the spin flip part of SOC
  vsocpm = 0.d0
  do ja2=1,nasusc2
    ja = iasusc2(ja2)
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      call projected_gf(ie,ia,ja,gfij,.true.,.true.)
      do j2=1,nlmsba(ja)
      do j1=1,nlmsba(ja)
        i3 = i2lmsb(:,j2,ia)
        jb2 = i3(1); jlm2 = i3(2); js2 = i3(3)
        i3 = i2lmsb(:,j1,ia)
        jb1 = i3(1); jlm1 = i3(2); js1 = i3(3)
        do i2=1,nlmsba(ia)
        do i1=1,nlmsba(ia)
          i3 = i2lmsb(:,i2,ia)
          ib2 = i3(1); ilm2 = i3(2); is2 = i3(3)
          i3 = i2lmsb(:,i1,ia)
          ib1 = i3(1); ilm1 = i3(2); is1 = i3(3)
          vsocpm(ib1,jb2,ilm1,jlm2,ia2,ja2) = vsocpm(ib1,jb2,ilm1,jlm2,ia2,ja2)  &
                                            + sum(magdir1(:,ia)*pauli(js2,is1,1:3))*vsocb(i1,i2,ia)*gfij(i2,j1)*vsocb(j1,j2,ja)
        end do
        end do
      end do
      end do
    end do
  end do
! All done!
  end subroutine build_vsocpm
