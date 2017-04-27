  subroutine build_vsocpm2(ie,vsocb,magdir1,vsocpm)
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
  complex(kind=c8b), intent(out) :: vsocpm(nbmax,nbmax,lmmax,lmmax,nasusc,nasusc)
! ----------------------------------------------------------------------
  integer(kind=i4b) :: i3(3)
  integer(kind=i4b) :: ia, ia2, i1, ib1, ilm1, is1, i2, ib2, ilm2, is2
  integer(kind=i4b) :: ja, ja2, j1, jb1, jlm1, js1, j2, jb2, jlm2, js2
  complex(kind=c8b) :: gfij(nlmsb,nlmsb)


! here the GFs don't include the effects of the spin flip part of SOC
  vsocpm = 0.d0
  do ja=1,nasusc
!    ja = iasusc2(ja2)
    if (isoc(ja) /= 0 .or. ibfield(ja) /= 0) then
      do ia=1,nasusc
!        ia = iasusc2(ia2)
        if (isoc(ia) /= 0 .or. ibfield(ia) /= 0) then
          call projected_gf(ie,ia,ja,gfij,.true.,.true.)
          do j2=1,nlmsba(ja)
          do j1=1,nlmsba(ja)
            i3 = i2lmsb(:,j2,ja)
            jb2 = i3(1); jlm2 = i3(2); js2 = i3(3)
            i3 = i2lmsb(:,j1,ja)
            jb1 = i3(1); jlm1 = i3(2); js1 = i3(3)
            do i2=1,nlmsba(ia)
            do i1=1,nlmsba(ia)
              i3 = i2lmsb(:,i2,ia)
              ib2 = i3(1); ilm2 = i3(2); is2 = i3(3)
              i3 = i2lmsb(:,i1,ia)
              ib1 = i3(1); ilm1 = i3(2); is1 = i3(3)
              vsocpm(ib1,jb2,ilm1,jlm2,ia,ja) = vsocpm(ib1,jb2,ilm1,jlm2,ia,ja)  &
                                              + sum(magdir1(:,ia)*pauli(js2,is1,1:3))*vsocb(i1,i2,ia)*gfij(i2,j1)*vsocb(j1,j2,ja)
            end do
            end do
          end do
          end do
        end if
      end do
    end if
  end do
! All done!
  end subroutine build_vsocpm2
