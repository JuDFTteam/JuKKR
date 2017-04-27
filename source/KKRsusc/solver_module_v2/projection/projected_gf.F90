  subroutine projected_gf(ie,ia,ja,gf,onsite,struct)
! Block of the projected GF
! This version does not interpolate yet
  use global

  implicit none

! which energy
  integer(kind=i4b), intent(in)  :: ie
! which atoms
  integer(kind=i4b), intent(in)  :: ia, ja
! the block of the GF in the projection basis
  complex(kind=c8b), intent(out) :: gf(nlmsb,nlmsb)
! whether to include the onsite and the structural parts
  logical,           intent(in)  :: onsite, struct
! -------------------------------------------------------------------
  complex(kind=c8b) :: left(nlms), right(nlms)
  complex(kind=c8b) :: gs(nlms,nlms), block(4,4)
! decoding
  integer(kind=i4b) :: j, jlms, l
  integer(kind=i4b) :: i, ilms, k

  gf = 0.d0
! fetch the structural GF block if needed
  if (struct) then
    do jlms=1,nlms
      j = alms2i(jlms,ja)
      do ilms=1,nlms
        i = alms2i(ilms,ia)
        gs(ilms,jlms) = gstruct(i,j,ie)
      end do
    end do
  end if
! put together the projected GF
  do j=1,nlmsba(ja)
    if (struct) left = pzl(j,:,ja,ie)
    do i=1,nlmsba(ia)
!   ****    add the onsite part    ****
      if (ia == ja .and. onsite) then
        gf(i,j) = gf(i,j) + gfpq(i,j,ia,ie)
      end if
!   ****  add the structural part  ****
      if (struct) then
        right = pzr(i,:,ia,ie)
        do l=1,nlms
          do k=1,nlms
            gf(i,j) = gf(i,j) + right(k)*gs(k,l)*left(l)
          end do
        end do
      end if
!   ***********************************
    end do
  end do
! All done!
  end subroutine projected_gf 
