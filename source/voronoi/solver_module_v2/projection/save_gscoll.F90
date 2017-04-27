  subroutine save_gscoll(gmat,is,ie)
! Puts the blocks of the structural GF in RAM
! Collinear spins here
  use global

  implicit none

  complex(kind=c8b), intent(in) :: gmat(lmmax,lmmax,nasusc,nasusc)
  integer(kind=i4b), intent(in) :: is, ie
! -----------------------------------------------------------------
  integer(kind=i4b) :: ja, ia, ilm, jlm, i, j

  if (is == 1 .and. ia == 1 .and. ja == 1) gstruct(:,:,ie) = 0.d0
  do ja=1,nasusc
  do jlm=1,lmmax
    j = lms2i(jlm,is)
    j = alms2i(j,ja)
    do ia=1,nasusc
    do ilm=1,lmmax
      i = lms2i(ilm,is)
      i = alms2i(i,ia)
      gstruct(i,j,ie) = gmat(ilm,jlm,ia,ja)
!      write(iodb,'(4i4,6es12.4)') i2alms(:,i), i2alms(:,j), gmat(ilm,jlm,ia,ja)
    end do
    end do
  end do
  end do
! Make it structurally symmetric: Gij(L1,s1,L2,s2) = Gji(L2,s2,L1,s1)
  return
  if (is == nsmax) then
  do j=1,nalms
    do i=j+1,nalms
      gstruct(i,j,ie) = 0.5d0*(gstruct(i,j,ie) + gstruct(j,i,ie))
      gstruct(j,i,ie) = gstruct(i,j,ie)
    end do
  end do
  end if
! All done
  end subroutine save_gscoll
