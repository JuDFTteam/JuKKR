  subroutine save_tmcoll(tmat,is,ie)
! Puts the blocks of the t-matrix in RAM
! Collinear spins here
  use global

  implicit none

  complex(kind=c8b), intent(in) :: tmat(lmmax,lmmax,nasusc)
  integer(kind=i4b), intent(in) :: is, ie
! -----------------------------------------------------------------
  integer(kind=i4b) :: ia, ilm, jlm, i, j

  if (is == 1) tmatrix(:,:,ia,ie) = 0.d0
  do ia=1,nasusc
    do jlm=1,lmmax
      j = lms2i(jlm,is)
      do ilm=1,lmmax
        i = lms2i(ilm,is)
        tmatrix(i,j,ia,ie) = tmat(ilm,jlm,ia)
!        tmatrix(i,j,ia,ie) = 0.5d0*(tmat(ilm,jlm,ia) + tmat(jlm,ilm,ia))
!        write(iodb,'(4i4,6es12.4)') i2lms(:,i), i2lms(:,j), tmat(ilm,jlm,ia)
      end do
    end do
  end do
! All done
  end subroutine save_tmcoll
