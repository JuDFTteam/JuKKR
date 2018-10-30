  subroutine save_tmati(tmat,is,ie,ia,lmsize)
! Puts the t-matrices in RAM
! Collinear spins here
! lmsize = 2*lmmax if new solver is used
  use global, only: i4b, r8b, c8b, nesusc, nasusc, nlms, lmmax, lms2i, i2lms_new, tmatrix

  implicit none

  complex(kind=c8b), intent(in) :: tmat(lmsize,lmsize)
  integer(kind=i4b), intent(in) :: is, ie, ia, lmsize
! -----------------------------------------------------------------
  integer(kind=i4b) :: ilm, jlm, i, j
  integer(kind=i4b) :: ilmn, jlmn, isn ! New indices
   
! Check wich solver is used ?
  if (lmsize == 2*lmmax) then
    tmatrix(:,:,ia,ie) = 0.d0
    do jlmn=1, lmsize
      jlm  = i2lms_new(1,jlmn)
      isn  = i2lms_new(2,jlmn) 
      j = lms2i(jlm,isn) 
      do ilmn=1, lmsize
        ilm  = i2lms_new(1,ilmn)
        isn  = i2lms_new(2,ilmn) 
        i    = lms2i(ilm,isn)
        tmatrix(i,j,ia,ie) = tmat(ilmn,jlmn)
!       tmatrix(i,j,ia,ie) = 0.5d0*(tmat(ilmn,jlmn) + tmat(jlmn,ilmn))
!       write(iodb,'(4i4,6es12.4)') i2alms(:,i), i2alms(:,j), tmat(ilmn,jlmn)
      end do
    end do 
  else ! Old solver
    if (is == 1) tmatrix(:,:,ia,ie) = 0.d0
    do jlm=1,lmsize
      j = lms2i(jlm,is)
      do ilm=1,lmsize
        i = lms2i(ilm,is)
        tmatrix(i,j,ia,ie) = tmat(ilm,jlm)
!       tmatrix(i,j,ia,ie) = 0.5d0*(tmat(ilm,jlm) + tmat(jlm,ilm))
!       write(iodb,'(4i4,6es12.4)') i2alms(:,i), i2alms(:,j), tmat(ilm,jlm)
      end do
    end do
  endif ! lmsize
! All done
  end subroutine save_tmati
