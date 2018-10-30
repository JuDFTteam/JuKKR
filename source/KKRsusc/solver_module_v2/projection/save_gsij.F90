  subroutine save_gsij(gmat,is,ie,ia,ja,lmsize)
! Puts the blocks of the structural GF in RAM
! gmat is used for temporary storage of free space GF => can be destroyed
! Collinear spins here
! lmsize = 2*lmmax if new solver is used
  use global, only: i4b, r8b, c8b, lfreegf, ri, nesusc, eksusc, nalms, lmmax, lms2i, alms2i, i2lms_new, gstruct 

  implicit none

  complex(kind=c8b), intent(inout) :: gmat(lmsize,lmsize)
  integer(kind=i4b), intent(in)    :: is, ie, ia, ja, lmsize
! -----------------------------------------------------------------
  integer(kind=i4b) :: ilm, jlm, ilms, jlms, i, j  
  integer(kind=i4b) :: ilmn, jlmn, isn ! New indices
 
! Check which solver is used
  if (lmsize == 2*lmmax) then 
    if (ia == 1 .and. ja == 1) gstruct(:,:,ie) = 0.d0
    do jlmn=1, lmsize
      jlm  = i2lms_new(1,jlmn)
      isn  = i2lms_new(2,jlmn) 
      jlms = lms2i(jlm,isn)
      j    = alms2i(jlms,ja)
      do ilmn=1, lmsize
        ilm  = i2lms_new(1,ilmn)
        isn  = i2lms_new(2,ilmn) 
        ilms = lms2i(ilm,isn)
        i    = alms2i(ilms,ia)
        gstruct(i,j,ie) = gmat(ilmn,jlmn)
!       gstruct(i,j,ie) = gstruct(i,j,ie) + 0.5d0*gmat(ilm,jlm)
!       gstruct(j,i,ie) = gstruct(j,i,ie) + 0.5d0*gmat(ilm,jlm)
!       write(iodb,'(4i4,6es12.4)') i2alms(:,i), i2alms(:,j), gmat(ilm,jlm)
      enddo
    enddo  
  else ! Old_solver
    if (is == 1 .and. ia == 1 .and. ja == 1) gstruct(:,:,ie) = 0.d0
!   Structural GF replaced with free space GF
    if (lfreegf) then
      call free_gf(eksusc(ie),ri(:,ja)-ri(:,ia),gmat)
    end if
!   Structural GF placed in memory
    do jlm=1,lmsize
      jlms = lms2i(jlm,is)
      j = alms2i(jlms,ja)
      do ilm=1,lmsize
        ilms = lms2i(ilm,is)
        i = alms2i(ilms,ia)
        gstruct(i,j,ie) = gmat(ilm,jlm)
!       gstruct(i,j,ie) = gstruct(i,j,ie) + 0.5d0*gmat(ilm,jlm)
!       gstruct(j,i,ie) = gstruct(j,i,ie) + 0.5d0*gmat(ilm,jlm)
!       write(iodb,'(4i4,6es12.4)') i2alms(:,i), i2alms(:,j), gmat(ilm,jlm)
      end do
    end do
  end if ! lmsize
! All done
  end subroutine save_gsij
