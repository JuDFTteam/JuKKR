  subroutine finit_gfpartsfit()
! Allocates arrays for different parts of projected GF
  use global

  implicit none

! -----------------------------------------------------------------------
!                             Pointers               
! -----------------------------------------------------------------------
  deallocate(lmsb2i,i2lmsb,nlmsba)
! -----------------------------------------------------------------------
!   Storage for coefficients in full projection basis
! -----------------------------------------------------------------------
  deallocate(overlap,pzl,pzr,gfpq)
  if (isra == 1) deallocate(fzl,fzr,gfps,gffq,gffs)
  deallocate(vlmsbgf)
  if (lfit .and. ifit == 2)  deallocate(gffit)
! -----------------------------------------------------------------------
!                           overlaps_gf              
! -----------------------------------------------------------------------
  deallocate(nalmsbgf,i2almsbgf,almsbgf2i,nalmsbgrad,i2almsbgrad,almsbgrad2i)
! All done!
  end subroutine finit_gfpartsfit
