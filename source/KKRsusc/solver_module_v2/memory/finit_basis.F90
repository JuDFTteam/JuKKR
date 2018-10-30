  subroutine finit_basis()
! Dellocate the storage allocated in init_basis
  use global

  implicit none

! -----------------------------------------------------------------------
!                  Basis and projection coefficients
! -----------------------------------------------------------------------
  deallocate(phiref)
  deallocate(pzc,pqc,noregcoeffs,noirrcoeffs)
  deallocate(noregcoeffs_soc,noirrcoeffs_soc)
  if (isra == 1) deallocate(fzc,psc,fqc,fsc)
  if (lsusc) deallocate(mtotsusc,mxcsusc,msocsusc)
!  deallocate(torque)
! -----------------------------------------------------------------------
! All done!
  end subroutine finit_basis
