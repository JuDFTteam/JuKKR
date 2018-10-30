  subroutine finit_arrays()
! Dellocate the storage allocated in init_arrays
  use global

  implicit none

! -----------------------------------------------------------------------
!                           Energy mesh          
! -----------------------------------------------------------------------
  deallocate(ensusc,esusc,eksusc,desusc,wsusc,escf,ekscf,descf)
! -----------------------------------------------------------------------
!                      Groundstate quantities 
! -----------------------------------------------------------------------
  deallocate(nrpts,nrpts0,nrpts1,rmesh,rsmesh,drmesh,drproj,normesh)
  deallocate(zat,ri,magdir,newdir,iarot,spinproj,vr,br,nrc,mrc,nrv)
  deallocate(old_rho2ns,gs_qlm,gs_mlm,new_rho2ns,rhomat,kxclm,ebandv,etorque,eldau,vshift,nlmpot,i2lmpot,lorb,rho_lm)
! -----------------------------------------------------------------------
!              Storage for t-matrices and structural GF
! -----------------------------------------------------------------------
  deallocate(lm2i,i2lm,lms2i,i2lms,lms2i_new,i2lms_new,alms2i,i2alms)
  deallocate(tmatrix,gstruct)
! -----------------------------------------------------------------------
!                              Basis
! -----------------------------------------------------------------------

  if(allocated(gradnorm)) deallocate(gradnorm)
  if(allocated(gradbasis_lm)) deallocate(gradbasis_lm)
  if(allocated(grad_mass)) deallocate(grad_mass)

! -----------------------------------------------------------------------
!                              Basis
! -----------------------------------------------------------------------

  if(lsusc) then
    if (allocated(i2almsbden)) deallocate(i2almsbden)
    if (allocated(suscbasis)) deallocate(suscbasis)
    if (allocated(iwsusc2)) deallocate(iwsusc2)
    if (allocated(nalmsbden)) deallocate(nalmsbden)
    if (allocated(dengaunt)) deallocate(dengaunt)
    if (allocated(vlmsbden)) deallocate(vlmsbden)
    if (allocated(suscnorm)) deallocate(suscnorm)
    if (allocated(kssusc)) deallocate(kssusc)
    if (allocated(kssusc0)) deallocate(kssusc0)
    if (allocated(kssusc1)) deallocate(kssusc1)
    if (allocated(kssusc2)) deallocate(kssusc2)
    if (allocated(kernel)) deallocate(kernel)
    if (allocated(denominator)) deallocate(denominator)
    if(lcurrcorr) then
      if (allocated(kscurrcorr)) deallocate(kscurrcorr)
      if (allocated(kscurrcorr0)) deallocate(kscurrcorr0)
      if (allocated(kscurrcorr1)) deallocate(kscurrcorr1)
      if (allocated(kscurrcorr2)) deallocate(kscurrcorr2)
    end if
  end if

! -----------------------------------------------------------------------
! All done!
  end subroutine finit_arrays
