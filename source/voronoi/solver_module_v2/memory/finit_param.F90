  subroutine finit_param()
! Deallocates arrays allocated by init_param
  use global

  implicit none

  deallocate(iasusc,iasusc2,igroup,isoc,socscaling,ildau,ueff,jeff,ibfield,blen,bdir,bconlen,bcondir,ijij)
  deallocate(inobxc,iadmat,ildmat,isusc,ikha,ikxc,iwsusc,issusc,ewsusc,nowfns,nobasis,ircutat,npanat)
! Flag these things as uninitialized
  noparameters = .true.
! All done!
  end subroutine finit_param
