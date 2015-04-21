program kkrcode

  use mod_main0
  use mod_main1a
  use mod_main1b
  use mod_main1c
  use mod_main2
  use mod_types

  implicit none
  
    
  ! KKR program, first do main0, where everything is read in and initialized.
  call main0()
  
  call init_tgmat(t_inc,t_tgmat)
  
  ! Then start scf iterations and do all steps of the KKR formalism until convergence
  do while ( (type0%i_iteration.lt.type0%N_iteration) .and. (type0%N_iteration.ne.-1) )
  
    call main1a()
   
    call main1b()
   
    call main1c()
   
    call main2()
    
  end do
  
  
  
end program
