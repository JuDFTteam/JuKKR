program FMSM

!###########   INPUT    ###############

! read input files



!###########   SET GEOMERTY    ###############



!###########   Read/Calculate host greenfunction properties   ###############
! T_E, t_ref, G_ref , G_ref T_E G_ref



Do iter = 1,maxiter
  
  call routines
  
END DO !iter = 1,maxiter

!###########   POST ITERATION CALCULATIONS    ###############

call DOS
call Jij


end program FMSM