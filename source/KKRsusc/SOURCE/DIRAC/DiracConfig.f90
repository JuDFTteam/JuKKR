module DiracConfig

! For calculations with B=0 and a spherical potential the calculation can be
! made significantly faster by using analytical expressions for the potential
! expansion and setting some (many) coefficients directly to 0 instead of 
! calculating them (and finding out they are 0).
! 0=calculate everything, 1=accelarate by setting the respective coefficients to 0
integer, parameter           :: spherical_only = 0

! Values for the v coefficients will be assumed to be zero if they are smaller
! than v_eps. If no such cutoff shall be done, set v_eps to 0. 
double precision, parameter  :: v_eps = 0
double precision, parameter  :: w_eps = 0

! define whether details of the current calculation should be written on the screen
! verbose = 0: little/no output
! verbose = 1: standard output
! verbose = 2: full output 
integer, parameter           :: verbose = 1
end module DiracConfig