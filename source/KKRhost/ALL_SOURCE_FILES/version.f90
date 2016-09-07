module mod_version

implicit none
private
public version

character(len=*), dimension(4), parameter :: version=(/ 'v2.0-38-g6593f48 ', 'mpi ', '-O2 -r8 -traceback -module ./OBJ ', '-L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core ' /)

end module mod_version
