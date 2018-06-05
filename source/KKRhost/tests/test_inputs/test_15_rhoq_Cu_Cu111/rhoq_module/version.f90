module mod_version

implicit none
private
public version

character(len=*), dimension(4), parameter :: version=(/ 'WARNING_no_current_version_info:v1.0-44-ge085edf ', 'mpi ', '-O3 -qopenmp -xhost -g -traceback ', '-mkl=sequential -liomp5 ' /)

end module mod_version
