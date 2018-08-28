module mod_version

implicit none
private
public version

character(len=*), dimension(4), parameter :: version=(/ 'v2.3-61-g4b60676 ', 'mpi ', '-O3 -openmp -xhost -g -traceback ', '-mkl=sequential -liomp5 ' /)

end module mod_version
