module Statistics_mod
  implicit none
  private
  
  public :: eval, add, allreduce

  interface add
    module procedure reduce_i8, reduce_i4, reduce_r8, reduce_r4
  endinterface

  interface allreduce
    module procedure allreduce_i8, allreduce_r8
  endinterface
  
  interface eval
    module procedure eval_int, eval_flt
  endinterface

  contains

  subroutine reduce_i8(x, xmom, xinf)
    integer(kind=8), intent(in) :: x
    integer(kind=8), intent(inout) :: xmom(0:) ! [size(x), sum(x), sum(x*x), ...]
    integer(kind=8), intent(inout) :: xinf(0:1) ! [maxval(x), -minval(x)]

    xmom = xmom + [1_8, x, x*x, x*x*x]
    xinf = max(xinf, [x, -x])
  endsubroutine ! reduce

  subroutine reduce_i4(x, xmom, xinf)
    integer(kind=4), intent(in) :: x
    integer(kind=8), intent(inout) :: xmom(0:), xinf(0:1)
    call reduce_i8(int(x, kind=8), xmom, xinf)
  endsubroutine ! reduce

  subroutine reduce_r8(x, xmom, xinf)
    real(kind=8), intent(in) :: x
    real(kind=8), intent(inout) :: xmom(0:) ! [size(x), sum(x), sum(x*x), ...]
    real(kind=8), intent(inout) :: xinf(0:1) ! [maxval(x), -minval(x)]

    xmom = xmom + [1.d0, x, x*x, x*x*x]
    xinf = max(xinf, [x, -x])
  endsubroutine ! reduce

  subroutine reduce_r4(x, xmom, xinf)
    real(kind=4), intent(in) :: x
    real(kind=8), intent(inout) :: xmom(0:), xinf(0:1)
    call reduce_r8(real(x, kind=8), xmom, xinf)
  endsubroutine ! reduce
  
  
  integer function allreduce_i8(xmom, xinf, communicator) result(ierr)
    include 'mpif.h'
    integer(kind=8), intent(inout) :: xmom(0:,:) ! [size(x), sum(x), sum(x*x), ...]
    integer(kind=8), intent(inout) :: xinf(0:,:) ! [maxval(x), -minval(x)]
    integer, intent(in) :: communicator   
    integer(kind=8), allocatable :: tmom(:,:), tinf(:,:)
    
    allocate(tmom(size(xmom, 1),size(xmom,2)), tinf(size(xinf, 1),size(xinf,2))) 
    tmom = xmom ; tinf = xinf ! copy 
    
    call MPI_Allreduce(xmom, tmom, size(xmom), MPI_INTEGER8, MPI_SUM, communicator, ierr)
    call MPI_Allreduce(xinf, tinf, size(xinf), MPI_INTEGER8, MPI_MAX, communicator, ierr)

    deallocate(tmom, tinf)
  endfunction ! allreduce
  
  integer function allreduce_r8(xmom, xinf, communicator) result(ierr)
    include 'mpif.h'
    real(kind=8), intent(inout) :: xmom(0:,:) ! [size(x), sum(x), sum(x*x), ...]
    real(kind=8), intent(inout) :: xinf(0:,:) ! [maxval(x), -minval(x)]
    integer, intent(in) :: communicator   
    real(kind=8), allocatable :: tmom(:,:), tinf(:,:)
    
    allocate(tmom(size(xmom, 1),size(xmom,2)), tinf(size(xinf, 1),size(xinf,2))) 
    tmom = xmom ; tinf = xinf ! copy 
    
    call MPI_Allreduce(xmom, tmom, size(xmom), MPI_REAL8, MPI_SUM, communicator, ierr)
    call MPI_Allreduce(xinf, tinf, size(xinf), MPI_REAL8, MPI_MAX, communicator, ierr)

    deallocate(tmom, tinf)
  endfunction ! allreduce
  
  character(len=64) function eval_int(xmom, xinf, ndigits) result(str)
    integer(kind=8), intent(in) :: xmom(0:) ! [size(x), sum(x), sum(x*x), ...]
    integer(kind=8), intent(in) :: xinf(0:1) ! [maxval(x), -minval(x)]
    integer, intent(in), optional :: ndigits ! digits behind the floating point
    
    integer :: ios, ndig, minv, maxv
    character(len=32) :: frmt
    double precision :: invN, mean, var, dev
    
    invN = 1.d0/dble(max(1, xmom(0)))
    mean = xmom(1)*invN ! divide by the number of samples
    
    ndig = 1 ; if (present(ndigits)) ndig = max(0, ndigits)
    frmt = '(9(i0,a))' ; if (ndig > 0) write(unit=frmt, fmt="(9(a,i0))", iostat=ios) "(2(f0.",ndig,",a),9(i0,a))"
    
    var  = 0 ; if (ubound(xmom, 1) > 1) &
    var  = xmom(2)*invN - mean*mean ! variance
    dev  = sqrt(max(0.d0, var)) ! compute sigma as sqrt(variance)
    minv = -xinf(1)
    maxv =  xinf(0)
    write(unit=str, fmt=frmt, iostat=ios) mean," +/-",dev,"  [",minv,", ",maxv,"]"
  endfunction ! eval

  character(len=96) function eval_flt(xmom, xinf, ndigits) result(str)
    real(kind=8), intent(in) :: xmom(0:) ! [size(x), sum(x), sum(x*x), ...]
    real(kind=8), intent(in) :: xinf(0:1) ! [maxval(x), -minval(x)]
    integer, intent(in), optional :: ndigits ! significant digits in total
    
    integer :: ios, ndig, nd10
    character(len=32) :: frmt
    real(kind=8) :: invN, mean, var, dev, minv, maxv
    
    invN = 1.d0/dble(max(1.d0, xmom(0)))
    mean = xmom(1)*invN ! divide by the number of samples
    
    ndig = 1 ; if (present(ndigits)) ndig = max(1, ndigits)
    nd10 = ceiling(log10(mean))
    frmt = '(9(f0.1,a))' !; if (mean > 0) write(unit=frmt, fmt="(9(a,i0))", iostat=ios) "(2(f0.",ndig,",a),9(i0,a))" ! ToDo

    var  = 0 ; if (ubound(xmom, 1) > 1) &
    var  = xmom(2)*invN - mean*mean ! variance
    dev  = sqrt(max(0.d0, var)) ! compute sigma as sqrt(variance)
    minv = -xinf(1)
    maxv =  xinf(0)
    write(unit=str, fmt=frmt, iostat=ios) mean," +/-",dev,"  [",minv,", ",maxv,"]"
  endfunction ! eval

endmodule ! Statistics_mod
