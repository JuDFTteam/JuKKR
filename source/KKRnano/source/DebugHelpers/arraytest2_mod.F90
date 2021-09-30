module arraytest2_mod
  implicit none
  private
  public :: testarray

  interface testarray
    module procedure dtest1d
    module procedure dtest2d
    module procedure dtest3d
    module procedure dtest4d

    module procedure ztest1d
    module procedure ztest2d
    module procedure ztest3d
    module procedure ztest4d
  endinterface

  contains

    character(len=80) function dtest4d(nr, msg, array)
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double precision, intent(in) :: array(:,:,:,:)

      dtest4d = doubleprectest(nr, msg, array, size(array))
    endfunction

    character(len=80) function dtest3d(nr, msg, array)
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double precision, intent(in) :: array(:,:,:)

      dtest3d = doubleprectest(nr, msg, array, size(array))
    endfunction

    character(len=80) function dtest2d(nr, msg, array)
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double precision, intent(in) :: array(:,:)

      dtest2d = doubleprectest(nr, msg, array, size(array))
    endfunction

    character(len=80) function dtest1d(nr, msg, array)
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double precision, intent(in) :: array(:)

      dtest1d = doubleprectest(nr, msg, array, size(array))
    endfunction

    character(len=80) function ztest4d(nr, msg, array)
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double complex, intent(in) :: array(:,:,:,:)

      ztest4d = doublecomplextest(nr, msg, array, size(array))
    endfunction

    character(len=80) function ztest3d(nr, msg, array)
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double complex, intent(in) :: array(:,:,:)

      ztest3d = doublecomplextest(nr, msg, array, size(array))
    endfunction

    character(len=80) function ztest2d(nr, msg, array)
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double complex, intent(in) :: array(:,:)

      ztest2d = doublecomplextest(nr, msg, array, size(array))
    endfunction

    character(len=80) function ztest1d(nr, msg, array)
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double complex, intent(in) :: array(:)

      ztest1d = doublecomplextest(nr, msg, array, size(array))
    endfunction

!=================== Helper routines ==========================================

   character(len=80) function doubleprectest(nr, msg, array, length) result(str)
     integer, intent(in) :: nr
     character(len=*), intent(in) :: msg
     double precision, intent(in) :: array(*)
     integer, intent(in) :: length

     double precision, external :: DNRM2 ! from LAPACK

     ! print norm and average
     write(unit=str, fmt='(a7,i4,x,a16,x,e16.9,x,e16.9)') &
       "DEBUG: ", nr, msg, DNRM2(length, array, 1), sum(array(1:length))/length
   endfunction

   character(len=80) function doublecomplextest(nr, msg, array, length) result(str)
     integer, intent(in) :: nr
     character(len=*), intent(in) :: msg
     double complex, intent(in) :: array(*)
     integer, intent(in) :: length

     double precision, external :: DZNRM2 ! from LAPACK
 
     ! print norm and average
     write(unit=str, fmt='(a7,i4,x,a16,x,e12.5,x,e12.5,x,e12.5)') &
       "DEBUG: ", nr, msg, DZNRM2(length, array, 1), sum(array(1:length))/length
   endfunction

endmodule
