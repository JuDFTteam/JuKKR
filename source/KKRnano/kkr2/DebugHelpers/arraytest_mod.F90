module arraytest_mod

  interface testarray
    module procedure dtest1d
    module procedure dtest2d
    module procedure dtest3d
    module procedure dtest4d

    module procedure ztest1d
    module procedure ztest2d
    module procedure ztest3d
    module procedure ztest4d

    !module procedure itest1d

    ! repeat until 4d
  end interface testarray

  private :: doubleprectest
  private :: doublecomplextest

  contains

    subroutine dtest4d(msg, array)
      implicit none
      character(len=*), intent(in) :: msg
      double precision, dimension(:,:,:,:), intent(in) :: array

      call doubleprectest(msg, array, size(array))

    end subroutine

    subroutine dtest3d(msg, array)
      implicit none
      character(len=*), intent(in) :: msg
      double precision, dimension(:,:,:), intent(in) :: array

      call doubleprectest(msg, array, size(array))

    end subroutine

    subroutine dtest2d(msg, array)
      implicit none
      character(len=*), intent(in) :: msg
      double precision, dimension(:,:), intent(in) :: array

      call doubleprectest(msg, array, size(array))

    end subroutine

    subroutine dtest1d(msg, array)
      implicit none
      character(len=*), intent(in) :: msg
      double precision, dimension(:), intent(in) :: array

      call doubleprectest(msg, array, size(array))

    end subroutine

    subroutine ztest4d(msg, array)
      implicit none
      character(len=*), intent(in) :: msg
      double complex, dimension(:,:,:,:), intent(in) :: array

      call doublecomplextest(msg, array, size(array))

    end subroutine

    subroutine ztest3d(msg, array)
      implicit none
      character(len=*), intent(in) :: msg
      double complex, dimension(:,:,:), intent(in) :: array

      call doublecomplextest(msg, array, size(array))

    end subroutine

    subroutine ztest2d(msg, array)
      implicit none
      character(len=*), intent(in) :: msg
      double complex, dimension(:,:), intent(in) :: array

      call doublecomplextest(msg, array, size(array))

    end subroutine

    subroutine ztest1d(msg, array)
      implicit none
      character(len=*), intent(in) :: msg
      double complex, dimension(:), intent(in) :: array

      call doublecomplextest(msg, array, size(array))

    end subroutine

   integer function test_getmyrank()
     implicit none
     include 'mpif.h'

     integer :: ierr

     call MPI_COMM_RANK(MPI_COMM_WORLD,test_getmyrank,ierr)
   end function

!=================== Helper routines ==========================================

   subroutine doubleprectest(msg, array, length)
     implicit none
     character(len=*), intent(in) :: msg
     double precision, dimension(*), intent(in) :: array
     integer, intent(in) :: length

     double precision, external :: DNRM2
     double precision :: asum
     integer :: ii

     asum = 0.0d0
     do ii = 1, length
       asum = asum + array(ii)
     end do

     ! print norm and average
     write(*,*) "DEBUG: ", msg, DNRM2(length, array, 1), &
                           asum / length        

   end subroutine

   subroutine doublecomplextest(msg, array, length)
     implicit none
     character(len=*), intent(in) :: msg
     double complex, dimension(*), intent(in) :: array
     integer, intent(in) :: length

     double precision, external :: DZNRM2
     double complex :: asum
     integer :: ii

     asum = dcmplx(0.0d0, 0.0d0)
     do ii = 1, length
       asum = asum + array(ii)
     end do
  
     ! print norm and average
     write(*,*) "DEBUG: ", msg, DZNRM2(length, array, 1), &
                           asum / length        

   end subroutine

end module

! !Example of use:

!#include 'test_macros.h'
!
!program test_it
!  USE_ARRAYTEST_MOD
!  implicit none
!  include 'mpif.h'
!
!  double precision :: darr(100)
!  double complex :: zarr(100)
!  double precision :: dmat(50,50)
!  double complex :: uninit(20,20)
!
!  integer :: ierr
!
!  call MPI_INIT(IERR)
!
!  darr = 1.0d0
!  dmat = 1.0d0
!
!  zarr = dcmplx(0.0d0, 2.0d0)
!
!  TESTARRAY(0, darr)
!  TESTARRAY(1, zarr)
!  TESTARRAY(2, dmat)
!  TESTARRAY(3, uninit)
!
!  call MPI_FINALIZE(IERR)
!
!end program
