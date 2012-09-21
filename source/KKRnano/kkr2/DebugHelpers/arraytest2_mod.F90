module arraytest2_mod

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

    character(len=80) function dtest4d(nr, msg, array)
      implicit none
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double precision, dimension(:,:,:,:), intent(in) :: array

      dtest4d = doubleprectest(nr, msg, array, size(array))

    end function

    character(len=80) function dtest3d(nr, msg, array)
      implicit none
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double precision, dimension(:,:,:), intent(in) :: array

      dtest3d = doubleprectest(nr, msg, array, size(array))

    end function

    character(len=80) function dtest2d(nr, msg, array)
      implicit none
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double precision, dimension(:,:), intent(in) :: array

      dtest2d = doubleprectest(nr, msg, array, size(array))

    end function

    character(len=80) function dtest1d(nr, msg, array)
      implicit none
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double precision, dimension(:), intent(in) :: array

      dtest1d = doubleprectest(nr, msg, array, size(array))

    end function

    character(len=80) function ztest4d(nr, msg, array)
      implicit none
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double complex, dimension(:,:,:,:), intent(in) :: array

      ztest4d = doublecomplextest(nr, msg, array, size(array))

    end function

    character(len=80) function ztest3d(nr, msg, array)
      implicit none
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double complex, dimension(:,:,:), intent(in) :: array

      ztest3d = doublecomplextest(nr, msg, array, size(array))

    end function

    character(len=80) function ztest2d(nr, msg, array)
      implicit none
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double complex, dimension(:,:), intent(in) :: array

      ztest2d = doublecomplextest(nr, msg, array, size(array))

    end function

    character(len=80) function ztest1d(nr, msg, array)
      implicit none
      integer, intent(in) :: nr
      character(len=*), intent(in) :: msg
      double complex, dimension(:), intent(in) :: array

      ztest1d = doublecomplextest(nr, msg, array, size(array))

    end function

!=================== Helper routines ==========================================

   character(len=80) function doubleprectest(nr, msg, array, length)
     implicit none
     integer, intent(in) :: nr
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
     write(doubleprectest,'(A7,I4,X,A16,X,E16.9,X,E16.9)') "DEBUG: ", nr, &
                              msg, DNRM2(length, array, 1), asum / length
     !write (doubleprectest,*) "Just a test"

   end function

   character(len=80) function doublecomplextest(nr, msg, array, length)
     implicit none
     integer, intent(in) :: nr
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
     write(doublecomplextest,'(A7,I4,X,A16,X,E12.5,X,E12.5,X,E12.5)') &
            "DEBUG: ", nr, msg, DZNRM2(length, array, 1), asum / length

     !write(doublecomplextest,*) "Just a test"

   end function

end module

