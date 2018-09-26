!------------------------------------------------------------------------------------
!> Summary: Short explanation module
!> Author: People who wrote it
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!>
!> One can write Latex comments like this \(i\hbar\frac{\partial \psi}{\partial t}=-\mathcal{H}\psi\)
!> or add labeled equations using the standard latex way
!> \begin{equation}
!> \mathbf{A} \mathbf{v}= \eta\mathbf{v}
!> \end{equation}
!> **FORd** also accepts markdown style so you can _write with style_ 
!> 
!> **IMPORTANT**
!> The JM-KKR follows the coding conventions noted in this example, one that is
!> not obvious is that **each level of indentation consists of two spaces**. Please keep this:
!> _do it for the children_.
!> So please keep the conventions.
! These are special boxes for ford, notice that this comment does not appear in the html file.
! These boxes contain important information and should be added when necessary. ALWAYS remember to close the box
! BEFORE oppening a new one or they will be nested.
!------------------------------------------------------------------------------------
!> @note Notes on the code
!> @endnote
!> @todo things that must be checked
!> @endtodo
!> @warning Important precautions
!> @endwarning
!> @bug If nasty things are found
!> @endbug
!------------------------------------------------------------------------------------
module test_docu

   implicit none

   integer :: some_global !! Explanation of the variable
   real(kind=dp), dimension(:,:), allocatable :: some_global_array !! Explanation of the array

contains

   !-------------------------------------------------------------------------------
   !> Summary: Short explanation of the subroutine
   !> Author: Who wrote this subroutine
   !> Category: TAGS for the code they must be written as TAG1, TAG2, ..., TAGN
   !> Deprecated: False ! This needs to be set to True for deprecated subroutines
   !> A More detailed explanation with the math, concepts, etc necessary to understand the routine
   !-------------------------------------------------------------------------------
   !> @note Notes on the code
   !> @endnote
   !> @todo things that must be checked
   !> @endtodo
   !> @warning Important precautions
   !> @endwarning
   !> @bug If nasty things are found
   !> @endbug
   !-------------------------------------------------------------------------------
   subroutine name_sub(dummy_in,dummy_out)

      implicit none
      integer, intent(in) :: dummy_in !! Short description of variables
      real(kind=dp), dimension(:,:), intent(out) :: dummy_out !! Short description of variables

      integer :: ii, N ! If one does not wish to track them, normal comments still work

      do ii=1, N
         dummy_out(ii,1)=dummy_in+1
         !! One can also have comments in the body of the code, best avoided unless necessary
      enddo

   end subroutine name_sub

   !-------------------------------------------------------------------------------
   !> Summary: Short explanation of the function
   !> Author: Who wrote this function
   !> Category: TAGS for the code they must be written as TAG1, TAG2, ..., TAGN
   !> Deprecated: False ! This needs to be set to True for deprecated functions
   !> A More detailed explanation with the math, concepts, etc necessary to understand the function
   !-------------------------------------------------------------------------------
   !> @note Notes on the code
   !> @endnote
   !> @todo things that must be checked
   !> @endtodo
   !> @warning Important precautions
   !> @endwarning
   !> @bug If nasty things are found
   !> @endbug
   !-------------------------------------------------------------------------------
   real(dp) function name_fun(in_fun1,in_fun2) result(output)

      implicit none

      integer, intent(in) :: in_fun1 !! Short description of variable
      integer, intent(in) :: in_fun2 !! Short description of variable

   end function name_fun

end module test_docu
