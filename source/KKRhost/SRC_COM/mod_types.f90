module mod_types

implicit none

      
   type :: type_main0
      
      integer :: i_iteration = -1 
      integer :: N_iteration = -1
      integer :: mit_bry = 1
      
   end type type_main0

!    type :: type_tgmatices
! 
!       integer ::
!    
! 
!    end type type_tgmatices


!    type :: type_lloyd
! 
!       integer ::
!    
! 
!    end type type_lloyd


   type (type_main0), save :: type0
!    type (type_lloyd), save :: t_lloyd
!    type (type_tgmatices), save :: t_tgmat


! contains
! 
!    subroutine bcast_
! 
!    end subroutine


end module