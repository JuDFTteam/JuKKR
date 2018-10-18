!-------------------------------------------------------------------------------
!> Summary: Test whether a complex or real variable is NaN
!> Author: 
!-------------------------------------------------------------------------------
!> @note
!> If the variable is real this routine can be replace by the command `ISNAN`
!> @endnote
!-------------------------------------------------------------------------------
module mod_checknan

      interface checknan
      module procedure  checknan_dim0_complex, checknan_dim0_real, &
                        checknan_dim1_complex, checknan_dim1_real, &
                        checknan_dim2_complex, checknan_dim2_real, &
                        checknan_dim3_complex, checknan_dim3_real, &
                        checknan_dim4_complex, checknan_dim4_real
      end interface

contains

!-------------------------------------------------------------------------------
!> Summary: Test whether a complex variable is NaN
!> Author: 
!> Category: KKRimp, numerical-tools
!> Deprecated: False
!-------------------------------------------------------------------------------
subroutine checknan_dim0_complex(array,ierror)
implicit none
!interface
double complex :: array
integer        :: ierror
ierror=0
  IF (dimag(array) .NE. dimag(array)) ierror=1
  IF (dreal(array) .NE. dreal(array)) ierror=1
end subroutine checknan_dim0_complex

!-------------------------------------------------------------------------------
!> Summary: Test whether a real variable is NaN
!> Author: 
!> Category: KKRimp, numerical-tools
!> Deprecated: False
!-------------------------------------------------------------------------------
subroutine checknan_dim0_real(array,ierror)
implicit none
!interface
double precision  :: array
integer           :: ierror
ierror=0
  IF (array .NE. array) ierror=1
end subroutine checknan_dim0_real


!-------------------------------------------------------------------------------
!> Summary: Test whether a complex 1-dim array is NaN
!> Author: 
!> Category: KKRimp, numerical-tools
!> Deprecated: False
!-------------------------------------------------------------------------------
subroutine checknan_dim1_complex(array,ierror)
implicit none
!interface
double complex :: array(:)
integer        :: ierror
!local
integer        :: ival1

ierror=0
do ival1= lbound(array,1),ubound(array,1)
  IF (dimag(array(ival1)) .NE. dimag(array(ival1))) ierror=1
  IF (dreal(array(ival1)) .NE. dreal(array(ival1))) ierror=1
end do
end subroutine checknan_dim1_complex

!-------------------------------------------------------------------------------
!> Summary: Test whether a real 1-dim array is NaN
!> Author: 
!> Category: KKRimp, numerical-tools
!> Deprecated: False
!-------------------------------------------------------------------------------
subroutine checknan_dim1_real(array,ierror)
implicit none
!interface
double precision  :: array(:)
integer           :: ierror
!local
integer           :: ival1

ierror=0
do ival1= lbound(array,1),ubound(array,1)
  IF (array(ival1) .NE. array(ival1)) ierror=1
end do
end subroutine checknan_dim1_real


!-------------------------------------------------------------------------------
!> Summary: Test whether a complex 2-dim array is NaN
!> Author: 
!> Category: KKRimp, numerical-tools
!> Deprecated: False
!-------------------------------------------------------------------------------
subroutine checknan_dim2_complex(array,ierror)
implicit none
!interface
double complex :: array(:,:)
integer        :: ierror
!local
integer        :: ival1,ival2

ierror=0
do ival2= lbound(array,2),ubound(array,2)
  do ival1= lbound(array,1),ubound(array,1)
    IF (dimag(array(ival1,ival2)) .NE. dimag(array(ival1,ival2))) ierror=1
    IF (dreal(array(ival1,ival2)) .NE. dreal(array(ival1,ival2))) ierror=1
  end do
end do
end subroutine checknan_dim2_complex

!-------------------------------------------------------------------------------
!> Summary: Test whether a real 2-dim array is NaN
!> Author: 
!> Category: KKRimp, numerical-tools
!> Deprecated: False
!-------------------------------------------------------------------------------
subroutine checknan_dim2_real(array,ierror)
implicit none
!interface
double precision :: array(:,:)
integer          :: ierror
!local
integer          :: ival1,ival2

ierror=0
do ival2= lbound(array,2),ubound(array,2)
  do ival1= lbound(array,1),ubound(array,1)
    IF ((array(ival1,ival2)) .NE. (array(ival1,ival2))) ierror=1
  end do
end do
end subroutine checknan_dim2_real

!-------------------------------------------------------------------------------
!> Summary: Test whether a complex 3-dim array is NaN
!> Author: 
!> Category: KKRimp, numerical-tools
!> Deprecated: False
!-------------------------------------------------------------------------------
subroutine checknan_dim3_complex(array,ierror)
implicit none
!interface
double complex :: array(:,:,:)
integer        :: ierror
!local
integer        :: ival1,ival2,ival3

ierror=0
do ival3= lbound(array,3),ubound(array,3)
  do ival2= lbound(array,2),ubound(array,2)
    do ival1= lbound(array,1),ubound(array,1)
      IF (dimag(array(ival1,ival2,ival3)) .NE. dimag(array(ival1,ival2,ival3))) ierror=1
      IF (dreal(array(ival1,ival2,ival3)) .NE. dreal(array(ival1,ival2,ival3))) ierror=1
    end do
  end do
end do
end subroutine checknan_dim3_complex

!-------------------------------------------------------------------------------
!> Summary: Test whether a real 3-dim array is NaN
!> Author: 
!> Category: KKRimp, numerical-tools
!> Deprecated: False
!-------------------------------------------------------------------------------
subroutine checknan_dim3_real(array,ierror)
implicit none
!interface
double precision :: array(:,:,:)
integer          :: ierror
!local
integer          :: ival1,ival2,ival3

ierror=0
do ival3= lbound(array,3),ubound(array,3)
  do ival2= lbound(array,2),ubound(array,2)
    do ival1= lbound(array,1),ubound(array,1)
      IF (array(ival1,ival2,ival3) .NE. array(ival1,ival2,ival3)) ierror=1
    end do
  end do
end do
end subroutine checknan_dim3_real

!-------------------------------------------------------------------------------
!> Summary: Test whether a complex 4-dim array is NaN
!> Author: 
!> Category: KKRimp, numerical-tools
!> Deprecated: False
!-------------------------------------------------------------------------------
subroutine checknan_dim4_complex(array,ierror)
implicit none
!interface
double complex :: array(:,:,:,:)
integer        :: ierror
!local
integer        :: ival1,ival2,ival3,ival4

ierror=0
do ival4= lbound(array,4),ubound(array,4)
  do ival3= lbound(array,3),ubound(array,3)
    do ival2= lbound(array,2),ubound(array,2)
      do ival1= lbound(array,1),ubound(array,1)
        IF (dimag(array(ival1,ival2,ival3,ival4)) .NE. dimag(array(ival1,ival2,ival3,ival4))) ierror=1
        IF (dreal(array(ival1,ival2,ival3,ival4)) .NE. dreal(array(ival1,ival2,ival3,ival4))) ierror=1
      end do
    end do
  end do
end do
end subroutine checknan_dim4_complex

!-------------------------------------------------------------------------------
!> Summary: Test whether a real 4-dim array is NaN
!> Author: 
!> Category: KKRimp, numerical-tools
!> Deprecated: False
!-------------------------------------------------------------------------------
subroutine checknan_dim4_real(array,ierror)
implicit none
!interface
double precision :: array(:,:,:,:)
integer        :: ierror
!local
integer        :: ival1,ival2,ival3,ival4

ierror=0
do ival4= lbound(array,4),ubound(array,4)
  do ival3= lbound(array,3),ubound(array,3)
    do ival2= lbound(array,2),ubound(array,2)
      do ival1= lbound(array,1),ubound(array,1)
        IF (array(ival1,ival2,ival3,ival4) .NE. array(ival1,ival2,ival3,ival4)) ierror=1
      end do
    end do
  end do
end do
end subroutine checknan_dim4_real



end module mod_checknan
