module cellproperties
!-------------------------------------------------------------------------------
!> Summary: Reads the 'kkrflex_hoststructure' file 
!> Author:
!> Category: KKRimp, input-output
!>           
!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> Summary: Reads the 'kkrflex_hoststructure' file 
!> Author:
!> Category: KKRimp, input-output
!>           
!-------------------------------------------------------------------------------
subroutine cellproperties_readhoststrucuture
use mod_version_info
implicit none
integer                       :: nbasis
double precision,allocatable  :: rbasis(:,:)
double precision              :: bravais(3,3)



integer                       :: ios,iline,iatom,idim1
character(len=200)            :: string1
open(unit=2382289,file='kkrflex_hoststructure', status='old', iostat=ios)
call version_check_header(2382289)
if (ios/=0) stop '[Error] could not open file kkrflex_hoststructure'
ios=0
do while (ios/=-1)
  iline=iline+1
  read(unit=2382289,fmt='(A)', iostat=ios) string1
  if (ios==-1) cycle
  if (string1(1:1) == "!" .or. string1(1:1) == "#") cycle
  if (len(trim(string1)) == 0) cycle

  if (trim(string1)=='[basis]') then 
    string1 = this_readline(2382289,ios)
    if (ios/=0) stop '[readhoststrucuture] Error'
    read(string1,*) nbasis 
    allocate(rbasis(3,nbasis))
    do iatom = 1 , nbasis
      string1 = this_readline(2382289,ios)
      if (ios/=0) stop '[readhoststrucuture] Error'
      read(string1,*) (rbasis(idim1,iatom),idim1=1,3) 
    end do !nbasis
  elseif (trim(string1)=='[bravias]') then
    do iatom=1,3
      string1 = this_readline(2382289,ios)
      read(string1,*) (bravais(idim1,iatom),idim1=1,3) 
    end do
  end if
end do !ios
end subroutine cellproperties_readhoststrucuture

function this_readline(ifile,ios)
!--------------------------------------------------------
!--  reads the next line in file unit ifile            --
!--------------------------------------------------------
!--  files starting with a dash (#) are treated as     --
!--  a comment !!!                                     --
!--  OUTPUT: next line which is not commented out      --
!--          error variable IOS (should be zero)       --
!--------------------------------------------------------
! input variables
  implicit none
integer,intent(in)               :: ifile
integer,intent(out)              :: ios
! local variables
character(len=200)  ::this_readline
do
  read(unit=ifile,fmt='(A)', iostat=ios) this_readline
  if (ios/=0 .or. this_readline(1:1)/='#') exit
end do
end function this_readline



end module cellproperties

