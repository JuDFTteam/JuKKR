!------------------------------------------------------------------------------------
!> Summary: Definitions of the compilation versions
!> Author: 
!> Category: version-control, deprecated, KKRhost
!> Deprecated: True
!> Definitions of the compilation versions 
!------------------------------------------------------------------------------------
!> @note Jonathan Chico: Maybe this is deprecated?
!------------------------------------------------------------------------------------
module mod_version

implicit none
private
public version1, version2, version3, version4

#if defined(COMPVER1)
character(len=*), parameter :: version1=COMPVER1
#else
character(len=*), parameter :: version1='xxx-version-unknown-xxx'
#endif
#if defined(COMPVER2)
character(len=*), parameter :: version2=COMPVER2
#else
character(len=*), parameter :: version2=''
#endif
#if defined(COMPVER3)
character(len=*), parameter :: version3=COMPVER3
#else
character(len=*), parameter :: version3=''
#endif
#if defined(COMPVER4)
character(len=*), parameter :: version4=COMPVER4
#else
character(len=*), parameter :: version4=''
#endif

end module mod_version
