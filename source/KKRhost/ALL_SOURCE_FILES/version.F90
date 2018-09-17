module mod_version

implicit none
private
public version1, version2, version3, version4

#if defined(compver1)
character(len=*), parameter :: version1=compver1
#else
character(len=*), parameter :: version1='xxx-version-unknown-xxxv2.3-168-g168b19b '
#endif
#if defined(compver2)
character(len=*), parameter :: version2=compver2
#else
character(len=*), parameter :: version2=''
#endif
#if defined(compver3)
character(len=*), parameter :: version3=compver3
#else
character(len=*), parameter :: version3=''
#endif
#if defined(compver4)
character(len=*), parameter :: version4=compver4
#else
character(len=*), parameter :: version4=''
#endif

end module mod_version
