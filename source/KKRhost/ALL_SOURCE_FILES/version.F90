module mod_version

implicit none
private
public version1, version2, version3, version4

#if defined(COMPVER1)
character(len=*), parameter :: version1=COMPVER1
#else
character(len=*), parameter :: version1='XXX-version-unknown-XXXv2.3-168-g168b19b '
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
