! ************************************************************************
SUBROUTINE spin_orbit_one_l(lmax,l_s)
! ************************************************************************
!      in this subroutine the matrix L*S is calculated for the basis of
!      real spherical harmonics

!      schematically it has the form
!      (  -L_z    L_+  )
!      (  L_-     L_z  )
IMPLICIT NONE

integer, intent(in)         ::    lmax
double complex, intent(out) ::    L_S((2*lmax+1)*2,(2*lmax+1)*2)

!  local variables 
integer                     ::    i1,i2,i1l 
double complex              ::    icompl
double complex,allocatable  ::    l_min(:,:)
double complex,allocatable  ::    l_up(:,:)
double precision            ::    lfac


icompl=(0D0,1D0)


allocate(l_min(-lmax:lmax,-lmax:lmax))
allocate(l_up(-lmax:lmax,-lmax:lmax))

!  initialize the matrix

DO i1=1,(2*lmax+1)*2
  DO i2=1,(2*lmax+1)*2
    l_s(i2,i1)=0D0
  END DO
END DO

DO i1=-lmax,lmax
  DO i2=-lmax,lmax
    l_min(i2,i1)=0D0
    l_up(i2,i1)=0D0
  END DO
END DO

!  fill the second and the forth quadrant with L_z
! (-L_z,respectively)


DO i1=1,2*lmax+1
  i1l=i1-lmax-1       ! the value of m (varies from -l to +l)
  i2=2*lmax+1-(i1-1)
  
!         L_S(i2,i1)=icompl*i1l
  l_s(i2,i1)=-icompl*i1l
  
END DO

DO i1=2*lmax+2,(2*lmax+1)*2
  i1l=i1-lmax-1-(2*lmax+1)       ! the value of m (varies from -l to +l)
  i2=(2*lmax+1)*2-(i1-(2*lmax+2))
  
!         L_S(i2,i1)=-icompl*i1l
  l_s(i2,i1)=icompl*i1l
  
END DO


!  implement now L_- in the third quadrant

IF (lmax>0) THEN
  
  lfac=SQRT(lmax*(lmax+1D0))/SQRT(2D0)
  l_min(0,-1)=-icompl*lfac
!         l_min(0,-1)=icompl*lfac
  l_min(0,1)=lfac
  l_min(-1,0)=icompl*lfac
  l_min(1,0)=-lfac
  
  IF (lmax > 1) THEN
    
    DO i1=2,lmax
      
      lfac=0.5D0*SQRT(lmax*(lmax+1D0)-i1*(i1-1D0))
      l_min(-i1,-i1+1)=-lfac
      l_min(-i1,i1-1)=icompl*lfac
      l_min(i1,-i1+1)=-icompl*lfac
      l_min(i1,i1-1)=-lfac
      
      lfac=0.5D0*SQRT(lmax*(lmax+1D0)-(i1-1)*(i1))
      l_min(-i1+1,-i1)=lfac
      l_min(-i1+1,i1)=icompl*lfac
      l_min(i1-1,-i1)=-icompl*lfac
      l_min(i1-1,i1)=lfac
      
    END DO
    
  END IF
END IF


DO i1=-lmax,lmax
  DO i2=-lmax,lmax
    l_s(i2+3*lmax+2,i1+lmax+1)=l_min(i1,i2)
  END DO
END DO


!  implement now L_+ in the   quadrant

IF (lmax>0) THEN
  
  lfac=SQRT(lmax*(lmax+1D0))/SQRT(2D0)
  l_up(0,-1)=-icompl*lfac
  l_up(0,1)=-lfac
  l_up(-1,0)=icompl*lfac
  l_up(1,0)=lfac
  
  IF (lmax > 1) THEN
    
    DO i1=2,lmax
      
      lfac=0.5D0*SQRT(lmax*(lmax+1D0)-i1*(i1-1D0))
      l_up(-i1,-i1+1)=lfac
      l_up(-i1,i1-1)=icompl*lfac
      l_up(i1,-i1+1)=-icompl*lfac
      l_up(i1,i1-1)=lfac
      
      lfac=0.5D0*SQRT(lmax*(lmax+1D0)-(i1-1)*(i1))
      l_up(-i1+1,-i1)=-lfac
      l_up(-i1+1,i1)=icompl*lfac
      l_up(i1-1,-i1)=-icompl*lfac
      l_up(i1-1,i1)=-lfac
      
    END DO
    
  END IF
END IF


DO i1=-lmax,lmax
  DO i2=-lmax,lmax
    l_s(i2+lmax+1,i1+3*lmax+2)=l_up(i1,i2)
  END DO
END DO



deallocate(l_min)
deallocate(l_up)


END SUBROUTINE spin_orbit_one_l
