! ************************************************************************
SUBROUTINE spin_orbit_compl(lmax,lmmaxd,l_s)
! ************************************************************************
!      in this subroutine the matrix L*S is calculated for the basis of
!      real spherical harmonics

IMPLICIT NONE

integer, intent(in)     ::     lmax,LMMAXD
double complex, intent(out)    ::     L_S(LMMAXD*2,LMMAXD*2)

! local variables 
integer                 ::     RL,lm1,lm2
double complex                 ::     icompl
double complex,allocatable     ::     LS_L(:,:)

icompl=(0D0,1D0)


CALL cinit((2*lmmaxd)**2,l_s)

DO rl=0,lmax
  
  allocate(ls_l((2*rl+1)*2,(2*rl+1)*2))
  CALL cinit(((2*rl+1)*2)**2,ls_l)
  
  
  CALL spin_orbit_one_l(rl,ls_l)
  
  DO lm1=1,(2*rl+1)*2
    
    IF (lm1 <= 2*rl+1 ) THEN
      DO lm2=1,(2*rl+1)
        l_s(rl**2+lm1,rl**2+lm2)=0.5D0*ls_l(lm1,lm2)
      END DO
      DO lm2=(2*rl+1)+1,(2*rl+1)*2
        l_s(rl**2+lm1,lmmaxd+rl**2-(2*rl+1)+lm2)= 0.5D0*ls_l(lm1,lm2)
      END DO
    ELSE
      DO lm2=1,(2*rl+1)
        l_s(lmmaxd+rl**2-(2*rl+1)+lm1,rl**2+lm2)= 0.5D0*ls_l(lm1,lm2)
      END DO
      DO lm2=(2*rl+1)+1,(2*rl+1)*2
        l_s(lmmaxd+rl**2-(2*rl+1)+lm1,lmmaxd+rl**2-(2*rl+1)+lm2)=  &
            0.5D0*ls_l(lm1,lm2)
      END DO
    END IF
    
  END DO    !lm1
  
  deallocate(ls_l)
  
  
END DO     !rl=0,lmax


END SUBROUTINE spin_orbit_compl
