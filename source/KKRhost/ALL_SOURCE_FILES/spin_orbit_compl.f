C ************************************************************************
      SUBROUTINE SPIN_ORBIT_COMPL(LMAX,LMMAXD,L_S)

      implicit none
C ************************************************************************
c      in this subroutine the matrix L*S is calculated for the basis of
c      real spherical harmonics 

c

       integer, intent(in)     ::     lmax,LMMAXD
       double complex, intent(out)    ::     L_S(LMMAXD*2,LMMAXD*2)

c  local variables 
       integer                 ::     i1,i2,i1l,RL,lm1,lm2
       double complex                 ::     icompl
       double complex,allocatable     ::     LS_L(:,:)


       icompl=(0d0,1d0) 


       CALL CINIT((2*LMMAXD)**2,L_S)

       do RL=0,lmax

         ALLOCATE(LS_L((2*RL+1)*2,(2*RL+1)*2))
         CALL CINIT(((2*RL+1)*2)**2,LS_L)


         CALL SPIN_ORBIT_ONE_L(RL,LS_L)

         do lm1=1,(2*RL+1)*2

           IF (lm1 <= 2*RL+1 ) THEN
             do lm2=1,(2*RL+1)
              L_S(RL**2+Lm1,RL**2+lm2)=0.5d0*LS_L(lm1,lm2)
             end do 
             do lm2=(2*RL+1)+1,(2*RL+1)*2
               L_S(RL**2+Lm1,LMMAXD+RL**2-(2*RL+1)+lm2)=
     +                                     0.5d0*LS_L(lm1,lm2)
             end do 
          ELSE
             do lm2=1,(2*RL+1)
               L_S(LMMAXD+RL**2-(2*RL+1)+Lm1,RL**2+lm2)=
     +                                     0.5d0*LS_L(lm1,lm2)
             end do 
             do lm2=(2*RL+1)+1,(2*RL+1)*2
               L_S(LMMAXD+RL**2-(2*RL+1)+Lm1,LMMAXD+RL**2-(2*RL+1)+lm2)=
     +                                     0.5d0*LS_L(lm1,lm2)
             end do
           END IF

         end do    !lm1

         DEALLOCATE(LS_L)


       end do     !rl=0,lmax


       END
