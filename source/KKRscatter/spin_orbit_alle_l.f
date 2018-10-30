C ************************************************************************
      SUBROUTINE SPIN_ORBIT_COMPL(LMAX,LMMAXD,L_S)

      implicit none
C ************************************************************************
c      in this subroutine the matrix L*S is calculated for the basis of
c      real spherical harmonics 

c

       integer, intent(in)     ::     lmax,LMMAXD

       complex, intent(out)    ::     L_S(LMMAXD*2,LMMAXD*2)
c       complex                  ::     L_S(LMMAXD*2,LMMAXD*2)

c  local variables 
       integer                 ::     i1,i2,i1l,RL,lm1,lm2
       complex                 ::     icompl
       complex,allocatable     ::     LS_L(:,:)


       icompl=(0d0,1d0) 

c       write(6,*) "lmax:",lmax
c       write(6,*) "LMMAXD",LMMAXD

c       allocate(LS_L((2*LMAX+1)*2,(2*LMAX+1)*2))
       CALL CINIT((2*LMMAXD)**2,L_S)
c       write(6,*) "after CINIT"

       do RL=0,lmax
c         write(6,*) "rl:",rl

         ALLOCATE(LS_L((2*RL+1)*2,(2*RL+1)*2))
         CALL CINIT(((2*RL+1)*2)**2,LS_L)

         CALL SPIN_ORBIT_ONE_L(RL,LS_L)
c        write(6,*) "SPIN_ORBIT_ONE_L"
c        do lm1=1,(2*RL+1)*2
c          do lm2=1,2*(2*RL+1)
c            write(6,*) lm1,lm2,LS_L(lm1,lm2)
c          end do
c        end do

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
c        write(6,*) "SPIN_ORBIT_GESAMT"
c        do lm1=1,LMMAXD*2
c          do lm2=1,LMMAXD*2
c            write(16+rl,*) lm1,lm2,L_S(lm1,lm2)
c          end do
c        end do


       end do     !rl=0,lmax

c        write(16,*) "SPIN_ORBIT_GESAMT"
c        do lm1=1,LMMAXD*2
c          do lm2=1,LMMAXD*2
c            write(16,*) lm1,lm2,L_S(lm1,lm2)
c          end do
c        end do

c       deallocate(LS_L)
c       write(6,*) "end of spin_orbit alle l"

       END
