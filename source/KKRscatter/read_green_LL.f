      SUBROUTINE READ_GREEN_LL(GLL_EXTRA,LMCL)

c     reads the Green functions for three complex energies and extrapolate to
c     the required Green function for real energies

      implicit none

c   input-output variables      
      integer,intent(in)   ::    LMCL
      complex,intent(out)  ::    GLL_EXTRA(LMCL,LMCL)


c   local variables
      integer              ::    LM,EN,LM1,LM2,status_open
      complex,allocatable  ::    GLL_3(:,:),
     +                           dG_ENERGY(:),d2G_ENERGY(:)
      complex              ::    ENERGY(3),dE1,dE2

       

      ALLOCATE(GLL_3(LMCL**2,3))
      ALLOCATE(dG_ENERGY(LMCL**2))
      ALLOCATE(d2G_ENERGY(LMCL**2))

      open(unit=12, file='GMATLL_GES', form='formatted')

      DO EN=1,3
 
        READ(12,"(2(e17.9,X))") ENERGY(EN)

        DO LM = 1,LMCL**2
          READ(12,"((2I5),(2e17.9))") LM1,LM2,GLL_3(LM,EN)
        END DO
      END DO

      CLOSE(12)

      dE1=ENERGY(2)-ENERGY(1)
      dE2=ENERGY(3)-ENERGY(2)


      DO LM=1,LMCL**2
        dG_ENERGY(LM)=0.5d0/dE1*(GLL_3(LM,3)-GLL_3(LM,1))
        d2G_ENERGY(LM)=(GLL_3(LM,3)-2*GLL_3(LM,2)+GLL_3(LM,1))/dE1**2
      END DO 

      LM=0
      DO LM1=1,LMCL
        DO LM2=1,LMCL
          LM=LM+1
          GLL_EXTRA(LM2,LM1)=GLL_3(LM,2)-(0d0,1d0)*aimag(ENERGY(2))*
     +                                               dG_ENERGY(LM)
     +       +0.5d0*((0d0,1d0)*aimag(ENERGY(2)))**2*d2G_ENERGY(LM)
        END DO
      END DO
      DEALLOCATE(GLL_3)
      DEALLOCATE(dG_ENERGY)
      DEALLOCATE(d2G_ENERGY)

 
      END SUBROUTINE READ_GREEN_LL
