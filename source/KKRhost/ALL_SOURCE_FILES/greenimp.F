       SUBROUTINE GREENIMP(NATOMIMP,DTMTRX,E)

#ifdef CPP_MPI
       use mpi
       use mod_mympi, only: myrank, master, nranks,
     +                      distribute_linear_on_tasks
#else
       use mod_mympi, only: myrank, master, nranks
#endif
       use mod_version_info
       use mod_wunfiles, only: t_params
       IMPLICIT NONE
c-----------------------------------------------------------------
c calculate impurity GF by solving dyson equation 
c N. H. Long, Juelich, 05.2013
c-----------------------------------------------------------------
       INCLUDE 'inc.p'
       INTEGER LMAXSQ
       PARAMETER (LMAXSQ=(LMAXD+1)**2)
       INTEGER LMMAXSO
       PARAMETER (LMMAXSO=(KORBIT+1)*LMAXSQ)
       DOUBLE COMPLEX E,E1
       DOUBLE COMPLEX CONE,CZERO
       PARAMETER (CONE=(1d0,0d0),CZERO=(0d0,0d0))
       INTEGER NATOMIMP,NDIM,INFO
       INTEGER I,J,LM1,LM2,ILM,JLM,ILM1,JLM1
       INTEGER IPVT1(NATOMIMP*LMMAXSO)
       DOUBLE COMPLEX DTMTRX(LMMAXSO*NATOMIMP,LMMAXSO*NATOMIMP)
       DOUBLE COMPLEX, ALLOCATABLE :: GIMP(:,:),GI(:,:)

c read in GF of the host
       ALLOCATE(GIMP(NATOMIMP*LMMAXSO,NATOMIMP*LMMAXSO))
       GIMP=CZERO
       WRITE(6,*) 'read in Green for host'
       READ(60,'(2e17.9)') E1
c       READ(60,*) E1   
       DO J=1,NATOMIMP
        DO LM2=1,LMMAXSO
         JLM=(J-1)*LMMAXSO+LM2
         DO I=1,NATOMIMP
          DO LM1=1,LMMAXSO
           ILM=(I-1)*LMMAXSO+LM1
c           READ(60,*) JLM1,ILM1,GIMP(ILM,JLM)  
            READ(60,'((2I5),(2e17.9))') JLM1,ILM1,GIMP(ILM,JLM)
          ENDDO
         ENDDO
        ENDDO
       ENDDO  

c calculate impurity GF
       ALLOCATE(GI(NATOMIMP*LMMAXSO,NATOMIMP*LMMAXSO))
       GI=CZERO
       NDIM=NATOMIMP*LMMAXSO
c -G_host * delta t
       CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,GIMP,NDIM,
     +            DTMTRX,NDIM,CZERO,GI,NDIM)
       DO I=1,NDIM
        GI(I,I)=CONE+GI(I,I)
       ENDDO

c solve linear equation
       CALL ZGETRF(NDIM,NDIM,GI,NDIM,IPVT1,INFO)
       CALL ZGETRS('N',NDIM,NDIM,GI,NDIM,IPVT1,GIMP,NDIM,INFO)

c write down to the file GMATLL_GES
       WRITE(59,'(2(e17.9,X))') E
       DO LM1=1,NDIM
        DO LM2=1,NDIM
         WRITE(59,'((2I5),(2e17.9))') LM2,LM1,GIMP(LM2,LM1)
        ENDDO
       ENDDO
       DEALLOCATE(GIMP)
       DEALLOCATE(GI)
       END
