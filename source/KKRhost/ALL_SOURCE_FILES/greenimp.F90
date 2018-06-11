SUBROUTINE greenimp(natomimp,dtmtrx,e)
#ifdef CPP_MPI
use mpi
use mod_mympi, only: myrank, master, nranks, distribute_linear_on_tasks
#ELSE
use mod_mympi, only: myrank, master, nranks
#ENDIF
use mod_version_info
use mod_wunfiles, only: t_params

       IMPLICIT NONE
!-----------------------------------------------------------------
! calculate impurity GF by solving dyson equation 
! N. H. Long, Juelich, 05.2013
!-----------------------------------------------------------------
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

! read in GF of the host
allocate(gimp(natomimp*lmmaxso,natomimp*lmmaxso))

gimp=czero
WRITE(6,*) 'read in Green for host'
READ(60,'(2e17.9)') e1
!       READ(60,*) E1
DO j=1,natomimp
  DO lm2=1,lmmaxso
    jlm=(j-1)*lmmaxso+lm2
    DO i=1,natomimp
      DO lm1=1,lmmaxso
        ilm=(i-1)*lmmaxso+lm1
!           READ(60,*) JLM1,ILM1,GIMP(ILM,JLM)
        READ(60,'((2I5),(2e17.9))') jlm1,ilm1,gimp(ilm,jlm)
      END DO
    END DO
  END DO
END DO

! calculate impurity GF
allocate(gi(natomimp*lmmaxso,natomimp*lmmaxso))
gi=czero
ndim=natomimp*lmmaxso
! -G_host * delta t
CALL zgemm('N','N',ndim,ndim,ndim,-cone,gimp,ndim, dtmtrx,ndim,czero,gi,ndim)
DO i=1,ndim
  gi(i,i)=cone+gi(i,i)
END DO

! solve linear equation
CALL zgetrf(ndim,ndim,gi,ndim,ipvt1,info)
CALL zgetrs('N',ndim,ndim,gi,ndim,ipvt1,gimp,ndim,info)

! write down to the file GMATLL_GES
WRITE(59,'(2(e17.9,1X))') e
DO lm1=1,ndim
  DO lm2=1,ndim
    WRITE(59,'((2I5),(2e17.9))') lm2,lm1,gimp(lm2,lm1)
  END DO
END DO
deallocate(gimp)
deallocate(gi)
END SUBROUTINE greenimp
