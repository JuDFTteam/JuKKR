        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar 23 23:37:03 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GATHER_TGMAT_MPI__genmod
          INTERFACE 
            SUBROUTINE GATHER_TGMAT_MPI(LMAXD,IELAST,EZ,WZ,MY_RANK,     &
     &MPI_SIZE,MPI_IEBOUNDS)
              INTEGER(KIND=4), INTENT(IN) :: MPI_SIZE
              INTEGER(KIND=4), INTENT(IN) :: IELAST
              INTEGER(KIND=4), INTENT(IN) :: LMAXD
              COMPLEX(KIND=8), INTENT(IN) :: EZ(IELAST)
              COMPLEX(KIND=8), INTENT(IN) :: WZ(IELAST)
              INTEGER(KIND=4), INTENT(IN) :: MY_RANK
              INTEGER(KIND=4), INTENT(IN) :: MPI_IEBOUNDS(2,0:MPI_SIZE-1&
     &)
            END SUBROUTINE GATHER_TGMAT_MPI
          END INTERFACE 
        END MODULE GATHER_TGMAT_MPI__genmod
