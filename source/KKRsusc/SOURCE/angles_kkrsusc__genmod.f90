        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar 23 23:37:01 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ANGLES_KKRSUSC__genmod
          INTERFACE 
            SUBROUTINE ANGLES_KKRSUSC(NATOM,THETA,PHI)
              INTEGER(KIND=4) :: NATOM
              REAL(KIND=8), INTENT(IN) :: THETA(NATOM)
              REAL(KIND=8), INTENT(IN) :: PHI(NATOM)
            END SUBROUTINE ANGLES_KKRSUSC
          END INTERFACE 
        END MODULE ANGLES_KKRSUSC__genmod
