        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar 23 23:37:02 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SCATT_SOL_ROTATION__genmod
          INTERFACE 
            SUBROUTINE SCATT_SOL_ROTATION(IE,NATOM,LMMAXD,THETA,PHI)
              USE GLOBAL
              INTEGER(KIND=4), INTENT(IN) :: LMMAXD
              INTEGER(KIND=4), INTENT(IN) :: NATOM
              INTEGER(KIND=4), INTENT(IN) :: IE
              REAL(KIND=8), INTENT(IN) :: THETA(NATOM)
              REAL(KIND=8), INTENT(IN) :: PHI(NATOM)
            END SUBROUTINE SCATT_SOL_ROTATION
          END INTERFACE 
        END MODULE SCATT_SOL_ROTATION__genmod
