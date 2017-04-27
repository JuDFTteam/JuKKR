        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar 23 23:37:04 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE KKRSUSC_PREPARE__genmod
          INTERFACE 
            SUBROUTINE KKRSUSC_PREPARE(NATOM,NSPIN,LMAXD,ZATOM,NRMAXD,  &
     &VPOT,DENSITY,CELL,CONFIG,ITSCF,INPSUSC,MY_RANK)
              USE TYPE_INPSUSC
              USE TYPE_CONFIG
              USE TYPE_CELL
              USE TYPE_DENSITY
              INTEGER(KIND=4), INTENT(IN) :: NRMAXD
              INTEGER(KIND=4), INTENT(IN) :: LMAXD
              INTEGER(KIND=4), INTENT(IN) :: NSPIN
              INTEGER(KIND=4), INTENT(IN) :: NATOM
              REAL(KIND=8), INTENT(IN) :: ZATOM(NATOM)
              REAL(KIND=8), INTENT(IN) :: VPOT(NRMAXD,(2*LMAXD+1)**2,   &
     &NSPIN,NATOM)
              TYPE (DENSITY_TYPE), INTENT(IN) :: DENSITY(NATOM)
              TYPE (CELL_TYPE), INTENT(IN) :: CELL(NATOM)
              TYPE (CONFIG_TYPE), INTENT(IN) :: CONFIG
              INTEGER(KIND=4), INTENT(IN) :: ITSCF
              TYPE (INPSUSC_TYPE), INTENT(OUT) :: INPSUSC
              INTEGER(KIND=4), INTENT(IN) :: MY_RANK
            END SUBROUTINE KKRSUSC_PREPARE
          END INTERFACE 
        END MODULE KKRSUSC_PREPARE__genmod
