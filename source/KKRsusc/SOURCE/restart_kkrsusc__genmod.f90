        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar 23 23:37:04 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RESTART_KKRSUSC__genmod
          INTERFACE 
            SUBROUTINE RESTART_KKRSUSC(ITC,LMAXD,CONFIG,IELAST,EZ,WZ,   &
     &MY_RANK)
              USE TYPE_CONFIG
              INTEGER(KIND=4), INTENT(IN) :: IELAST
              INTEGER(KIND=4), INTENT(IN) :: ITC
              INTEGER(KIND=4), INTENT(IN) :: LMAXD
              TYPE (CONFIG_TYPE) :: CONFIG
              COMPLEX(KIND=8), INTENT(IN) :: EZ(IELAST)
              COMPLEX(KIND=8), INTENT(IN) :: WZ(IELAST)
              INTEGER(KIND=4) :: MY_RANK
            END SUBROUTINE RESTART_KKRSUSC
          END INTERFACE 
        END MODULE RESTART_KKRSUSC__genmod
