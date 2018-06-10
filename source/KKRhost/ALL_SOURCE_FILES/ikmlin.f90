SUBROUTINE ikmlin(iprint,nsollm,ikm1lin,ikm2lin,nlmax,nmuemax,  &
        linmax,nl)
!   ********************************************************************
!   *                                                                  *
!   * SETUP TABLE OF INDICES    IKM(INT)                               *
!   *                                                                  *
!   *  IKM IS STANDARD INDEX IN  (KAPPA,MUE)-REPRESENTATION            *
!   *  IKM = 2*L*(J+1/2) + J + MUE + 1                                 *
!   *                                                                  *
!   *  INT NUMBERS LINEARLY ONLY NON-VANISHING ELEMENTS OF M-SS        *
!   *  USED TO CALCULATE DOS ...                                       *
!   *                                                                  *
!   ********************************************************************
use mod_types, only: t_inc
IMPLICIT NONE


! Dummy arguments
INTEGER IPRINT,LINMAX,NL,NLMAX,NMUEMAX
INTEGER IKM1LIN(LINMAX),IKM2LIN(LINMAX),NSOLLM(NLMAX,NMUEMAX)

! Local variables
INTEGER I,IL,IMUE,K1,K2,KAP(2),L,LIN,MUEM05,NSOL
INTEGER IKAPMUE

lin = 0

DO il = 1,nl
  l = il - 1
  muem05 = -il - 1
  kap(1) = -l - 1
  kap(2) = +l
  
  DO imue = 1,2*il
    muem05 = muem05 + 1
    nsol = nsollm(il,imue)
    
    DO k2 = 1,nsol
      DO k1 = 1,nsol
        lin = lin + 1
        ikm1lin(lin) = ikapmue(kap(k1),muem05)
        ikm2lin(lin) = ikapmue(kap(k2),muem05)
      END DO
    END DO
    
  END DO
END DO

IF ( iprint < 2 ) RETURN
IF(t_inc%i_write>0) THEN
  WRITE (1337,FMT='('' INT='',I3,''  IKM=('',I3,'','',I3,'')'')')  &
      (i,ikm1lin(i),ikm2lin(i),i=1,lin)
END IF
END SUBROUTINE ikmlin
