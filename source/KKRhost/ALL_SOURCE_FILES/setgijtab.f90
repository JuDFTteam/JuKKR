SUBROUTINE setgijtab(linterface,icc,naez,iqat,rbasis,bravais,  &
    natomimp,atomimp,rclsimp,nofgij,ijtabcalc,  &
    iofgij,jofgij,nqcalc,iqcalc,natomimpd, ijtabcalc_i)
! **********************************************************************
! * Task-specific settings of Gij elements that need to be calculated  *
! * Subroutine (called for ICC=-1) sets up the arrays                  *
! * NATOMIMP    : number of different sites i,j = 1,NATOMIMP           *
! * RCLSIMP     : site coordinates                                     *
! * ATOMIMP     : index of the corresponding site in the unit cell     *
! * IJTABCALC   : flag specifying wehter pair (I,J) needs to be        *
! *               calculated - linear pointer (I-1)*NATOMIMP + J = 1/0 *
! *               for YES/NO                                           *
! * NOFGIJ      : number of all (I,J) pairs - sum of all non-zero I,J  *
! * IOFGIJ      : I index in the list 1..NATOMIMP for pair I,J         *
! * JOFGIJ      : J index                                              *
! **********************************************************************
IMPLICIT NONE

!Scalar arguments
INTEGER ICC,NAEZ,NATOMIMP,NATOMIMPD,NOFGIJ,NQCALC
LOGICAL LINTERFACE
 
!Array arguments
INTEGER ATOMIMP(*),IJTABCALC(*),IJTABCALC_I(*),IOFGIJ(*),IQAT(*), &
        IQCALC(*),JOFGIJ(*)
DOUBLE PRECISION BRAVAIS(3,3),RBASIS(3,*),RCLSIMP(3,*)

!Local scalars
INTEGER I,IDO,II,J,JJ,NN
LOGICAL OPT

!External subroutines
EXTERNAL GIJCOND,GIJXCPL,OPT


! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,'(79(1H=),/,15X,A)') 'SETGIJTAB: setting task-specific Gij pairs'
WRITE (1337,'(79(1H=),/)')
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

ido = 0
! ======================================================================
IF ( opt('CONDUCT ') ) CALL gijcond(ido,naez,rbasis,iqat,natomimp,  &
    rclsimp,atomimp,ijtabcalc,natomimpd)
! ======================================================================
IF ( opt('XCPL    ') ) CALL gijxcpl(ido,naez,rbasis,bravais,  &
    linterface,nqcalc,iqcalc,natomimp,rclsimp,atomimp,ijtabcalc,  &
    ijtabcalc_i,natomimpd)
! ======================================================================
IF ( ido == 0 ) THEN
  icc = 0
  WRITE (6,99002)
  RETURN
END IF
! ======================================================================
nofgij = 0
DO i = 1,natomimp
  nn = (i-1)*natomimp
  DO j = 1,natomimp
    IF ( ijtabcalc(nn+j) > 0 ) THEN
      nofgij = nofgij + 1
      IF ( nofgij > natomimpd*natomimpd ) THEN
        WRITE (6,99001) 'NATOMIMPD',nofgij/natomimp
        STOP
      END IF
      iofgij(nofgij) = i
      jofgij(nofgij) = j
    END IF
  END DO
END DO
IF ( nofgij == 0 ) THEN
  icc = 0
  WRITE (6,99002)
  RETURN
END IF

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,99003) natomimp,nofgij
WRITE (1337,99004)
WRITE (1337,99005)
WRITE (1337,99004)
DO i = 1,nofgij
  ii = iofgij(i)
  jj = jofgij(i)
  WRITE (1337,99006) i,ii,atomimp(ii),(rclsimp(j,ii),j=1,3),jj,  &
      atomimp(jj),(rclsimp(j,jj),j=1,3)
END DO
WRITE (1337,99004)
WRITE (1337,*)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

99001 FORMAT (6X,'brahim ERROR: please increase the global parameter'  &
    ,/,6X,a,' to a value >=',i5,/)
99002 FORMAT (6X,'WARNING: Subroutine entered with invalid task ',  &
    'specification',/,6X,  &
    '         ICC will be set to 0 - no Gij calculated - ', 'input check? ',/)
99003 FORMAT (6X,'Number of different sites (NATOMIMP) :',i4,/,6X,  &
    'Number of pairs set       (NOFGIJ)   :',i4)
99004 FORMAT (8X,71('-'))
99005 FORMAT (9X,'pair|',' I  IQ           position',9X,  &
    'J  JQ           position')
99006 FORMAT (9X,i3,' |',2(i3,1X),3F8.4,1X,2(i3,1X),3F8.4)
99007 FORMAT (i5,2(i5,1X),3F10.6,1X,2(i5,1X),3F10.6)
END SUBROUTINE setgijtab
