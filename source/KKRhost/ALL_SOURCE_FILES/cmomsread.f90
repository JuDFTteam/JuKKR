SUBROUTINE cmomsread(nlbasis,nrbasis,naez,cmomhost,vacflag,kaoez,  &
    natypd,nembd1,lmpotd)
! **********************************************************************
! *                                                                    *
! * This subroutine reads in the CMOMHOST from the decimation files    *
! * Note that they are needed only in case of an SCF decimation cal-   *
! * culation (SCFSTEPS > 1 )                                           *
! *                                                                    *
! * The t-matrices are writen out in kloopz1  (option 'deci-out')      *
! *                                                                    *
! * This subroutine must be called after the t-matrices for all the    *
! * energies are read in (see < decimaread > )                         *
! * It returns the CMOMHOST array. First the left bulk (unit 37) then  *
! * the right bulk (unit 38) are indexed.                              *
! * CMOMHOST(*,NEMBD1) =                                               *
! *                CMOMHOST(*,1..NLBASIS,NLBASIS+1..NLBASIS+NRBASIS)   *
! * Condider this mapping for further use.                             *
! *                                                                    *
! *                                                  29.10.99          *
! *                                                  05.06.04          *
! **********************************************************************
      IMPLICIT NONE
!..
!.. Arguments
      INTEGER LMPOTD,NATYPD,NEMBD1
      INTEGER NAEZ,NLBASIS,NRBASIS
      DOUBLE PRECISION CMOMHOST(LMPOTD,NEMBD1)
      INTEGER KAOEZ(NATYPD,*)
      LOGICAL VACFLAG(2)
!..
!.. Local variables ..
      DOUBLE PRECISION C00(LMPOTD)
      CHARACTER*5 CHHOST(2)
      INTEGER IH,IH1,IHL,IHOST,LM,LMPOTL,NAEZL,NATHOST
!..
!.. Data statements
      DATA CHHOST/'LEFT ','RIGHT'/
!..
WRITE (1337,'(5X,A,/,8X,30(1H-),/,8X,3A6,A10,/,8X,30(1H-))')  &
    'Reading in host charge moments ( SCFSTEPS > 1 )',  &
    ' HOST ','  IBAS','  ATOM','   CMOM(1)'
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
DO ihost = 1,2
  nathost = nlbasis
  IF ( ihost == 2 ) nathost = nrbasis
  WRITE (1337,'(8X,A5,1X,$)') chhost(ihost)
! ----------------------------------------------------------------------
  IF ( vacflag(ihost) ) THEN
    DO ih = 1,nlbasis
      DO lm = 1,lmpotd
        cmomhost(lm,(ihost-1)*nlbasis+ih) = 0.d0
! mapping the CMOMHOST array, ordering is important
      END DO
    END DO
    WRITE (1337,'(A)') ' Vacuum setting    0.000'
    IF ( ihost == 1 ) THEN
      WRITE (1337,'(14X,24(1H-))')
    ELSE
      WRITE (1337,'(8X,30(1H-))')
    END IF
! ----------------------------------------------------------------------
  ELSE
    READ (36+ihost,99002) naezl,lmpotl
! ......................................................................
    IF ( naezl /= nathost ) THEN
      WRITE (6,'(/,5X,2A)') 'ERROR: ', 'host not compatible with your input.'
      WRITE (6,'(/,12X,A,I3,A,I3)') 'Charge moments tabulated for',naezl,  &
          ' host atoms, input NBASIS =',nathost
      STOP '       < CMOMSREAD > '
    END IF
! ......................................................................
    DO ih = 1,naezl
      READ (36+ihost,*) ihl
      IF ( ihl /= ih ) THEN
        WRITE (6,'(/,5X,2A,/)') 'ERROR reading host file',  &
            ' basis indexing wrong'
        STOP '       < CMOMSREAD > '
      END IF
      READ (36+ihost,99001) (c00(lm),lm=1,lmpotl)
      ih1 = kaoez(1,naez+(ihost-1)*nlbasis+ih)
      
      DO lm = 1,lmpotl
        cmomhost(lm,(ihost-1)*nlbasis+ih) = c00(lm)
! mapping the CMOMHOST array, ordering is important
      END DO
      
      IF ( ih == 1 ) THEN
        WRITE (1337,'(1X,2I6,D12.4)') ih,ih1,c00(1)
      ELSE
        WRITE (1337,'(14X,2I6,D12.4)') ih,ih1,c00(1)
      END IF
    END DO
! ......................................................................
    IF ( ihost == 1 ) THEN
      WRITE (1337,'(14X,24(1H-))')
    ELSE
      WRITE (1337,'(8X,30(1H-))')
    END IF
  END IF
! ----------------------------------------------------------------------
END DO
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
WRITE (1337,*)

99001 FORMAT (4D22.14)
99002 FORMAT (5X,2I6)
END SUBROUTINE cmomsread
