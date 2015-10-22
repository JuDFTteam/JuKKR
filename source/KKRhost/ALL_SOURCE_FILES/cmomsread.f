C*==cmomsread.f    processed by SPAG 6.05Rc at 19:12 on  5 Jun 2004
      SUBROUTINE CMOMSREAD(NLBASIS,NRBASIS,NAEZ,CMOMHOST,VACFLAG,KAOEZ,
     &                     NATYPD,NEMBD1,LMPOTD)
C **********************************************************************
C *                                                                    *
C * This subroutine reads in the CMOMHOST from the decimation files    *
C * Note that they are needed only in case of an SCF decimation cal-   *
C * culation (SCFSTEPS > 1 )                                           *
C *                                                                    *
C * The t-matrices are writen out in kloopz1  (option 'deci-out')      *
C *                                                                    *
C * This subroutine must be called after the t-matrices for all the    *
C * energies are read in (see < decimaread > )                         *
C * It returns the CMOMHOST array. First the left bulk (unit 37) then  *
C * the right bulk (unit 38) are indexed.                              *
C * CMOMHOST(*,NEMBD1) =                                               *
C *                CMOMHOST(*,1..NLBASIS,NLBASIS+1..NLBASIS+NRBASIS)   *
C * Condider this mapping for further use.                             *
C *                                                                    *
C *                                                  29.10.99          *
C *                                                  05.06.04          *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Arguments
      INTEGER LMPOTD,NATYPD,NEMBD1
      INTEGER NAEZ,NLBASIS,NRBASIS
      DOUBLE PRECISION CMOMHOST(LMPOTD,NEMBD1)
      INTEGER KAOEZ(NATYPD,*)
      LOGICAL VACFLAG(2)
C     ..
C     .. Local variables ..
      DOUBLE PRECISION C00(LMPOTD)
      CHARACTER*5 CHHOST(2)
      INTEGER IH,IH1,IHL,IHOST,LM,LMPOTL,NAEZL,NATHOST
C     ..
C     .. Data statements
      DATA CHHOST/'LEFT ','RIGHT'/
C     ..
      WRITE (1337,'(5X,A,/,8X,30(1H-),/,8X,3A6,A10,/,8X,30(1H-))')
     &        'Reading in host charge moments ( SCFSTEPS > 1 )',
     &       ' HOST ','  IBAS','  ATOM','   CMOM(1)'
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
      DO IHOST = 1,2
         NATHOST = NLBASIS
         IF ( IHOST.EQ.2 ) NATHOST = NRBASIS
         WRITE (1337,'(8X,A5,1X,$)') CHHOST(IHOST)
C ----------------------------------------------------------------------
         IF ( VACFLAG(IHOST) ) THEN
            DO IH = 1,NLBASIS
               DO LM = 1,LMPOTD
                  CMOMHOST(LM,(IHOST-1)*NLBASIS+IH) = 0.D0
C mapping the CMOMHOST array, ordering is important
               END DO
            END DO
            WRITE (1337,'(A)') ' Vacuum setting    0.000'
            IF ( IHOST.EQ.1 ) THEN
               WRITE (1337,'(14X,24(1H-))')
            ELSE
               WRITE (1337,'(8X,30(1H-))')
            END IF
C ----------------------------------------------------------------------
         ELSE
            READ (36+IHOST,99002) NAEZL,LMPOTL
C ......................................................................
            IF ( NAEZL.NE.NATHOST ) THEN
               WRITE (6,'(/,5X,2A)') 'ERROR: ',
     &                            'host not compatible with your input.'
               WRITE (6,'(/,12X,A,I3,A,I3)')
     &                 'Charge moments tabulated for',NAEZL,
     &                ' host atoms, input NBASIS =',NATHOST
               STOP '       < CMOMSREAD > '
            END IF
C ......................................................................
            DO IH = 1,NAEZL
               READ (36+IHOST,*) IHL
               IF ( IHL.NE.IH ) THEN
                  WRITE (6,'(/,5X,2A,/)') 'ERROR reading host file',
     &                   ' basis indexing wrong'
                  STOP '       < CMOMSREAD > '
               END IF
               READ (36+IHOST,99001) (C00(LM),LM=1,LMPOTL)
               IH1 = KAOEZ(1,NAEZ+(IHOST-1)*NLBASIS+IH)
C
               DO LM = 1,LMPOTL
                  CMOMHOST(LM,(IHOST-1)*NLBASIS+IH) = C00(LM)
C mapping the CMOMHOST array, ordering is important
               END DO
C
               IF ( IH.EQ.1 ) THEN
                  WRITE (1337,'(1X,2I6,D12.4)') IH,IH1,C00(1)
               ELSE
                  WRITE (1337,'(14X,2I6,D12.4)') IH,IH1,C00(1)
               END IF
            END DO
C ......................................................................
            IF ( IHOST.EQ.1 ) THEN
               WRITE (1337,'(14X,24(1H-))')
            ELSE
               WRITE (1337,'(8X,30(1H-))')
            END IF
         END IF
C ----------------------------------------------------------------------
      END DO
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
      WRITE (1337,*)
C
99001 FORMAT (4D22.14)
99002 FORMAT (5X,2I6)
      END
