      SUBROUTINE MIXLDAU(
     >     MMAXD,NSPIND,NATYPD,NATYP,NSPIN,LOPT,WLDAUOLD,
     X     WLDAU)
      implicit none
c Input:
      INTEGER NATYPD,NSPIND,MMAXD
      INTEGER LOPT(NATYPD)
      INTEGER NATYP,NSPIN
      DOUBLE PRECISION WLDAUOLD(MMAXD,MMAXD,NSPIND,NATYPD)
c Input/Output:
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
c Inside:
      INTEGER IAT,IS,M1,M2,MMAX
      INTEGER IER
      DOUBLE PRECISION XMIX,XMIX2,RMSERR
      CHARACTER*256 UIO ! NCOLIO=256


      EXTERNAL IOINPUT


c First calculate rms error in interaction matrix

      DO IAT = 1,NATYP
         RMSERR = 0.D0
         IF (LOPT(IAT).GE.0) THEN
            MMAX = 2*LOPT(IAT) + 1
            DO IS = 1,NSPIN
               DO M2 = 1,MMAX
                  DO M1 = 1,MMAX
                     RMSERR = RMSERR + 
     &               (WLDAU(M1,M2,IS,IAT) - WLDAUOLD(M1,M2,IS,IAT))**2
                  ENDDO
               ENDDO
            ENDDO
            RMSERR = DSQRT(RMSERR)
            WRITE(1337,9000) IAT,RMSERR
 9000       FORMAT('LDA+U interaction matrix rms error for atom',
     &           I6,' = ',E10.2)
         ENDIF
      ENDDO

c Now mix old/new interaction matrices
      IER = 0
      CALL IOINPUT('MIXLDAU         ',UIO,1,7,IER)
      IF ( IER.NE.0 ) THEN 
         WRITE(*,*) 'MIXLDAU not found, setting to 1.'
         RETURN
      ELSE
         READ (UNIT=UIO,FMT=*) XMIX
         WRITE(1337,*) 'Using MIXLDAU = ',XMIX
      ENDIF

      XMIX2 = 1.D0 - XMIX

      DO IAT = 1,NATYP
         IF (LOPT(IAT).GE.0) THEN
            DO IS = 1,NSPIN
               DO M2 = 1,MMAXD
                  DO M1 = 1,MMAXD
                     WLDAU(M1,M2,IS,IAT) = XMIX * WLDAU(M1,M2,IS,IAT)
     &                               + XMIX2 * WLDAUOLD(M1,M2,IS,IAT)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO

      END
