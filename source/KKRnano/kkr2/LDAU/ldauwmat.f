C*==wmatldau.f    processed by SPAG 6.05Rc at 16:58 on 22 Dec 2004
      SUBROUTINE LDAUWMAT(I1,NSPIN,ITER,MIXING,DMATLDAU,NLDAU,LLDAU,
     &                    ULDAU,JLDAU,UMLDAU,WMLDAU,
     &                    EULDAU,EDCLDAU,
C                         new input parameters after inc.p removal
     &                    lmaxd)
C **********************************************************************
C *                                                                    *
C * Calculation of Coulomb interaction potential in LDA+U              *
C * non-relativistic case -- otherwise matrices DENMAT and VLDAU must  *
C *                          have double dimension                     *
C *                                                                    *
C * Uses the Coulomb matrix U (array ULDAU), the density matrix n      *
C * (array DENMAT) and the occupation numbers dentot (total) and n_s   *
C * (array DENTOTS) (per spin).                                        *
C *                                                                    *
C * The expression evaluated (array VLDAU) is                          *
C *                                                                    *
C *       V_{m1,s,m2,s'} =                                             *
C * delta_{ss'} Sum_{s'',m3,m4} U_{m1,m2,m3,m4} n_{m3,s'',m4,s''}      *
C * - Sum_{m3,m4} U_{m1,m4,m3,m2} n_{m3,s',m4,s}                       *
C * - [Ueff (dentot-1/2) - Jeff (n_s - 1/2)] delta_{ss'} delta_{m1,m2} *
C *                                                                    *
C *                  ph. mavropoulos, h.ebert munich/juelich 2002-2004 *
C **********************************************************************
      IMPLICIT NONE

      INTEGER lmaxd

C PARAMETER definitions
C
C     INTEGER            LMAXD1
C     PARAMETER          (LMAXD1 = LMAXD + 1 )
C     INTEGER            MMAXD
C     PARAMETER          (MMAXD  = 2*LMAXD + 1 )

      DOUBLE COMPLEX     CZERO
      PARAMETER          (CZERO=(0.0D0,0.0D0))
C
C Global arguments
C
      INTEGER            I1,NSPIN,ITER,NLDAU
      INTEGER            LLDAU(LMAXD + 1)
C     DOUBLE PRECISION   UMLDAU(MMAXD,MMAXD,MMAXD,MMAXD,LMAXD1)
      DOUBLE PRECISION   UMLDAU(2*LMAXD + 1,2*LMAXD + 1,
     &                          2*LMAXD + 1,2*LMAXD + 1,LMAXD + 1)
      DOUBLE PRECISION   ULDAU(LMAXD + 1)
      DOUBLE PRECISION   JLDAU(LMAXD + 1)
      DOUBLE PRECISION   EDCLDAU
      DOUBLE PRECISION   EULDAU
      DOUBLE PRECISION   EDC(LMAXD + 1)
      DOUBLE PRECISION   EU(LMAXD + 1)

      DOUBLE PRECISION   WMLDAU(2*LMAXD + 1,2*LMAXD + 1,LMAXD + 1,NSPIN)
      DOUBLE PRECISION   MIXING
C     DOUBLE COMPLEX     DMATLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
C     DOUBLE COMPLEX     DDUMMY(MMAXD,MMAXD,NSPIND)

      DOUBLE COMPLEX     DMATLDAU(2*LMAXD + 1,2*LMAXD + 1,NSPIN,LMAXD+1)
      DOUBLE COMPLEX     DDUMMY(2*LMAXD + 1,2*LMAXD + 1,NSPIN)
C
C Local variables
C
      DOUBLE COMPLEX     CSUM,CSUM2
C     DOUBLE COMPLEX     VMLDAU(MMAXD,MMAXD,NSPIND)
C     DOUBLE PRECISION   DENMAT(MMAXD,MMAXD,NSPIND)

      DOUBLE COMPLEX     VMLDAU(2*LMAXD + 1, 2*LMAXD + 1, NSPIN)
      DOUBLE PRECISION   DENMAT(2*LMAXD + 1, 2*LMAXD + 1, NSPIN)

      DOUBLE PRECISION   DENTOT,
     &                   DENTOTS(NSPIN),
     &                   ALPHA
      INTEGER            ISPIN,JSPIN,M1,M2,M3,M4,MM,MMAX,ILDAU
      INTEGER            LRECLDAU
      LOGICAL            SLOW,FREEZE

      INTEGER            LMAXD1
      INTEGER            MMAXD
      INTEGER            NSPIND

      LMAXD1 = LMAXD + 1
      MMAXD  = 2*LMAXD + 1
      NSPIND = NSPIN

      DMATLDAU = CZERO
      EULDAU = 0.0D0
      EDCLDAU = 0.0D0

      WRITE (6,'(/,79(1H#),/,16X,A,/,79(1H#))') 
     &     'LDA+U: Calculating interaction potential VLDAU'
C
C ..
C
      SLOW = .FALSE.
      INQUIRE(FILE='SLOW',EXIST=SLOW)
C
      FREEZE = .FALSE.
      INQUIRE(FILE='FREEZE',EXIST=FREEZE)
C
      ALPHA = 0.0D0
C
      IF (MOD(ITER-1,10).EQ.0) THEN
C      IF (MOD(ITER-1,10).EQ.0) THEN
        IF (.NOT.SLOW.AND.ITER.LE.1) THEN
          ALPHA = 1.0D0
        ELSE
          ALPHA = MIXING
        ENDIF
      ENDIF
C
      IF (FREEZE) ALPHA = 0.0D0
C
      WRITE(6,*) 'LDA+U: MIXING=',ALPHA,' for iteration',ITER
C
C=======================================================================
C=======================================================================
      DO ILDAU=1,NLDAU
      IF (LLDAU(ILDAU).GE.0) THEN
C=======================================================================
C
        EDC(ILDAU) = 0.0D0
        EU(ILDAU)  = 0.0D0
C
        MMAX = 2*LLDAU(ILDAU) + 1
        WRITE (6,99001) I1,LLDAU(ILDAU)
C
        DO M3=1,MMAXD
          DO M4=1,MMAXD
            DO ISPIN = 1,NSPIN
              DDUMMY(M3,M4,ISPIN) = DMATLDAU(M3,M4,ISPIN,ILDAU)
            ENDDO
          ENDDO
        ENDDO
C ..
C
        CALL RINIT(MMAXD*MMAXD*NSPIND,DENMAT(1,1,1))
C
C Result is in real Ylm basis. 
C It must be converted to complex Ylm basis:
C
C ----------------------------------------------------------------------
C -> Convert DENMATC and DENMAT to complex spherical harmonics.
C
        DO ISPIN = 1,NSPIN
          CALL RCLM(1,LLDAU(ILDAU),LMAXD,DDUMMY(1,1,ISPIN))
        ENDDO
C
C ----------------------------------------------------------------------
        DENTOT = 0.D0
C ----------------------------------------------------------------------
C
        DO ISPIN = 1,NSPIN
C
C -> DENMAT is real: (imag(denmatc))
C
          DO M2 = 1,MMAX
            DO M1 = 1,MMAX
              DENMAT(M1,M2,ISPIN) = DIMAG(DDUMMY(M1,M2,ISPIN))
            ENDDO
          ENDDO
C
C 2.  Calculate total occupation numbers:
C ntot_s = Sum_m n_{m,s,m,s}, ntot = n_1 + n_2
C
          DENTOTS(ISPIN) = 0.D0
          DO MM = 1,MMAX
            DENTOTS(ISPIN) = DENTOTS(ISPIN) + DENMAT(MM,MM,ISPIN)
          ENDDO
          DENTOT = DENTOT + DENTOTS(ISPIN)
C
        ENDDO
C
C ----------------------------------------------------------------------
C
C In paramagnetic case the spin degeneracy has been accounted
C for by the weight DF in tmatrho.
C
C ----------------------------------------------------------------------
        CALL CINIT(MMAXD*MMAXD*NSPIND,VMLDAU(1,1,1))
C
        DO ISPIN = 1,NSPIN
C
C 3.  Use density matrix and Coulomb matrix ULDAU to calculate the
C interaction potential VLDAU
C 3a. First part (always diagonal in spin).
C
          DO M2 = 1,MMAX
            DO M1 = 1,MMAX
              CSUM = CZERO
              DO M4 = 1,MMAX
                DO M3 = 1,MMAX
                  CSUM2 = CZERO
                  DO JSPIN = 1,NSPIN
                    CSUM2 = CSUM2 + DENMAT(M3,M4,JSPIN)
                  ENDDO
                  CSUM = CSUM + UMLDAU(M1,M2,M3,M4,ILDAU)*CSUM2
                ENDDO
              ENDDO
              VMLDAU(M1,M2,ISPIN) = VMLDAU(M1,M2,ISPIN) + CSUM
            ENDDO
          ENDDO
C
C 3b. Second part (in fully rel. case not diagonal in spin; then this
C loop must be changed accordingly).
C
          DO M2 = 1,MMAX
            DO M1 = 1,MMAX
              CSUM = CZERO
              DO M4 = 1,MMAX
                DO M3 = 1,MMAX
                  CSUM = CSUM 
     &                 - UMLDAU(M1,M4,M3,M2,ILDAU)*DENMAT(M3,M4,ISPIN)
                ENDDO
              ENDDO
              VMLDAU(M1,M2,ISPIN) = VMLDAU(M1,M2,ISPIN) + CSUM
            ENDDO
          ENDDO
C
C 3c. Third part (always spin- and m-diagonal).
C
          DO M1 = 1,MMAX
            VMLDAU(M1,M1,ISPIN) = VMLDAU(M1,M1,ISPIN) 
     &                      - ULDAU(ILDAU)*(DENTOT-0.5D0)
     &                      + JLDAU(ILDAU)*(DENTOTS(ISPIN)-0.5D0)
          ENDDO
C
C 4. Calculate total-energy corrections EU and EDC (double-counting).
C Then the correction is EU-EDC.
C Note: EU,EDC initialised outside the routine
C
C Here VLDAU is assumed spin-diagonal (contrary to the spin-orbit case).
C
          DO M2 = 1,MMAX
            DO M1 = 1,MMAX
              EU(ILDAU) = EU(ILDAU)
     +           + DENMAT(M1,M2,ISPIN) * DREAL(VMLDAU(M1,M2,ISPIN))
            ENDDO
          ENDDO
          EDC(ILDAU) = EDC(ILDAU)
     &        - JLDAU(ILDAU)*DENTOTS(ISPIN)*(DENTOTS(ISPIN)-1.D0)
        ENDDO

C ----------------------------------------------------------------------
!           IF ( IPRINT.GT.0 ) WRITE (6,99002)
!      &         'Interaction potential in COMPLEX basis:'
C ----------------------------------------------------------------------
        DO ISPIN = 1,NSPIN
C
C 5.  Transform VLDAU into real spherical harmonics basis
C
          CALL RCLM(2,LLDAU(ILDAU),LMAXD,VMLDAU(1,1,ISPIN))
C
C Copy transformed VLDAU to real WLDAU
C
C Apply damping to the interaction matrix WLDAU ? Here not.
C
C
          DO M2 = 1,MMAX
            DO M1 = 1,MMAX
              WMLDAU(M1,M2,ILDAU,ISPIN) =
     +           ALPHA     * DREAL(VMLDAU(M1,M2,ISPIN))
     +         + (1.0D0-ALPHA) * WMLDAU(M1,M2,ILDAU,ISPIN)
            ENDDO
          ENDDO
C
        ENDDO

C ----------------------------------------------------------------------
        WRITE (6,99002) 'Interaction potential (real):'
        DO ISPIN = 1,NSPIN
          WRITE(6,99003) ISPIN
          CALL RWRITE(WMLDAU(1,1,ILDAU,ISPIN),MMAXD,MMAX,6)
        ENDDO
        WRITE(6,*)

C ----------------------------------------------------------------------
C
C Corrections in total energy:
C
        EU(ILDAU)  =0.5D0* EU(ILDAU)
        EDC(ILDAU) =0.5D0*(ULDAU(ILDAU)*DENTOT*(DENTOT-1.D0)-EDC(ILDAU))
C
C sum up cntribution of all orbitals treated by LDA+U
C
        EULDAU  = EULDAU  + EU(ILDAU)
        EDCLDAU = EDCLDAU + EDC(ILDAU)
C
C -> Write out corrections on energy:
C    E[LDA+U] = E[LDA] + EU - EDC
C
        IF (ILDAU.EQ.NLDAU) THEN
          WRITE(6,99002) 'Corrections to the total energy:'
          WRITE(6,*)
          WRITE(6,99004) 'EU  =',EULDAU
          WRITE(6,99004) 'Edc =',EDCLDAU
          WRITE(6,99006) 'E[LDA+U] = E[LDA] + EU - Edc'
        ENDIF
C
C=======================================================================
      ENDIF
      ENDDO
C=======================================================================
C=======================================================================
C
C
C
C update file 'wldau.unf':
      LRECLDAU = 4*(1+LMAXD1)                   ! NLDAU & LLDAU
     +         + 8*2*LMAXD1                     ! ULDAU & JLDAU
     +         + 8*MMAXD*MMAXD*NSPIND*LMAXD1    ! WMLDAU
C
      OPEN (65,ACCESS='direct',RECL=LRECLDAU,FILE='wldau.unf',
     +      FORM='unformatted')
      WRITE (65,REC=I1) NLDAU,LLDAU,ULDAU,JLDAU,WMLDAU
      CLOSE(65)
C
C
C ..
C
C
99001 FORMAT(/,6X,65(1H=),/,6X,'Atom :',I3,' (l =',I2,')',/,6X,18(1H=))
99002 FORMAT(8X,'* ',A)
99003 FORMAT(/,15X,'> ISPIN =',I1)
99004 FORMAT(10X,A,F10.6)
99005 FORMAT(10X,21(1H-),/,10X,A,F10.6,/,10X,60(1H-),/)
99006 FORMAT(27X,A,/)
C
      END
C
C
C
C
C
C*==rwrite.f    processed by SPAG 6.05Rc at 16:58 on 22 Dec 2004
      SUBROUTINE RWRITE(Z,MMAXD,MMAX,IFILE)
      IMPLICIT NONE
      INTEGER IFILE,MMAX,MMAXD
      DOUBLE PRECISION Z(MMAXD,MMAXD)
      INTEGER M1,M2
C
      WRITE(IFILE,99000)
      DO M2 = 1,MMAX
         WRITE (IFILE,99001) (Z(M1,M2),M1=1,MIN(MMAX,7))
      END DO
      WRITE(IFILE,99000)
99000 FORMAT(10X,60(1H-))
99001 FORMAT(10X,7F10.6)
      END
