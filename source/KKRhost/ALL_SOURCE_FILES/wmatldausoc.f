      SUBROUTINE WMATLDAUSOC(NTLDAU,ITLDAU,NSPIN,DENMATN,LOPT,
     &                    UEFF,JEFF,ULDAU,WLDAU,EU,EDC,MMAXD)
C **********************************************************************
C *                                                                    *
C * Calculation of Coulomb interaction potential in LDA+U              *
C * relativistic + SOC (new solver)                                    *
C *                                                                    *
C * The expression evaluated (array VLDAU) is                          *
C *                                                                    *
C *       V_{m1,s,m2,s'} =                                             *
C * delta_{ss'} Sum_{s'',m3,m4} U_{m1,m2,m3,m4} n_{m3,s'',m4,s''}      *
C * - Sum_{m3,m4} U_{m1,m4,m3,m2} n_{m3,s',m4,s}                       *
C * - [Ueff (dentot-1/2) - Jeff (n_s - 1/2)] delta_{ss'} delta_{m1,m2} *
C *                                                                    *
C * details see H. Ebert at al., Sol. Stat. Comm. 127 (2003) 443        *
C *                                                                    *
C *                  n.long,  April 2016, Juelich                      *
C **********************************************************************
      IMPLICIT NONE
      INCLUDE 'inc.p'
C
C PARAMETER definitions
C
      DOUBLE COMPLEX CZERO
      PARAMETER (CZERO=(0.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER NTLDAU,NSPIN,MMAXD
      INTEGER ITLDAU(NATYPD),LOPT(NATYPD)
      DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD),
     &                 UEFF(NATYPD),JEFF(NATYPD),EDC(NATYPD),EU(NATYPD),
     &                 WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
      DOUBLE COMPLEX DENMATN(MMAXD,MMAXD,2,2,NATYPD)
C
C Local variables
C
      DOUBLE COMPLEX CSUM,CSUM2,VLDAU(MMAXD,MMAXD,2,2)
      DOUBLE PRECISION DENMAT(MMAXD,MMAXD,2,2),DENTOT,
     &                 DENTOTS(NSPIND),FACTOR
      INTEGER I1,IT,IS,JS,M1,M2,M3,M4,MM,MMAX
      INTEGER IPRINT
      CHARACTER*15 STR15
C     ..
      DATA IPRINT /1/
      DATA FACTOR /1.D0/  ! if this is 1. then: n*(n-1) in Edc and potential
                          ! if this is 0. then: n**2    in Edc and potential
C    ..



      WRITE (1337,'(/,79(1H#),/,16X,A,/,79(1H#))') 
     &     'LDA+U: Calculating interaction potential VLDAU'
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      DO IT = 1,NTLDAU
         I1 = ITLDAU(IT)
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         IF ( LOPT(I1).GE.0 ) THEN
            CALL RINIT(MMAXD*MMAXD*2*2,DENMAT(1,1,1,1))
            MMAX = 2*LOPT(I1) + 1
            WRITE (1337,99001) I1,LOPT(I1)
C
C Result is in real Ylm basis. 
C It must be converted to complex Ylm basis:
C
            IF ( IPRINT.GT.1 ) 
     &           WRITE (1337,99002) 'Occupation matrix in REAL basis:'
C ----------------------------------------------------------------------
            DO IS = 1,NSPIN
              IF ( IPRINT.GT.1 ) THEN
                 WRITE (STR15,'(4X,"> ",A,I1)') 'ISPIN = ',IS
                 CALL CMATSTR(STR15,15,DENMATN(1,1,IS,IS,I1),
     &                                    MMAXD,MMAX,0,0,0,1d-8,1337)
              END IF
C
C -> Convert DENMATC and DENMAT to complex spherical harmonics.
C
             DO JS = 1,NSPIN
               CALL RCLM(1,LOPT(I1),LMAXD,DENMATN(1,1,JS,IS,I1))
             END DO ! JS
            END DO ! IS
C ----------------------------------------------------------------------
            IF ( IPRINT.GT.1 ) 
     &        WRITE (1337,99002) 'Occupation matrix in COMPLEX basis:'
            DENTOT = 0.D0
C ----------------------------------------------------------------------
            DO IS = 1,NSPIN
               IF ( IPRINT.GT.1 ) THEN
                  WRITE (STR15,'(4X,"> ",A,I1)') 'ISPIN = ',IS
                  CALL CMATSTR(STR15,15,DENMATN(1,1,IS,IS,I1),
     &                                      MMAXD,MMAX,0,0,0,1d-8,1337)
               END IF
C
C -> DENMAT is real: (imag(denmatc))
C
             DO JS = 1,NSPIN
               DO M2 = 1,MMAX
                  DO M1 = 1,MMAX
                     DENMAT(M1,M2,JS,IS) = (DENMATN(M1,M2,JS,IS,I1))
                  END DO
                ENDDO
               ENDDO ! JS
              ENDDO ! IS
C
C 2.  Calculate total occupation numbers:
C ntot_s = Sum_m n_{m,s,m,s}, ntot = n_1 + n_2
C
             DO IS=1,NSPIN
              DENTOTS(IS) = 0.D0
               DO JS=1, NSPIN
                DO MM = 1,MMAX
                 DENTOTS(IS) = DENTOTS(IS) + DENMAT(MM,MM,JS,IS)
                END DO
               ENDDO ! JS
              DENTOT = DENTOT + DENTOTS(IS)
             END DO ! IS
C ----------------------------------------------------------------------
            IF ( IPRINT.GT.0 ) THEN
               WRITE (1337,99002) 'Occupation matrix (real):'
               DO IS=1,NSPIN
                  WRITE(1337,99003) IS
                  CALL RWRITE(DENMAT(1,1,IS,IS),MMAXD,MMAX,1337)
                  WRITE(1337,99004) 'Trace     =',DENTOTS(IS)
               END DO
                  WRITE(1337,99005) 'Spins sum =',DENTOT
            END IF
C
C ----------------------------------------------------------------------
            CALL CINIT(MMAXD*MMAXD*2*2,VLDAU(1,1,1,1))
            DO IS = 1,NSPIN
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
                           DO JS = 1,NSPIN
                              CSUM2 = CSUM2 + DENMAT(M3,M4,JS,JS)
                           END DO
                           CSUM = CSUM + ULDAU(M1,M2,M3,M4,I1)*CSUM2
                        END DO
                     END DO
                     VLDAU(M1,M2,IS,IS) = VLDAU(M1,M2,IS,IS) + CSUM
                  END DO
               END DO
C
C 3b. Second part
C
             DO JS = 1,NSPIN
               DO M2 = 1,MMAX
                  DO M1 = 1,MMAX
                     CSUM = CZERO
                     DO M4 = 1,MMAX
                        DO M3 = 1,MMAX
                           CSUM = CSUM 
     &                          - ULDAU(M1,M4,M3,M2,I1)*
     &                            DENMAT(M3,M4,JS,IS)
                        END DO
                     END DO
                     VLDAU(M1,M2,JS,IS) = VLDAU(M1,M2,JS,IS) + CSUM
                  END DO
               END DO
              ENDDO ! JS
C
C 3c. Third part (always spin- and m-diagonal).
C
                DO M1 = 1,MMAX
                  VLDAU(M1,M1,IS,IS) = VLDAU(M1,M1,IS,IS) 
     &                          - UEFF(I1)*(DENTOT-0.5D0*FACTOR)
     &                          + JEFF(I1)*(DENTOTS(IS)-0.5D0*FACTOR   )
                END DO
            END DO ! IS
         
C
C 4. Calculate total-energy corrections EU and EDC (double-counting).
C Then the correction is EU - EDC.
C L[LDA+U]=E[LDA]+E[U]-E[DC]
C Note: EU,EDC initialised outside the routine

! Calculate EDC
            DO IS=1,NSPIN
             EDC(I1) = EDC(I1) + 
     &                   JEFF(I1)*DENTOTS(IS)*(DENTOTS(IS)-FACTOR)
            ENDDO

            EDC(I1) = 0.5D0*(UEFF(I1)*DENTOT*(DENTOT-1.D0)-EDC(I1))

! Calculate EU
            DO IS = 1,NSPIN
             DO JS = 1,NSPIN
              DO M4 = 1,MMAX
               DO M3 = 1,MMAX
                DO M2 = 1,MMAX
                 DO M1 = 1,MMAX
                  EU(I1) = EU(I1) + DENMAT(M1,M2,IS,IS)*
     &                       ULDAU(M1,M2,M3,M4,I1)*DENMAT(M3,M4,JS,JS)
                 ENDDO
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ENDDO

            DO IS = 1,NSPIN
             DO JS = 1,NSPIN
              DO M4 = 1,MMAX
               DO M3 = 1,MMAX
                DO M2 = 1,MMAX
                 DO M1 = 1,MMAX
                  EU(I1) = EU(I1) - DENMAT(M1,M2,IS,JS)*
     &                        ULDAU(M1,M4,M3,M2,I1)*DENMAT(M3,M4,JS,IS)
                 ENDDO
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ENDDO

            EU(I1) = 0.5D0*EU(I1)



C ----------------------------------------------------------------------
            IF ( IPRINT.GT.0 ) WRITE (1337,99002)
     &           'Interaction potential in COMPLEX basis:'
C ----------------------------------------------------------------------
            DO IS = 1,NSPIN
             WLDAU(:,:,IS,I1)=0d0
               IF ( IPRINT.GT.0 ) THEN
                  WRITE (STR15,'(4X,"> ",A,I1)') 'ISPIN = ',IS
                  CALL CMATSTR(STR15,15,VLDAU(1,1,IS,IS),
     &                         MMAXD,MMAX,0,0,0,1d-8,1337)
               END IF
C
C 5.  Transform VLDAU into real spherical harmonics basis
C
              DO JS=1,NSPIN
               CALL RCLM(2,LOPT(I1),LMAXD,VLDAU(1,1,JS,IS))
C
C Copy transformed VLDAU to real WLDAU
C
C Apply damping to the interaction matrix WLDAU ? Here not.
C
C
               DO M2 = 1,MMAX
                  DO M1 = 1,MMAX
                     WLDAU(M1,M2,IS,I1) = 
     +                 WLDAU(M1,M2,IS,I1)+DREAL(VLDAU(M1,M2,JS,IS))
                  END DO
               END DO
             ENDDO ! JS
            END DO ! IS
C ----------------------------------------------------------------------
            IF ( IPRINT.GT.0 ) THEN
               WRITE (1337,99002) 'Interaction potential in REAL basis:'
               DO IS = 1,NSPIN
                  WRITE (STR15,'(4X,"> ",A,I1)') 'ISPIN = ',IS
                  CALL CMATSTR(STR15,15,VLDAU(1,1,IS,IS),
     &                 MMAXD,MMAX,0,0,0,1d-8,1337)
               END DO
            END IF
C ----------------------------------------------------------------------
            WRITE (1337,99002) 'Interaction potential (real):'
            DO IS = 1,NSPIN
               WRITE(1337,99003) IS
               CALL RWRITE(WLDAU(1,1,IS,I1),MMAXD,MMAX,1337)
            END DO
            WRITE(1337,*)
C ----------------------------------------------------------------------
C
C Corrections in total energy:
C
C
C -> Write out corrections on energy:
C    E[LDA+U] = E[LDA] + EU - EDC
C
            WRITE(1337,99002) 'Corrections to the total energy:'
            WRITE(1337,*)
            WRITE(1337,99004) 'EU  =',EU(I1)
            WRITE(1337,99004) 'Edc =',EDC(I1)
            WRITE(1337,99006) 'E[LDA+U] = E[LDA] + EU - Edc'
         END IF
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      END DO                    ! I1 = 1,NTLDAU
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
99001 FORMAT(/,6X,65(1H=),/,6X,'Atom :',I3,' (l =',I2,')',/,6X,18(1H=))
99002 FORMAT(8X,'* ',A)
99003 FORMAT(/,15X,'> ISPIN =',I1)
99004 FORMAT(10X,A,F10.6)
99005 FORMAT(10X,21(1H-),/,10X,A,F10.6,/,10X,60(1H-),/)
99006 FORMAT(27X,A,/)
      END
