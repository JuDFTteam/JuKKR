SUBROUTINE wmatldau(ntldau,itldau,nspin,denmatc,lopt,  &
    ueff,jeff,uldau,wldau,eu,edc,mmaxd,npotd)
! **********************************************************************
! *                                                                    *
! * Calculation of Coulomb interaction potential in LDA+U              *
! * non-relativistic case -- otherwise matrices DENMAT and VLDAU must  *
! *                          have double dimension                     *
! *                                                                    *
! * Uses the Coulomb matrix U (array ULDAU), the density matrix n      *
! * (array DENMAT) and the occupation numbers dentot (total) and n_s   *
! * (array DENTOTS) (per spin).                                        *
! *                                                                    *
! * The expression evaluated (array VLDAU) is                          *
! *                                                                    *
! *       V_{m1,s,m2,s'} =                                             *
! * delta_{ss'} Sum_{s'',m3,m4} U_{m1,m2,m3,m4} n_{m3,s'',m4,s''}      *
! * - Sum_{m3,m4} U_{m1,m4,m3,m2} n_{m3,s',m4,s}                       *
! * - [Ueff (dentot-1/2) - Jeff (n_s - 1/2)] delta_{ss'} delta_{m1,m2} *
! *                                                                    *
! *                  ph. mavropoulos, h.ebert munich/juelich 2002-2004 *
! **********************************************************************
IMPLICIT NONE
INCLUDE 'inc.p'

! PARAMETER definitions

DOUBLE COMPLEX CZERO
PARAMETER (CZERO=(0.0D0,0.0D0))

! Dummy arguments

INTEGER NTLDAU,NSPIN,MMAXD,NPOTD
INTEGER ITLDAU(NATYPD),LOPT(NATYPD)
DOUBLE PRECISION &
                 UEFF(NATYPD),JEFF(NATYPD),EDC(NATYPD),EU(NATYPD), &
                 WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
DOUBLE PRECISION, allocatable :: ULDAU(:,:,:,:,:) 
DOUBLE COMPLEX DENMATC(MMAXD,MMAXD,NPOTD)

! Local variables
DOUBLE COMPLEX CSUM,CSUM2,VLDAU(MMAXD,MMAXD,NSPIND)
DOUBLE PRECISION DENMAT(MMAXD,MMAXD,NSPIND),DENTOT, &
                 DENTOTS(NSPIND)
INTEGER I1,IT,IPOT,IS,JS,M1,M2,M3,M4,MM,MMAX
INTEGER IPRINT
CHARACTER*15 STR15
!..
DATA IPRINT /1/
!    ..

WRITE (1337,'(/,79(1H#),/,16X,A,/,79(1H#))')  &
    'LDA+U: Calculating interaction potential VLDAU'

allocate( uldau(mmaxd,mmaxd,mmaxd,mmaxd,natypd) )

! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
DO it = 1,ntldau
  i1 = itldau(it)
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
  IF ( lopt(i1) >= 0 ) THEN
    CALL rinit(mmaxd*mmaxd*nspind,denmat(1,1,1))
    mmax = 2*lopt(i1) + 1
    WRITE (1337,99001) i1,lopt(i1)
    
! Result is in real Ylm basis.
! It must be converted to complex Ylm basis:
    
    IF ( iprint > 1 ) WRITE (1337,99002) 'Occupation matrix in REAL basis:'
! ----------------------------------------------------------------------
    DO is = 1,nspin
      ipot = (i1-1)*nspin + is
      IF ( iprint > 1 ) THEN
        WRITE (str15,'(4X,"> ",A,I1)') 'ISPIN = ',is
        CALL cmatstr(str15,15,denmatc(1,1,ipot), mmaxd,mmax,0,0,0,1D-8,6)
      END IF
      
! -> Convert DENMATC and DENMAT to complex spherical harmonics.
      
      CALL rclm(1,lopt(i1),lmaxd,denmatc(1,1,ipot))
    END DO
! ----------------------------------------------------------------------
    IF ( iprint > 1 )  &
        WRITE (1337,99002) 'Occupation matrix in COMPLEX basis:'
    dentot = 0.d0
! ----------------------------------------------------------------------
    DO is = 1,nspin
      ipot = (i1-1)*nspin + is
      IF ( iprint > 1 ) THEN
        WRITE (str15,'(4X,"> ",A,I1)') 'ISPIN = ',is
        CALL cmatstr(str15,15,denmatc(1,1,ipot), mmaxd,mmax,0,0,0,1D-8,6)
      END IF
      
! -> DENMAT is real: (imag(denmatc))
      
      DO m2 = 1,mmax
        DO m1 = 1,mmax
          denmat(m1,m2,is) = DIMAG(denmatc(m1,m2,ipot))
        END DO
      END DO
      
! 2.  Calculate total occupation numbers:
! ntot_s = Sum_m n_{m,s,m,s}, ntot = n_1 + n_2
      
      dentots(is) = 0.d0
      DO mm = 1,mmax
        dentots(is) = dentots(is) + denmat(mm,mm,is)
      END DO
      dentot = dentot + dentots(is)
    END DO
! ----------------------------------------------------------------------
    IF ( iprint > 0 ) THEN
      WRITE (1337,99002) 'Occupation matrix (real):'
      DO is = 1,nspin
        WRITE(1337,99003) is
        CALL rwrite(denmat(1,1,is),mmaxd,mmax,1337)
        WRITE(1337,99004) 'Trace     =',dentots(is)
      END DO
      WRITE(1337,99005) 'Spins sum =',dentot
    END IF
    
! In paramagnetic case the spin degeneracy has been accounted
! for by the weight DF in tmatrho.
    
! ----------------------------------------------------------------------
    CALL cinit(mmaxd*mmaxd*nspind,vldau(1,1,1))
    DO is = 1,nspin
      
! 3.  Use density matrix and Coulomb matrix ULDAU to calculate the
! interaction potential VLDAU
! 3a. First part (always diagonal in spin).
      
      DO m2 = 1,mmax
        DO m1 = 1,mmax
          csum = czero
          DO m4 = 1,mmax
            DO m3 = 1,mmax
              csum2 = czero
              DO js = 1,nspin
                csum2 = csum2 + denmat(m3,m4,js)
              END DO
              csum = csum + uldau(m1,m2,m3,m4,i1)*csum2
            END DO
          END DO
          vldau(m1,m2,is) = vldau(m1,m2,is) + csum
        END DO
      END DO
      
! 3b. Second part (in fully rel. case not diagonal in spin; then this
! loop must be changed accordingly).
      
      DO m2 = 1,mmax
        DO m1 = 1,mmax
          csum = czero
          DO m4 = 1,mmax
            DO m3 = 1,mmax
              csum = csum - uldau(m1,m4,m3,m2,i1)*denmat(m3,m4,is)
            END DO
          END DO
          vldau(m1,m2,is) = vldau(m1,m2,is) + csum
        END DO
      END DO
      
! 3c. Third part (always spin- and m-diagonal).
      
      DO m1 = 1,mmax
        vldau(m1,m1,is) = vldau(m1,m1,is) - ueff(i1)*(dentot-0.5D0)  &
            + jeff(i1)*(dentots(is)-0.5D0)
      END DO
      
! 4. Calculate total-energy corrections EU and EDC (double-counting).
! Then the correction is EU-EDC.
! Note: EU,EDC initialised outside the routine
      
! Here VLDAU is assumed spin-diagonal (contrary to the spin-orbit case).
      
      DO m2 = 1,mmax
        DO m1 = 1,mmax
          eu(i1) = eu(i1) + denmat(m1,m2,is) * dreal(vldau(m1,m2,is))
        END DO
      END DO
      edc(i1) = edc(i1) + jeff(i1)*dentots(is)*(dentots(is)-1.d0)
    END DO
! ----------------------------------------------------------------------
    IF ( iprint > 0 ) WRITE (1337,99002)  &
        'Interaction potential in COMPLEX basis:'
! ----------------------------------------------------------------------
    DO is = 1,nspin
      IF ( iprint > 0 ) THEN
        WRITE (str15,'(4X,"> ",A,I1)') 'ISPIN = ',is
        CALL cmatstr(str15,15,vldau(1,1,is), mmaxd,mmax,0,0,0,1D-8,6)
      END IF
      
! 5.  Transform VLDAU into real spherical harmonics basis
      
      CALL rclm(2,lopt(i1),lmaxd,vldau(1,1,is))
      
! Copy transformed VLDAU to real WLDAU
      
! Apply damping to the interaction matrix WLDAU ? Here not.
      
      
      DO m2 = 1,mmax
        DO m1 = 1,mmax
          wldau(m1,m2,is,i1) = dreal(vldau(m1,m2,is))
        END DO
      END DO
      
    END DO
! ----------------------------------------------------------------------
    IF ( iprint > 0 ) THEN
      WRITE (1337,99002) 'Interaction potential in REAL basis:'
      DO is = 1,nspin
        WRITE (str15,'(4X,"> ",A,I1)') 'ISPIN = ',is
        CALL cmatstr(str15,15,vldau(1,1,is), mmaxd,mmax,0,0,0,1D-8,6)
      END DO
    END IF
! ----------------------------------------------------------------------
    WRITE (1337,99002) 'Interaction potential (real):'
    DO is = 1,nspin
      WRITE(1337,99003) is
      CALL rwrite(wldau(1,1,is,i1),mmaxd,mmax,1337)
    END DO
    WRITE(1337,*)
! ----------------------------------------------------------------------
    
! Corrections in total energy:
    
    eu(i1) = 0.5D0*eu(i1)
    edc(i1) = 0.5D0*(ueff(i1)*dentot*(dentot-1.d0)-edc(i1))
    
! -> Write out corrections on energy:
!    E[LDA+U] = E[LDA] + EU - EDC
    
    WRITE(1337,99002) 'Corrections to the total energy:'
    WRITE(1337,*)
    WRITE(1337,99004) 'EU  =',eu(i1)
    WRITE(1337,99004) 'Edc =',edc(i1)
    WRITE(1337,99006) 'E[LDA+U] = E[LDA] + EU - Edc'
  END IF
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
END DO                    ! I1 = 1,NTLDAU
! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
99001 FORMAT(/,6X,65(1H=),/,6X,'Atom :',i3,' (l =',i2,')',/,6X,18(1H=))
99002 FORMAT(8X,'* ',a)
99003 FORMAT(/,15X,'> ISPIN =',i1)
99004 FORMAT(10X,a,f10.6)
99005 FORMAT(10X,21(1H-),/,10X,a,f10.6,/,10X,60(1H-),/)
99006 FORMAT(27X,a,/)
END SUBROUTINE wmatldau
!*==rwrite.f    processed by SPAG 6.05Rc at 16:58 on 22 Dec 2004

SUBROUTINE rwrite(z,mmaxd,mmax,ifile)

DOUBLE PRECISION, INTENT(IN OUT)         :: z(mmaxd,mmaxd)
INTEGER, INTENT(IN OUT)                  :: mmaxd
INTEGER, INTENT(IN)                      :: mmax
INTEGER, INTENT(IN OUT)                  :: ifile
IMPLICIT NONE


INTEGER :: m1,m2

WRITE(ifile,99000)

DO m2 = 1,mmax
  WRITE (ifile,99001) (z(m1,m2),m1=1,MIN(mmax,7))
END DO
WRITE(ifile,99000)
99000 FORMAT(10X,60(1H-))
99001 FORMAT(10X,7F10.6)
END SUBROUTINE rwrite
