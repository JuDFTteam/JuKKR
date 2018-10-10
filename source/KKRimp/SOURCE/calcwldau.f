      MODULE MOD_CALCWLDAU
!-------------------------------------------------------------------------------
!> Summary: Calculation of Coulomb interaction potential in LDA+U in the
!> non-relativistic case
!> Author:
!> Category: KKRimp, single-site, electrostatics, potential, lda+u 
!>           
!-------------------------------------------------------------------------------
      CONTAINS

      SUBROUTINE CALCWLDAU(
     >     NSPIN,NATOM,LMAXD,IRMD,LMAXATOM,DENSITY,STRMIX,
     X     LDAU)

C **********************************************************************
C
C Calculation of Coulomb interaction potential in LDA+U
C Non-relativistic case. (Otherwise matrices DENMAT and VLDAU must
C have double dimension)
C
C Uses the Coulomb matrix U (array ULDAU), the density matrix n
C (array DENMATC) and the occupation numbers dentot (total) and n_s (array
C DENTOTS) (per spin).
C
C The expression evaluated (array VLDAU) is
C
C V_{m1,s,m2,s'} = 
C   delta_{ss'} Sum_{s'',m3,m4} U_{m1,m2,m3,m4} n_{m3,s'',m4,s''}
C - Sum_{m3,m4} U_{m1,m4,m3,m2} n_{m3,s',m4,s}
C - [Ueff (dentot-1/2) - Jeff (n_s - 1/2)] delta_{ss'} delta_{m1,m2}
C
C This is expressed in complex spherical harmonics basis. Then it is
C transformed into the real spherical harmonics basis via subr. rclm.
C
!    Old method:
C The density matrix n is calculated using the Green function and the
C reference functions Phi as
C
C n_{m,s,m',s'} = -1/pi Int dE 
C [ Sum_{L''L'''} (Phi_L,R_L'') G_{L''L'''}(E) (R_L''',Phi_L') +
C   Sum_L'' (Phi_L,R_L'')(H_L'',Phi_L')                          ]
C
C
!     New method:
! The density matrix is obtained by the Green-function matrix elements
! integrated in energy up to EF: n_{m,m',s} = (G-G^{+})_{mm',ss}  / 2i    
!
C
C **********************************************************************
C
C PARAMETER definitions
C
      USE TYPE_LDAU
      USE NRTYPE
      USE TYPE_DENSITY
      USE MOD_RWLDAUPOT
      USE MOD_RCLM
      implicit none
      INTEGER LMAXD,MMAXD,IRMD
      DOUBLE COMPLEX CZERO,CONE,CI
      PARAMETER (CZERO=(0.0D0,0.0D0),CONE=(1.D0,0.D0),CI=(0.D0,1.D0))
      REAL*8 DZERO
      PARAMETER (DZERO=0.0D0)
C
C Dummy arguments
C
      INTEGER NATOM,NSPIN
      INTEGER LMAXATOM(:)
      TYPE(LDAU_TYPE),ALLOCATABLE :: LDAU(:) ! lda+u variables dimension: (NATOM)
      TYPE(DENSITY_TYPE) :: DENSITY(:)

C
C Local variables
C

      REAL*8,ALLOCATABLE  ::  WLDAU_OLD(:,:,:),RMS(:)   ! for mixing purposes
      REAL*8 DENTOT,DENTOTS(2),RMSTOT
      COMPLEX*16 CSUM,CSUM2
      COMPLEX*16,ALLOCATABLE :: VLDAU(:,:)
      REAL*8 XMIX,STRMIX  ! mixing factor, input straight-mixing factor
      INTEGER IRUNLDAU, IAT
     &     ,L1,M1,M2,M3,M4,MM,IS,JS,MMAX
     &     ,LM1,LMLO,LMHI


!      SAVE 

!      EXTERNAL CINIT,RCLM,RWLDAUPOT

      MMAXD = 2*LMAXD + 1
C--------------------------------------------------------------------
      ALLOCATE( WLDAU_OLD(MMAXD,MMAXD,NSPIN) )
      ALLOCATE( VLDAU(MMAXD,MMAXD) )
      ALLOCATE( RMS(NATOM) )

      RMSTOT = 0.D0


      DO 200 IAT = 1,NATOM

      LDAU(IAT)%EDC = 0.D0
      LDAU(IAT)%EU = 0.D0
 
      L1 = LDAU(IAT)%LOPT

      IF (L1.GE.0) THEN         

         MMAX = 2*L1 + 1

! Save old potential for mixing later
         WLDAU_OLD(1:MMAX,1:MMAX,1:NSPIN) =  
     &           LDAU(IAT)%WLDAU(1:MMAX,1:MMAX,1:NSPIN)


          WRITE(6,*) 'Atom Number:',IAT

! Calculate density matrix per spin direction for this atom.
! Remember, GFINT is always allocated as 2x2 matrix in spin-space
          LMLO = L1**2 + 1
          LMHI = (L1 + 1)**2
          DO M1 = 1,MMAX
             LM1 = LMLO - 1 + M1
             LDAU(IAT)%DENMATC(1:MMAX,M1,1) = 
     &               DENSITY(IAT)%GFINT(LMLO:LMHI,LM1) - 
     &       DCONJG( DENSITY(IAT)%GFINT(LM1,LMLO:LMHI) )      ! ( G-G^{+} ), down-down
          ENDDO
          IF (NSPIN.EQ.2) THEN
             LMLO = (LMAXATOM(IAT) + 1)**2 + L1**2 + 1
             LMHI = (LMAXATOM(IAT) + 1)**2 + (L1 + 1)**2
             DO M1 = 1,MMAX
                LM1 = LMLO - 1 + M1
                LDAU(IAT)%DENMATC(1:MMAX,M1,2) = 
     &                       DENSITY(IAT)%GFINT(LMLO:LMHI,LM1) - 
     &               DCONJG( DENSITY(IAT)%GFINT(LM1,LMLO:LMHI) ) ! ( G-G^{+} ) , up-up
             ENDDO
          ENDIF

          LDAU(IAT)%DENMATC(1:MMAX,1:MMAX,1:NSPIN) = 
     &         LDAU(IAT)%DENMATC(1:MMAX,1:MMAX,1:NSPIN) / 2.D0 / CI   ! ( G-G^{+} ) / 2i
          ! Factor (-1/pi) has been included in integr. weight

C Result is in real Ylm basis. It must be converted to complex Ylm basis:
          WRITE(6,*) 'Occupation matrix (complex) in real basis:'
          DO IS = 1,NSPIN
              WRITE(6,*) 'IS=',IS
              CALL ZWRITE(LDAU(IAT)%DENMATC(:,:,IS),MMAX,MMAX,6)
          ENDDO
C Convert DENMATC to complex spherical harmonics.
          DO IS = 1,NSPIN
             CALL RCLM(1,L1,LMAXD,LDAU(IAT)%DENMATC(1,1,IS))
          ENDDO
          WRITE(6,*) 'Occupation matrix (complex) in complex basis:'
          DO IS = 1,NSPIN
              WRITE(6,*) 'IS=',IS
              CALL ZWRITE(LDAU(IAT)%DENMATC(:,:,IS),MMAX,MMAX,6)
          ENDDO



C 2.  Calculate total occupation numbers: 
C ntot_s = Sum_m n_{m,s,m,s}, ntot = n_1 + n_2
C

          DENTOT = 0.D0
          DO IS = 1,NSPIN
              DENTOTS(IS) = 0.D0
              DO MM = 1,MMAX
                 DENTOTS(IS) = DENTOTS(IS) + LDAU(IAT)%DENMATC(MM,MM,IS)
              ENDDO
              DENTOT = DENTOT + DENTOTS(IS)
              DENTOTS(IS) = DENTOTS(IS) / (3-NSPIN) ! If nspin=1, denmatc is spin-summed
          ENDDO                                     ! but dentots should be per spin.
          WRITE(6,FMT='(A7,F8.5,A11,2F8.5)') 
     &         'DENTOT=',DENTOT,'  DENTOTS(IS)=',(DENTOTS(IS),IS=1,2)
C In paramagnetic case the spin degeneracy has been accounted
C for by the weight DF.


          DO 100 IS = 1,NSPIN

          VLDAU(1:MMAXD,1:MMAXD) = CZERO

C 3.  Use density matrix and Coulomb matrix ULDAU to calculate the
C interaction potential VLDAU
C 3a. First part (always diagonal in spin).
          DO M1 = 1,MMAX
              DO M2 = 1,MMAX
                  CSUM = CZERO
                  DO M4 = 1,MMAX
                      DO M3 = 1,MMAX
                          CSUM2 = CZERO
                          DO JS = 1,NSPIN
                             CSUM2 = CSUM2 + LDAU(IAT)%DENMATC(M3,M4,JS)
                          ENDDO
                          CSUM = CSUM 
     &                         + LDAU(IAT)%ULDAU(M1,M2,M3,M4) * CSUM2
                      ENDDO
                  ENDDO
                  VLDAU(M1,M2) = VLDAU(M1,M2) + CSUM
              ENDDO
          ENDDO
C 3b. Second part (in fully rel. case not diagonal in spin; then this
C loop must be changed accordingly).
          DO M2 = 1,MMAX
              DO M1 = 1,MMAX
                  CSUM = CZERO
                  DO M4 = 1,MMAX
                     DO M3 = 1,MMAX
                        CSUM = CSUM - LDAU(IAT)%ULDAU(M1,M4,M3,M2) 
     &                   * LDAU(IAT)%DENMATC(M3,M4,IS) / DFLOAT(3-NSPIN) ! If nspin=1, denmatc is spin-summed
                     ENDDO                                         ! but here we need it per spin.
                  ENDDO
                  VLDAU(M1,M2) = VLDAU(M1,M2) + CSUM
              ENDDO
          ENDDO
C 3c. Third part (always spin- and m-diagonal).
          DO M1 = 1,MMAX
              VLDAU(M1,M1) = VLDAU(M1,M1)
     &             - LDAU(IAT)%UEFF*(DENTOT-0.5D0) 
     &             + LDAU(IAT)%JEFF*(DENTOTS(IS)-0.5D0)
          ENDDO

C 4. Calculate total-energy corrections EU and EDC (double-counting).
C Then the correction is EU-EDC.
C Here VLDAU is assumed spin-diagonal (contrary to the spin-orbit case).
          DO M2 = 1,MMAX
              DO M1 = 1,MMAX
                  LDAU(IAT)%EU = LDAU(IAT)%EU 
     &               + LDAU(IAT)%DENMATC(M1,M2,IS) * DREAL(VLDAU(M1,M2))
              ENDDO
          ENDDO

          LDAU(IAT)%EDC = LDAU(IAT)%EDC - 
     &                    DENTOTS(IS) * (DENTOTS(IS) - 1.D0)



          WRITE(6,*) 'Interaction potential (complex) for spin:',IS
          WRITE(6,*) '(in complex basis)'
          CALL ZWRITE (VLDAU,MMAXD,MMAX,6)

C 5.  Transform VLDAU into real spherical harmonics basis
          CALL RCLM(2,L1,LMAXD,VLDAU)

          WRITE(6,*) 'Interaction potential (complex) for spin:',IS
          WRITE(6,*) '(in real basis)'
          CALL ZWRITE (VLDAU,MMAXD,MMAX,6)

C Copy transformed VLDAU
          DO M2 = 1,MMAX
              DO M1 = 1,MMAX
                  LDAU(IAT)%WLDAU(M1,M2,IS) = DREAL(VLDAU(M1,M2))
              ENDDO
          ENDDO

 100      ENDDO                 ! IS = 1,NSPIN


C Corrections in total energy:
          LDAU(IAT)%EU = 0.5D0 * LDAU(IAT)%EU 
          LDAU(IAT)%EDC = 0.5D0 
     &     * (LDAU(IAT)%UEFF * DENTOT * (DENTOT - 1.D0) - LDAU(IAT)%EDC)


C Write out corrections on energy:
      WRITE(*,*) 'Correction to energy for atom ',IAT
      WRITE(*,*) 'EU = ',LDAU(IAT)%EU,' EDC= ',LDAU(IAT)%EDC
      WRITE(*,*) 'Total energy in LDA+U should be calculated as:'
      WRITE(*,*) 'E[LDA+U] = E[LDA] + EU - EDC'
      WRITE(*,*) 'only for atoms treated with LDA+U.'

C Calculate rms error for wldau
      RMS(IAT) = 0.D0
      DO IS = 1,NSPIN
         DO M2 = 1,MMAX
            DO M1 = 1,MMAX
               RMS(IAT) = RMS(IAT) + 
     &              (LDAU(IAT)%WLDAU(M1,M2,IS) - WLDAU_OLD(M1,M2,IS))**2
            ENDDO
         ENDDO
      ENDDO
      RMSTOT = RMSTOT + RMS(IAT)

C Apply damping to the interaction matrix WLDAU (mixing):
      XMIX =STRMIX
      write(*,*) 'Mixing old-new wldau with xmix=',xmix
      CALL WMIX(XMIX,LDAU(IAT)%WLDAU(:,:,:),WLDAU_OLD(:,:,:),MMAX,NSPIN)
      
      ENDIF                     !  (L1.GE.0) 
 200  ENDDO                     ! IAT = 1,NATOM

C 6.  Write result on disk (to be used in the next iteration)

      IRUNLDAU = 1
      CALL RWLDAUPOT(.TRUE.,4,NATOM,NSPIN,LMAXD,IRMD,IRUNLDAU,LDAU)


C 7. Write out rms error
      DO IAT = 1,NATOM
         IF (LDAU(IAT)%LOPT.GE.0) 
     &    WRITE(*,*) 'RMS-ERROR for LDA+U pot. atom',IAT,DSQRT(RMS(IAT))
      ENDDO
      WRITE(*,*) 'TOTAL RMS-ERROR for LDA+U:',DSQRT(RMSTOT)


      DEALLOCATE( WLDAU_OLD )
      DEALLOCATE( VLDAU )
      DEALLOCATE( RMS )

      RETURN

      END SUBROUTINE CALCWLDAU

!*********************************************************************
      SUBROUTINE RWRITE(Z,MMAXD,MMAX,IFILE)
      implicit none
      INTEGER MMAXD,MMAX,M1,M2,IFILE
      REAL*8 Z(MMAXD,MMAXD)
      DO M2=1,MMAX
          WRITE(IFILE,9000) (Z(M1,M2),M1=1,MMAX)
      ENDDO
 9000 FORMAT(20F8.4)
      RETURN
      END SUBROUTINE RWRITE
!*********************************************************************
      SUBROUTINE ZWRITE(Z,MMAXD,MMAX,IFILE)
      implicit none
      INTEGER MMAXD,MMAX,M1,M2,IFILE
      COMPLEX*16 Z(:,:)
      DO M2=1,MMAX
          WRITE(IFILE,9000) (Z(M1,M2),M1=1,MMAX)
      ENDDO
 9000 FORMAT(20F12.8)
      RETURN
      END SUBROUTINE ZWRITE


!*********************************************************************

      SUBROUTINE WMIX(XMIX,WLDAU,WLDAU_OLD,MMAXD,NSPIND)
c Mix old and new potential. Linear mixing with factor xmix.
      implicit none
      INTEGER MMAXD,NSPIND
      INTEGER M1,M2,IS!,IAT
      REAL*8 WLDAU_OLD(:,:,:) ! for mixing purposes
      REAL*8 WLDAU(:,:,:)
      REAL*8 XMIX,ONEMXMIX

      WRITE(*,*) 'Mixing old and new WLDAU with mixing coef.=',XMIX

      ONEMXMIX = 1.D0 - XMIX

!     DO IAT = 0,NATLDAUD
         DO IS = 1,NSPIND
            DO M2 = 1,MMAXD
               DO M1 = 1,MMAXD
                  WLDAU(M1,M2,IS) = XMIX * WLDAU(M1,M2,IS) +
     &                 ONEMXMIX * WLDAU_OLD(M1,M2,IS)
               ENDDO
            ENDDO
         ENDDO
!     ENDDO
      RETURN
      END SUBROUTINE WMIX

      END MODULE MOD_CALCWLDAU
