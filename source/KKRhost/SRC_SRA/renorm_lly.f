      SUBROUTINE RENORM_LLY(                                       ! LLY Lloyd
     >     CDOS_LLY,IELAST,NSPIN,NATYP,CDEN,LMAXP1,CONC,
     >     IESTART,IEEND,WEZ,IRCUT,IPAN,
     X     RHO2NS)

! Renormalize the valence charge according to Lloyd's formula.
! Find renormalization constant per energy, then renormalize 
! charge/atom/energy, then integrate over energies to find
! the renormalized charge/atom. Use it to renormalize the density.
! Phivos Mavropoulos, July 2014
      implicit none
      INCLUDE 'inc.p'
      INTEGER LMAXD1
      PARAMETER (LMAXD1=LMAXD+1)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER NPOTD
      PARAMETER (NPOTD= (2*(KREL+KORBIT) + 
     +           (1-(KREL+KORBIT))*NSPIND)*NATYPD)
! Input:
      INTEGER LMAXP1,NATYP,NSPIN,IRMAX         ! LMAX+1 for FP or LMAX for ASA, number of atoms,radial pts
      INTEGER IESTART,IEEND,IELAST             ! Starting and ending energy point for renormalization
      INTEGER IRCUT(0:IPAND,NATYPD),IPAN(NATYPD) ! Mesh info
      REAL*8 CONC(NATYPD)                      ! Concentration (for cpa)
      COMPLEX*16 CDEN(0:LMAXD1,IEMXD,NPOTD)     ! Non-renormalized density per atom (density=-cden/pi)
      COMPLEX*16 CDOS_LLY(IEMXD,NSPIND)        ! DOS according to Lloyd's formula
      COMPLEX*16 WEZ(IEMXD)
! Input/Output:
      REAL*8 RHO2NS(IRMD,LMPOTD,NATYPD,2)
! Internal:
      INTEGER LL,IE,I1,ISPIN,IPOT,SPINDEGEN,IRC1,SIGNSP
      REAL*8 RENORM_AT(NATYPD,2)                   ! 1: charge renormalization per atom (energy-integrated)
                                                   ! 2: same for spin moment
      COMPLEX*16 CDOS_LOC                          ! Density from local summation and from Lloyd's formula
      COMPLEX*16 CREN(IEMXD,2)                     ! Renormalization constant for charge and spin density
      REAL*8 CHARGE(NATYPD,2),CHARGE_LLY(NATYPD,2) ! Atomic charge per spin (local summation and renormalized)
      COMPLEX*16 CHADD(NATYPD),CDOS_ADD            ! Integration step for charge/atom/spin
      REAL*8 SUM0(2),SUM1(2)
      COMPLEX*16 CZERO
      REAL*8 PI

      CZERO = (0.D0,0.D0)
      PI = 4.D0 * DATAN(1.D0)

      SPINDEGEN = 3 - NSPIN ! Spin degeneracy, 2 if nspin=1, 1 if nspin=2

      CREN(:,:) = CZERO
      RENORM_AT(:,:) = 1.D0
      CHARGE_LLY(:,:) = 0.D0
      CHARGE(:,:) = 0.D0

! First find renormalization factor per energy and atomic charges
      DO IE = IESTART,IEEND
         DO ISPIN = 1,NSPIN
            CDOS_LOC = CZERO
            CHADD(1:NATYPD) = CZERO
            DO I1 = 1,NATYP
               IPOT = (I1-1) * NSPIN + ISPIN
               CDOS_ADD = CZERO
               DO LL = 0,LMAXP1
                  CDOS_ADD = CDOS_ADD + CONC(I1) * CDEN(LL,IE,IPOT)    ! Factor 1/pi included in WEZ 
               END DO
               CDOS_LOC = CDOS_LOC + CDOS_ADD
               CHADD(I1) = WEZ(IE) * CDOS_ADD                          ! Complex charge
               CHARGE(I1,ISPIN) = CHARGE(I1,ISPIN) + DIMAG( CHADD(I1) )
            ENDDO
            ! Now the locally-summed charge/energy is in cdos_loc, charge/energy/atom in chadd

            ! Renormalization factor per energy:
            CDOS_LOC = -CDOS_LOC / PI
            CREN(IE,ISPIN) = CDOS_LLY(IE,ISPIN) / CDOS_LOC

            ! Apply to DOS of each atom:
            CHARGE_LLY(1:NATYP,ISPIN) = CHARGE_LLY(1:NATYP,ISPIN) +  
     &               DIMAG( CHADD(1:NATYP) * CREN(IE,ISPIN) )

         ENDDO ! ISPIN = 1,NSPIN

         IF (LNC.OR.NSPIN.EQ.1) CREN(IE,2) = CREN(IE,1) ! If spins are coupled, then only charge density 
                                                        ! is given by Lloyd's formula
      ENDDO ! IE = IESTART,IEEND


! Now apply renormalization to energy-integrated charge density & spin density
      DO I1 = 1,NATYP

         RENORM_AT(I1,1) = ( CHARGE_LLY(I1,1) + CHARGE_LLY(I1,2) ) /      ! Charge renormalization
     &                         ( CHARGE(I1,1) + CHARGE(I1,2) ) 
         IF (NSPIN.EQ.2)
     &       RENORM_AT(I1,2) = ( CHARGE_LLY(I1,2) - CHARGE_LLY(I1,1) ) /  ! Spin renormalization
     &                             ( CHARGE(I1,1) + CHARGE(I1,2) ) 
         IF (LNC.OR.NSPIN.EQ.1) RENORM_AT(I1,2) = RENORM_AT(I1,2)


         ! Renormalize charge density (ispin=1,2 means here charge,spin density)
         IRC1 = IRCUT(IPAN(I1),I1) ! Index of outmost radial point
         DO ISPIN = 1,NSPIN
            RHO2NS(1:IRC1,1:LMPOTD,I1,ISPIN) = 
     &           RHO2NS(1:IRC1,1:LMPOTD,I1,ISPIN) * RENORM_AT(I1,ISPIN)
         ENDDO

      ENDDO

! Write out renormalization factors
      WRITE(*,*) 'Information on renormalization by Lloyds formula'
      WRITE(*,*)'RENORM_LLY: Complex renormalization factor per energy:'
      WRITE(*,FMT='(A5,2A32)') 'IE',
     &               'Spin 1 (down)           ','Spin 2 (up)           '
      DO IE = IESTART,IEEND
         WRITE(*,FMT='(I5,4F16.12)') IE,(CREN(IE,ISPIN),ISPIN=1,NSPIN)
      ENDDO
      WRITE(*,*) 'RENORM_LLY: renormalization factor per atom:'
      WRITE(*,FMT='(A5,2A16)') 'IAT','Charge','Moment'
      DO I1 = 1,NATYP
         WRITE(*,FMT='(I5,2F16.12)') 
     &        I1,(RENORM_AT(I1,ISPIN),ISPIN=1,NSPIN)
      ENDDO
      WRITE(*,*) 'RENORM_LLY: Renormalized charge per atom:'
      WRITE(*,FMT='(A5,2A16)') 'IAT','Spin up','Spin down'
      DO I1 = 1,NATYP
         WRITE(*,FMT='(I5,2F16.12)') 
     &        I1,(CHARGE_LLY(I1,ISPIN),ISPIN=1,NSPIN)
      ENDDO
      SUM0(:) = 0.D0
      SUM1(:) = 0.D0
      DO ISPIN = 1,NSPIN
         SIGNSP = 2*ISPIN - 3       ! -1,+1 for spin down,up (ispin=1,2)
         IF (NSPIN.EQ.1) SIGNSP = 1
         DO I1 = 1,NATYP
            SUM0(ISPIN) = SUM0(ISPIN) + 
     &                    SIGNSP * CONC(I1) * CHARGE(I1,ISPIN)
            SUM1(ISPIN) = SUM1(ISPIN) + 
     &                    SIGNSP * CONC(I1) * CHARGE_LLY(I1,ISPIN)
         ENDDO
      ENDDO
      WRITE(*,FMT='(A45,2F14.10)') 
     &              'RENORM_LLY: Locally summed charge and moment:', 
     &               (SUM0(ISPIN),ISPIN=1,NSPIN) 
      WRITE(*,FMT='(A45,2F14.10)') 
     &              'RENORM_LLY: Renormalized charge and moment:  ',
     &               (SUM1(ISPIN),ISPIN=1,NSPIN) 
      WRITE(*,FMT='(A50,2F14.10)') 
     &            'RENORM_LLY: Renormalization factor of total charge:',
     &             SUM1(1)/SUM0(1)



      END SUBROUTINE RENORM_LLY
