      SUBROUTINE RENORM_LLY(                                       ! LLY Lloyd
     >     CDOS_LLY,IELAST,NSPIN,NATYP,CDEN,LMAXP1,CONC,
     >     IESTART,IEEND,WEZ,IRCUT,IPAN,EZ,
     >     RHO2NS,CDOS1,CDOS2,R2NEF,DENEF,DENEFAT,ESPV)

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
      COMPLEX*16 CDEN(0:LMAXD1,IEMXD,NPOTD)    ! Non-renormalized density per atom (density=-cden/pi)
      COMPLEX*16 CDOS_LLY(IEMXD,NSPIND),
     +           CDOS1(IEMXD),CDOS2(NATYPD)    ! DOS according to Lloyd's formula
      COMPLEX*16 WEZ(IEMXD),EZ(IEMXD)
! Input/Output:
      DOUBLE PRECISION RHO2NS(IRMD,LMPOTD,NATYPD,2)
      DOUBLE PRECISION R2NEF(IRMD,LMPOTD,NATYPD,2)
      DOUBLE PRECISION DENEF,DENEFAT(NATYPD)
      DOUBLE PRECISION ESPV(0:LMAXD1,NPOTD)
! Internal:
      INTEGER LL,IE,I1,ISPIN,IPOT,SPINDEGEN,IRC1,SIGNSP,IDIM
      REAL*8 RENORM_AT(NATYPD,2)                   ! 1: charge renormalization per atom (energy-integrated)
                                                   ! 2: same for spin moment
      COMPLEX*16 CDOS_LOC(IEMXD,NSPIND)            ! Density from local summation and from Lloyd's formula
      REAL*8 CREN(IEMXD,2)                     ! Renormalization constant for charge and spin density
      REAL*8 CHARGE(NATYPD,2),CHARGE_LLY(NATYPD,2) ! Atomic charge per spin (local summation and renormalized)
      COMPLEX*16 CHADD(IEMXD,NATYPD,NSPIND),CDOS_ADD     ! Integration step for charge/atom/spin
      COMPLEX*16 QLLY(2),QSTAR(2)
      REAL*8 SUM0(2),SUM1(2)
      COMPLEX*16 CZERO
      REAL*8 PI
      LOGICAL OPT,TEST
      EXTERNAL OPT,TEST

      CZERO = (0.D0,0.D0)
      PI = 4.D0 * DATAN(1.D0)

      SPINDEGEN = 3 - NSPIN ! Spin degeneracy, 2 if nspin=1, 1 if nspin=2

      CREN(:,:) = 0d0
      RENORM_AT(:,:) = 1.D0
      CHARGE_LLY(:,:) = 0.D0
      CHARGE(:,:) = 0.D0
      QLLY(:)=CZERO
      QSTAR(:)=CZERO

! First find renormalization factor per energy and atomic charges
c      DO IE=IESTART,IEEND
c       WRITE(65,'(2E17.9)') WEZ(IE)
c       DO ISPIN=1,NSPIN
c        DO I1=1,NATYP
c         DO LL=0,LMAXP1
c          IPOT=(I1-1)*NSPIN+ISPIN
c          WRITE(66,'(3I5,2E17.9)') IE,IPOT,LL,CDEN(LL,IE,IPOT)
c         ENDDO
c        ENDDO
c       ENDDO
c      ENDDO
      CDOS_LOC = CZERO
      CHADD=CZERO
      DO IE = IESTART,IEEND
         DO ISPIN = 1,NSPIN
            DO I1 = 1,NATYP
               IPOT = (I1-1) * NSPIN + ISPIN
               CDOS_ADD = CZERO
               DO LL = 0,LMAXP1
                  CDOS_ADD = CDOS_ADD + CONC(I1) * CDEN(LL,IE,IPOT)    ! Factor 1/pi included in WEZ 
               END DO
               CDOS_LOC(IE,ISPIN) = CDOS_LOC(IE,ISPIN) + CDOS_ADD
               CHADD(IE,I1,ISPIN) = WEZ(IE) * CDOS_ADD                 ! Complex charge
c        WRITE(67,'(3I5,4E17.9)') IE,I1,ISPIN,CDOS_ADD,CHADD(IE,I1,ISPIN)
               CHARGE(I1,ISPIN) = CHARGE(I1,ISPIN) +  
     &                            DIMAG(CHADD(IE,I1,ISPIN))/DBLE(NSPIN)
            ENDDO ! I1=1,NATYP
            CDOS_LOC(IE,ISPIN) = -CDOS_LOC(IE,ISPIN) / PI
           ENDDO ! ISPIN = 1,NSPIN
          ENDDO ! IE = IESTART,IEEND
! Now the locally-summed charge/energy is in cdos_loc, charge/energy/atom in chadd
          IF (.NOT.OPT('NEWSOSOL')) THEN
           DO IE=IESTART,IEEND
            DO ISPIN=1,NSPIN
            ! Renormalization factor per energy:
            CREN(IE,ISPIN) = DIMAG(CDOS_LLY(IE,ISPIN)*WEZ(IE))/
     +                        DIMAG(CDOS_LOC(IE,ISPIN)*WEZ(IE))
c            CREN(IE,ISPIN) = DIMAG((CDOS_LLY(IE,ISPIN)-CDOS1(IE))*
c     +                        WEZ(IE))/
c     +                        DIMAG(CDOS_LOC(IE,ISPIN)*WEZ(IE))

            ! Apply to DOS of each atom:
             DO I1=1,NATYPD
              CHARGE_LLY(I1,ISPIN) = CHARGE_LLY(I1,ISPIN) +
     &              CREN(IE,ISPIN)*DIMAG(CHADD(IE,I1,ISPIN))/DBLE(NSPIN)
             ENDDO
            ENDDO ! ISPIN = 1,NSPIN
c             WRITE(68,'(1I5,4E17.9)') IE,(CREN(IE,ISPIN),ISPIN=1,NSPIN)
           ENDDO ! IE = IESTART,IEEND
          ELSE 
           DO IE=IESTART,IEEND
            ! Renormalization factor per energy:
            CREN(IE,1) = DIMAG(CDOS_LLY(IE,1)*WEZ(IE))/
     +                    DIMAG((CDOS_LOC(IE,1)+CDOS_LOC(IE,2))*WEZ(IE))
c            WEZ(IE) = CREN(IE,1)*WEZ(IE)
c            CREN(IE,1) = DIMAG((CDOS_LLY(IE,1)-2d0*CDOS1(IE))*
c     +                        WEZ(IE))/
c     +                    DIMAG((CDOS_LOC(IE,1)+CDOS_LOC(IE,2))*WEZ(IE))
c             WRITE(68,'(1I5,4E17.9)') IE,CREN(IE,1)
            ! Apply to DOS of each atom:
            DO ISPIN=1,NSPIN
             DO I1=1,NATYPD
              CHARGE_LLY(I1,ISPIN) = CHARGE_LLY(I1,ISPIN) +  
     &               CREN(IE,1)*DIMAG(CHADD(IE,I1,ISPIN))/DBLE(NSPIN)
             ENDDO
            ENDDO
           ENDDO ! IE = IESTART,IEEND
          ENDIF
 
c add term from sum from l>lmax to infinity
           DO I1=1,NATYPD
            DO ISPIN=1,NSPIN
c             CHARGE_LLY(I1,ISPIN)=CHARGE_LLY(I1,ISPIN)-DIMAG(CDOS2(I1))
            ENDDO
           ENDDO

      IF (NSPIN.EQ.1.OR.OPT('NEWSOSOL')) CREN(:,2) = CREN(:,1) 

! Now apply renormalization to energy-integrated density
! If spins are coupled, then only charge density 
      IF (NSPIN.EQ.1) THEN

       DO I1 = 1,NATYP
        RENORM_AT(I1,1) = CHARGE_LLY(I1,1)/CHARGE(I1,1)  
        RENORM_AT(I1,2) = RENORM_AT(I1,1)  
        IRC1 = IRCUT(IPAN(I1),I1) ! Index of outmost radial point
        RHO2NS(1:IRC1,1:LMPOTD,I1,1) = 
     &           RHO2NS(1:IRC1,1:LMPOTD,I1,1)*RENORM_AT(I1,1)
        R2NEF(1:IRC1,1:LMPOTD,I1,1) = 
     &           R2NEF(1:IRC1,1:LMPOTD,I1,1)*RENORM_AT(I1,1)
       ENDDO
       
      ELSE

! First decouple charge and spin density to the density of each channels
        IDIM = IRMD*LMPOTD*NATYPD
        CALL DAXPY(IDIM,1.0D0,RHO2NS(1,1,1,1),1,RHO2NS(1,1,1,2),1)
        CALL DSCAL(IDIM,0.5D0,RHO2NS(1,1,1,2),1)
        CALL DAXPY(IDIM,-1.0D0,RHO2NS(1,1,1,2),1,RHO2NS(1,1,1,1),1)

        CALL DAXPY(IDIM,1.0D0,R2NEF(1,1,1,1),1,R2NEF(1,1,1,2),1)
        CALL DSCAL(IDIM,0.5D0,R2NEF(1,1,1,2),1)
        CALL DAXPY(IDIM,-1.0D0,R2NEF(1,1,1,2),1,R2NEF(1,1,1,1),1)
        DO I1 = 1,NATYP
         IRC1 = IRCUT(IPAN(I1),I1) ! Index of outmost radial point
         DO ISPIN=1,NSPIN
          RENORM_AT(I1,ISPIN) = CHARGE_LLY(I1,ISPIN)/CHARGE(I1,ISPIN)
          RHO2NS(1:IRC1,1:LMPOTD,I1,ISPIN) = 
     &           RHO2NS(1:IRC1,1:LMPOTD,I1,ISPIN) * RENORM_AT(I1,ISPIN)
          R2NEF(1:IRC1,1:LMPOTD,I1,ISPIN) = 
     &           R2NEF(1:IRC1,1:LMPOTD,I1,ISPIN)*RENORM_AT(I1,ISPIN)
         ENDDO
        ENDDO
! Second merge density of each channels to charge and spin density
        CALL DSCAL(IDIM,2.0D0,RHO2NS(1,1,1,1),1)
        CALL DAXPY(IDIM,-0.5D0,RHO2NS(1,1,1,1),1,RHO2NS(1,1,1,2),1)
        CALL DAXPY(IDIM,1.0D0,RHO2NS(1,1,1,2),1,RHO2NS(1,1,1,1),1)

        CALL DSCAL(IDIM,2.0D0,R2NEF(1,1,1,1),1)
        CALL DAXPY(IDIM,-0.5D0,R2NEF(1,1,1,1),1,R2NEF(1,1,1,2),1)
        CALL DAXPY(IDIM,1.0D0,R2NEF(1,1,1,2),1,R2NEF(1,1,1,1),1)
      ENDIF

! calculate density at Fermi level
      DENEF=0d0
      DO I1=1,NATYP
       DENEFAT(I1)=0d0
       DO ISPIN=1,NSPIN
        IPOT=(I1-1)*NSPIN+ISPIN
        DO LL=0,LMAXP1
         DENEF = DENEF - 2.0d0*CONC(I1)*CREN(IELAST,ISPIN)*
     &             DIMAG(CDEN(LL,IELAST,IPOT))/PI/DBLE(NSPIN)
         DENEFAT(I1) = DENEFAT(I1) - 2.0D0 * CREN(IELAST,ISPIN)*
     &              DIMAG(CDEN(LL,IELAST,IPOT))/PI/DBLE(NSPIN) 
         ESPV(LL,IPOT)=0d0
         DO IE=1,IELAST
          ESPV(LL,IPOT) = ESPV(LL,IPOT)+CREN(IE,ISPIN)*
     &               DIMAG(EZ(IE)*CDEN(LL,IE,IPOT)*WEZ(IE)/DBLE(NSPIN))
         ENDDO
        ENDDO
       ENDDO
      ENDDO
! Write out renormalization factors
      WRITE(*,*) 'Information on renormalization by Lloyds formula'
      WRITE(*,*)'RENORM_LLY: Complex renormalization factor per energy:'
      WRITE(*,FMT='(A5,2A32)') 'IE',
     &               'Spin 1 (down)           ','Spin 2 (up)           '
      DO IE = IESTART,IEEND
         WRITE(*,FMT='(I5,4F16.12)') IE,
     +                (CREN(IE,ISPIN),ISPIN=1,NSPIN)
      ENDDO
      WRITE(*,*) 'RENORM_LLY: renormalization factor per atom:'
      WRITE(*,FMT='(A5,2A16)') 'IAT','Spin down','Spin up'
      DO I1 = 1,NATYP
c         WRITE(*,FMT='(I5,2F16.12)') 
         WRITE(*,FMT='(I5,2E17.9)') 
     &        I1,(RENORM_AT(I1,ISPIN),ISPIN=1,NSPIN)
      ENDDO
      WRITE(*,*) 'RENORM_LLY: Renormalized charge per atom:'
      WRITE(*,FMT='(A5,2A16)') 'IAT','Spin down','Spin up'
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
      WRITE(*,FMT='(A45,2E17.9)') 
     &              'RENORM_LLY: Locally summed charge and moment:', 
     &               (SUM0(ISPIN),ISPIN=1,NSPIN) 
      WRITE(*,FMT='(A45,2E17.9)') 
     &              'RENORM_LLY: Renormalized charge and moment:  ',
     &               (SUM1(ISPIN),ISPIN=1,NSPIN) 
      WRITE(*,FMT='(A50,2E17.9)') 
     &            'RENORM_LLY: Renormalization factor of total charge:',
     &             SUM1(1)/SUM0(1)



      END SUBROUTINE RENORM_LLY
