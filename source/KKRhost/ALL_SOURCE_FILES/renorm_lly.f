      SUBROUTINE RENORM_LLY(                                       ! LLY Lloyd
     &     CDOS_LLY,IELAST,NSPIN,NATYP,CDEN,LMAXP1,CONC,
     &     IESTART,IEEND,WEZ,IRCUT,IPAN,EZ,ZAT,
     &     RHO2NS,R2NEF,DENEF,DENEFAT,ESPV)

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
      INTEGER LMAXP1,NATYP,NSPIN         ! LMAX+1 for FP or LMAX for ASA, number of atoms,radial pts
      INTEGER IESTART,IEEND,IELAST             ! Starting and ending energy point for renormalization
      INTEGER IRCUT(0:IPAND,NATYPD),IPAN(NATYPD) ! Mesh info
      DOUBLE PRECISION CONC(NATYPD)                      ! Concentration (for cpa)
      DOUBLE COMPLEX CDEN(0:LMAXD1,IELAST,NPOTD)    ! Non-renormalized density per atom (density=-cden/pi)
      DOUBLE COMPLEX CDOS_LLY(IEMXD,NSPIND)        ! DOS according to Lloyd's formula
      DOUBLE COMPLEX WEZ(IEMXD),EZ(IEMXD)
      DOUBLE PRECISION ZAT(NATYPD)
! Input/Output:
      DOUBLE PRECISION RHO2NS(IRMD,LMPOTD,NATYPD,2)
      DOUBLE PRECISION R2NEF(IRMD,LMPOTD,NATYPD,2)
      DOUBLE PRECISION DENEF,DENEFAT(NATYPD)
      DOUBLE PRECISION ESPV(0:LMAXD1,NPOTD)
! Internal:
      INTEGER LL,IE,I1,ISPIN,IPOT,SPINDEGEN,IRC1,SIGNSP,IDIM
      DOUBLE PRECISION RENORM_AT(NATYPD,2)                   ! 1: charge renormalization per atom (energy-integrated)
                                                             ! 2: same for spin moment
      DOUBLE COMPLEX CDOS_LOC(IEMXD,(1+KREL)*NSPIND)     !  Density from local summation
      DOUBLE COMPLEX CDOS_LOCVC(IEMXD,(1+KREL)*NSPIND)   !      and from Lloyd's formula
      DOUBLE PRECISION CREN(IEMXD,2)                     ! Renormalization constant for charge and spin density
      DOUBLE PRECISION CHARGE(NATYPD,2),CHARGE_LLY(NATYPD,2) ! Atomic charge per spin (local summation and renormalized)
      DOUBLE COMPLEX CHADD(IEMXD,NATYPD,NSPIND),CDOS_ADD     ! Integration step for charge/atom/spin
      DOUBLE COMPLEX QLLY(2),QSTAR(2)
      DOUBLE PRECISION SUM0(2),SUM1(2)
      DOUBLE COMPLEX CZERO
      DOUBLE PRECISION PI
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
      CDOS_LOC = CZERO
      CDOS_LOCVC = CZERO
      CHADD=CZERO
      DO IE = IESTART,IEEND
         DO ISPIN = 1,NSPIN
            DO I1 = 1,NATYP
               IPOT = (I1-1) * NSPIN + ISPIN
               CDOS_ADD = CZERO
               DO LL = 0,LMAXP1
                  CDOS_ADD = CDOS_ADD + CONC(I1) * CDEN(LL,IE,IPOT)    ! Factor 1/pi included in WEZ 
               END DO
               IF (ZAT(I1).GT.1D-06) THEN
               CDOS_LOC(IE,ISPIN) = CDOS_LOC(IE,ISPIN) + CDOS_ADD
               ELSE
               CDOS_LOCVC(IE,ISPIN) = CDOS_LOCVC(IE,ISPIN) + CDOS_ADD
               ENDIF
               CHADD(IE,I1,ISPIN) = WEZ(IE) * CDOS_ADD                 ! Complex charge
               CHARGE(I1,ISPIN) = CHARGE(I1,ISPIN) +  
     &                            DIMAG(CHADD(IE,I1,ISPIN))/DBLE(NSPIN)
            ENDDO ! I1=1,NATYP
            CDOS_LOC(IE,ISPIN) = -CDOS_LOC(IE,ISPIN) / PI
            CDOS_LOCVC(IE,ISPIN) = -CDOS_LOCVC(IE,ISPIN) / PI
           ENDDO ! ISPIN = 1,NSPIN
          ENDDO ! IE = IESTART,IEEND
! Now the locally-summed charge/energy is in cdos_loc, charge/energy/atom in chadd
          IF (.NOT.OPT('NEWSOSOL')) THEN
           DO IE=IESTART,IEEND
            DO ISPIN=1,NSPIN
            ! Renormalization factor per energy:
            CREN(IE,ISPIN) = DIMAG((CDOS_LLY(IE,ISPIN)-
     +                              CDOS_LOCVC(IE,ISPIN))*WEZ(IE))/
     +                        DIMAG(CDOS_LOC(IE,ISPIN)*WEZ(IE))
            ! Apply to DOS of each atom:
             DO I1=1,NATYPD
              IF (ZAT(I1).GT.1D-06) THEN
               CHARGE_LLY(I1,ISPIN) = CHARGE_LLY(I1,ISPIN) +
     &              CREN(IE,ISPIN)*DIMAG(CHADD(IE,I1,ISPIN))/DBLE(NSPIN)
              ELSE
               CHARGE_LLY(I1,ISPIN) = CHARGE_LLY(I1,ISPIN) +
     &                             DIMAG(CHADD(IE,I1,ISPIN))/DBLE(NSPIN)
              ENDIF
             ENDDO
            ENDDO ! ISPIN = 1,NSPIN
           ENDDO ! IE = IESTART,IEEND
          ELSE 
           DO IE=IESTART,IEEND
            ! Renormalization factor per energy:
            CREN(IE,1) = DIMAG((CDOS_LLY(IE,1)-CDOS_LOCVC(IE,1)-
     +                          CDOS_LOCVC(IE,2))*WEZ(IE))/
     +                    DIMAG((CDOS_LOC(IE,1)+CDOS_LOC(IE,2))*WEZ(IE))
            ! Apply to DOS of each atom:
            DO ISPIN=1,NSPIN
             DO I1=1,NATYPD
              IF (ZAT(I1).GT.1D-06) THEN
               CHARGE_LLY(I1,ISPIN) = CHARGE_LLY(I1,ISPIN) +  
     &               CREN(IE,1)*DIMAG(CHADD(IE,I1,ISPIN))/DBLE(NSPIN)
              ELSE
               CHARGE_LLY(I1,ISPIN) = CHARGE_LLY(I1,ISPIN) +  
     &                          DIMAG(CHADD(IE,I1,ISPIN))/DBLE(NSPIN)
              ENDIF
             ENDDO
            ENDDO
           ENDDO ! IE = IESTART,IEEND
          ENDIF
 
c add term from sum from l>lmax to infinity
!            DO I1=1,NATYPD
!             DO ISPIN=1,NSPIN
c             CHARGE_LLY(I1,ISPIN)=CHARGE_LLY(I1,ISPIN)-DIMAG(CDOS2(I1))
!             ENDDO
!            ENDDO

      IF (NSPIN.EQ.1.OR.OPT('NEWSOSOL')) CREN(:,2) = CREN(:,1) 

! Now apply renormalization to energy-integrated density
! If spins are coupled, then only charge density 
      IF (NSPIN.EQ.1) THEN
       DO I1 = 1,NATYP
        if (CHARGE(I1,1)>0) then
           RENORM_AT(I1,1) = CHARGE_LLY(I1,1)/CHARGE(I1,1)  
        else
           RENORM_AT(I1,1) = 1.0d0  
        end if
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
         IF (ZAT(I1).GT.1D-06) THEN
          DENEFAT(I1) = DENEFAT(I1) - 2.0D0*CONC(I1)*CREN(IELAST,ISPIN)*
     &              DIMAG(CDEN(LL,IELAST,IPOT))/PI/DBLE(NSPIN) 
         ELSE
          DENEFAT(I1) = DENEFAT(I1) - 2.0D0*CONC(I1)*
     &              DIMAG(CDEN(LL,IELAST,IPOT))/PI/DBLE(NSPIN) 
         ENDIF
         ESPV(LL,IPOT)=0d0
         IF (ZAT(I1).GT.1D-06) THEN
          DO IE=1,IELAST
           ESPV(LL,IPOT) = ESPV(LL,IPOT)+CREN(IE,ISPIN)*
     &               DIMAG(EZ(IE)*CDEN(LL,IE,IPOT)*WEZ(IE)/DBLE(NSPIN))
          ENDDO
         ELSE
          DO IE=1,IELAST
           ESPV(LL,IPOT) = ESPV(LL,IPOT)+
     &               DIMAG(EZ(IE)*CDEN(LL,IE,IPOT)*WEZ(IE)/DBLE(NSPIN))
          ENDDO
         ENDIF
        ENDDO ! LL
       ENDDO ! ISPIN
       DENEF = DENEF+DENEFAT(I1)
      ENDDO ! I1
! Write out renormalization factors
      WRITE(1337,*) 'Information on renormalization by Lloyds formula'
      WRITE(1337,*) 
     &        'RENORM_LLY: Complex renormalization factor per energy:'
      WRITE(1337,FMT='(A5,2A32)') 'IE',
     &               'Spin 1 (down)           ','Spin 2 (up)           '
      DO IE = IESTART,IEEND
         WRITE(1337,FMT='(I5,4F16.12)') IE,
     +                (CREN(IE,ISPIN),ISPIN=1,NSPIN)
      ENDDO
      WRITE(1337,*) 'RENORM_LLY: renormalization factor per atom:'
      WRITE(1337,FMT='(A5,2A16)') 'IAT','Spin down','Spin up'
      DO I1 = 1,NATYP
         WRITE(1337,FMT='(I5,2E17.9)') 
     &        I1,(RENORM_AT(I1,ISPIN),ISPIN=1,NSPIN)
      ENDDO
      WRITE(1337,*) 'RENORM_LLY: Renormalized charge per atom:'
      WRITE(1337,FMT='(A5,2A16)') 'IAT','Spin down','Spin up'
      DO I1 = 1,NATYP
         WRITE(1337,FMT='(I5,2F16.12)') 
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
      WRITE(1337,FMT='(A45,2E17.9)') 
     &              'RENORM_LLY: Locally summed charge and moment:', 
     &               (SUM0(ISPIN),ISPIN=1,NSPIN) 
      WRITE(1337,FMT='(A45,2E17.9)') 
     &              'RENORM_LLY: Renormalized charge and moment:  ',
     &               (SUM1(ISPIN),ISPIN=1,NSPIN) 
      WRITE(1337,FMT='(A50,2E17.9)') 
     &            'RENORM_LLY: Renormalization factor of total charge:',
     &             SUM1(1)/SUM0(1)



      END SUBROUTINE RENORM_LLY
