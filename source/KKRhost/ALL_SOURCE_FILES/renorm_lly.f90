SUBROUTINE renorm_lly( &                                      ! LLY Lloyd  &
        cdos_lly,ielast,nspin,natyp,cden,lmaxp1,conc,  &
        iestart,ieend,wez,ircut,ipan,ez,zat,  &
        rho2ns,r2nef,denef,denefat,espv)
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
      PARAMETER (NPOTD= (2*(KREL+KORBIT) + (1-(KREL+KORBIT))*NSPIND)*NATYPD)
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


czero = (0.d0,0.d0)
pi = 4.d0 * DATAN(1.d0)

spindegen = 3 - nspin ! Spin degeneracy, 2 if nspin=1, 1 if nspin=2

cren(:,:) = 0D0
renorm_at(:,:) = 1.d0
charge_lly(:,:) = 0.d0
charge(:,:) = 0.d0
qlly(:)=czero
qstar(:)=czero

! First find renormalization factor per energy and atomic charges
cdos_loc = czero
cdos_locvc = czero
chadd=czero
DO ie = iestart,ieend
  DO ispin = 1,nspin
    DO i1 = 1,natyp
      ipot = (i1-1) * nspin + ispin
      cdos_add = czero
      DO ll = 0,lmaxp1
        cdos_add = cdos_add + conc(i1) * cden(ll,ie,ipot)    ! Factor 1/pi included in Wez
      END DO
      IF (zat(i1) > 1D-06) THEN
        cdos_loc(ie,ispin) = cdos_loc(ie,ispin) + cdos_add
      ELSE
        cdos_locvc(ie,ispin) = cdos_locvc(ie,ispin) + cdos_add
      END IF
      chadd(ie,i1,ispin) = wez(ie) * cdos_add                 ! Complex charge
      charge(i1,ispin) = charge(i1,ispin) +  &
          DIMAG(chadd(ie,i1,ispin))/DBLE(nspin)
    END DO ! I1=1,NATYP
    cdos_loc(ie,ispin) = -cdos_loc(ie,ispin) / pi
    cdos_locvc(ie,ispin) = -cdos_locvc(ie,ispin) / pi
  END DO ! ISPIN = 1,NSPIN
END DO ! IE = IESTART,IEEND
! Now the locally-summed charge/energy is in cdos_loc, charge/energy/atom in chadd
IF (.NOT.opt('NEWSOSOL')) THEN
  DO ie=iestart,ieend
    DO ispin=1,nspin
! Renormalization factor per energy:
      cren(ie,ispin) = DIMAG((cdos_lly(ie,ispin)-  &
          cdos_locvc(ie,ispin))*wez(ie))/ DIMAG(cdos_loc(ie,ispin)*wez(ie))
! Apply to DOS of each atom:
      DO i1=1,natypd
        IF (zat(i1) > 1D-06) THEN
          charge_lly(i1,ispin) = charge_lly(i1,ispin) +  &
              cren(ie,ispin)*DIMAG(chadd(ie,i1,ispin))/DBLE(nspin)
        ELSE
          charge_lly(i1,ispin) = charge_lly(i1,ispin) +  &
              DIMAG(chadd(ie,i1,ispin))/DBLE(nspin)
        END IF
      END DO
    END DO ! ISPIN = 1,NSPIN
  END DO ! IE = IESTART,IEEND
ELSE
  DO ie=iestart,ieend
! Renormalization factor per energy:
    cren(ie,1) = DIMAG((cdos_lly(ie,1)-cdos_locvc(ie,1)-  &
        cdos_locvc(ie,2))*wez(ie))/ DIMAG((cdos_loc(ie,1)+cdos_loc(ie,2))*wez(ie))
! Apply to DOS of each atom:
    DO ispin=1,nspin
      DO i1=1,natypd
        IF (zat(i1) > 1D-06) THEN
          charge_lly(i1,ispin) = charge_lly(i1,ispin) +  &
              cren(ie,1)*DIMAG(chadd(ie,i1,ispin))/DBLE(nspin)
        ELSE
          charge_lly(i1,ispin) = charge_lly(i1,ispin) +  &
              DIMAG(chadd(ie,i1,ispin))/DBLE(nspin)
        END IF
      END DO
    END DO
  END DO ! IE = IESTART,IEEND
END IF

! add term from sum from l>lmax to infinity
!            DO I1=1,NATYPD
!             DO ISPIN=1,NSPIN
!             CHARGE_LLY(I1,ISPIN)=CHARGE_LLY(I1,ISPIN)-DIMAG(CDOS2(I1))
!             ENDDO
!            ENDDO

IF (nspin == 1.OR.opt('NEWSOSOL')) cren(:,2) = cren(:,1)

! Now apply renormalization to energy-integrated density
! If spins are coupled, then only charge density
IF (nspin == 1) THEN
  DO i1 = 1,natyp
    IF (charge(i1,1)>0) THEN
      renorm_at(i1,1) = charge_lly(i1,1)/charge(i1,1)
    ELSE
      renorm_at(i1,1) = 1.0D0
    END IF
    renorm_at(i1,2) = renorm_at(i1,1)
    irc1 = ircut(ipan(i1),i1) ! Index of outmost radial point
    rho2ns(1:irc1,1:lmpotd,i1,1) =  &
        rho2ns(1:irc1,1:lmpotd,i1,1)*renorm_at(i1,1)
    r2nef(1:irc1,1:lmpotd,i1,1) = r2nef(1:irc1,1:lmpotd,i1,1)*renorm_at(i1,1)
  END DO
  
ELSE
  
! First decouple charge and spin density to the density of each channels
  idim = irmd*lmpotd*natypd
  CALL daxpy(idim,1.0D0,rho2ns(1,1,1,1),1,rho2ns(1,1,1,2),1)
  CALL dscal(idim,0.5D0,rho2ns(1,1,1,2),1)
  CALL daxpy(idim,-1.0D0,rho2ns(1,1,1,2),1,rho2ns(1,1,1,1),1)
  
  CALL daxpy(idim,1.0D0,r2nef(1,1,1,1),1,r2nef(1,1,1,2),1)
  CALL dscal(idim,0.5D0,r2nef(1,1,1,2),1)
  CALL daxpy(idim,-1.0D0,r2nef(1,1,1,2),1,r2nef(1,1,1,1),1)
  DO i1 = 1,natyp
    irc1 = ircut(ipan(i1),i1) ! Index of outmost radial point
    DO ispin=1,nspin
      renorm_at(i1,ispin) = charge_lly(i1,ispin)/charge(i1,ispin)
      rho2ns(1:irc1,1:lmpotd,i1,ispin) =  &
          rho2ns(1:irc1,1:lmpotd,i1,ispin) * renorm_at(i1,ispin)
      r2nef(1:irc1,1:lmpotd,i1,ispin) =  &
          r2nef(1:irc1,1:lmpotd,i1,ispin)*renorm_at(i1,ispin)
    END DO
  END DO
! Second merge density of each channels to charge and spin density
  CALL dscal(idim,2.0D0,rho2ns(1,1,1,1),1)
  CALL daxpy(idim,-0.5D0,rho2ns(1,1,1,1),1,rho2ns(1,1,1,2),1)
  CALL daxpy(idim,1.0D0,rho2ns(1,1,1,2),1,rho2ns(1,1,1,1),1)
  
  CALL dscal(idim,2.0D0,r2nef(1,1,1,1),1)
  CALL daxpy(idim,-0.5D0,r2nef(1,1,1,1),1,r2nef(1,1,1,2),1)
  CALL daxpy(idim,1.0D0,r2nef(1,1,1,2),1,r2nef(1,1,1,1),1)
END IF

! calculate density at Fermi level
denef=0D0
DO i1=1,natyp
  denefat(i1)=0D0
  DO ispin=1,nspin
    ipot=(i1-1)*nspin+ispin
    DO ll=0,lmaxp1
      IF (zat(i1) > 1D-06) THEN
        denefat(i1) = denefat(i1) - 2.0D0*conc(i1)*cren(ielast,ispin)*  &
            DIMAG(cden(ll,ielast,ipot))/pi/DBLE(nspin)
      ELSE
        denefat(i1) = denefat(i1) - 2.0D0*conc(i1)*  &
            DIMAG(cden(ll,ielast,ipot))/pi/DBLE(nspin)
      END IF
      espv(ll,ipot)=0D0
      IF (zat(i1) > 1D-06) THEN
        DO ie=1,ielast
          espv(ll,ipot) = espv(ll,ipot)+cren(ie,ispin)*  &
              DIMAG(ez(ie)*cden(ll,ie,ipot)*wez(ie)/DBLE(nspin))
        END DO
      ELSE
        DO ie=1,ielast
          espv(ll,ipot) = espv(ll,ipot)+  &
              DIMAG(ez(ie)*cden(ll,ie,ipot)*wez(ie)/DBLE(nspin))
        END DO
      END IF
    END DO ! LL
  END DO ! ISPIN
  denef = denef+denefat(i1)
END DO ! I1
! Write out renormalization factors
WRITE(1337,*) 'Information on renormalization by Lloyds formula'
WRITE(1337,*) 'RENORM_LLY: Complex renormalization factor per energy:'
WRITE(1337,FMT='(A5,2A32)') 'IE',  &
    'Spin 1 (down)           ','Spin 2 (up)           '
DO ie = iestart,ieend
  WRITE(1337,FMT='(I5,4F16.12)') ie, (cren(ie,ispin),ispin=1,nspin)
END DO
WRITE(1337,*) 'RENORM_LLY: renormalization factor per atom:'
WRITE(1337,FMT='(A5,2A16)') 'IAT','Spin down','Spin up'
DO i1 = 1,natyp
  WRITE(1337,FMT='(I5,2E17.9)') i1,(renorm_at(i1,ispin),ispin=1,nspin)
END DO
WRITE(1337,*) 'RENORM_LLY: Renormalized charge per atom:'
WRITE(1337,FMT='(A5,2A16)') 'IAT','Spin down','Spin up'
DO i1 = 1,natyp
  WRITE(1337,FMT='(I5,2F16.12)') i1,(charge_lly(i1,ispin),ispin=1,nspin)
END DO
sum0(:) = 0.d0
sum1(:) = 0.d0
DO ispin = 1,nspin
  signsp = 2*ispin - 3       ! -1,+1 for spin down,up (ispin=1,2)
  IF (nspin == 1) signsp = 1
  DO i1 = 1,natyp
    sum0(ispin) = sum0(ispin) + signsp * conc(i1) * charge(i1,ispin)
    sum1(ispin) = sum1(ispin) + signsp * conc(i1) * charge_lly(i1,ispin)
    
  END DO
END DO
WRITE(1337,FMT='(A45,2E17.9)')  &
    'RENORM_LLY: Locally summed charge and moment:', (sum0(ispin),ispin=1,nspin)
WRITE(1337,FMT='(A45,2E17.9)')  &
    'RENORM_LLY: Renormalized charge and moment:  ', (sum1(ispin),ispin=1,nspin)
WRITE(1337,FMT='(A50,2E17.9)')  &
    'RENORM_LLY: Renormalization factor of total charge:', sum1(1)/sum0(1)



END SUBROUTINE renorm_lly
