SUBROUTINE decitset(alat,bravsys,ez,ielast,  &
        nlbasis,nrbasis,fileleft,fileright,  &
        ins,kvrel,krel,nspin,kmrot,  &
        vref,rmtref,nref,refpot,  &
        lefttinv,righttinv,vacflag,  &
        nembd1,iemxd,irmd,ipand,  &
        lmaxd,lmgf0d,lmmaxd,lm2d,nspind)
! **********************************************************************
! *                                                                    *
! * This subroutine is thought as an alternative to the < decimaread > *
! * which requires an a priori calculated set of single-site matrices  *
! * over a fixed energy mesh. It is using the potential as written out *
! * in < outpothost > routine and determines the matrix                *
! *                                                                    *
! *            /         \-1   /             \-1                       *
! *            | Delta t |   = |  t    - t   |                         *
! *            \         /     \  sys    ref /                         *
! *                                                                    *
! * for the left and the right host, using the energy mesh as read in  *
! * from the input-file                                                *
! *                                                                    *
! *                                        v.popescu - munich, Dec 04  *
! *                                                                    *
! * Notes: - no charge moments are calculated -- thus this option CAN  *
! *          NOT be used in SCF calculations                           *
! *        - non-spherical case not implemented, neither LDA+U (al-    *
! *          though the interface to regsol in decitmat is supplied)   *
! *        - CPA case not implemented - requires BZ integration        *
! *                                                                    *
! **********************************************************************
      IMPLICIT NONE
!..
!.. Scalars arguments ..
      INTEGER IEMXD,NEMBD1,LMMAXD,IPAND,NSPIND,IRMD,LMAXD,LM2D,LMGF0D
      INTEGER IELAST,KMROT,NLBASIS,NRBASIS,NREF,INS,KVREL,KREL,NSPIN
      DOUBLE PRECISION ALAT
      CHARACTER*40 FILELEFT,FILERIGHT
!..
!.. Array arguments ..
      INTEGER REFPOT(NEMBD1)
      DOUBLE PRECISION VREF(*),RMTREF(*),BRAVSYS(3,3)
      DOUBLE COMPLEX EZ(IEMXD)
      DOUBLE COMPLEX LEFTTINV(LMMAXD,LMMAXD,NEMBD1,NSPIND,IEMXD), &
                     RIGHTTINV(LMMAXD,LMMAXD,NEMBD1,NSPIND,IEMXD)
      LOGICAL VACFLAG(2)
!..
!.. Local scalars ..
      INTEGER IHOST,I,LL,MM,LNGSTRING,NQHOST,ILHOST
      INTEGER NHOST
      INTEGER NQ,NT,IQOFF,ITOFF,IE,IH,IQH,IOQ,INFO
      INTEGER IPOT,I1,ISPIN,NSRA,LM1,LM2,IRC1,IREF
      INTEGER NTLEFT,NTRIGHT,NTHOST
!.. LDA+U
      INTEGER IDOLDAU,LOPT
      DOUBLE PRECISION WLDAUAV
!..
      DOUBLE PRECISION EFERMI,RIRC
      DOUBLE COMPLEX ERYD,CARG,CFCTOR
      LOGICAL TEST
      CHARACTER*40 FILEHOST
      CHARACTER*10 SOLVER
!..
!.. Local arrays
      INTEGER KRELH(2),NSPINH(2),INSH(2),IPVT(LMMAXD)
      INTEGER NOQ(NEMBD1),KAOEZ(NEMBD1,NEMBD1),INHOST(2)
      DOUBLE PRECISION BRAVAIS(3,3,2),RBASIS(3,NEMBD1),QMTET(NEMBD1), &
                       QMPHI(NEMBD1)
      CHARACTER*5 CHHOST(2)
      CHARACTER*9 TXTS(2)
!..
!.. Allocatable local arrays 
      INTEGER NTMAX
      DOUBLE PRECISION ZAT(:),RWS(:),RMT(:),CONC(:)
      DOUBLE PRECISION RR(:,:),DRDI(:,:),VISP(:,:),DROR(:,:)
      DOUBLE PRECISION SOCSCL(:,:),CSCL(:,:)
      INTEGER IRWS(:),IPAN(:),IQAT(:,:),IRCUT(:,:),LOFLM(:)
      DOUBLE COMPLEX TREFLL(:,:,:),TMATLL(:,:),DHMAT(:,:,:)
      DOUBLE COMPLEX DTREFLL(:,:,:) ! LLY Lloyd
      DOUBLE COMPLEX ALPHAREF(:,:),DALPHAREF(:,:) ! LLY Lloyd Alpha matrix and deriv.
      DOUBLE PRECISION VTREL(:,:),BTREL(:,:),R2DRDIREL(:,:)
      INTEGER ZREL(:)
      ALLOCATABLE ZAT,RWS,RMT,CONC,RR,DRDI,VISP,DROR,SOCSCL,CSCL
      ALLOCATABLE VTREL,BTREL,R2DRDIREL
      ALLOCATABLE IRWS,IPAN,IQAT,IRCUT,LOFLM,ZREL
      ALLOCATABLE TREFLL,TMATLL,DHMAT
      ALLOCATABLE DTREFLL,ALPHAREF,DALPHAREF ! LLY
!.. 
!.. External subroutines
      EXTERNAL CALCTREF13,CHANGEREP,CINIT,CMATSTR,DECIPOTBAS, &
               DECIPOTHEAD,DECITMAT,ZAXPY,ZCOPY,ZGETRF,ZGETRI
!..
!.. External Functions ..
      EXTERNAL LNGSTRING,TEST
!..
!.. Data statements
      DATA CHHOST/'LEFT ','RIGHT'/
      DATA TXTS /'spin   UP','spin DOWN'/
! ......................................................................

cfctor = alat/(8.d0*ATAN(1.0D0))           ! = ALAT/(2*PI)

idoldau = 0
lopt = -1
wldauav = 0D0
allocate(loflm(lm2d),stat=i1)
IF ( i1 /= 0 ) STOP '    Allocate LOFLM'
WRITE (6,'(5X,A,/,8X,65(1H-))') 'Reading in host potentials'
vacflag(1) = .false.
vacflag(2) = .false.
nsra = 1
IF ( kvrel >= 1 ) nsra = 2
i = 1
DO ll = 0,2*lmaxd
  DO mm = -ll,ll
    loflm(i) = ll
    i = i + 1
  END DO
END DO
ntleft = 0
ntright = 0
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
nhost = 0
DO ihost = 1,2
  filehost = fileleft
  nqhost = nlbasis
  IF ( ihost == 2 ) THEN
    filehost = fileright
    nqhost = nrbasis
  END IF
  ilhost = lngstring(filehost,40)
  CALL decipothead(ihost,filehost,ilhost,nqhost,  &
      vacflag,alat,bravsys,nq,nt,bravais(1,1,ihost),  &
      efermi,insh(ihost),krelh(ihost),nspinh(ihost), ins,krel,nspin,kmrot)
  
  IF ( .NOT.vacflag(ihost) ) THEN
    nhost = nhost + 1
    inhost(nhost) = ihost
    IF ( ihost == 1 ) THEN
      ntleft = nt
    ELSE
      ntright = nt
    END IF
  END IF
  
  
END DO

IF ( ntleft+ntright <= 0 ) THEN
  WRITE (6,'(8X,"Vacuum will be considered on both sides",/, 8X,65(1H-))')
  RETURN
END IF

ntmax = ntleft+ntright
allocate(zat(ntmax),rws(ntmax),rmt(ntmax),conc(ntmax),stat=i1)
IF ( i1 /= 0 ) STOP '    Allocate ZAT/RWS/RMT/CONC'
allocate(rr(irmd,ntmax),drdi(irmd,ntmax),stat=i1)
IF ( i1 /= 0 ) STOP '    Allocate RR/DRDI'
allocate(visp(irmd,ntmax*nspind),stat=i1)
IF ( i1 /= 0 ) STOP '    Allocate VISP'
allocate(irws(ntmax),ipan(ntmax),iqat(nembd1,ntmax),stat=i1)
IF ( i1 /= 0 ) STOP '    Allocate IRWS/IPAN/IQAT'
allocate(ircut(0:ipand,ntmax),stat=i1)
IF ( i1 /= 0 ) STOP '    Allocate IRCUT'
allocate(socscl(krel*lmaxd+1,krel*ntmax+(1-krel)),stat=i1)
IF ( i1 /= 0 ) STOP '    Allocate SOCSCL'
allocate(cscl(krel*lmaxd+1,krel*ntmax+(1-krel)),stat=i1)
IF ( i1 /= 0 ) STOP '    Allocate CSCL'
allocate(vtrel(irmd*krel+(1-krel),ntmax),stat=i1)
IF ( i1 /= 0 ) STOP '    Allocate VTREL'
allocate(btrel(irmd*krel+(1-krel),ntmax),stat=i1)
IF ( i1 /= 0 ) STOP '    Allocate BTREL'

! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
DO ihost = 1,2
  WRITE (6,'(8X,A5," side host: ",$)') chhost(ihost)
  iqoff = 0
  itoff = 0
  nqhost = nlbasis
  nthost = ntleft
  filehost = fileleft
  IF ( ihost == 2 ) THEN
    nqhost = nrbasis
    nthost = ntright
    iqoff = nlbasis
    itoff = ntleft
    filehost = fileright
  END IF
  ilhost = lngstring(filehost,40)
  
  IF ( filehost(1:7) == 'vacuum' ) THEN
    WRITE (6,'(A,/,8X,65(1H-))') 'VACUUM will be used'
  ELSE
    WRITE (6,'(A,/)') filehost(1:ilhost)
    WRITE (6,99005) krelh(ihost),nspinh(ihost),insh(ihost),  &
        kmrot,nqhost,alat,efermi
    WRITE (6,99006) ((bravais(ll,mm,ihost),mm=1,3),ll=1,3)
    CALL decipotbas(ihost,iqoff,itoff,nqhost,nthost,  &
        rbasis,qmtet,qmphi,noq,kaoez, zat,iqat,conc,irws,ipan,ircut,  &
        rr,drdi,visp,nspinh(ihost),krelh(ihost),  &
        solver,socscl,cscl,vtrel,btrel, irmd,ipand,nembd1,ntmax,nspind,lmaxd)
    WRITE (6,'(8X,65(1H-))')
  END IF
END DO
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

allocate(dror(irmd,ntmax),stat=i1)
IF ( i1 /= 0 ) STOP '    Allocate DROR'
allocate(r2drdirel(irmd*krel+(1-krel),ntmax),stat=i1)
IF ( i1 /= 0 ) STOP '    Allocate R2DRDIREL'
allocate(zrel(ntmax),stat=i1)
IF ( i1 /= 0 ) STOP '    Allocate ZREL'

WRITE (6,'(/,5X,A,/)') 'Calculating host (Delta_t)^(-1) matrices'
IF ( krel == 0 ) THEN
  DO i = 1,ntleft+ntright
    irc1 = ircut(ipan(i),i)
    DO i1 = 2,irc1
      dror(i1,i) = drdi(i1,i)/rr(i1,i)
    END DO
  END DO
ELSE
  DO i = 1,ntleft+ntright
    irc1 = ircut(ipan(i),i)
    DO i1 = 1,irc1
      r2drdirel(i1,i) = rr(i1,i)*rr(i1,i)*drdi(i1,i)
    END DO
    zrel(i) = nint(zat(i))
  END DO
END IF

! ******************************************************* energy loop IE
allocate(trefll(lmmaxd,lmmaxd,nref),stat=i1)
IF ( i1 /= 0 ) STOP '    Allocate TREFLL'
allocate(dtrefll(lmmaxd,lmmaxd,nref),stat=i1)     ! LLY
IF ( i1 /= 0 ) STOP '    Allocate DTREFLL'        ! LLY
allocate(tmatll(lmmaxd,lmmaxd),dhmat(lmmaxd,lmmaxd,2),stat=i1)
IF ( i1 /= 0 ) STOP '    Allocate TMATLL/DHMAT'
allocate( alpharef(0:lmaxd,nref),dalpharef(0:lmaxd,nref) ,stat=i1) ! LLY Lloyd Alpha matrix AND deriv.
IF ( i1 /= 0 ) STOP '    Allocate ALPHAREF/DALPHAREF'

DO ie = 1,ielast
  eryd = ez(ie)
  
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
! -> set up t matrices for the reference system
  
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  IF ( krel == 0 ) THEN
    DO i1 = 1,nref
      CALL calctref13(eryd,vref(i1),rmtref(i1),lmaxd,ih,  &
          trefll(1,1,i1),dtrefll(1,1,i1),  &
          alpharef(0,i1),dalpharef(0,i1),lmaxd+1,lmgf0d)
    END DO
  END IF
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
  DO ilhost = 1,nhost
    ihost = inhost(ilhost)
    iqoff = 0
    itoff = 0
    nqhost = nlbasis
    IF ( ihost == 2 ) THEN
      nqhost = nrbasis
      iqoff = nlbasis
      itoff = ntleft
    END IF
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ sites in host
    DO ih = 1,nqhost
      
! -> assign Delta_t = -t_ref
      
!  Note: REFPOT(1) = REFPOT(NAEZ) (i.e., of the 2D system)
      
      iqh = iqoff + ih
      iref = refpot(iqh+1)
      DO lm2 = 1,lmmaxd
        DO lm1 = 1,lmmaxd
          dhmat(lm1,lm2,1) = -trefll(lm1,lm2,iref)
        END DO
      END DO
      
      IF ( nspinh(ihost) > 1 ) THEN
        DO lm2 = 1,lmmaxd
          CALL zcopy(lmmaxd,dhmat(1,lm2,1),1, dhmat(1,lm2,2),1)
        END DO
      END IF
! ====================================================== spins and atoms
      DO ispin = 1,nspinh(ihost)
! ----------------------------------------------------------------------
        DO ioq = 1,noq(iqh)
          i1 = kaoez(ioq,iqh) + itoff
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -> calculate t_sys for the atom I1 located on site IH
          
          ipot = (i1-1) * nspinh(ihost) + ispin
          irc1 = ircut(ipan(i1),i1)
          rirc = rr(irc1,i1)
          
          CALL decitmat(eryd,zat(i1),ipan(i1),rr(1,i1),  &
              dror(1,i1),visp(1,ipot),ircut(0,i1),  &
              rirc,krel,nsra,ins,tmatll,loflm, idoldau,lopt,wldauav,  &
              solver,socscl(1,krel*i1+(1-krel)),  &
              cscl(1,krel*i1+(1-krel)),zrel(i1), vtrel(1,i1),btrel(1,i1),  &
              drdi(1,i1),r2drdirel(1,i1), ipand,irmd,lmaxd,lmaxd+1,lm2d,lmmaxd)
          
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tmat calculated
          
! -> Delta_t = Delta_t + CONC(I1)*t_mat(I1)
          
          carg = conc(i1)
          DO lm2 = 1,lmmaxd
            CALL zaxpy(lmmaxd,carg,tmatll(1,lm2),1, dhmat(1,lm2,ispin),1)
          END DO
        END DO
! ----------------------------------------------------------------------
! tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
        IF ( test('tmat    ') ) THEN
          WRITE (1337,*)
          WRITE (1337,99002) '      ---> Delta_t  matrix for site: ',iqh
          IF ( krel == 0 ) WRITE (1337,99003) txts(ispin)
          WRITE (1337,99004) ', energy: ',eryd
          CALL cmatstr(' ',1,dhmat(1,1,ispin),lmmaxd,lmmaxd,  &
              2*krel+1,2*krel+1,0,1D-8,6)
          WRITE (1337,*)
        END IF
! tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
        
! --> inversion
        
        CALL zgetrf(lmmaxd,lmmaxd,dhmat(1,1,ispin),lmmaxd, ipvt,info)
        CALL zgetri(lmmaxd,dhmat(1,1,ispin),lmmaxd,ipvt,  &
            tmatll,lmmaxd*lmmaxd,info)
      END DO
! ======================================================================
      
! --> scaling the host t-matrices to p.u.
      
      IF ( ihost == 1 ) THEN
        DO ispin = 1,nspinh(ihost)
          DO lm2 = 1,lmmaxd
            DO lm1 = 1,lmmaxd
              lefttinv(lm1,lm2,ih,ispin,ie) = cfctor * dhmat(lm1,lm2,ispin)
            END DO
          END DO
        END DO
      ELSE
        DO ispin = 1,nspinh(ihost)
          DO lm2 = 1,lmmaxd
            DO lm1 = 1,lmmaxd
              righttinv(lm1,lm2,ih,ispin,ie) = cfctor * dhmat(lm1,lm2,ispin)
            END DO
          END DO
        END DO
      END IF
    END DO
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  END DO
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END DO
! **********************************************************************
deallocate(zat,rws,rmt,conc,rr,drdi,visp,stat=i1)
IF ( i1 /= 0 ) STOP '   Deallocate ZAT/RWS/RMT/.../VISP'
deallocate(irws,ipan,iqat,ircut,loflm,stat=i1)
IF ( i1 /= 0 ) STOP '   Deallocate IRWS/IPAN/IQAT/IRCUT/LOFLM'
deallocate(trefll,tmatll,dhmat,stat=i1)
IF ( i1 /= 0 ) STOP '   Deallocate TREFLL/TMATLL/DHMAT'
deallocate(socscl,cscl,vtrel,btrel,stat=i1)
IF ( i1 /= 0 ) STOP '   Deallocate SOCSCL/CSCL/VTREL/BTREL'
deallocate(alpharef,dalpharef ,stat=i1)
IF ( i1 /= 0 ) STOP '   Deallocate ALPHAREF/DALPHAREF'
IF ( krel == 0 ) THEN
  deallocate(dror,stat=i1)
  IF ( i1 /= 0 ) STOP '   Deallocate DROR'
ELSE
  deallocate(r2drdirel,zrel,stat=i1)
  IF ( i1 /= 0 ) STOP '   Deallocate R2DRDIREL/ZREL'
END IF
99002 FORMAT (a,i3,$)
99003 FORMAT (', ',a,$)
99004 FORMAT (a,2F10.6)
99005 FORMAT(10X,'KREL= ',i1,' NSPIN= ',i1,' INS= ',i1,' KMROT= ',i1,/,  &
    10X,'NAEZ=',i3,' ALAT= ',f9.6,' EFERMI= ',f9.6)
99006 FORMAT(10X,'BRAVAIS '/10X,3F8.4/10X,3F8.4/10X,3F8.4/10X,'RBASIS')
END SUBROUTINE decitset
