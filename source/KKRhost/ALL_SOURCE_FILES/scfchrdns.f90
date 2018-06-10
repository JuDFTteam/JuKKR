SUBROUTINE scfchrdns(nfilcbwf,r2drdi,jws,imt,shftef,totdos,muespn,  &
    mueorb,irel,iprint,nt,nl,nkm,eryd,we,efermi,  &
    iecurr,netab,dos,smt,omt,hff,dosi,smti,omti,  &
    hffi,dosm,dosl0,dosint,smtm,smtl0,smtint,  &
    omtm,omtl0,omtint,hffm,hffl0,hffint,bcor,  &
    bcors,dzz,dzj,szz,szj,ozz,ozj,bzz,bzj,  &
    ozzs,ozjs,omtls0,tautlin,nvaltot,txtt,  &
    conc,nat,rhochr,rhospn,rhoorb,qel,gdia,gmdia,  &
    goff,ntmax,nlmax,nmuemax,linmax,nrmax,nmmax, nkmmax,eband,ebandt)
!   ********************************************************************
!   *                                                                  *
!   * SUBROUTINE TO CALCULATE THE  CHARGE, SPIN  AND  ORBITAL DENSITY  *
!   *                  WITHIN AN ATOMIC CELL                           *
!   *                                                                  *
!   * 12/03/96 HE                                                      *
!   ********************************************************************
IMPLICIT NONE

! PARAMETER definitions
INTEGER NTMAXCHK
PARAMETER (NTMAXCHK=10)
DOUBLE COMPLEX C0
PARAMETER (C0=(0.0D0,0.0D0))
DOUBLE PRECISION PI
PARAMETER (PI=3.141592653589793238462643D0)

! Dummy arguments
DOUBLE COMPLEX EBAND,ERYD,WE
DOUBLE PRECISION EFERMI,MUEORB,MUESPN,SHFTEF,TOTDOS,NVALTOT
INTEGER IECURR,IPRINT,IREL,LINMAX,NETAB,NFILCBWF,NKM,NKMMAX,NL, &
        NLMAX,NMMAX,NMUEMAX,NRMAX,NT,NTMAX
DOUBLE PRECISION BCOR(NTMAX),BCORS(NTMAX),CONC(NTMAX),DOS(NTMAX), &
       DOSI(NTMAX),GDIA(NKMMAX),GMDIA(NKMMAX),GOFF(NKMMAX), &
       HFF(NTMAX),HFFI(NTMAX),OMT(NTMAX),OMTI(NTMAX),QEL(NTMAX), &
       R2DRDI(NRMAX,NMMAX),RHOCHR(NRMAX,NTMAX),RHOORB(NRMAX,NTMAX) &
       ,RHOSPN(NRMAX,NTMAX),SMT(NTMAX),SMTI(NTMAX)
DOUBLE COMPLEX BZJ(LINMAX,NTMAX),BZZ(LINMAX,NTMAX), &
           DOSINT(NLMAX,NTMAX) &
           ,DOSL0(NLMAX,NTMAX),DOSM(NMUEMAX),DZJ(LINMAX,NTMAX), &
           DZZ(LINMAX,NTMAX),EBANDT(NTMAX),HFFINT(NLMAX,NTMAX), &
           HFFL0(NLMAX,NTMAX),HFFM(NMUEMAX),OMTINT(NLMAX,NTMAX), &
           OMTL0(NLMAX,NTMAX),OMTM(NMUEMAX),OZJ(LINMAX,NTMAX), &
           OZZ(LINMAX,NTMAX),SMTINT(NLMAX,NTMAX), &
           SMTL0(NLMAX,NTMAX),SMTM(NMUEMAX),SZJ(LINMAX,NTMAX), &
           SZZ(LINMAX,NTMAX),TAUTLIN(LINMAX,NTMAX)
DOUBLE COMPLEX OZZS(LINMAX,NTMAX,2),OZJS(LINMAX,NTMAX,2), &
           OMTLS0(NLMAX,NTMAX,2)
INTEGER IMT(NTMAX),JWS(NMMAX),NAT(NTMAX)
CHARACTER (len=4) :: TXTT(NTMAX)

! Local variables
DOUBLE PRECISION AUX,BDUM(3),CFF(2,2),CFG(2,2),CGF(2,2),CGG(2,2), &
       CHKO(NTMAXCHK),CHKQ(NTMAXCHK),CHKS(NTMAXCHK),DEFERMI,DQ,MJ, &
       MJMAX,MJMIN,R1M(2,2),RINT(NRMAX),TOTNOS
DOUBLE PRECISION DBLE,DSQRT
DOUBLE COMPLEX DOSL,HFFL,JF(NRMAX,2,2),JG(NRMAX,2,2),OMTL,SMTL, &
           WDS,WOF,WOG,WSF,WSG,WT,ZF(NRMAX,2,2),ZFJF,ZFZF, &
           ZG(NRMAX,2,2),ZGJG,ZGZG
DOUBLE COMPLEX OMTLS(2),OMTMS(NMUEMAX,2)
INTEGER I,IFLAG,IKM1,IKM2,IL,IM,IS,IT,JJ,JTOP,K1,K2,KA, &
        KAP1,KAP2,KB,L,LIN,LMAX,MM,MUE,NSOL
INTEGER IKAPMUE, IMJ
INTEGER NINT

SAVE CHKO,CHKQ,CHKS

DATA R1M/1.0D0,0.0D0,0.0D0,1.0D0/

IF ( iecurr == 1 ) THEN
  
  DO it = 1,nt
    DO il = 1,nl
      dosint(il,it) = c0
      smtint(il,it) = c0
      omtint(il,it) = c0
      hffint(il,it) = c0
    END DO
    ebandt(it) = c0
  END DO
  
! ----------------------------- account for spin degeneracy for IREL <=1
  IF ( irel <= 1 ) THEN
    DO it = 1,nt
      im = imt(it)
      jtop = jws(im)
      DO i = 1,jtop
        rhochr(i,it) = rhochr(i,it)/2.0D0
        rhospn(i,it) = 0.0D0
        rhoorb(i,it) = 0.0D0
      END DO
    END DO
  END IF
  
  IF ( iprint > 0 ) THEN
    IF ( nt > ntmaxchk ) STOP '<SCFCHRDNS> NT > NTMAXCHK'
    DO it = 1,nt
      im = imt(it)
      jtop = jws(im)
      DO i = 1,jtop
        rint(i) = rhochr(i,it)*r2drdi(i,im)
      END DO
      CALL rintsimp(rint,jtop,chkq(it))
      DO i = 1,jtop
        rint(i) = rhospn(i,it)*r2drdi(i,im)
      END DO
      CALL rintsimp(rint,jtop,chks(it))
      DO i = 1,jtop
        rint(i) = rhoorb(i,it)*r2drdi(i,im)
      END DO
      CALL rintsimp(rint,jtop,chko(it))
    END DO
  END IF
END IF

100  CONTINUE
totnos = 0.0D0
totdos = 0.0D0
muespn = 0.0D0
mueorb = 0.0D0


! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
DO it = 1,nt
  im = imt(it)
  jtop = jws(im)
  
  lmax = nl - 1
  lin = 0
  
  dos(it) = 0.0D0
  smt(it) = 0.0D0
  omt(it) = 0.0D0
  hff(it) = 0.0D0
  dosi(it) = 0.0D0
  smti(it) = 0.0D0
  omti(it) = 0.0D0
  hffi(it) = 0.0D0
  
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
  DO l = 0,lmax
    il = l + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    kap1 = -l - 1
    kap2 = l
    IF ( l == 0 ) kap2 = kap1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    dosl = 0.0D0
    smtl = 0.0D0
    omtl = 0.0D0
    hffl = 0.0D0
    omtls(1) = omtl
    omtls(2) = omtl
    
    IF ( irel > 1 ) THEN
      mjmax = DBLE(l) + 0.5D0
    ELSE
      mjmax = DBLE(l)
    END IF
    mjmin = -mjmax
    mue = 0
    
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!DO MJ = MJMIN,MJMAX,1.0D0
    DO imj = nint(mjmin),nint(mjmax),1
      mj = DBLE(imj)
      mue = mue + 1
      dosm(mue) = 0.0D0
      smtm(mue) = 0.0D0
      omtm(mue) = 0.0D0
      hffm(mue) = 0.0D0
      
      DO is = 1,2
        omtms(mue,is) = 0.0D0
      END DO
      
      IF ( irel <= 1 ) THEN
        nsol = 1
!           no coupling for:  abs(mue)= j   +  j=l+1/2 == kap=-l-1
      ELSE IF ( ABS(mj) > DBLE(l) ) THEN
        nsol = 1
      ELSE
        nsol = 2
      END IF
!------------------------------------------------------------------------
      ikm1 = ikapmue(kap1,nint(mj-0.5D0))
      ikm2 = ikapmue(kap2,nint(mj-0.5D0))
!------------------------------------------------------------------------
      IF ( irel <= 1 ) THEN
        ikm1 = il
        ikm2 = il
        IF ( nkm /= nl**2 ) WRITE (1337,99001) nkm
      END IF
      
      
!   COEFFICIENTS TO CALCULATE THE SPIN  MAGNETISATION
      
      cgg(1,1) = gdia(ikm1)
      cgg(1,2) = goff(ikm1)
      cgg(2,1) = goff(ikm1)
      cgg(2,2) = gdia(ikm2)
      CALL rinit(4,cgf)
      cgf(1,1) = gmdia(ikm1)
      cgf(2,2) = gmdia(ikm2)
      
!   COEFFICIENTS TO CALCULATE THE ORBITAL MAGNETISATION
      
      cfg(1,1) = mj*(kap1+1.0D0)/(kap1+0.5D0)
      cfg(2,2) = mj*(kap2+1.0D0)/(kap2+0.5D0)
      cfg(1,2) = 0.5D0*DSQRT(1.0D0-(mj/(kap1+0.5D0))**2)
      cfg(2,1) = cfg(1,2)
      CALL rinit(4,cff)
      cff(1,1) = mj*(-kap1+1.0D0)/(-kap1+0.5D0)
      cff(2,2) = mj*(-kap2+1.0D0)/(-kap2+0.5D0)
      
! -------------------------------------------------- read wave functions
      
      CALL readwfun(nfilcbwf,it,l,mj,nsol,'REG','IRR',  &
          ikm1,kap1,ikm2,kap2,nt,nkm,zg,zf,jg,jf, jtop,nrmax)
      
      DO k1 = 1,nsol
        DO k2 = 1,nsol
          lin = lin + 1
          wt = -tautlin(lin,it)/pi
          dosm(mue) = dosm(mue) + wt*dzz(lin,it)
          smtm(mue) = smtm(mue) + wt*szz(lin,it)
          omtm(mue) = omtm(mue) + wt*ozz(lin,it)
          hffm(mue) = hffm(mue) + wt*bzz(lin,it)
          
          DO is = 1,2
            omtms(mue,is) = omtms(mue,is) + wt*ozzs(lin,it,is)
          END DO
          
          DO ka = 1,nsol
            DO kb = 1,nsol
              wds = we*wt*r1m(ka,kb)
              wsg = we*wt*cgg(ka,kb)
              wsf = we*wt*cgf(ka,kb)
              wog = we*wt*cfg(ka,kb)
              wof = we*wt*cff(ka,kb)
              DO i = 1,jtop
                zgzg = zg(i,ka,k1)*zg(i,kb,k2)
                zfzf = zf(i,ka,k1)*zf(i,kb,k2)
                rhochr(i,it) = rhochr(i,it) + DIMAG(wds*zgzg+wds*zfzf)
                rhospn(i,it) = rhospn(i,it) + DIMAG(wsg*zgzg-wsf*zfzf)
                rhoorb(i,it) = rhoorb(i,it) + DIMAG(wog*zgzg-wof*zfzf)
              END DO
            END DO
          END DO
          
!    NO IRREGULAR CONTRIBUTIONS TO THE BACKSCATTERING TERMS
          IF ( k1 == k2 ) THEN
            dosm(mue) = dosm(mue) + dzj(lin,it)/pi
            smtm(mue) = smtm(mue) + szj(lin,it)/pi
            omtm(mue) = omtm(mue) + ozj(lin,it)/pi
            hffm(mue) = hffm(mue) + bzj(lin,it)/pi
            
            DO is = 1,2
              omtms(mue,is) = omtms(mue,is) + ozjs(lin,it,is)/pi
            END DO
            
            DO ka = 1,nsol
              DO kb = 1,nsol
                wds = we*r1m(ka,kb)/pi
                wsg = we*cgg(ka,kb)/pi
                wsf = we*cgf(ka,kb)/pi
                wog = we*cfg(ka,kb)/pi
                wof = we*cff(ka,kb)/pi
                DO i = 1,jtop
                  zgjg = zg(i,ka,k1)*jg(i,kb,k2)
                  zfjf = zf(i,ka,k1)*jf(i,kb,k2)
                  rhochr(i,it) = rhochr(i,it) + DIMAG(wds*zgjg+wds*zfjf)
                  rhospn(i,it) = rhospn(i,it) + DIMAG(wsg*zgjg-wsf*zfjf)
                  rhoorb(i,it) = rhoorb(i,it) + DIMAG(wog*zgjg-wof*zfjf)
                END DO
              END DO
            END DO
          END IF
        END DO
      END DO
      
      
      dosl = dosl + dosm(mue)
      smtl = smtl + smtm(mue)
      omtl = omtl + omtm(mue)
      hffl = hffl + hffm(mue)
      
      DO is = 1,2
        omtls(is) = omtls(is) + omtms(mue,is)
      END DO
      
    END DO
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    
    ebandt(it) = ebandt(it) + we*dosl*eryd
    
    dosint(il,it) = dosint(il,it) + we*dosl
    smtint(il,it) = smtint(il,it) + we*smtl
    omtint(il,it) = omtint(il,it) + we*omtl
    hffint(il,it) = hffint(il,it) + we*hffl
    
    dosl0(il,it) = dosl
    smtl0(il,it) = smtl
    omtl0(il,it) = omtl
    hffl0(il,it) = hffl
    
    DO is = 1,2
      omtls0(il,it,is) = omtls(is)
    END DO
! ----------------------------------------------------------------------
    IF ( (iprint > 0) .AND. (iecurr == netab) ) THEN
      
      IF ( irel > 1 ) THEN
        jj = 2*l + 2
        
        WRITE (1337,99005) iecurr,eryd,l,it,txtt(it),  &
            'CRYSTAL TERMS       ',dosint(il,it), smtint(il,it),omtint(il,it),  &
            (hffint(il,it)*1D-3),dosl,smtl,omtl, (hffl*1D-6),  &
            ((-jj-1+2*mue),dosm(mue),smtm(mue),  &
            omtm(mue),(hffm(mue)*1D-6),mue=1,jj)
        
      ELSE
        jj = 2*l + 1
        
        WRITE (1337,99014) iecurr,eryd,l,it,txtt(it),  &
            'CRYSTAL TERMS       ',dosint(il,it), smtint(il,it),omtint(il,it),  &
            (hffint(il,it)*1D-3),dosl,smtl,omtl, (hffl*1D-6),  &
            ((-l-1+mm),dosm(mm),smtm(mm),omtm(mm), (hffm(mm)*1D-6),mm=1,jj)
        
      END IF
    END IF
! ---------------------------------------------------------------------
    
    
    dos(it) = dos(it) + DIMAG(dosl)
    smt(it) = smt(it) + DIMAG(smtl)
    omt(it) = omt(it) + DIMAG(omtl)
    hff(it) = hff(it) + DIMAG(hffl)
    dosi(it) = dosi(it) + DIMAG(dosint(il,it))
    smti(it) = smti(it) + DIMAG(smtint(il,it))
    omti(it) = omti(it) + DIMAG(omtint(il,it))
    hffi(it) = hffi(it) + DIMAG(hffint(il,it))
  END DO
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
  
  
  totnos = totnos + dosi(it)*conc(it)*nat(it)
  totdos = totdos + dos(it)*conc(it)*nat(it)
  muespn = muespn + smti(it)*conc(it)*nat(it)
  mueorb = mueorb + omti(it)*conc(it)*nat(it)
  
  IF ( iprint > 0 ) THEN
    
    WRITE (1337,99006) iecurr,eryd,it,txtt(it),dosi(it),  &
        smti(it),omti(it),(hffi(it)*1D-3),dos(it), smt(it),omt(it),(hff(it)*1D-3)
    
    IF ( it < nt ) THEN
      WRITE (1337,'(1X,79(''-''))')
    ELSE IF ( (iprint > 0) .OR. (iecurr == netab) ) THEN
      WRITE (1337,99010) totdos,totnos,muespn,mueorb
    ELSE
      WRITE (1337,'('' '',79(''=''))')
    END IF
  END IF
  
END DO
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

IF ( iecurr == netab ) THEN
  
  eband = c0
  DO it=1,nt
    eband = eband + ebandt(it)*conc(it)*nat(it)
  END DO
  
  IF ( irel > 1 ) THEN
    dq = nvaltot - totnos
  ELSE
    dq = nvaltot/2.0D0 - totnos
  END IF
  
  defermi = dq/totdos
  
  IF ( ABS(dq) > 1D-06 ) THEN
    we = defermi
    efermi = efermi + defermi
    shftef = defermi
    
    WRITE (1337,'(/)')
    WRITE (1337,99012) (txtt(it),conc(it),it=1,nt)
    WRITE (1337,99013) dq,defermi,efermi
    
    GO TO 100
  END IF
  
  IF ( iprint > 0 ) THEN
    iflag = 0
    DO it = 1,nt
      im = imt(it)
      jtop = jws(im)
      DO i = 1,jtop
        rint(i) = rhochr(i,it)*r2drdi(i,im)
      END DO
      CALL rintsimp(rint,jtop,aux)
      chkq(it) = aux - chkq(it)
      IF ( ABS(chkq(it)-dosi(it)) > 1.0D-8 ) THEN
        iflag = 1
        WRITE (1337,99004) it,'Q',dosi(it),chkq(it),dosi(it) /chkq(it)
      END IF
      DO i = 1,jtop
        rint(i) = rhospn(i,it)*r2drdi(i,im)
      END DO
      CALL rintsimp(rint,jtop,aux)
      chks(it) = aux - chks(it)
      IF ( ABS(chks(it)-smti(it)) > 1.0D-8 ) THEN
        iflag = 1
        WRITE (1337,99004) it,'S',smti(it),chks(it),smti(it) /chks(it)
      END IF
      DO i = 1,jtop
        rint(i) = rhoorb(i,it)*r2drdi(i,im)
      END DO
      CALL rintsimp(rint,jtop,aux)
      chko(it) = aux - chko(it)
      IF ( ABS(chko(it)-omti(it)) > 1.0D-8 ) THEN
        iflag = 1
        WRITE (1337,99004) it,'O',omti(it),chko(it),omti(it) /chko(it)
      END IF
    END DO
    
    IF ( iflag == 0 ) THEN
      WRITE (1337,99002)
    ELSE
      WRITE (1337,99003)
    END IF
  END IF
  
  
  DO it = 1,nt
    
    WRITE (1337,99006) (iecurr+1),efermi,0.0D0,it,txtt(it)
    
    bdum(1) = bcors(it)*1D-3
    bdum(2) = (bcor(it)-bcors(it))*1D-3
    bdum(3) = bcor(it)*1D-3
    
    WRITE (1337,99007) (DIMAG(dosl0(il,it)),DIMAG(dosint(il,it))  &
        ,DIMAG(smtl0(il,it)),DIMAG(smtint(il,it)),  &
        DIMAG(omtl0(il,it)),DIMAG(omtint(il,it)),  &
        DIMAG(hffint(il,it))*1D-3,bdum(il),il=1, MIN(3,nl))
    IF ( nl > 3 ) WRITE (1337,99008) (DIMAG(dosl0(il,it)),DIMAG(dosint(il,  &
        it)),DIMAG(smtl0(il,it)), DIMAG(smtint(il,it)),  &
        DIMAG(omtl0(il,it)), DIMAG(omtint(il,it)),  &
        DIMAG(hffint(il,it))*1D-3,il=4,nl)
    
    WRITE (1337,99009) dos(it),dosi(it),smt(it),smti(it),omt(it)  &
        ,omti(it),(hffi(it)*1D-3), ((hffi(it)+bcor(it))*1D-3)
    
    IF ( it < nt ) THEN
      WRITE (1337,'(1X,79(''-''))')
    ELSE
      WRITE (1337,99011) totdos,totnos,muespn,mueorb, DIMAG(eband)
    END IF
    
    im = imt(it)
    jtop = jws(im)
    
    DO i = 1,jtop
      rint(i) = rhochr(i,it)*r2drdi(i,im)
    END DO
    
    CALL rintsimp(rint,jtop,qel(it))
    
! ----------------------------- account for spin degeneracy for IREL <=1
    IF ( irel <= 1 ) THEN
      qel(it) = qel(it)*2.0D0
      DO i = 1,jtop
        rhochr(i,it) = rhochr(i,it)*2.0D0
        rhospn(i,it) = 0.0D0
        rhoorb(i,it) = 0.0D0
      END DO
    END IF
    
  END DO
  
END IF

99001 FORMAT ('warning in <SCFCHRDNS>:  IREL<=1 and  NL**2 <> NKM=',i5)
99002 FORMAT (/,10X,'integrals in <SCFCHRDNS> agree within 1D-8',/)
99003 FORMAT (/,10X,'... integrals in <SCFCHRDNS>  NOT OK')
99004 FORMAT (' IT ',i3,2X,a,2X,f20.10,/,12X,4F20.10)
99005 FORMAT (/,i4,' E=',2F7.4,3X,'L=',i2,3X,'IT=',i2,2X,a,2X,a20,/,  &
    15X,'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |',  &
    '   B_tot   [kG]',/,' INT(DE)  ',2F8.3,2X,2F8.3,2X,2F8.3,  &
    f10.1,f8.1,/,' SUM(MJ)  ',2F8.3,2X,2F8.3,2X,2F8.3,f10.1,  &
    f8.1,20(:,/,' MJ= ',i2,'/2 ',2F8.3,2X,2F8.3,2X,2F8.3, f10.1,f8.1))
99006 FORMAT (/,i4,' E=',2F7.4,10X,'IT=',i2,2X,a,:,/,15X,  &
    'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |',  &
    '   B_tot   [kG]',/,' INT(DE) crystal  ',f8.3,10X,f8.3,  &
    10X,f8.3,10X,f8.1,/,' TOTAL   crystal  ',f8.3,10X,f8.3, 10X,f8.3,10X,f8.1)
99007 FORMAT ('         DOS      NOS     P_spin   m_spin',  &
    '    P_orb    m_orb    B_val      B_core',/,'  s ',2F9.4,  &
    f10.4,f9.4,f10.5,f9.5,f8.2,' s  ',f8.2,:,/,'  p ',2F9.4,  &
    f10.4,f9.4,f10.5,f9.5,f8.2,' ns ',f8.2,:,/,'  d ',2F9.4,  &
    f10.4,f9.4,f10.5,f9.5,f8.2,' cor',f8.2)
99008 FORMAT ('  f ',2F9.4,f10.4,f9.4,f10.5,f9.5,f8.2,:,/,'  g ',2F9.4,  &
    f10.4,f9.4,f10.5,f9.5,f8.2)
99009 FORMAT (' sum',2F9.4,f10.4,f9.4,f10.5,f9.5,f8.2,' v+c',f8.2)
99010 FORMAT (' ',79('-'),/,' TDOS/NOS ',2F8.3,' MUE-SPIN:',f8.3,  &
    '  MUE-ORB:',f8.3)
99011 FORMAT (' ',79('-'),/,' TOT',2F9.4,10X,f9.4,10X,f9.5,/,' E_band',  &
    f17.6,' [Ry]',/,' ',79('='))
99012 FORMAT ((' ',79('*'),/),/,' KKR-run for: ',15(a,f5.2))
99013 FORMAT (/,' results extrapolated to corrected FERMI - ENERGY:',/,  &
    ' CHARGE MISFIT     ',f9.5,/,' E_F CORRECTION    ',f9.5,/,  &
    ' NEW FERMI ENERGY  ',f9.5,/)
99014 FORMAT (/,i4,' E=',2F7.4,3X,'L=',i2,3X,'IT=',i2,2X,a,2X,a20,/,  &
    15X,'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |',  &
    '   B_tot   [kG]',/,' INT(DE)  ',2F8.3,2X,2F8.3,2X,2F8.3,  &
    f10.1,f8.1,/,' SUM(ML)  ',2F8.3,2X,2F8.3,2X,2F8.3,f10.1,  &
    f8.1,20(:,/,' ML= ',i2,'   ',2F8.3,2X,2F8.3,2X,2F8.3, f10.1,f8.1))
END SUBROUTINE scfchrdns
