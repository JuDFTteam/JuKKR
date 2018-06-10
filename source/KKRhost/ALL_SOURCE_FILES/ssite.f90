SUBROUTINE ssite(iwrregwf,iwrirrwf,nfilcbwf,calcint,getirrsol,  &
    soctl,ctl,eryd,p,ihyper,iprint,ikm1lin,ikm2lin,  &
    nlq,nkmq,nlinq,nt,nkm,iqat,tsst,msst,tsstlin,dzz,  &
    dzj,szz,szj,ozz,ozj,bzz,bzj,qzz,qzj,tzz,tzj,vt,  &
    bt,at,z,nucleus,r,drdi,r2drdi,jws,imt,  &
    ameopo,lopt,solver,cgc,ozzs,ozjs,nlmax,nqmax,  &
    linmax,nrmax,nmmax,ntmax,nkmmax,nkmpmax,nlamax)
!   ********************************************************************
!   *                                                                  *
!   * ASSIGN QUANTUM NUMBERS AND CALL ROUTINE TO SOLVE                 *
!   * 8 COUPLED PARTIAL DIFFERENTIAL RADIAL DIRAC EQUATIONS. THE       *
!   * RESULTING WAVEFUNCTIONS ARE USED TO CALCULATE T-MATRICES IN      *
!   * THE KAPPA-MU REPRESENTATION                                      *
!   *                                                                  *
!   * + CALCULATION OF THE RADIAL INTEGRALS                            *
!   *   [ G1*G2 + F1*F2 ] R**2 DR                                      *
!   *                                                                  *
!   * FOR IHYPER <> 0 :                                                *
!   * CALCULATION OF THE HYPERFINE MATRIX ELEMENTS                     *
!   *                                                                  *
!   * RYD-UNITS USED THROUGHOUT                                        *
!   *                                                                  *
!   * NOTE: to save storage force  JG/JF  and  PR/QR  to share the     *
!   *       same storage by corresponding argument list in CALL ....   *
!   *                                                                  *
!   * 28/10/94  HE  tidy up,  P,Q used in <DIRAC> instead of g,f       *
!   * 05/10/96  HE  workspace for wavefunctions and matrices           *
!   *               is allocated dynamically !!!                       *
!   * 07/02/05  VP  few changes connected to the calculation of orbital*
!   *               polarisation                                       *
!   ********************************************************************


IMPLICIT NONE

!PARAMETER definitions
COMPLEX*16 CI
PARAMETER ( CI =(0.0D0, 1.0D0) )

REAL*8 F1,E0,A0,CAUTOG

! conversion factor for hyperfine fields from A.U. to GAUSS
!                                 electron charge     in esu
!                                 Bohr-radius         in cm

PARAMETER (F1=1.0D0,E0=1.6021892D-19*2.997930D+09, &
           A0=0.52917706D-08,CAUTOG=E0/(A0*A0))

! Dummy arguments
LOGICAL CALCINT,GETIRRSOL
COMPLEX*16 ERYD,P
INTEGER IHYPER,IPRINT,IWRIRRWF,IWRREGWF,LINMAX,NFILCBWF,NKM, &
        NKMMAX,NLAMAX,NLMAX,NMMAX,NQMAX,NRMAX,NT,NTMAX,NUCLEUS
INTEGER NKMPMAX
CHARACTER*10 SOLVER
REAL*8 AMEOPO(NKMMAX,NKMMAX,NLAMAX,3),AT(NRMAX,NLAMAX,3,NTMAX), &
       BT(NRMAX,NTMAX),CTL(NTMAX,NLMAX),DRDI(NRMAX,NMMAX), &
       R(NRMAX,NMMAX),R2DRDI(NRMAX,NMMAX),SOCTL(NTMAX,NLMAX), &
       VT(NRMAX,NTMAX)
REAL*8 CGC(NKMPMAX,2)
COMPLEX*16 OZJS(LINMAX,NTMAX,2),OZZS(LINMAX,NTMAX,2)
COMPLEX*16 BZJ(LINMAX,NTMAX),BZZ(LINMAX,NTMAX),DZJ(LINMAX,NTMAX), &
           DZZ(LINMAX,NTMAX),MSST(NKMMAX,NKMMAX,NTMAX), &
           OZJ(LINMAX,NTMAX),OZZ(LINMAX,NTMAX),QZJ(LINMAX,NTMAX), &
           QZZ(LINMAX,NTMAX),SZJ(LINMAX,NTMAX),SZZ(LINMAX,NTMAX), &
           TSST(NKMMAX,NKMMAX,NTMAX),TSSTLIN(LINMAX,NTMAX), &
           TZJ(LINMAX,NTMAX),TZZ(LINMAX,NTMAX)
INTEGER IKM1LIN(LINMAX),IKM2LIN(LINMAX),IMT(NTMAX), &
        IQAT(NQMAX,NTMAX),JWS(NMMAX),LOPT(NTMAX),NKMQ(NQMAX), &
        NLINQ(NQMAX),NLQ(NQMAX),Z(NTMAX),MUEM05

! Local variables
COMPLEX*16 A(2,2),ARG,B1,B2,CINT(NRMAX),CRSQ,CSUM,DET,DXP(2,2), &
           F11,F12,F21,F22,G11,G11P,G12,G12P,G21,G21P,G22,G22P, &
           GAM(2,2),GAMINV(2,2),HL,HLB1,HLB2,JF(NRMAX,2,2), &
           JG(NRMAX,2,2),JL,JLB1,JLB2,JLP,MAUX(NKMMAX,NKMMAX), &
           MSST2(2,2),NL,NLB1,NLB2,NLP,NORM,PFAC,PI(2,2,NRMAX), &
           PR(2,2,NRMAX),QI(2,2,NRMAX),QR(2,2,NRMAX),RMEHF(2,2), &
           RMEHF1(2,2),RMEHF2(2,2),SIG(2,2),TSST2(2,2),XSST2(2,2), &
           ZF(NRMAX,2,2),ZFJF(2,2),ZFZF(2,2),ZG(NRMAX,2,2), &
           ZGJG(2,2),ZGZG(2,2)
REAL*8 AP(2,2,NRMAX),AQ(2,2,NRMAX),C,CFF(2,2),CFG(2,2),CG1,CG2, &
       CG4,CG5,CG8,CGF(2,2),CGG(2,2),CH(2,2),CSQR,CTF(2,2), &
       CTG(2,2),DOVR(NRMAX),DROVRN(NRMAX),MJ,R1M(2,2),RKD(2,2), &
       RNUC,SK1,SK2,TDIA1,TDIA2,TOFF
REAL*8 COG(2,2,2),COF(2,2,2)
COMPLEX*16 CDJLZDZ,CDNLZDZ,CJLZ,CNLZ
DOUBLE PRECISION DBLE,DSQRT
INTEGER I,I1,I2,I3,I5,IKM1,IKM2,IL,IM,IN,INFO,IPIV(NKMMAX),IQ, &
        ISK1,ISK2,IT,J,JLIM,JTOP,K,K1,K2,KAP1,KAP2,KC,L,L1,LB1, &
        LB2,LIN,N,NSOL,IMKM1,IMKM2,IS, IMJ
INTEGER IKAPMUE
INTEGER ISIGN,NINT
REAL*8 RNUCTAB
LOGICAL WRONSKI

DATA R1M/1.0D0,0.0D0,0.0D0,1.0D0/
DATA RKD/1.0D0,0.0D0,0.0D0, - 1.0D0/
!     DATA RKD / 1.0D0, 0.0D0, 0.0D0, 1.0D0 /

CALL cinit(ntmax*linmax,dzz)
CALL cinit(ntmax*linmax,dzj)
CALL cinit(ntmax*linmax,szz)
CALL cinit(ntmax*linmax,szj)
CALL cinit(ntmax*linmax,ozz)
CALL cinit(ntmax*linmax,ozj)
CALL cinit(ntmax*linmax,bzz)
CALL cinit(ntmax*linmax,bzj)
CALL cinit(ntmax*linmax,qzz)
CALL cinit(ntmax*linmax,qzj)
CALL cinit(ntmax*linmax,tzz)
CALL cinit(ntmax*linmax,tzj)

wronski = .true.
wronski = .false.
!------------------------------------------------------------------------

c = ctl(1,1)
csqr = c*c

!     calculate relativistic momentum

p = SQRT(eryd*(1.0D0+eryd/csqr))

pfac = p/(1.0D0+eryd/csqr)


! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
DO it = 1,nt
  
  IF ( iprint > 0 ) WRITE (1337,'(A,I3,A,10F10.4)') ' SOLVER ',  &
      it,solver,(soctl(it,il),il=1,nlmax)
  
  iq = iqat(1,it)
  im = imt(it)
  jtop = jws(im)
  IF ( nucleus /= 0 ) THEN
    rnuc = rnuctab(z(it))
    in = 1
    DO WHILE ( r(in,im) < rnuc )
      in = in + 1
    END DO
!        RLIM=R(IN,IM)
    jlim = in
    IF ( MOD(jlim,2) == 0 ) jlim = jlim - 1
    rnuc = r(jlim,im)
  END IF
  DO i = 1,jtop
    dovr(i) = drdi(i,im)/r(i,im)
    IF ( nucleus /= 0 ) drovrn(i) = (r(i,im)/rnuc)**3*drdi(i,im)
  END DO
  
  lin = 0
  
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
  DO l = 0,(nlq(iq)-1)
    il = l + 1
    c = ctl(it,il)
    csqr = c*c
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    kap1 = -l - 1
    kap2 = l
    IF ( l == 0 ) kap2 = kap1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    isk1 = ISIGN(1,kap1)
    isk2 = ISIGN(1,kap2)
    sk1 = DBLE(isk1)
    sk2 = DBLE(isk2)
    l1 = l
    lb1 = l - isk1
    lb2 = l - isk2
    
    arg = p*r(jtop,im)
    jl = cjlz(l1,arg)
    jlb1 = cjlz(lb1,arg)
    jlb2 = cjlz(lb2,arg)
    nl = cnlz(l1,arg)
    nlb1 = cnlz(lb1,arg)
    nlb2 = cnlz(lb2,arg)
    hl = jl + ci*nl
    hlb1 = jlb1 + ci*nlb1
    hlb2 = jlb2 + ci*nlb2
    
    IF ( solver(1:7) == 'ABM-SOC' ) THEN
      jlp = cdjlzdz(l1,arg,1)*p
      nlp = cdnlzdz(l1,arg,1)*p
    END IF
    
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!DO MJ = -(DBLE(L)+0.5D0), + (DBLE(L)+0.5D0),1.0D0
    DO imj = 1, 2*l+1
      mj = -(DBLE(l)+0.5D0) + DBLE(imj-1)
      
!------------------------------------------------------------------------
! NO COUPLING FOR:  ABS(MUE)= J   +  J=L+1/2 == KAP=-L-1
      IF ( ABS(mj) >= DBLE(l) ) THEN
        nsol = 1
      ELSE
        nsol = 2
      END IF
!------------------------------------------------------------------------
      muem05 = nint(mj-0.5D0)
      ikm1 = ikapmue(kap1,muem05)
      ikm2 = ikapmue(kap2,muem05)
      imkm1 = ikapmue(-kap1,muem05)
      imkm2 = ikapmue(-kap2,muem05)
!------------------------------------------------------------------------
      i5 = nrmax*2*2
      CALL cinit(i5,zg)
      CALL cinit(i5,zf)
      CALL cinit(i5,jg)
      CALL cinit(i5,jf)
      CALL cinit(i5,pr)
      CALL cinit(i5,qr)
      CALL cinit(i5,pi)
      CALL cinit(i5,qi)
      
      IF ( solver(1:2) == 'BS' ) THEN
        CALL dirbs(getirrsol,ctl(it,il),eryd,l,mj,kap1,  &
            kap2,p,cg1,cg2,cg4,cg5,cg8,vt(1,it),  &
            bt(1,it),z(it),nucleus,r(1,im),drdi(1,im),  &
            dovr,jtop,pr,qr,pi,qi,zg,zf)
      ELSE IF ( solver == 'ABM-BI' ) THEN
        STOP ' < DIRABMBI > : Not implemented. Set SOLVER=BS in inputcard'
!                   CALL DIRABMBI(GETIRRSOL,CTL(IT,IL),IT,ERYD,L,MJ,KAP1,
!      &                          KAP2,P,CG1,CG2,CG4,CG5,CG8,AMEOPC,
!      &                          AMEOPO,VT(1,IT),BT(1,IT),AT,Z(IT),
!      &                          NUCLEUS,R(1,IM),DRDI(1,IM),DOVR,JTOP,PR,
!      &                          QR,PI,QI,ZG,ZF,AP,AQ,NTMAX,NLAMAX,
!      &                          NKMMAX,NRMAX)
      ELSE IF ( solver(1:6) == 'ABM-OP' ) THEN
        CALL dirabmop(getirrsol,ctl(it,il),it,eryd,l,mj,kap1,  &
            kap2,p,cg1,cg2,cg4,cg5,cg8,ameopo,  &
            vt(1,it),bt(1,it),at,z(it),nucleus,  &
            r(1,im),drdi(1,im),dovr,jtop,pr,qr,pi,  &
            qi,zg,zf,ap,aq,lopt(it),ntmax,nlamax, nkmmax,nrmax)
      ELSE IF ( solver == 'ABM-SOC   ' ) THEN
        CALL dirabmsoc(getirrsol,ctl(it,il),soctl(it,il),it,  &
            eryd,l,mj,kap1,kap2,p,cg1,cg2,cg4,cg5,  &
            cg8,vt(1,it),bt(1,it),z(it),nucleus,  &
            r(1,im),drdi(1,im),dovr,jtop,dxp,pr,qr, pi,qi,zg,zf,nrmax)
      ELSE IF ( solver == 'ABM-SOC-II' ) THEN
        CALL dirabmsoc2(getirrsol,ctl(it,il),soctl(it,il),it,  &
            eryd,l,mj,kap1,kap2,p,cg1,cg2,cg4,cg5,  &
            cg8,vt(1,it),bt(1,it),z(it),nucleus,  &
            r(1,im),drdi(1,im),dovr,jtop,dxp,pr, qr,pi,qi,zg,zf,nrmax)
      ELSE
        WRITE (6,*) 'No solver found for: ',solver
        STOP
      END IF
      
!  wavefunctions at the muffin-tin-radius
      
      n = jtop
      
      g11 = pr(1,1,n)/r(n,im)
      f11 = qr(1,1,n)/(r(n,im)*c)
      g21 = pr(2,1,n)/r(n,im)
      f21 = qr(2,1,n)/(r(n,im)*c)
      g22 = pr(2,2,n)/r(n,im)
      f22 = qr(2,2,n)/(r(n,im)*c)
      g12 = pr(1,2,n)/r(n,im)
      f12 = qr(1,2,n)/(r(n,im)*c)
      
      IF ( solver(1:7) == 'ABM-SOC' ) THEN
        g11p = (dxp(1,1)/drdi(n,im)-g11)/r(n,im)
        g21p = (dxp(2,1)/drdi(n,im)-g21)/r(n,im)
        g12p = (dxp(1,2)/drdi(n,im)-g12)/r(n,im)
        g22p = (dxp(2,2)/drdi(n,im)-g22)/r(n,im)
      END IF
      
! ------- the minor component for the soc-manipulated wf is meaningless
      
      IF ( solver(1:7) == 'ABM-SOC' ) THEN
        CALL cinit(2*2*nrmax,qr)
        CALL cinit(2*2*nrmax,qi)
      END IF
      
!      COSD= NL * C * F11 - PFAC * SK1 * NLB1 * G11
!      SIND= JL * C * F11 - PFAC * SK1 * JLB1 * G11
      
! -------------------------------------------------------------------
!       T-SS  CONSTRUCTED USING EXPRESSIONS FROM H.E. + B.L.G. (1988)
! -------------------------------------------------------------------
      
      nl = (hl-jl)/ci
      nlb1 = (hlb1-jlb1)/ci
      nlb2 = (hlb2-jlb2)/ci
      
      
      IF ( solver(1:7) == 'ABM-SOC' ) THEN
        gam(1,1) = jl*g11p - jlp*g11
        gam(2,1) = jl*g21p - jlp*g21
        gam(1,2) = jl*g12p - jlp*g12
        gam(2,2) = jl*g22p - jlp*g22
        
        sig(1,1) = nl*g11p - nlp*g11
        sig(2,1) = nl*g21p - nlp*g21
        sig(1,2) = nl*g12p - nlp*g12
        sig(2,2) = nl*g22p - nlp*g22
      ELSE
        gam(1,1) = jl*c*f11 - pfac*sk1*jlb1*g11
        gam(2,1) = jl*c*f21 - pfac*sk2*jlb2*g21
        gam(1,2) = jl*c*f12 - pfac*sk1*jlb1*g12
        gam(2,2) = jl*c*f22 - pfac*sk2*jlb2*g22
        
        sig(1,1) = nl*c*f11 - pfac*sk1*nlb1*g11
        sig(2,1) = nl*c*f21 - pfac*sk2*nlb2*g21
        sig(1,2) = nl*c*f12 - pfac*sk1*nlb1*g12
        sig(2,2) = nl*c*f22 - pfac*sk2*nlb2*g22
      END IF
      
      CALL zcopy(nsol*nsol,gam,1,gaminv,1)
      CALL zgetrf(nsol,nsol,gaminv,2,ipiv,info)
      CALL zgetri(nsol,gaminv,2,ipiv,maux,2*2,info)
      
      DO i2 = 1,nsol
        DO i1 = 1,nsol
          csum = 0.0D0
          DO i3 = 1,nsol
            csum = csum + sig(i1,i3)*gaminv(i3,i2)
          END DO
          xsst2(i1,i2) = p*csum
        END DO
      END DO
      
      DO i1 = 1,nsol
        DO i2 = 1,nsol
          msst2(i1,i2) = -xsst2(i1,i2)
        END DO
        msst2(i1,i1) = msst2(i1,i1) + ci*p
      END DO
      
      CALL zcopy(nsol*nsol,msst2,1,tsst2,1)
      
      CALL zgetrf(nsol,nsol,tsst2,2,ipiv,info)
      CALL zgetri(nsol,tsst2,2,ipiv,maux,2*2,info)
      
      IF ( iprint >= 3 ) WRITE (1337,99001) it,l,mj,  &
          ((tsst2(i1,i2),i2=1,nsol),i1=1,nsol)
!------------------------------------------------------------------------
      
!   COEFFICIENTS TO CALCULATE THE SPIN  MAGNETISATION
      
      cgg(1,1) = cg1
      cgg(1,2) = cg2
      cgg(2,1) = cg2
      cgg(2,2) = cg4
      CALL rinit(4,cgf)
      cgf(1,1) = cg5
      cgf(2,2) = cg8
      
      
!   COEFFICIENTS TO CALCULATE THE SPIN  DIPOLAR MOMENT TZ
      
      tdia1 = 2*mj/DBLE((2*l1+1)*(2*lb1+1))
      tdia2 = 2*mj/DBLE((2*l1+1)*(2*lb2+1))
      toff = -SQRT((l1+0.5D0)**2-mj**2)/DBLE(2*l1+1)
      
      ctg(1,1) = 0.5D0*(cg1-3.0D0*tdia1)
      ctg(1,2) = 0.5D0*(cg2-3.0D0*toff)
      ctg(2,1) = 0.5D0*(cg2-3.0D0*toff)
      ctg(2,2) = 0.5D0*(cg4-3.0D0*tdia2)
      CALL rinit(4,ctf)
      
      
!   COEFFICIENTS TO CALCULATE THE ORBITAL MAGNETISATION
      
      cfg(1,1) = mj*(kap1+1.0D0)/(kap1+0.5D0)
      cfg(2,2) = mj*(kap2+1.0D0)/(kap2+0.5D0)
      cfg(1,2) = 0.5D0*DSQRT(1.0D0-(mj/(kap1+0.5D0))**2)
      cfg(2,1) = cfg(1,2)
      CALL rinit(4,cff)
      cff(1,1) = mj*(-kap1+1.0D0)/(-kap1+0.5D0)
      cff(2,2) = mj*(-kap2+1.0D0)/(-kap2+0.5D0)
      
!-----------------------------------------------------------------------
!   COEFFICIENTS TO CALCULATE THE ORBITAL POLARISATION
      
      CALL rinit(4*2,cog)
      CALL rinit(4*2,cof)
      DO is = 1,2
        cog(1,1,is) = cgc(ikm1,is)*cgc(ikm1,is) * DBLE(muem05-is+2)
        cof(1,1,is) = cgc(imkm1,is)*cgc(imkm1,is) * DBLE(muem05-is+2)
      END DO
      
      IF ( nsol == 2 ) THEN
        DO is = 1,2
          cog(2,2,is) = cgc(ikm2,is)*cgc(ikm2,is) * DBLE(muem05-is+2)
          cof(2,2,is) = cgc(imkm2,is)*cgc(imkm2,is) * DBLE(muem05-is+2)
          
          cog(1,2,is) = cgc(ikm1,is)*cgc(ikm2,is) * DBLE(muem05-is+2)
          cog(2,1,is) = cog(1,2,is)
        END DO
      END IF
!-----------------------------------------------------------------------
!   ANGULAR HYPERFINE MATRIX ELEMENTS   SEE E.G.  E.M.ROSE
!        THE FACTOR  I  HAS BEEN OMITTED
      
      ch(1,1) = 4.0D0*kap1*mj/(4.0D0*kap1*kap1-1.0D0)
      ch(2,2) = 4.0D0*kap2*mj/(4.0D0*kap2*kap2-1.0D0)
      IF ( nsol == 2 ) THEN
        ch(1,2) = DSQRT(0.25D0-(mj/DBLE(kap1-kap2))**2)
        ch(2,1) = ch(1,2)
      END IF
!-----------------------------------------------------------------------
!ALCULATE RADIAL INTEGRALS  UP TO   OR RWS(JTOP=JWS)
      
      
      IF ( nsol == 1 ) THEN
!====================================================================
! NO COUPLING TO OTHER SCATTERING CHANNELS
! REGULAR PART    Z*Z    Z = (GRA,FRA)
        
        norm = (p*nl-jl*xsst2(1,1))/g11
        
        DO i = 1,jtop
          zg(i,1,1) = (pr(1,1,i)/r(i,im))*norm
          zf(i,1,1) = (qr(1,1,i)/r(i,im)/c)*norm
          jg(i,1,1) = pi(1,1,i)/r(i,im)
          jf(i,1,1) = qi(1,1,i)/r(i,im)/c
        END DO
        
        
        IF ( iwrregwf /= 0 ) WRITE (nfilcbwf,REC=ikm1+(it-1)*nkm) it,l,mj,  &
            nsol,'REG',kap1,ikm1, (zg(i,1,1),zf(i,1,1),i=1,jtop)
        
        IF ( iwrirrwf /= 0 )  &
            WRITE (nfilcbwf,REC=ikm1+(it-1+nt)*nkm) it,l,mj,  &
            nsol,'IRR',kap1,ikm1, (jg(i,1,1),jf(i,1,1),i=1,jtop)
        
!============================================== NO COUPLING = END ===
      ELSE
!====================================================================
! COUPLING OF TWO SCATTERING CHANNELS
!   Z(K1,K2):  INDEX 1: SPIN-ANGULAR CHARACTER
!              INDEX 2: BOUNDARY CONDITION
        
        
        det = g11*g22 - g12*g21
        
!OEFFICIENTS TO GET:   Z(K1,K1)  Z(K2,K1)
        b1 = p*nl - xsst2(1,1)*jl
        b2 = -xsst2(2,1)*jl
        a(1,1) = (g22*b1-g12*b2)/det
        a(2,1) = (g11*b2-g21*b1)/det
        
!OEFFICIENTS TO GET:   Z(K1,K2)  Z(K2,K2)
        b1 = -xsst2(1,2)*jl
        b2 = p*nl - xsst2(2,2)*jl
        a(1,2) = (g22*b1-g12*b2)/det
        a(2,2) = (g11*b2-g21*b1)/det
        
!ALCULATE FUNCTIONS: Z(K1,K1), Z(K2,K1), Z(K1,K2), Z(K2,K2)
        DO k = 1,nsol
          DO i = 1,jtop
            zg(i,1,k) = (pr(1,1,i)*a(1,k)+pr(1,2,i)*a(2,k)) /r(i,im)
            zf(i,1,k) = (qr(1,1,i)*a(1,k)+qr(1,2,i)*a(2,k)) /r(i,im)/c
            
            zg(i,2,k) = (pr(2,1,i)*a(1,k)+pr(2,2,i)*a(2,k)) /r(i,im)
            zf(i,2,k) = (qr(2,1,i)*a(1,k)+qr(2,2,i)*a(2,k)) /r(i,im)/c
          END DO
        END DO
        DO k = 1,nsol
          kc = 3 - k
          DO i = 1,jtop
            jg(i,k,k) = pi(k,k,i)/r(i,im)
            jf(i,k,k) = qi(k,k,i)/r(i,im)/c
            jg(i,kc,k) = pi(kc,k,i)/r(i,im)
            jf(i,kc,k) = qi(kc,k,i)/r(i,im)/c
          END DO
        END DO
        
!-----------------------------------------------------------------------
        IF ( iwrregwf /= 0 ) THEN
! solution 1
          WRITE (nfilcbwf,REC=ikm1+(it-1)*nkm) it,l,mj,nsol,  &
              'REG',kap1,ikm1, (zg(i,1,1),zf(i,1,1),i=1,jtop),kap2,ikm2,  &
              (zg(i,2,1),zf(i,2,1),i=1,jtop)
          
! solution 2
          WRITE (nfilcbwf,REC=ikm2+(it-1)*nkm) it,l,mj,nsol,  &
              'REG',kap2,ikm2, (zg(i,2,2),zf(i,2,2),i=1,jtop),kap1,ikm1,  &
              (zg(i,1,2),zf(i,1,2),i=1,jtop)
        END IF
        
        IF ( iwrirrwf /= 0 ) THEN
! solution 1
          WRITE (nfilcbwf,REC=ikm1+(it-1+nt)*nkm) it,l,mj,  &
              nsol,'IRR',kap1,ikm1, (jg(i,1,1),jf(i,1,1),i=1,jtop),kap2,ikm2,  &
              (jg(i,2,1),jf(i,2,1),i=1,jtop)
          
! solution 2
          WRITE (nfilcbwf,REC=ikm2+(it-1+nt)*nkm) it,l,mj,  &
              nsol,'IRR',kap2,ikm2, (jg(i,2,2),jf(i,2,2),i=1,jtop),kap1,ikm1,  &
              (jg(i,1,2),jf(i,1,2),i=1,jtop)
          
        END IF
!================================================= COUPLING = END ===
      END IF
      
      
      
!ALCULATE SUM OF INTEGRALS TO BE MULTIPLIED TO   TAU(K1,K2)
      DO k1 = 1,nsol
        DO k2 = 1,nsol
          
          lin = lin + 1
          tsstlin(lin,it) = tsst2(k1,k2)
          IF ( calcint ) THEN
! REGULAR PART    Z*Z
            
            CALL cintabr(zg(1,1,k1),zg(1,1,k2),zgzg,  &
                zf(1,1,k1),zf(1,1,k2),zfzf, r2drdi(1,im),nsol,nsol,jtop,nrmax)
            
            CALL sumupint(dzz(lin,it),f1,zgzg,r1m,f1,zfzf, r1m,nsol)
            CALL sumupint(szz(lin,it),f1,zgzg,cgg,-f1,zfzf, cgf,nsol)
            CALL sumupint(ozz(lin,it),f1,zgzg,cfg,-f1,zfzf, cff,nsol)
! ----------------------------------------------------------------------
            DO is = 1,2
              CALL sumupint(ozzs(lin,it,is),f1,zgzg, cog(1,1,is),-f1,zfzf,  &
                  cof(1,1,is),nsol)
            END DO
! ----------------------------------------------------------------------
            CALL sumupint(qzz(lin,it),f1,zgzg,rkd,f1,zfzf, rkd,nsol)
            CALL sumupint(tzz(lin,it),f1,zgzg,ctg,-f1,zfzf, ctf,nsol)
            
!       write(66,'(3I3,2e16.7)') it,nsol,lin,DZZ(LIN,IT)
!       write(66,'(4e16.7)') ((ZGZG(ii,jj),ii=1,nsol),jj=1,nsol)
!       write(66,'(4e16.7)') ((ZFZF(ii,jj),ii=1,nsol),jj=1,nsol)
            
            
!-----------------------------------------------------------------------
            IF ( ihyper == 1 ) THEN
              CALL cinthff(zg(1,1,k1),zf(1,1,k1),zg(1,1,k2)  &
                  ,zf(1,1,k2),rmehf,nsol,nsol, jtop,cint,r(1,im),drdi(1,im),  &
                  nrmax)
              
              IF ( nucleus /= 0 ) THEN
! calculates integrals inside nucleus but up to now only
! approximately because jlim is not the nuclear radius
! the same arguments are valid for the irregular parts below
                CALL cinthff(zg(1,1,k1),zf(1,1,k1),  &
                    zg(1,1,k2),zf(1,1,k2),rmehf1,nsol,nsol,  &
                    jlim,cint,r(1,im),drdi(1,im),nrmax)
                CALL cinthff(zg(1,1,k1),zf(1,1,k1),  &
                    zg(1,1,k2),zf(1,1,k2),rmehf2,nsol,nsol,  &
                    jlim,cint,r(1,im),drovrn,nrmax)
                DO i = 1,nsol
                  DO j = 1,nsol
                    rmehf(i,j) = rmehf(i,j) - rmehf1(i,j) + rmehf2(i,j)
                  END DO
                END DO
              END IF
!                !end of nucleus.eq.0
              CALL sumupint(bzz(lin,it),cautog,rmehf,ch, 0.0D0,rmehf,ch,nsol)
            END IF
!-----------------------------------------------------------------------
            
! IRREGULAR PART    Z*J
! THE  ENDING  A (B)  STANDS FOR THE DOMINATING (DIMINATED)
! SET OF SPIN-ANGULAR-CHAR:  I.E.  J==J(A,A)  FOR R>RMT
            
            IF ( k1 == k2 ) THEN
              
              CALL cintabr(zg(1,1,k1),jg(1,1,k1),zgjg,  &
                  zf(1,1,k1),jf(1,1,k1),zfjf, r2drdi(1,im),nsol,nsol,jtop,  &
                  nrmax)
              
              CALL sumupint(dzj(lin,it),f1,zgjg,r1m,f1, zfjf,r1m,nsol)
              CALL sumupint(szj(lin,it),f1,zgjg,cgg,-f1, zfjf,cgf,nsol)
              CALL sumupint(ozj(lin,it),f1,zgjg,cfg,-f1, zfjf,cff,nsol)
! ----------------------------------------------------------------------
              DO is = 1,2
                CALL sumupint(ozjs(lin,it,is),f1,zgjg,  &
                    cog(1,1,is),-f1,zfjf, cof(1,1,is),nsol)
              END DO
! ----------------------------------------------------------------------
              CALL sumupint(qzj(lin,it),f1,zgjg,rkd,f1, zfjf,rkd,nsol)
              CALL sumupint(tzj(lin,it),f1,zgjg,ctg,-f1, zfjf,ctf,nsol)
              
!-----------------------------------------------------------------------
              IF ( ihyper == 1 ) THEN
                CALL cinthff(zg(1,1,k1),zf(1,1,k1),  &
                    jg(1,1,k1),jf(1,1,k1),rmehf,nsol,nsol,  &
                    jtop,cint,r(1,im),drdi(1,im),nrmax)
                IF ( nucleus /= 0 ) THEN
! calculates integrals inside nucleus but up to now only
! approximately because jlim is not the nuclear radius
! the same arguments are valid for the irregular parts below
                  CALL cinthff(zg(1,1,k1),zf(1,1,k1),  &
                      jg(1,1,k1),jf(1,1,k1),rmehf1,nsol,  &
                      nsol,jlim,cint,r(1,im),drdi(1,im), nrmax)
                  CALL cinthff(zg(1,1,k1),zf(1,1,k1),  &
                      jg(1,1,k1),jf(1,1,k1),rmehf2,nsol,  &
                      nsol,jlim,cint,r(1,im),drovrn,nrmax)
                  DO i = 1,nsol
                    DO j = 1,nsol
                      rmehf(i,j) = rmehf(i,j) - rmehf1(i,j) + rmehf2(i,j)
                    END DO
                  END DO
                END IF
!                !end of nucleus.eq.0
                
                CALL sumupint(bzj(lin,it),cautog,rmehf,ch,  &
                    0.0D0,rmehf,ch,nsol)
              END IF
            END IF
!-----------------------------------------------------------------------
          END IF
!           ! OF IF (.CALCINT.)
        END DO
      END DO
      
!heck WRONSKI-relationship
      
      IF ( wronski ) THEN
        DO i = 1,jtop,40
          crsq = c*r(i,im)**2
          WRITE (1337,99002) it,l,nint(2*mj),i,r(i,im),  &
              1.0D0 - (zf(i,1,1)*jg(i,1,1) -zg(i,1,1)*jf(i,1,1)+zf(i,2,1)  &
              *jg(i,2,1)-zg(i,2,1)*jf(i,2,1)) *crsq
          IF ( nsol == 2 ) THEN
            WRITE (1337,99003) 1.0D0 - (zf(i,1,2)*jg(i,1,2)-zg(  &
                i,1,2)*jf(i,1,2)+zf(i,2,2) *jg(i,2,2)-zg(i,2,2)*jf(i,2,2))  &
                *crsq, 1.0D0 - (zf(i,1,2)*jg(i,1,1)  &
                -zg(i,1,2)*jf(i,1,1)+zf(i,2,2) *jg(i,2,1)-zg(i,2,2)*jf(i,2,1))  &
                *crsq, 1.0D0 - (zf(i,1,1)*jg(i,1,2)  &
                -zg(i,1,1)*jf(i,1,2)+zf(i,2,1) *jg(i,2,2)-zg(i,2,1)*jf(i,2,2))  &
                *crsq
          ELSE
            WRITE (1337,*)
          END IF
        END DO
      END IF
      
    END DO
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
  END DO
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
  
  CALL cinit(nkmmax*nkmmax,tsst(1,1,it))
  
  DO lin = 1,nlinq(iq)
    i1 = ikm1lin(lin)
    i2 = ikm2lin(lin)
    tsst(i1,i2,it) = tsstlin(lin,it)
  END DO
  
  DO j = 1,nkmq(iq)
    CALL zcopy(nkmq(iq),tsst(1,j,it),1,msst(1,j,it),1)
  END DO
  
  CALL zgetrf(nkmq(iq),nkmq(iq),msst(1,1,it),nkmmax,ipiv,info)
  CALL zgetri(nkmq(iq),msst(1,1,it),nkmmax,ipiv,maux, nkmmax*nkmmax,info)
  
  IF ( iprint >= 4 ) THEN
    DO lin = 1,nlinq(iq)
      i1 = ikm1lin(lin)
      i2 = ikm2lin(lin)
      WRITE (1337,99004) it,lin
      WRITE (1337,99004) it,i1,i2,' DZZ ', dzz(lin,it),dzj(lin,it)
      WRITE (1337,99004) it,i1,i2,' SZZ ', szz(lin,it),szj(lin,it)
      WRITE (1337,99004) it,i1,i2,' OZZ ', ozz(lin,it),ozj(lin,it)
    END DO
  END IF
END DO
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

99001 FORMAT ('  t-ss(TLM)',2I3,f5.1,2E14.5,2X,2E14.5,:,/,22X,2E14.5,2X,  &
    2E14.5)
99002 FORMAT (' IT=',i2,' L=',i2,' MJ=',i2,'/2',i7,f10.6,1(2X,2F9.6),$)
99003 FORMAT (3(2X,2F9.6))
99004 FORMAT (' IT=',i2,2I3,a,2X,2E14.5,2X,2E14.5)
END SUBROUTINE ssite
!*==readwfun.f    processed by SPAG 6.05Rc at 17:31 on 29 Apr 2001

SUBROUTINE readwfun(nfil,it,l,mj,nsol,sreg,sirr,ikm1,kap1,ikm2,  &
    kap2,nt,nkm,zg,zf,jg,jf,jtop,nrmax)
!   ********************************************************************
!   *                                                                  *
!   *  reread the wave functions written by  <SSITE>  or  <CORE>       *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

! Dummy arguments
INTEGER IKM1,IKM2,IT,JTOP,KAP1,KAP2,L,NFIL,NKM,NRMAX,NSOL,NT
REAL*8 MJ
CHARACTER*3 SIRR,SREG
COMPLEX*16 JF(NRMAX,2,2),JG(NRMAX,2,2),ZF(NRMAX,2,2),ZG(NRMAX,2,2)

! Local variables
INTEGER I,IFLAG,IKMIN(2),ITP,K,KAPIN(2),LP,NSOLP
REAL*8 MJP
CHARACTER*3 STRP

iflag = 0
!-----------------------------------------------------------------------
!                                    REGULAR wave function -- solution 1
IF ( sreg == 'REG' .OR. sreg == 'COR' ) THEN
  READ (nfil,REC=ikm1+(it-1)*nkm) itp,lp,mjp,nsolp,strp,  &
      (kapin(k),ikmin(k),(zg(i,k,1),zf(i,k,1),i=1,jtop),k=1, nsol)
  IF ( itp /= it .OR. lp /= l .OR. ABS(mjp-mj) > 0.001D0 .OR.  &
      nsolp /= nsol .OR. strp /= 'REG' ) iflag = iflag + 1
  IF ( kap1 /= kapin(1) ) iflag = iflag + 1
  IF ( ikm1 /= ikmin(1) ) iflag = iflag + 1
  IF ( nsol > 1 ) THEN
    IF ( kap2 /= kapin(2) ) iflag = iflag + 1
    IF ( ikm2 /= ikmin(2) ) iflag = iflag + 1
  END IF
END IF

!-----------------------------------------------------------------------
!                                  IRREGULAR wave function -- solution 1
IF ( sirr == 'IRR' ) THEN
  READ (nfil,REC=ikm1+(it-1+nt)*nkm) itp,lp,mjp,nsolp,strp,  &
      (kapin(k),ikmin(k),(jg(i,k,1),jf(i,k,1),i=1,jtop),k=1, nsol)
  IF ( itp /= it .OR. lp /= l .OR. ABS(mjp-mj) > 0.001D0 .OR.  &
      nsolp /= nsol .OR. strp /= 'IRR' ) iflag = iflag + 1
  IF ( kap1 /= kapin(1) ) iflag = iflag + 1
  IF ( ikm1 /= ikmin(1) ) iflag = iflag + 1
  IF ( nsol > 1 ) THEN
    IF ( kap2 /= kapin(2) ) iflag = iflag + 1
    IF ( ikm2 /= ikmin(2) ) iflag = iflag + 1
  END IF
END IF

IF ( nsol == 2 ) THEN
!-----------------------------------------------------------------------
!                                    REGULAR wave function -- solution 2
  IF ( sreg == 'REG' .OR. sreg == 'COR' ) THEN
    READ (nfil,REC=ikm2+(it-1)*nkm) itp,lp,mjp,nsolp,strp,  &
        (kapin(k),ikmin(k),(zg(i,k,2),zf(i,k,2),i=1,jtop),k=2, 1,-1)
    IF ( itp /= it .OR. lp /= l .OR. ABS(mjp-mj) > 0.001D0 .OR.  &
        nsolp /= nsol .OR. strp /= 'REG' ) iflag = iflag + 1
  END IF
  
!-----------------------------------------------------------------------
!                                  IRREGULAR wave function -- solution 2
  IF ( sirr == 'IRR' ) THEN
    READ (nfil,REC=ikm2+(it-1+nt)*nkm) itp,lp,mjp,nsolp,strp,  &
        (kapin(k),ikmin(k),(jg(i,k,2),jf(i,k,2),i=1,jtop),k=2, 1,-1)
    IF ( itp /= it .OR. lp /= l .OR. ABS(mjp-mj) > 0.001D0 .OR.  &
        nsolp /= nsol .OR. strp /= 'IRR' ) iflag = iflag + 1
  END IF
  
END IF

!-----------------------------------------------------------------------

IF ( iflag > 0 ) THEN
  WRITE (*,99001) iflag,it,l,mj,sreg,sirr
  STOP ' in <READWFUN>'
END IF
99001 FORMAT (//,1X,79('*'),/,10X,'error reading the wave functions',  &
    ' IFLAG = 1',/,10X,'for  IT=',i2,' L=',i2,' MJ=',f4.1, ' KEYS=',a,a)
END SUBROUTINE readwfun
