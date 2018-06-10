SUBROUTINE hffcore(rnuc,jtop,kap1,kap2,nsol,mj,gc,fc,nrc,shf,s,  &
        nmemax,nkmmax,r,drdi,sdia,smdia,soff,smoff,  &
        qdia,qoff,qmdia,qmoff,nucleus,jlim)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Calculates matrix elements of several hyperfine interaction
!     connected quantities in the core.
!     All the related arrays have a counting index as
!     the last index of the array indicates the corresponding physical
!     property.
!     Index-list
!     1      electron-Spin-electron-Spin Hyperfine field
!     2      nuclear-spin-electron-orbit hyperfine field
!     3      electron-spin-nulceus-spin-contact hyperfine field
!     4      expectation value of (1/r)^3
!     5      Total Hyperfine Field (see Rose (1961))
!     called by core
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IMPLICIT NONE


! PARAMETER definitions
DOUBLE PRECISION MB,A0,F1,F2

!BOHR-MAGNETON       IN ERG/GAUSS
!CONVERSION FACTOR FOR HYPERFINE FIELDS FROM A.U. TO GAUSS
!ELECTRON CHARGE     IN ESU

PARAMETER (MB=9.274078D-21,A0=0.52917706D-08,F1=1.0D0, &
           F2=2.0D0*MB/(A0*A0*A0))

! Dummy arguments
INTEGER JLIM,JTOP,KAP1,KAP2,NKMMAX,NMEMAX,NRC,NSOL,NUCLEUS,S
DOUBLE PRECISION MJ,RNUC
DOUBLE PRECISION DRDI(NRC),FC(2,2,NRC),GC(2,2,NRC),QDIA(NKMMAX), &
       QMDIA(NKMMAX),QMOFF(NKMMAX),QOFF(NKMMAX),R(NRC), &
       SDIA(NKMMAX),SHF(2,2,NMEMAX),SMDIA(NKMMAX),SMOFF(NKMMAX), &
       SOFF(NKMMAX)

! Local variables
DOUBLE PRECISION AME(2,2),CFF(2,2),CFG(2,2),CGF(2,2),CGG(2,2), &
       CQF(2,2), &
       CQG(2,2),CSF(2,2),CSG(2,2),DOVR(NRC),DROVRN(NRC), &
       DROVRN1(NRC),F(NRC,2),FF(2,2),FF1(2,2),FF2(2,2),FG(2,2), &
       FG1(2,2),FG2(2,2),G(NRC,2),GF(2,2),GF1(2,2),GF2(2,2), &
       GG(2,2),GG1(2,2),GG2(2,2)
DOUBLE PRECISION DBLE,DSQRT
INTEGER I,IKM1,IKM2,J,K,K1,K2,N
INTEGER IKAPMUE
INTEGER NINT

IF ( kap2 == 0 ) kap2 = kap1

CALL rinit(4,gg)
CALL rinit(4,ff)
CALL rinit(4,gg1)
CALL rinit(4,ff1)
CALL rinit(4,gg2)
CALL rinit(4,ff2)
CALL rinit(4,gf)
CALL rinit(4,fg)
CALL rinit(4,gf1)
CALL rinit(4,fg1)
CALL rinit(4,gf2)
CALL rinit(4,fg2)
CALL rinit(2*nrc,g)
CALL rinit(2*nrc,f)

DO k = 1,2
  DO n = 1,jtop
    g(n,k) = gc(k,s,n)
    f(n,k) = fc(k,s,n)
  END DO
END DO
!     prepare meshes for finite nucleus calculation
DO i = 1,jtop
  dovr(i) = drdi(i)/r(i)
  IF ( nucleus /= 0 ) THEN
    drovrn1(i) = (r(i)/rnuc)**3*drdi(i)
    drovrn(i) = drovrn1(i)/r(i)
  END IF
END DO
ikm1 = ikapmue(kap1,nint(mj-0.5D0))
ikm2 = ikapmue(kap2,nint(mj-0.5D0))
!     angular hyperfine matrix elements   see e.g.  E.M.Rose
!     the factor  i  has been omitted
CALL rinit(4,ame)
ame(1,1) = 4.0D0*kap1*mj/(4.0D0*kap1*kap1-1.0D0)
IF ( nsol == 2 ) THEN
  ame(2,2) = 4.0D0*kap2*mj/(4.0D0*kap2*kap2-1.0D0)
  ame(1,2) = DSQRT(0.25D0-(mj/DBLE(kap1-kap2))**2)
  ame(2,1) = ame(1,2)
END IF
!     coefficients for the spin-dipolar matrix elements
CALL rinit(4,csf)
CALL rinit(4,csg)
csg(1,1) = sdia(ikm1)
csf(1,1) = smdia(ikm1)
IF ( nsol == 2 ) THEN
  csg(2,2) = sdia(ikm2)
  csg(1,2) = soff(ikm1)
  csg(2,1) = csg(1,2)
  csf(2,2) = smdia(ikm2)
  csf(1,2) = smoff(ikm1)
  csf(2,1) = smoff(ikm1)
END IF
!     COEFFICIENTS FOR THE QUADRUPOLAR MATRIX ELEMENTS
cqg(1,1) = qdia(ikm1)
cqg(2,2) = qdia(ikm2)
cqg(1,2) = qoff(ikm1)
cqg(2,1) = cqg(1,2)
CALL rinit(4,cqf)
cqf(1,1) = qmdia(ikm1)
cqf(2,2) = qmdia(ikm2)
cqf(1,2) = qmoff(ikm1)
cqf(2,1) = cqf(1,2)
!     coefficients to calculate the spin-spin field
CALL rinit(4,cgg)
CALL rinit(4,cgf)
cgg(1,1) = -mj/(kap1+0.5D0)
cgf(1,1) = -mj/(-kap1+0.5D0)
IF ( nsol == 2 ) THEN
  cgg(1,2) = -DSQRT(1.0D0-(mj/(kap1+0.5D0))**2)
  cgg(2,1) = cgg(1,2)
  cgg(2,2) = -mj/(kap2+0.5D0)
  cgf(2,2) = -mj/(-kap2+0.5D0)
!     CGF(1,2) = -DSQRT( 1.0D0 - (MJ/(- KAP1+0.5D0))**2 )
!     CGF(2,1) = CGF(1,2)
END IF
!     coefficients to calculate the orbital field
CALL rinit(4,cfg)
CALL rinit(4,cff)
cfg(1,1) = mj*(kap1+1.0D0)/(kap1+0.5D0)
cff(1,1) = mj*(-kap1+1.0D0)/(-kap1+0.5D0)
IF ( nsol == 2 ) THEN
  cfg(2,2) = mj*(kap2+1.0D0)/(kap2+0.5D0)
  cfg(1,2) = 0.5D0*DSQRT(1.0D0-(mj/(kap1+0.5D0))**2)
  cfg(2,1) = cfg(1,2)
  cff(2,2) = mj*(-kap2+1.0D0)/(-kap2+0.5D0)
!     CFF(1,2) = 0.5D0 * DSQRT( 1.0D0 - (MJ/(- KAP1 + 0.5D0))**2 )
!     CFF(2,1) = CFF(1,2)
END IF
! Calculates integrals from 0.0 to jtop
CALL hffint(gg,g,g,dovr,r,0.0D0,nsol,jtop,nrc)
CALL hffint(ff,f,f,dovr,r,0.0D0,nsol,jtop,nrc)
CALL hffint(gf,g,f,drdi,r,0.0D0,nsol,jtop,nrc)
CALL hffint(fg,f,g,drdi,r,0.0D0,nsol,jtop,nrc)
CALL rsumupint(shf(1,1,5),f1,gg,cqg,f1,ff,cqf,nsol)
IF ( nucleus /= 0 ) THEN
!     calculates integrals inside nucleus at RNUC in order to get
!     contribution outside the nucleus
  CALL hffint(gg1,g,g,dovr,r,rnuc,nsol,jlim,nrc)
  CALL hffint(ff1,f,f,dovr,r,rnuc,nsol,jlim,nrc)
  CALL hffint(gf1,g,f,drdi,r,rnuc,nsol,jlim,nrc)
  CALL hffint(fg1,f,g,drdi,r,rnuc,nsol,jlim,nrc)
!     calculates contribution from RNUC to jtop
  DO i = 1,nsol
    DO j = 1,nsol
      gg(i,j) = gg(i,j) - gg1(i,j)
      ff(i,j) = ff(i,j) - ff1(i,j)
      gf(i,j) = gf(i,j) - gf1(i,j)
      fg(i,j) = fg(i,j) - fg1(i,j)
    END DO
  END DO
END IF                    !end of nucleus.eq.0
!     calculates B_sp which is zero inside the nucleus
CALL rsumupint(shf(1,1,1),f1,gg,csg,-f1,ff,csf,nsol)
!     calculates hyperfine integrals from 0.0 to RNUC which are added
!     external integrals
IF ( nucleus /= 0 ) THEN
  CALL hffint(gg2,g,g,drovrn,r,rnuc,nsol,jlim,nrc)
  CALL hffint(ff2,f,f,drovrn,r,rnuc,nsol,jlim,nrc)
  CALL hffint(gf2,g,f,drovrn1,r,rnuc,nsol,jlim,nrc)
  CALL hffint(fg2,f,g,drovrn1,r,rnuc,nsol,jlim,nrc)
  DO i = 1,nsol
    DO j = 1,nsol
      gg(i,j) = gg(i,j) + gg2(i,j)
      ff(i,j) = ff(i,j) + ff2(i,j)
      gf(i,j) = gf(i,j) + gf2(i,j)
      fg(i,j) = fg(i,j) + fg2(i,j)
    END DO
  END DO
END IF
!     calculates B_nseo and B_tot
CALL rsumupint(shf(1,1,2),f2,gg,cfg,-f2,ff,cff,nsol)
!      CALL RSUMUPINT(SHF(1,1,5),CAUTOG,GF,AME,CAUTOG,FG,AME,NSOL)
!     modifications for B_ssc which is zero outside the nucleus
IF ( nucleus /= 0 ) THEN
  DO i = 1,nsol
    DO j = 1,nsol
      gg(i,j) = gg2(i,j)
      ff(i,j) = ff2(i,j)
    END DO
  END DO
END IF
!     calculates B_ssc
CALL rsumupint(shf(1,1,3),f2,gg,cgg,-f2,ff,cgf,nsol)
!     for testing purposes write in 4 the sum of 1,2,3
DO k1 = 1,nsol
  DO k2 = 1,nsol
    shf(k1,k2,4) = shf(k1,k2,1) + shf(k1,k2,2) + shf(k1,k2,3)
  END DO
END DO

END SUBROUTINE hffcore

SUBROUTINE rsumupint(sum,vg,g,wg,vf,f,wf,n)
IMPLICIT NONE

! Dummy arguments
INTEGER N
DOUBLE PRECISION VF,VG
DOUBLE PRECISION F(N,N),G(N,N),SUM(N,N),WF(N,N),WG(N,N)

! Local variables
INTEGER I,J

DO j = 1,n
  DO i = 1,n
    sum(i,j) = vg*g(i,j)*wg(i,j) + vf*f(i,j)*wf(i,j)
  END DO
END DO
END SUBROUTINE rsumupint

SUBROUTINE hffint(gg,ga,gb,dr,r,rnuc,nsol,jtop,nrc)
!     Calculates Hyperfine integrals, extrapolates to zero and
!     intrapolates to exact nuclear radius RNUC
IMPLICIT NONE

! Dummy arguments
INTEGER JTOP,NRC,NSOL
DOUBLE PRECISION RNUC
DOUBLE PRECISION DR(NRC),GA(NRC,2),GB(NRC,2),GG(2,2),R(NRC)

! Local variables
INTEGER I,K1,K2
DOUBLE PRECISION X(5),Y(5),YI(NRC),ZI(NRC)
DOUBLE PRECISION YLAG

DO k1 = 1,nsol
  DO k2 = 1,nsol
    DO i = 1,jtop
      yi(i) = ga(i,k1)*gb(i,k2)*dr(i)
    END DO
    CALL rint4pts(yi,jtop,zi)
!     Intrapolation
    IF ( rnuc /= 0.0D0 ) THEN
      DO i = 1,5
        x(i) = r(jtop-5+i)
        y(i) = zi(jtop-5+i)
      END DO
      zi(jtop) = ylag(rnuc,x,y,0,4,5)
    END IF
!     Extrapolation
    x(1) = 1.0D0
    x(2) = 6.0D0
    x(3) = 11.0D0
    y(1) = zi(jtop) - zi(1)
    y(2) = zi(jtop) - zi(5)
    y(3) = zi(jtop) - zi(9)
    gg(k1,k2) = ylag(0.0D0,x,y,0,2,3)
  END DO
END DO
END SUBROUTINE hffint

