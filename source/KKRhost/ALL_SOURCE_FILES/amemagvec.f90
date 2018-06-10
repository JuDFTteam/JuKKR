SUBROUTINE amemagvec(irel,iprint,nkm,amemvec,ikmllim1,ikmllim2,  &
    imkmtab,cgc,nlmax,nkmmax,nkmpmax,nmvecmax)
!   ********************************************************************
!   *                                                                  *
!   *   calculate the angular matrix elements connected with           *
!   *                                                                  *
!   *   1:       < LAM | sigma(ipol) | LAM' >    spin moment           *
!   *   2:       < LAM |     l(ipol) | LAM' >    orbital moment        *
!   *   3:       < LAM |     T(ipol) | LAM' >    spin dipole moment    *
!   *   4:       < LAM |  B_hf(ipol) | LAM' >    hyperfine field       *
!   *                                                                  *
!   *   ipol= 1,2,3  ==  (+),(-),(z):                                  *
!   *                                                                  *
!   ********************************************************************
IMPLICIT REAL*8(A-H,O-Z)

! Dummy arguments

INTEGER IPRINT,IREL,NKM,NKMMAX,NKMPMAX,NLMAX,NMVECMAX
DOUBLE PRECISION AMEMVEC(NKMMAX,NKMMAX,3,NMVECMAX),CGC(NKMPMAX,2)
INTEGER IKMLLIM1(NKMMAX),IKMLLIM2(NKMMAX),IMKMTAB(NKMMAX)

! Local variables

CHARACTER (len=1) :: CHPOL(3)
DOUBLE PRECISION DBLE
INTEGER I,IKM1,IKM2,IMKM,IMV,IMVEC,IPOL,J1P05,J2P05,K,K1,K2,KAP1, &
        KAP2,L,L1,L2,LB1,LB2,M2,MSM05,MUE1M05,MUE2M05,NK,NMVEC,IXM
INTEGER IABS,NINT
CHARACTER (len=20) :: STR20
DOUBLE PRECISION SUM,XJ,XJM,XJP,XM,XYNORM
CHARACTER (len=4) :: TXTMVEC(4)

DATA CHPOL/'+','-','z'/
DATA TXTMVEC/'spin','orb ','T_z ','B_hf'/

nk = 2*nlmax - 1
!     XYNORM = DSQRT(2.0D0)
xynorm = 2.0D0

CALL rinit(nkmmax*nkmmax*3*nmvecmax,amemvec)

IF ( irel <= 1 ) RETURN

! ----------------------------------------------------------------------
!     find the bounding indices  IKMLLIM1  and  IKMLLIM2  for IKM-loops
!     assuming that the matrix elements are diagonal with respect to l
!     this does not hold for B_hf for which there are l-(l+/-2)-terms
! ----------------------------------------------------------------------

i = 0
DO k = 1,nk
  l = k/2
  xjm = l - 0.5D0
  xjp = l + 0.5D0
  IF ( MOD(k,2) == 1 ) THEN
    xj = l + 0.5D0
  ELSE
    xj = l - 0.5D0
  END IF
!DO XM = -XJ, + XJ
  DO ixm = 1, 2*nint(xj)+1
    xm = -xj + DBLE(ixm-1)
    i = i + 1
    ikmllim1(i) = nint(l*2*(xjm+0.5D0)+1)
    ikmllim2(i) = nint(l*2*(xjp+0.5D0)+2*xjp+1)
  END DO
END DO

! ----------------------------------------------------------------------

ikm1 = 0
DO k1 = 1,nk
  l1 = k1/2
  IF ( MOD(k1,2) == 0 ) THEN
    kap1 = l1
    lb1 = l1 - 1
  ELSE
    kap1 = -l1 - 1
    lb1 = l1 + 1
  END IF
  j1p05 = IABS(kap1)
  
  DO mue1m05 = -j1p05,j1p05 - 1
    ikm1 = ikm1 + 1
    imkm = lb1*2*j1p05 + j1p05 + mue1m05 + 1
    imkmtab(ikm1) = imkm
    
    ikm2 = 0
    DO k2 = 1,nk
      l2 = k2/2
      IF ( MOD(k2,2) == 0 ) THEN
        kap2 = l2
        lb2 = l2 - 1
      ELSE
        kap2 = -l2 - 1
        lb2 = l2 + 1
      END IF
      j2p05 = IABS(kap2)
      
      DO mue2m05 = -j2p05,j2p05 - 1
        ikm2 = ikm2 + 1
! ----------------------------------------------------------------------
        IF ( (mue1m05-mue2m05) == +1 ) THEN
          amemvec(ikm1,ikm2,1,1) = xynorm*cgc(ikm1,2) *cgc(ikm2,1)
          
          sum = 0D0
          DO msm05 = -1,0
            m2 = mue2m05 - msm05
            IF (ABS(m2) <= l2) sum = sum +  &
                cgc(ikm1,msm05+2) * cgc(ikm2,msm05+2)  &
                * SQRT(DBLE((l2-m2)*(l2+m2+1)))
          END DO
          amemvec(ikm1,ikm2,1,2) = sum
        END IF
        
        IF ( (mue1m05-mue2m05) == -1 ) THEN
          amemvec(ikm1,ikm2,2,1) = xynorm*cgc(ikm1,1) *cgc(ikm2,2)
          
          sum = 0D0
          DO msm05 = -1,0
            m2 = mue2m05 - msm05
            IF (ABS(m2) <= l2) sum = sum +  &
                cgc(ikm1,msm05+2) * cgc(ikm2,msm05+2)  &
                * SQRT(DBLE((l2+m2)*(l2-m2+1)))
            
          END DO
          amemvec(ikm1,ikm2,2,2) = sum
        END IF
        
        IF ( (mue1m05-mue2m05) == 0 ) THEN
          amemvec(ikm1,ikm2,3,1) = cgc(ikm1,2)*cgc(ikm2,2)  &
              - cgc(ikm1,1)*cgc(ikm2,1)
          
          sum = 0D0
          DO msm05 = -1,0
            m2 = mue2m05 - msm05
            sum = sum + cgc(ikm1,msm05+2)*cgc(ikm2,msm05+2) *m2
          END DO
          amemvec(ikm1,ikm2,3,2) = sum
        END IF
        
! ----------------------------------------------------------------------
      END DO
    END DO
  END DO
END DO

nmvec = 2
! ----------------------------------------------------------------------
IF ( iprint <= 90 ) RETURN

DO imvec = 1,nmvec
  imv = MIN(imvec,2)
  DO ipol = 1,3
    str20 = 'A  '//txtmvec(imv)//'  ('//chpol(ipol)//')'
    CALL rmatstr(str20,12,amemvec(1,1,ipol,imvec),nkm,nkmmax,3, 3,1D-8,6)
  END DO
END DO
END SUBROUTINE amemagvec
