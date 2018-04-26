      SUBROUTINE SPHERE_GGA(LMAX,YR,WTYR,RIJ,IJD,LMMAXD,THET,YLM,
     +                      DYLMT1,DYLMT2,DYLMF1,DYLMF2,DYLMTF)
      use Harmonics_mod, only: ymy
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     generate an angular mesh and spherical harmonics at those
c     mesh points. For an angular integration the weights are ge-
c     rated .
c
c     R. Zeller      Feb. 1996
c     New call to subr. ylmderiv for accurate derivatives of
c     spherical harmonics.
c     Phivos Mavropoulos, July 2007.
c-----------------------------------------------------------------------

C     .. Scalar Arguments ..
      INTEGER IJD,LMAX,LMMAXD
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION PI,R,R1,R2,R3
      INTEGER IJ,LM1
C     ..
C     .. External Subroutines ..
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DYLMF1(IJD,LMMAXD),DYLMF2(IJD,LMMAXD),
     +                 DYLMT1(IJD,LMMAXD),DYLMT2(IJD,LMMAXD),
     +                 DYLMTF(IJD,LMMAXD),RIJ(IJD,3),THET(IJD),
     +                 WTYR(IJD,*),YLM(IJD,LMMAXD),YR(IJD,*)
     &                ,DYDTH(LMMAXD),DYDFI(LMMAXD),D2YDTH2(LMMAXD)
     &                ,D2YDFI2(LMMAXD),D2YDTHDFI(LMMAXD)
C     ..
C     .. Local Arrays ..
C
      DOUBLE PRECISION WGHT,Y(1000)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ACOS,ATAN,COS,SIN,SQRT
C     ..
      PI = 4.D0*ATAN(1.D0)
C      WRITE (1337,*) 'SPHERE for GGA: read LEBEDEV mesh'
      IF (IJD.GT.1000) STOP 'SPHERE'
c
c
      DO 30 IJ = 1,IJD
         CALL LEBEDEV (IJ,R1,R2,R3,WGHT)
         RIJ(IJ,1) = R1
         RIJ(IJ,2) = R2
         RIJ(IJ,3) = R3


C
c For the needs of GGA PW91 as implemented here, ylm and derivatives
c come with a different sign convention compared to the usual in the
c program: sin(fi)**m --> -sin(fi)**m. Thus some signs change
c also in array ylm compared to array yr (below).
         CALL DERIVYLM(
     >        R1,R2,R3,LMAX,
     <        R,Y,DYDTH,DYDFI,D2YDTH2,D2YDFI2,D2YDTHDFI)

         THET(IJ) = DACOS(R3/R)


         DO 10 LM1 = 1, (LMAX+1)**2
            YLM(IJ,LM1) = Y(LM1)
            DYLMT1(IJ,LM1) = DYDTH(LM1)
            DYLMF1(IJ,LM1) = DYDFI(LM1)
            DYLMT2(IJ,LM1) = D2YDTH2(LM1)
            DYLMF2(IJ,LM1) = D2YDFI2(LM1)
            DYLMTF(IJ,LM1) = D2YDTHDFI(LM1)
 10      CONTINUE

c Call ymy to obtain sher. harmonics with usual convention
c
c---> multiply the spherical harmonics with the weights
c

         CALL YMY(R1,R2,R3,R,Y,LMAX)
         DO 20 LM1 = 1, (LMAX+1)**2
            YR(IJ,LM1) = Y(LM1)
            WTYR(IJ,LM1) = YR(IJ,LM1)*WGHT*PI*4.D0
 20      CONTINUE




   30 CONTINUE


      END

c-----------------------------------------------------------------------
      subroutine derivylm(
     > v1,v2,v3,lmax,
     < Rabs,YLM,dYdth,dYdfi,d2Ydth2,d2Ydfi2,d2Ydthdfi)
c Calculate the 1st and 2nd derivatives of real spherical harmonics
c with respect to theta, fi.
c
c Use recursion relations for the assoc. Legendre functions P[l,m] to generate
c the derivatives. These are (taken from Abramowitz and Stegun, 
c Handbook of Mathematical Functions, chapt. 8.):
c
c P[l,m+1] = (x**2-1)**(-1/2) ( (l-m)*x*P[l,m] - (l+m)*P[l-1,m] ) (8.5.1)
c
c (x**2-1)*dP[l,m]/dx = (l+m)*(l-m+1)*(x**2-1)**(1/2) P[l,m-1] - m*x*P[l,m]  (8.5.2)
c
c (x**2-1)*dP[l,m]/dx = l*x*P[l,m] - (l+m)*P[l-1,m]           (8.5.4)
c
c where x=cos(th), (x**2-1)**(1/2) = -sin(th), d/dth = -sin(th) d/dx.
c
c Adding (8.5.2)+(8.5.4) and using (8.5.1) results in:
c
c dP[l,m](cos(th)) / dth = (1/2) * ( -(l+m)*(l-m+1)*P[l,m-1] + P[l,m+1] )   (A)
c
c
c It is implied that P[l,m]=0 if m>l or m<-l. Also, the term (x**2-1)**(1/2)
c is ambiguous for real x, 0<x<1; here it is interpreted as 
c (x**2-1)**(1/2)=-sin(th), but (x**2-1)**(-1/2)=+1/sin(th) (in 8.5.1), 
c otherwise the result (A) (which is cross-checked and correct) does not follow.
c
c For the 2nd derivative apply (A) twice. Result:
c
c ddP[l,m](cos(th))/dth/dth = (1/4) * 
c                         (    (l+m)*(l-m+1)*(l+m-1)*(l-m+2) * P[l,m-2] 
c                           - ( (l-m)*(l+m+1)+(l+m)*(l-m+1) )* P[l,m]   
c                           +                                  P[l,m+2]  )   (B)
c
c The fi-derivatives act on cos(fi),sin(fi) and are trivial.
c
c For the associated Legendre functions use the recursion formulas:
c
c (l-m+1)*P[l+1,m] = (2l+1)*cos(th)*P[l,m] - (l+m)*P[l-1,m]   (8.5.3)
c
c P[l+1,m] = P[l-1,m] - (2*l+1)*sin(th)*P[l,m-1]              (8.5.5)
c
c ( with x=cos(th) ).
c
c Recursion algorithm for the calculation of P[l,m] and calculation of Ylm 
c taken over from subr. ymy of KKR program 
c (implemented there by M. Weinert, B. Drittler).
c
c For m<0, use P[l,-m] = P[l,m] (l-m)!/(l+m)!     (C)
c
c Taking into account the lm-prefactors of the spherical harmonics, 
c we construct and use the functions
c 
c Q[l,m] = sqrt((2*l+1)/(4*pi)) * sqrt((l-m)!/(l+m)!) * P[l,m]
c
c whence (A) and (B) become
c
c dQ[l,m]/dth = (1/2) *
c     ( -sqrt((l+m)*(l-m+1))*Q[l,m-1] + sqrt((l+m+1)*(l-m))*Q[l,m+1] )   (A1)
c
c ddQ[l,m]/dth/dth = (1/4) *
c     (       sqrt((l+m)*(l+m-1)*(l-m+1)*(l-m+2)) * Q[l,m-2]
c         +   ((l-m)*(l+m+1)+(l+m)*(l-m+1))       * Q[l,m]
c         +   sqrt((l-m)*(l-m-1)*(l+m+1)*(l+m+2)) * Q[l,m+2]      )      (A2)
c
c
c Note on sign convension:
c
c For the needs of GGA PW91 as implemented here, ylm and derivatives
c come with a different sign convention compared to the usual in the
c program: sin(fi)**m --> (-1)**m * sin(fi)**m. Thus some signs change.

c
c Ph.Mavropoulos, Juelich, July 2007
c
      implicit none
c Parameters:
      integer lmaxd,l4maxd
      parameter (lmaxd=4,l4maxd=4*lmaxd)
c Input:
      integer lmax  ! up to which l to calculate
      real*8 v1,v2,v3 ! vector where Ylm etc are calculated (not necessarily normalized)
c Output:
c Y[l,m], dY/dth, dY/dfi, d(dY/dth)/dth, d(dY/dfi)/dfi, d(dY/dth)/dfi
      real*8 Ylm(*),dYdth(*),dYdfi(*),d2Ydth2(*),d2Ydfi2(*),d2Ydthdfi(*)
      real*8 Rabs  ! Norm of input vector (V1,V2,V3)
c Inside:
      real*8 cth,sth,cfi,sfi  ! cos and sin of th and fi
      real*8 pi,fpi,rtwo   ! pi (what else?), 4*pi, sqrt(2)
      real*8 fac      ! factor in construction of polynomials.
      real*8 Plm(0:l4maxd,0:l4maxd)  ! Legendre polynomials
      real*8 Qlm((l4maxd+1)**2)      ! Ylm/cos(m*fi) (m>0) and Ylm/sin(m*fi) (m<0)
      real*8 cmfi(0:l4maxd),smfi(0:l4maxd) ! cos(m*fi) and sin(m*fi)
      real*8 xy,xyz,sgm,sgmm,fi
      real*8 aux
      real*8 tiny
      parameter(tiny=1.d-20)  ! if th < tiny set th=0
      real*8 tt,aa,cd  ! factors in calcul. of Ylm
      integer ll,mm,ii   ! l and m indexes
      integer lmmax      ! (lmax+1)**2, total number of spher. harmonics.
      integer imm,ipm,lpm,lmm,lpmp1,lmmp1 ! i-m,i+m,l+m,l-m,l+m+1,l-m-1

      pi = 4.d0*datan(1.d0)
      fpi = 4.d0*pi
      rtwo = dsqrt(2.d0)
      lmmax = (lmax+1)**2

      if (lmax.gt.l4maxd) stop 'derivylm: lmax out of range.'

c
c--->    calculate sin and cos of theta and phi
c
      xy = v1**2 + v2**2
      xyz = xy + v3**2
c
      Rabs = dsqrt(xyz)
      if (xyz.le.0.0d0) stop 'derivylm: v=0.'

      if (xy.gt.tiny*xyz) then
         xy = dsqrt(xy)
         xyz = dsqrt(xyz)
         cth = v3/xyz
         sth = xy/xyz
         cfi = v1/xy
         sfi = v2/xy

      else
           
         sth = 0.0d0
         cth = 1.0d0
         if (v3.lt.0) cth = -1.0d0
         cfi = 1.0d0
         sfi = 0.0d0
      end if

c First calculate Legendre functions. Use recursion formulas (8.5.3,8.5.5).
c Following taken from KKR program (routine ymy, by M.Weinert).
      fac = 1.0d0
      do 20 mm = 0,lmax - 1
         fac = - dfloat(2*mm-1) * fac
         Plm(mm,mm) = fac
         Plm(mm+1,mm) = dfloat(2*mm+1) * cth * fac
c
c--->    recurse upward in l
c
         do 10 ll = mm + 2,lmax
            Plm(ll,mm) 
     &      = (   dfloat(2*ll-1) * cth * Plm(ll-1,mm)
     &          - dfloat(ll+mm-1) * Plm(ll-2,mm)       ) / dfloat(ll-mm)
 10      continue
         fac = fac*sth
 20   continue
      Plm(lmax,lmax) = - (2*lmax - 1) * fac


c Next calculate Ylm and derivatives.
c
c--->    determine powers of sin and cos of phi
c
      smfi(0) = 0.0d0
      smfi(1) = sfi
      cmfi(0) = 1.0d0
      cmfi(1) = cfi
      do 30 mm = 2,lmax
         smfi(mm) = 2.d0*cfi*smfi(mm-1) - smfi(mm-2)
         cmfi(mm) = 2.d0*cfi*cmfi(mm-1) - cmfi(mm-2)
 30   continue

c For the needs of GGA PW91 as implemented here, ylm and derivatives
c come with a different sign convention compared to the usual in the
c program: sin(fi)**m --> (-1)**m * sin(fi)**m. Thus some signs change.
c This is taken care of here:
      fi = datan2(v2,v1)
!      if (fi.lt.0.d0) then
!         do mm = 1,lmax
!            smfi(mm) = -smfi(mm) 
!         enddo
!      endif
c
c
c
c--->    multiply in the normalization factors;
c        calculate Ylm and derivatives with respect to fi.
c
      ii = 0
      do 50 ll = 0,lmax
         ii = ii + ll + 1
         aa = dsqrt(dfloat(2*ll+1)/fpi)
         cd = 1.d0
         Ylm(ii) = aa * Plm(ll,0)
         dYdfi(ii) = 0.d0
         d2Ydfi2(ii) = 0.d0

         Qlm(ii) = rtwo * aa * Plm(ll,0)


         sgm = -rtwo   ! updated to (-1)**m * rtwo
         sgmm = -1     ! updated to (-1)**m
         do 40 mm = 1,ll
            ipm = ii + mm
            imm = ii - mm
            tt = dfloat((ll+1-mm)*(ll+mm))
            cd = cd / tt
            tt = aa * dsqrt(cd)

            Qlm(ipm) = sgm * tt * Plm(ll,mm)
            Qlm(imm) = sgmm * Qlm(ipm)

            Ylm(ipm) =        Qlm(ipm) * cmfi(mm)
            Ylm(imm) = sgmm * Qlm(imm) * smfi(mm)

            dYdfi(ipm) = -dfloat(mm) *  Ylm(imm)
            dYdfi(imm) =  dfloat(mm) *  Ylm(ipm)
            d2Ydfi2(ipm) = -dfloat(mm*mm) * Ylm(ipm)
            d2Ydfi2(imm) = -dfloat(mm*mm) * Ylm(imm)

            sgm = -sgm
            sgmm = -sgmm

 40      continue
         ii = ii + ll
 50   continue


c Derivatives with respect to th
      call rinit(lmmax,dYdth)
      call rinit(lmmax,d2Ydth2)
      call rinit(lmmax,d2Ydthdfi)
c The l=0 derivatives are zero (established by initialization above).
c Start with l=1.      
      do 70 ll = 1,lmax
         ii = ll*(ll+1) + 1       ! (position of m=0 harmonic in array)
         aa = dfloat(ll*(ll+1))

c Take special care of m=0 harmonic due to 1/sqrt(2)
         dYdth(ii) = -dsqrt(aa) * Qlm(ii+1) / rtwo

         aux = -2.d0 * aa * Qlm(ii)
         if (ll.gt.1) aux = aux + 
     &    (Qlm(ii-2)+Qlm(ii+2)) * dsqrt(dfloat((ll-1)*ll*(ll+1)*(ll+2)))
         d2Ydth2(ii) = 0.25d0 * aux / rtwo
         
         do 60 mm = 1,ll
            ipm = ii + mm
            imm = ii - mm

            lpm = ll + mm
            lmm = ll - mm
            lpmp1 = ll + mm + 1
            lmmp1 = ll - mm + 1
c Apply Eq. (A1)            
            aux = Qlm(ipm-1) * dsqrt(dfloat(lpm*lmmp1))  
            if (mm.lt.ll) 
     &      aux = aux - Qlm(ipm+1) * dsqrt(dfloat(lpmp1*lmm))
            aux = 0.5d0 * aux

            dYdth(ipm) = aux * cmfi(mm)
            dYdth(imm) = aux * smfi(mm)

            d2Ydthdfi(ipm) = -dfloat(mm) * aux * smfi(mm)
            d2Ydthdfi(imm) =  dfloat(mm) * aux * cmfi(mm)

c Apply Eq. (B1)
            aux = -Qlm(ipm) 
     &       * dfloat( lmm*lpmp1 + lpm*lmmp1 )
            if (mm.lt.ll-1) aux = aux + Qlm(ipm+2) 
     &       * dsqrt(dfloat( (lmm-1)*lmm*lpmp1*(lpm+2) ))
            aux = aux + Qlm(ipm-2)
     &       * dsqrt(dfloat( lmmp1*(lmm+2)*(lpm-1)*lpm ))
            aux =  0.25d0 * aux

            d2Ydth2(ipm) = aux * cmfi(mm)
            d2Ydth2(imm) = aux * smfi(mm)

 60      continue
 70   continue


      end

c-----------------------------------------------------------------------
      SUBROUTINE RINIT(N,A)
C **********************************************************************
C * Setting the first N values of a double precision array A to zero   *
C **********************************************************************
C     ..
C     .. Arguments ..
      INTEGER N
      DOUBLE PRECISION A(*)
C     ..
C     .. Locals ..
      INTEGER I,M,MP1
      DOUBLE PRECISION DZERO
C     ..
      DATA DZERO / 0.0D0 /
C     ..
C     ..
      M = MOD(N,5)
      IF ( M.NE.0 ) THEN
         DO I = 1,M
            A(I) = DZERO
         END DO
         IF ( N.LT.5 ) RETURN
      END IF
      MP1 = M + 1
      DO I = MP1,N,5
        A(I  ) = DZERO
        A(I+1) = DZERO
        A(I+2) = DZERO
        A(I+3) = DZERO
        A(I+4) = DZERO
      END DO
C
      END
