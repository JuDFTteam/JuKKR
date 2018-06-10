SUBROUTINE cylm02(lmax,cosx,fai,lpot2p,lmmaxd,thet,ylm,dylmt1,  &
        dylmt2,dylmf1,dylmf2,dylmtf)
!.....------------------------------------------------------------------
!     preparation of cylm0(=ylm(ip,i)), cylmt1(=dylm/dtheta),
!     cylmt2(=d2ylm/dt2),
!     cylmf1, cylmf2 are for fai.
!     cylmtf=d2ylm/dfdt
!     i=1,2,....,(lmax+1)**2
!.....------------------------------------------------------------------
implicit none

!.. Parameters ..
INTEGER IJD
PARAMETER (IJD=434)
!..
!.. Scalar Arguments ..
INTEGER LMAX,LMMAXD,LPOT2P
!..
!.. Array Arguments ..
DOUBLE PRECISION COSX(IJD),DYLMF1(IJD,LMMAXD),DYLMF2(IJD,LMMAXD), &
                 DYLMT1(IJD,LMMAXD),DYLMT2(IJD,LMMAXD), &
                 DYLMTF(IJD,LMMAXD),FAI(IJD),THET(IJD), &
                 YLM(IJD,LMMAXD)
!..
!.. Local Scalars ..
DOUBLE COMPLEX CI,EM1F,EM2F,EP1F,EP2F
DOUBLE PRECISION AAA,CCC,DI,FI,ONE,PI,SSS
INTEGER I,IP,L,LLMAX,LM,LM1,LM1M,LM2,LMM,LMM1,LMM1M,LMM2,M,MM
!..
!.. Local Arrays ..
DOUBLE COMPLEX CYLM0(LMMAXD),CYLMF1(LMMAXD),CYLMF2(LMMAXD), &
               CYLMT1(LMMAXD),CYLMT2(LMMAXD),CYLMTF(LMMAXD)
DOUBLE PRECISION BB1(LMMAXD),YL(LPOT2P)
!..
!.. External Subroutines ..
EXTERNAL SPHER,TRAREA
!..
!.. Intrinsic Functions ..
INTRINSIC ACOS,ATAN,CMPLX,CONJG,COS,DBLE,SIN,SQRT
!..

ci = CMPLX(0.d0,1.d0)
one = 1.d0
pi = 4.d0*ATAN(one)
llmax = (lmax+1)**2

DO  ip = 1,ijd
  
  thet(ip) = ACOS(cosx(ip))
  fi = fai(ip)
  di = 2*fai(ip)
  ep1f = CMPLX(COS(fi),SIN(fi))
  em1f = CONJG(ep1f)
  ep2f = CMPLX(COS(di),SIN(di))
  em2f = CONJG(ep2f)
  
  DO  l = 0,lmax
    
    CALL spher(yl,l,cosx(ip))
    DO  m = -l,l
      mm = l + m + 1
      i = (l+1)**2 - l + m
      aaa = m*fai(ip)
      ccc = COS(aaa)
      sss = SIN(aaa)
      cylm0(i) = yl(mm)*CMPLX(ccc,sss)
    END DO
    
    DO  m = -l,l
      i = (l+1)**2 - l + m
      cylmt1(i) = 0.d0
      cylmt2(i) = 0.d0
      cylmtf(i) = 0.d0
    END DO
    
    DO  m = -l,l
      i = (l+1)**2 - l + m
      
      lmm1m = l - m - 1
      lmm = l - m
      lmm1 = l - m + 1
      lmm2 = l - m + 2
      lm1m = l + m - 1
      lm = l + m
      lm1 = l + m + 1
      lm2 = l + m + 2
      
      cylmt2(i) = cylmt2(i) - (lmm*lm1+lmm1*lm)/4.d0*cylm0(i)
      
      IF (m+2 <= l) cylmt2(i) = cylmt2(i) +  &
          SQRT(DBLE(lmm1m*lmm*lm1*lm2))/4* cylm0(i+2)*em2f
      
      IF (m+1 <= l) cylmt1(i) = cylmt1(i) +  &
          SQRT(DBLE(lmm*lm1))/2*cylm0(i+1)* em1f
      
      IF (m-1 >= -l) cylmt1(i) = cylmt1(i) -  &
          SQRT(DBLE(lm*lmm1))/2*cylm0(i-1)* ep1f
      
      IF (m-2 >= -l) cylmt2(i) = cylmt2(i) +  &
          SQRT(DBLE(lmm1*lmm2*lm1m*lm))/4* cylm0(i-2)*ep2f
      
    END DO
    
    DO  m = -l,l
      i = (l+1)**2 - l + m
      cylmf1(i) = ci*m*cylm0(i)
      cylmf2(i) = -m*m*cylm0(i)
      cylmtf(i) = ci*m*cylmt1(i)
    END DO
    
  END DO
  
!        calculate real spherical harmonics differenciated
  
  
!        write(6,9005) (cylm0(i),i=1,5)
!9005 format(1x,' cylm0',4f10.5)
  CALL trarea(cylm0,bb1,lmax)
  
  DO  m = 1,llmax
    ylm(ip,m) = bb1(m)
  END DO
  
!        write(6,9006) (ylm(ip,i),i=1,5)
!9006 format(1x,' ylm',10f10.5)
  
  
  CALL trarea(cylmt1,bb1,lmax)
  DO  m = 1,llmax
    dylmt1(ip,m) = bb1(m)
  END DO
  
  CALL trarea(cylmt2,bb1,lmax)
  DO  m = 1,llmax
    dylmt2(ip,m) = bb1(m)
  END DO
  
  CALL trarea(cylmf1,bb1,lmax)
  DO  m = 1,llmax
    dylmf1(ip,m) = bb1(m)
  END DO
  
  CALL trarea(cylmf2,bb1,lmax)
  DO  m = 1,llmax
    dylmf2(ip,m) = bb1(m)
  END DO
  
  CALL trarea(cylmtf,bb1,lmax)
  DO  m = 1,llmax
    dylmtf(ip,m) = bb1(m)
  END DO
  
END DO
RETURN
END SUBROUTINE cylm02
