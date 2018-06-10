SUBROUTINE cinthff(ag,af,bg,bf,rmehf,nka,nkb,jtop,fx,r,drdi,nrmax)
!   ********************************************************************
!   *                                                                  *
!   *  routine to calculate the radial hyperfine matrixelement         *
!   *  by extrapolating the lower integration boundary to r -> 0       *
!   *                                                                  *
!   ********************************************************************
IMPLICIT NONE

!Dummy arguments
INTEGER JTOP,NKA,NKB,NRMAX
DOUBLE COMPLEX AF(NRMAX,NKA),AG(NRMAX,NKA),BF(NRMAX,NKB), &
               BG(NRMAX,NKB),FX(JTOP),RMEHF(2,2)
DOUBLE PRECISION DRDI(NRMAX),R(NRMAX)

! Local variables
DOUBLE PRECISION AI,AR,BI,BR,DELTA,HIBAR,R1BAR,RELDIF, &
                 RLIM1,RLIM2,S,SI,SR,SX,SXI,SXR,SXX,YI,YR

DOUBLE COMPLEX HI,HIF,Z(NRMAX)
INTEGER I,IMAX,IMIN,KA,KB,N

DO kb = 1,nkb
  DO ka = 1,nka
    
    DO i = 1,jtop
      fx(i) = (ag(i,ka)*bf(i,kb)+af(i,ka)*bg(i,kb))*drdi(i)
    END DO
    
    CALL cint4pts(fx,jtop,z)
    
    rlim1 = 0.5D-5
    rlim2 = 0.5D-4
!            RLIM1 = 1D-4
!            RLIM2 = 5D-4
    IF ( r(1) > rlim1 ) THEN
      rlim1 = r(1)
      rlim2 = r(20)
    END IF
    
    imin = 0
    imax = 0
    DO i = 1,jtop
      IF ( r(i) <= rlim1 ) imin = i
      IF ( r(i) >= rlim2 ) THEN
        imax = i
        GO TO 20
      END IF
    END DO
    20         CONTINUE
    IF ( imin == 0 ) STOP '<CINTHFF> IMIN = 0'
    IF ( imin > jtop ) STOP '<CINTHFF> IMIN > JTOP'
    IF ( imax == 0 ) STOP '<CINTHFF> IMAX = 0'
    IF ( imax > jtop ) STOP '<CINTHFF> IMAX > JTOP'
    
    n = 0
    s = 0.0D0
    sx = 0.0D0
    sxx = 0.0D0
    sr = 0.0D0
    si = 0.0D0
    sxr = 0.0D0
    sxi = 0.0D0
    
    DO i = imin,imax
      yr = dreal(z(jtop)-z(i))
      yi = DIMAG(z(jtop)-z(i))
      n = n + 1
      s = s + 1.0D0
      sx = sx + r(i)
      sxx = sxx + r(i)**2
      sr = sr + yr
      si = si + yi
      sxr = sxr + yr*r(i)
      sxi = sxi + yi*r(i)
    END DO
    
    delta = s*sxx - sx*sx
    ar = (sxx*sr-sx*sxr)/delta
    ai = (sxx*si-sx*sxi)/delta
    br = (s*sxr-sx*sr)/delta
    bi = (s*sxi-sx*si)/delta
    
    hibar = 0.0D0
    r1bar = 0.0D0
    DO i = imin,imax
      hibar = hibar + DBLE(z(jtop)-z(i))
      r1bar = r1bar + r(i)
    END DO
    hibar = hibar/DBLE(imax-imin+1)
    r1bar = r1bar/DBLE(imax-imin+1)
    
    sxx = 0.0D0
    sxr = 0.0D0
    sxi = 0.0D0
    reldif = 0.0D0
    DO i = imin,imax
      hi = z(jtop) - z(i)
      hif = DCMPLX(ar,ai) + DCMPLX(br,bi)*r(i)
      IF ( ABS(hif) /= 0.0D0 ) reldif = MAX(reldif,ABS(1.0D0-hi/hif))
    END DO
    
    rmehf(ka,kb) = DCMPLX(ar,ai)
    
  END DO
END DO

END SUBROUTINE cinthff
