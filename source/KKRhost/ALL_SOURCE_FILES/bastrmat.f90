SUBROUTINE bastrmat(lmax,cgc,rc,crel,rrel,nkmmax,nkmpmax)
!   ********************************************************************
!   *                                                                  *
!   *    INITIALIZE TRANSFORMATION MATRIX THAT TAKES MATRICES FROM     *
!   *    RELATIVISTIC  TO  REAL SPERICAL HARM.  REPRESENTATION         *
!   *                                                                  *
!   *    this is a special version of <STRSMAT> passing the            *
!   *    full BASis TRansformation MATrices  RC, CREL and RREL         *
!   *                                                                  *
!   * 13/01/98  HE                                                     *
!   ********************************************************************
IMPLICIT none

! PARAMETER definitions
COMPLEX*16 CI,C1,C0
PARAMETER (CI=(0.0D0,1.0D0),C1=(1.0D0,0.0D0),C0=(0.0D0,0.0D0))

! Dummy arguments
INTEGER LMAX,NKMMAX,NKMPMAX
REAL*8 CGC(NKMPMAX,2)
COMPLEX*16 CREL(NKMMAX,NKMMAX),RC(NKMMAX,NKMMAX), &
           RREL(NKMMAX,NKMMAX)

! Local variables
INTEGER I,IKM,J,JP05,K,L,LM,LNR,M,MUEM05,MUEP05,NK,NKM,NLM
REAL*8 W

nk = 2*(lmax+1) + 1
nlm = (lmax+1)**2
nkm = 2*nlm
!     ===================================================
!     INDEXING:
!     IKM  = L*2*(J+1/2) + J + MUE + 1
!     LM   = L*(L+1)     +     M   + 1
!     ===================================================

! ----------------------------------------------------------------------
! CREL  transforms from  COMPLEX (L,M,S)  to  (KAP,MUE) - representation
!                 |LAM> = sum[LC] |LC> * CREL(LC,LAM)
! ----------------------------------------------------------------------
CALL cinit(nkmmax*nkmmax,crel)

lm = 0
DO lnr = 0,lmax
  DO m = -lnr,lnr
    lm = lm + 1
    
    ikm = 0
    DO k = 1,nk
      l = k/2
      IF ( 2*l == k ) THEN
        jp05 = l
      ELSE
        jp05 = l + 1
      END IF
      
      DO muem05 = -jp05,(jp05-1)
        muep05 = muem05 + 1
        ikm = ikm + 1
        
        IF ( l == lnr ) THEN
          IF ( muep05 == m ) crel(lm,ikm) = cgc(ikm,1)
          IF ( muem05 == m ) crel(lm+nlm,ikm) = cgc(ikm,2)
        END IF
        
      END DO
    END DO
    
  END DO
END DO

! ----------------------------------------------------------------------
!    RC  transforms from  REAL to  COMPLEX (L,M,S) - representation
!                 |LC> = sum[LR] |LR> * RC(LR,LC)
! ----------------------------------------------------------------------
CALL cinit(nkmmax*nkmmax,rc)

w = 1.0D0/SQRT(2.0D0)

DO l = 0,lmax
  DO m = -l,l
    i = l*(l+1) + m + 1
    j = l*(l+1) - m + 1
    
    IF ( m < 0 ) THEN
      rc(i,i) = -ci*w
      rc(j,i) = w
      rc(i+nlm,i+nlm) = -ci*w
      rc(j+nlm,i+nlm) = w
    END IF
    IF ( m == 0 ) THEN
      rc(i,i) = c1
      rc(i+nlm,i+nlm) = c1
    END IF
    IF ( m > 0 ) THEN
      rc(i,i) = w*(-1.0D0)**m
      rc(j,i) = ci*w*(-1.0D0)**m
      rc(i+nlm,i+nlm) = w*(-1.0D0)**m
      rc(j+nlm,i+nlm) = ci*w*(-1.0D0)**m
    END IF
  END DO
END DO

! ----------------------------------------------------------------------
! RREL  transforms from   REAL (L,M,S)  to  (KAP,MUE) - representation
!                 |LAM> = sum[LR] |LR> * RREL(LR,LAM)
! ----------------------------------------------------------------------

CALL zgemm('N','N',nkm,nkm,nkm,c1,rc,nkmmax,crel,nkmmax,c0,rrel, nkmmax)

END SUBROUTINE bastrmat

