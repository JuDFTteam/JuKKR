SUBROUTINE strsmat(lmax,cgc,srrel,nrrel,irrel,nkmmax,nkmpmax)
!   ********************************************************************
!   *                                                                  *
!   *    INITIALIZE TRANSFORMATION MATRIX THAT TAKES MATRICES FROM     *
!   *    RELATIVISTIC  TO  REAL SPERICAL HARM.  REPRESENTATION         *
!   *                                                                  *
!   *    ONLY THE NON-0 ELEMENTS OF THE MATRIX ARE STORED              *
!   *                                                                  *
!   * 25/10/95  HE  proper convention of trans. matrix introduced      *
!   ********************************************************************

IMPLICIT NONE

! PARAMETER definitions
DOUBLE COMPLEX CI,C1,C0
PARAMETER (CI=(0.0D0,1.0D0),C1=(1.0D0,0.0D0),C0=(0.0D0,0.0D0))

! Dummy arguments
INTEGER LMAX,NKMMAX,NKMPMAX
DOUBLE PRECISION CGC(NKMPMAX,2)
INTEGER IRREL(2,2,NKMMAX),NRREL(2,NKMMAX)
DOUBLE COMPLEX SRREL(2,2,NKMMAX)

! Local variables
DOUBLE COMPLEX CREL(NKMMAX,NKMMAX),RC(NKMMAX,NKMMAX), &
           RREL(NKMMAX,NKMMAX)
INTEGER I,IKM,J,JP05,K,L,LAM,LM,LNR,LR,M,MUEM05,MUEP05,NK,NKM,NLM, &
        NS1,NS2
DOUBLE PRECISION W

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

!     ---------------------------------------------------
!     store the elements of  RREL
!     ---------------------------------------------------
DO lam = 1,nkm
  ns1 = 0
  ns2 = 0
  
  DO lr = 1,2*nlm
!             IF ( CDABS(RREL(LR,LAM)).GT.1D-6 ) THEN
    IF ( CDABS(rrel(lr,lam)) > 1D-4 ) THEN
      IF ( lr <= nlm ) THEN
        ns1 = ns1 + 1
        IF ( ns1 > 2 ) STOP ' IN <STRSMAT>   NS1 > 2'
        srrel(ns1,1,lam) = rrel(lr,lam)
        irrel(ns1,1,lam) = lr
      ELSE
        ns2 = ns2 + 1
        IF ( ns2 > 2 ) STOP ' IN <STRSMAT>   NS2 > 2'
        srrel(ns2,2,lam) = rrel(lr,lam)
        irrel(ns2,2,lam) = lr - nlm
      END IF
    END IF
  END DO
  
  nrrel(1,lam) = ns1
  nrrel(2,lam) = ns2
END DO

END SUBROUTINE strsmat
