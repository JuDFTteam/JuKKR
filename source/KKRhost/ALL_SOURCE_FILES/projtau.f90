SUBROUTINE projtau(icpaflag,cpachng,kmrot,wrtau,wrtaumq,ifiltau,  &
    eryd,nt,nq,nkmq,msst,mssq,nlinq,iqat,conc,tauq,  &
    taut,tautlin,ikm1lin,ikm2lin,drotq,ntmax,nqmax, nkmmax,linmax)
!   ********************************************************************
!   *                                                                  *
!   *   calculate the component projected TAU - matrices               *
!   *                                                                  *
!   *      TAU(IT) =  TAU(IQ) * ( 1 + (m(t)-m(c))*TAU(IQ) )**(-1)      *
!   *                                                                  *
!   *   NOTE: it is assumed that all equivalent sites  IQ  have the    *
!   *   same TAU-matrix  TAUQ(IQ). To get  TAU(IT)  the first site IQ  *
!   *   occupied by type IT is taken to be representative for          *
!   *   all other (NAT(IT)-1) sites occupied by IT                     *
!   *                                                                  *
!   *   allows an atom type IT to have different orientation of        *
!   *   its moment on different but equivalent sites  IQ               *
!   *                                                                  *
!   * 01/11/00                                                         *
!   ********************************************************************
IMPLICIT NONE

! PARAMETER definitions
DOUBLE PRECISION TOL
PARAMETER (TOL=1.0D-6)
DOUBLE COMPLEX C0,C1
PARAMETER (C0=(0.0D0,0.0D0),C1=(1.0D0,0.0D0))

! Dummy arguments
DOUBLE PRECISION CPACHNG
DOUBLE COMPLEX ERYD
INTEGER ICPAFLAG,IFILTAU,KMROT,LINMAX,NKMMAX,NQ,NQMAX,NT,NTMAX
LOGICAL WRTAU,WRTAUMQ
DOUBLE PRECISION CONC(NTMAX)
DOUBLE COMPLEX DROTQ(NKMMAX,NKMMAX,NQMAX), &
           MSSQ(NKMMAX,NKMMAX,NQMAX), &
           MSST(NKMMAX,NKMMAX,NTMAX),TAUQ(NKMMAX,NKMMAX,NQMAX), &
           TAUT(NKMMAX,NKMMAX,NTMAX),TAUTLIN(LINMAX,NTMAX)
INTEGER IKM1LIN(LINMAX),IKM2LIN(LINMAX),IQAT(NTMAX), &
        NKMQ(NQMAX),NLINQ(NQMAX)

! Local variables
DOUBLE PRECISION CPAC
DOUBLE COMPLEX DMAMC(NKMMAX,NKMMAX),DMATTG(NKMMAX,NKMMAX), &
           DTILTG(NKMMAX,NKMMAX),RMSS,RTAU,W1(NKMMAX,NKMMAX)
INTEGER I,ICPAF,IQ,IT,J,LIN,M,N

DO it = 1,nt
  
! ---------- pick first site IQ occupied by type IT to be representative
! ----------- all other (NAT(IT)-1) occupied sites have to be equivalent
  
  iq = iqat(it)
  m = nkmmax
  n = nkmq(iq)
  
  IF ( conc(it) < 0.995 ) THEN
    
! ------------------------- rotate the single site m-matrix if necessary
    IF ( kmrot /= 0 ) THEN
      
      CALL rotate(msst(1,1,it),'L->G',w1,n,drotq(1,1,iq),m)
      
      CALL getdmat(tauq(1,1,iq),dmattg,dtiltg,dmamc,n, mssq(1,1,iq),w1,m)
      
    ELSE
      
      CALL getdmat(tauq(1,1,iq),dmattg,dtiltg,dmamc,n,  &
          mssq(1,1,iq),msst(1,1,it),m)
      
    END IF
    
!     -------------------------------------------
!              TAU(t) = TAU * D~(t)
!     ----------------------------------------
    CALL zgemm('N','N',n,n,n,c1,tauq(1,1,iq),m,dtiltg,m,c0, taut(1,1,it),m)
    
    icpaf = icpaflag
    cpac = cpachng
    
  ELSE
    
!     CONC > 0.995:  COPY TAU TO TAUTLIN
    
    DO j = 1,n
      CALL zcopy(n,tauq(1,j,iq),1,taut(1,j,it),1)
    END DO
    
    icpaf = 0
    cpac = 0D0
    
  END IF
  
!     -------------------------------------------
!            rotate  TAU(t)  if required
!     -------------------------------------------
  
  IF ( kmrot /= 0 ) THEN
    
    DO j = 1,n
      CALL zcopy(n,taut(1,j,it),1,w1(1,j),1)
    END DO
    
    CALL rotate(w1,'G->L',taut(1,1,it),n,drotq(1,1,iq),m)
    
  END IF
  
!     -------------------------------------------
!        STORE TAU(t) IN LINEAR ARRAY TAUTLIN
!     -------------------------------------------
  
  DO lin = 1,nlinq(iq)
    tautlin(lin,it) = taut(ikm1lin(lin),ikm2lin(lin),it)
  END DO
  
  IF ( wrtau ) THEN
    WRITE (ifiltau,99001) eryd,it,iq,icpaf,cpac
    DO i = 1,n
      DO j = 1,n
        IF ( i == j ) THEN
          WRITE (ifiltau,99003) i,j,taut(i,j,it)
        ELSE
          IF ( CDABS(taut(i,j,it)/taut(i,i,it)) > tol )  &
              WRITE (ifiltau,99003) i,j,taut(i,j,it)
        END IF
      END DO
    END DO
  END IF
  
END DO
!================================================================= IT ==

IF ( wrtaumq ) THEN
  DO iq = 1,nq
    WRITE (ifiltau,99002) eryd,iq,icpaflag,cpachng
    DO i = 1,n
      DO j = 1,n
        IF ( i == j ) THEN
          WRITE (ifiltau,99003) i,j,tauq(i,j,iq),mssq(i,j,iq)
        ELSE
          rtau = tauq(i,j,iq)/tauq(i,i,iq)
          rmss = mssq(i,j,iq)/mssq(i,i,iq)
          IF ( (CDABS(rtau) > tol) .OR. (CDABS(rmss) > tol)  &
              ) WRITE (ifiltau,99003) i,j,tauq(i,j,iq), mssq(i,j,iq)
        END IF
        
      END DO
    END DO
    
  END DO
END IF
!--------------------------------------------------------- FORMAT IFMT=2
99001 FORMAT (/,80('*')/,2F21.15,' RYD   TAU FOR IT=',i2,'  IQ=',i2,:,  &
    '  CPA:',i2,f15.6)
99002 FORMAT (/,80('*')/,2F21.15,' RYD   TAU-C M-C  FOR IQ=',i2,:,  &
    '  CPA:',i2,f15.6)
99003 FORMAT (2I5,1P,4E22.14)

END SUBROUTINE projtau
