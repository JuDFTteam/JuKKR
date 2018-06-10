SUBROUTINE cpamillsx(itcpa,cpaerr,cpacorr,cpachng,iprint,icpa,nq,  &
    nkmq,noq,itoq,conc,mssq,msst,tauq,dmssq, kmrot,drotq,ntmax,nqmax,nkmmax)
!   ********************************************************************
!   *                                                                  *
!   *   perform  CPA-iteration according    MILLS's  algorithm         *
!   *                                                                  *
!   *   the CPA - iteration step for site IQ is omitted if             *
!   *   ICPA(IQ) = 0 ( set in <INITALL> )                              *
!   *                                                                  *
!   *   only the projection matrix  DTILT(G) is needed                 *
!   *   this is set up with respect to the global frame                *
!   *   for this reason MSST has to be rotated prior calling <GETDMAT> *
!   *                                                                  *
!   *   allows an atom type IT to have different orientation of        *
!   *   its moment on different but equivalent sites  IQ               *
!   *                                                                  *
!   * 15/12/03                                                         *
!   ********************************************************************
use mod_types, only: t_inc
IMPLICIT COMPLEX*16(A-H,O-Z)

! PARAMETER definitions
DOUBLE PRECISION TOL,SCLSTD
PARAMETER (TOL=10D0,SCLSTD=1D0)
DOUBLE COMPLEX C0,C1
PARAMETER (C0=(0.0D0,0.0D0),C1=(1.0D0,0.0D0))

! Dummy arguments
DOUBLE PRECISION CPACHNG,CPACORR,CPAERR
INTEGER IPRINT,ITCPA,KMROT,NKMMAX,NQ,NQMAX,NTMAX
DOUBLE PRECISION CONC(NTMAX)
DOUBLE COMPLEX DROTQ(NKMMAX,NKMMAX,NQMAX), &
     MSSQ(NKMMAX,NKMMAX,NQMAX),DMSSQ(NKMMAX,NKMMAX,NQMAX), &
     MSST(NKMMAX,NKMMAX,NTMAX),TAUQ(NKMMAX,NKMMAX,NQMAX)
INTEGER ICPA(NQMAX),ITOQ(NTMAX,NQMAX),NKMQ(NQMAX),NOQ(NQMAX)

! Local variables
LOGICAL CHECK
DOUBLE PRECISION CPACHNGL,CPACORRL,SCL
DOUBLE COMPLEX CSUM
DOUBLE COMPLEX DMAMC(NKMMAX,NKMMAX),DMATTG(NKMMAX,NKMMAX), &
     DQ(NKMMAX,NQMAX), &
     DTILTG(NKMMAX,NKMMAX),ERR(NKMMAX,NKMMAX), &
     W1(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX)
INTEGER I,ICPARUN,INFO,IO
INTEGER IPIV(NKMMAX),IQ,IT,IW,IW0,J,M,N
DOUBLE PRECISION P1,P2
SAVE CPACHNGL,CPACORRL,SCL

DATA ICPARUN/0/

cpaerr = 0.0D0
cpacorr = 0.0D0
cpachng = 0.0D0
check = .true.
check = .false.

IF ( itcpa == 1 ) THEN
  scl = sclstd
  cpachngl = 1D+20
  cpacorrl = 1D+20
  icparun = icparun + 1
END IF

DO iq = 1,nq
  IF ( icpa(iq) /= 0 ) THEN
    
    m = nkmmax
    n = nkmq(iq)
    
    DO j = 1,n
      dq(j,iq) = -c1
      CALL cinit(n,ERR(1,j))
    END DO
    
!================================================================= IT ==
    DO io = 1,noq(iq)
      it = itoq(io,iq)
      
! ------------------------- rotate the single site m-matrix if necessary
      IF ( kmrot /= 0 ) THEN
        
        CALL rotate(msst(1,1,it),'L->G',w1,n,drotq(1,1,iq),m)
        
        CALL getdmat(tauq(1,1,iq),dmattg,dtiltg,dmamc,n, mssq(1,1,iq),w1,m)
        
      ELSE
        
        CALL getdmat(tauq(1,1,iq),dmattg,dtiltg,dmamc,n,  &
            mssq(1,1,iq),msst(1,1,it),m)
        
      END IF
      
      DO i = 1,n
        dq(i,iq) = dq(i,iq) + conc(it)*dtiltg(i,i)
      END DO
      
!     -------------------------------------------
!            - E[a] = D~[a] * ( m[a] - m[c] )
!     -------------------------------------------
      CALL zgemm('N','N',n,n,n,c1,dtiltg(1,1),m,dmamc,m,c0,w1, m)
      
!     -------------------------------------------
!            E = SUM[a]  c[a] *  E[a]
!     -------------------------------------------
      DO j = 1,n
        DO i = 1,n
          ERR(i,j) = ERR(i,j) - conc(it)*w1(i,j)
        END DO
      END DO
      
    END DO
!================================================================= IT ==
    
!     -------------------------------------------
!                   E * TAU
!     -------------------------------------------
    
    CALL zgemm('N','N',n,n,n,c1,ERR,m,tauq(1,1,iq),m,c0,w2,m)
    
!     -------------------------------------------
!                1 + E * TAU
!     -------------------------------------------
    DO i = 1,n
      w2(i,i) = c1 + w2(i,i)
    END DO
    
!     -------------------------------------------
!               ( 1 + E * TAU )**(-1)
!     -------------------------------------------
    
    CALL zgetrf(n,n,w2,m,ipiv,info)
    CALL zgetri(n,w2,m,ipiv,w1,m*m,info)
    
!     -------------------------------------------
!           ( 1 + E * TAU )**(-1) * E
!     -------------------------------------------
    
    CALL zgemm('N','N',n,n,n,c1,w2,m,ERR,m,c0,w1,m)
    
!     -------------------------------------------
!           APPLY CORRECTION  TO  MEFF
!     m{n+1} = m{n} -  ( 1 + E * TAU )**(-1) * E
!     -------------------------------------------
    DO j = 1,n
      cpaerr = cpaerr + ABS(dreal(dq(j,iq))) + ABS(DIMAG(dq(j,iq)))
      cpacorr = cpacorr + ABS(dreal(w1(j,j))) + ABS(DIMAG(w1(j,j)))
      cpachng = MAX(cpachng,ABS(w1(j,j)/mssq(j,j,iq)))
    END DO
    cpachng = scl*cpachng
    cpacorr = scl*cpacorr
    
    IF ( cpachng > tol*cpachngl .OR. cpacorr > cpacorrl ) THEN
      WRITE (*,*) '############### CPA step back'
      
      p1 = 0.5D0   ! P1 = 0.05D0
      p2 = 1D0 - p1
      cpachng = p1*cpachngl
      cpacorr = p1*cpacorrl
      cpachng = cpachngl
      cpacorr = cpacorrl
      scl = p1
      DO j = 1,n
        DO i = 1,n
          mssq(i,j,iq) = mssq(i,j,iq) + p2*dmssq(i,j,iq)
          dmssq(i,j,iq) = dmssq(i,j,iq)*p1
        END DO
      END DO
    ELSE
      DO j = 1,n
        DO i = 1,n
          w1(i,j) = scl*w1(i,j)
          mssq(i,j,iq) = mssq(i,j,iq) - w1(i,j)
          dmssq(i,j,iq) = w1(i,j)
        END DO
      END DO
    END IF
    
    cpaerr = cpachng
    
    IF ( iprint >= 2 .OR. check ) THEN
      csum = c0
      DO i = 1,n
        csum = csum + mssq(i,i,iq)
      END DO
      IF(t_inc%i_write>0) THEN
        WRITE (1337,99001) iq,cpaerr,cpacorr,csum
      END IF
    END IF
!-----------------------------------------------------------------------
    IF ( check ) THEN
      IF ( icparun == 2 ) STOP 'CPA-iter written to for...'
      iw0 = 100*iq
      DO i = 1,n
        iw = iw0 + i
        WRITE (iw,'(I4,4E14.4)') itcpa,mssq(i,i,iq),w1(i,i)
      END DO
    END IF
!-----------------------------------------------------------------------
  END IF
  
END DO
!================================================================= IQ ==

cpachngl = cpachng
cpacorrl = cpacorr

99001 FORMAT (' CPA:  IQ',i3,'  ERR',f12.5,'  CORR',f13.5,'  M',  &
    18(1X,2(1PE14.6)))
END SUBROUTINE cpamillsx
