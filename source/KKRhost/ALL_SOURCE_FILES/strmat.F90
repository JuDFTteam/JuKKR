SUBROUTINE strmat(alat,lpot,naez,ngmax,nrmax,nsg,nsr,nshlg,nshlr,  &
    gn,rm,qi0,smat,vol,iprint,lassld,lmxspd,naezd)
! **********************************************************************
! *                                                                    *
! *  calculation of lattice sums for l .le. 2*lpot :                   *
! *                                                                    *
! *                   ylm( q(i) - q(j) - rm )                          *
! *        sum      ===========================                        *
! *                 | q(i) - q(j) - rm |**(l+1)                        *
! *                                                                    *
! *         - summed over all lattice vectors rm  -                    *
! *                                                                    *
! *  ylm       : real spherical harmic to given l,m                    *
! *  q(i),q(j) : basis vectors of the unit cell                        *
! *                                                                    *
! *  in the case of i = j, rm = 0 is omitted.                          *
! *                                                                    *
! *  the ewald method is used to perform the lattice summations        *
! *  the splitting parameter lamda is set equal sqrt(pi)/alat          *
! *  (alat is the lattice constant) .                                  *
! *                                                                    *
! *  if the contribution of the last shell of the direct and the       *
! *  reciprocal lattice is greater than 1.0e-8 a message is written    *
! *                                                                    *
! *                                    b.drittler may 1989             *
! *                                                                    *
! *  Dimension of arrays GN,RM changed from (4,*) to (3,*), the 4th    *
! *  one not being used (see also lattice3d)     v.popescu May 2004    *
! *                                                                    *
! **********************************************************************
#ifdef CPP_HYBRID
#define CPP_OMPSTUFF
#endif
#ifdef CPP_OMP
#define CPP_OMPSTUFF
#endif
#ifdef CPP_OMPSTUFF
use omp_lib
#endif

use constants
      Use mod_datatypes, Only: dp

      IMPLICIT NONE
!..
!.. Parameters ..
      real (kind=dp) BOUND
      PARAMETER ( BOUND=1D-8 )
!..
!.. Scalar arguments ..
      real (kind=dp) ALAT,VOL
      INTEGER IPRINT,LPOT,NAEZ,NGMAX,NRMAX,NSHLG,NSHLR
      INTEGER LASSLD,LMXSPD,NAEZD
!..
!.. Array arguments ..
      real (kind=dp)  GN(3,*),QI0(3,*),RM(3,*),SMAT(LMXSPD,NAEZD,*)
      INTEGER NSG(*),NSR(*)
!..
!.. Local scalars ..
      complex (kind=dp) BFAC
      real (kind=dp) ALPHA,BETA,DQ1,DQ2,DQ3,DQDOTG,EXPBSQ,FPI, &
                       G1,G2,G3,GA,LAMDA,R,R1,R2,R3,RFAC,S
      INTEGER I,I1,I2,IT,L,LM,LMX,LMXSP,LFMT,M,NGE,NGS,NRE,NRS,NSTART
      CHARACTER*80 FMT
!..
!.. Local arrays ..
      complex (kind=dp) STEST(LMXSPD)
      real (kind=dp) G(0:LASSLD),YLM(LMXSPD),QI(3,NAEZD)
!..
!.. External subroutines ..
      EXTERNAL GAMFC,YMY
!..
!.. Intrinsic functions ..
      INTRINSIC ATAN,ABS,DBLE,EXP,SQRT
!     ..................................................................

lmx = 2*lpot
lmxsp = (lmx+1)*(lmx+1)
fpi = 4.0D0*pi

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,'(5X,2A,/)') '< STRMAT > : ', 'calculating lattice sums'
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT


! --> choose proper splitting parameter

lamda = SQRT(pi)/alat

! --> loop over atoms per unit cell -- scale basis atoms with alat

DO i2 = 1,naez
  DO i1 = 1,3
    qi(i1,i2) = qi0(i1,i2)*alat
  END DO
END DO

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
IF ( iprint > 2 ) THEN
  WRITE (1337,99002)
  DO i2 = 1,naez
    WRITE (1337,99003) i2,(qi0(i,i2),i=1,3)
  END DO
  WRITE (1337,99004)
endif
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! **********************************************************************
#ifdef CPP_OMPSTUFF
!$omp parallel do default(shared) private(DQ1, DQ2, DQ3, STEST)
!$omp& private(LM, NSTART, IT, NRS, NGS, NRE, NGE, I, R1, R2)
!$omp& private(R3 , R, YLM, ALPHA, G, L, RFAC, M, G1, G2)
!$omp& private(G3, GA, BETA, EXPBSQ, DQDOTG, BFAC, S, I1, I2)
#endif
DO i1 = 1,naez
  DO i2 = 1,naez
!======================================================================
    dq1 = qi(1,i1) - qi(1,i2)
    dq2 = qi(2,i1) - qi(2,i2)
    dq3 = qi(3,i1) - qi(3,i2)
    
    stest(1) = -SQRT(fpi)/vol/(4D0*lamda*lamda)
    DO lm = 2,lmxsp
      stest(lm) = 0.0D0
    END DO
    
! --> exclude the origine and add correction if i1.eq.i2
    
    IF ( i1 == i2 ) THEN
      stest(1) = stest(1) - lamda/pi
      nstart = 2
    ELSE
      nstart = 1
    endif
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
! --> loop first over n-1 shells of real and reciprocal lattice - then
!     add the contribution of the last shells to see convergence
    
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO it = 1,2
      IF ( it == 1 ) THEN
        nrs = nstart
        ngs = 2
        nre = nrmax - nsr(nshlr)
        nge = ngmax - nsg(nshlg)
      ELSE
        nrs = nre + 1
        ngs = nge + 1
        nre = nrmax
        nge = ngmax
      endif
      
! --> sum over real lattice
      
! ---------------------------------------------------------------------
      DO i = nrs,nre
        r1 = dq1 - rm(1,i)
        r2 = dq2 - rm(2,i)
        r3 = dq3 - rm(3,i)
        
        CALL ymy(r1,r2,r3,r,ylm,lmx)
        alpha = lamda*r
        CALL gamfc(alpha,g,lmx,r)
        
!                   IF (R>0.3D0 .or. .not. OPT('VIRATOMS') ) THEN         ! added by bauer VIRTAOM
        DO l = 0,lmx
          rfac = g(l)/SQRT(pi)
          DO m = -l,l
            lm = l*(l+1) + m + 1
            stest(lm) = stest(lm) + ylm(lm)*rfac
          END DO
        END DO
!                   ELSE
!                   write(*,*) 'omitting ',R
!                   endif !(R>0.3) THEN         ! added by bauer VIRTAOM
        
      END DO
! ---------------------------------------------------------------------
      
! --> sum over reciprocal lattice
      
! ---------------------------------------------------------------------
      DO i = ngs,nge
        g1 = gn(1,i)
        g2 = gn(2,i)
        g3 = gn(3,i)
        
        CALL ymy(g1,g2,g3,ga,ylm,lmx)
        beta = ga/lamda
        expbsq = EXP(beta*beta/4.0D0)
        dqdotg = dq1*g1 + dq2*g2 + dq3*g3
        
        bfac = fpi*EXP(ci*dqdotg)/(ga*ga*expbsq*vol)
        
        DO l = 0,lmx
          DO m = -l,l
            lm = l*(l+1) + m + 1
            stest(lm) = stest(lm) + ylm(lm)*bfac
          END DO
          bfac = bfac*ga/ci/real(2*l+1, kind=dp)
        END DO
      END DO
! ---------------------------------------------------------------------
      IF ( it == 1 ) THEN
        DO lm = 1,lmxsp
          IF ( ABS(AIMAG(stest(lm))) > bound ) THEN
            WRITE (6,*) ' ERROR: Imaginary contribution',  &
                ' to REAL lattice sum'
            STOP
          endif
          smat(lm,i1,i2) = real(stest(lm), kind=dp)
          stest(lm) = 0.0D0
        END DO
      ELSE
        
! --> test convergence
        
#ifdef CPP_OMPSTUFF
!$omp critical
#endif
      DO lm = 1,lmxsp
        s = real(stest(lm), kind=dp)
        smat(lm,i1,i2) = smat(lm,i1,i2) + s
        IF ( ABS(s) > bound ) WRITE (6,FMT=99001) i1,i2, lm,ABS(s)
      END DO
#ifdef CPP_OMPSTUFF
!$omp end critical
#endif
  endif
! ---------------------------------------------------------------------
END DO
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END DO
END DO
#ifdef CPP_OMPSTUFF
!$omp end parallel do
#endif
! **********************************************************************

IF ( iprint < 2 ) RETURN

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,99005) MIN(6,naez)
FMT = ' '
lfmt = 0
DO i = 1,MIN(6,naez)
  FMT = FMT(1:lfmt)//'------------'
  lfmt = lfmt + 12
END DO
WRITE (1337,'(7X,A)') FMT(1:lfmt)
DO i1 = 1,MIN(6,naez)
  WRITE (1337,99006) (smat(1,i1,i2),i2=1,MIN(6,naez))
END DO
WRITE (1337,'(7X,A,/)') FMT
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

99001 FORMAT (5X,'WARNING : Convergence of SMAT(',i2,',',i2,') ',  &
    ' for LMXSP =',i3,' is ',1P,d9.2,' > 1D-8',/,15X,  &
    'You should use more lattice vectors (RMAX/GMAX)')
99002 FORMAT (12X,47('-'),/,16X,'      Positions of atomic sites',/,16X,  &
    '    in CARTESIAN coordinates (a.u.)',/,12X,47('-'),/,15X,  &
    'IQ       x           y           z     ',/,12X,47('-'))
99003 FORMAT (13X,i5,3F12.6)
99004 FORMAT (12X,47('-'),/)
99005 FORMAT (8X,'Lattice sum (LMXSP = 1) up to NAEZ =',i2)
99006 FORMAT (7X,6(d12.4))
END SUBROUTINE strmat
