! 02.08.95 *************************************************************
SUBROUTINE rrgen (bv1,lsurf,rr,nr,nrd)
! **********************************************************************
! *                                                                    *
! * generates a number of real space vectors to construct the          *
! * clusters representing the local surrounding of the atoms in        *
! * routine CLSGEN99                                                   *
! *                                                                    *
! **********************************************************************
      IMPLICIT NONE
!..
!.. Scalar arguments ..
      LOGICAL LSURF
      INTEGER NR,NRD
!    ..
!    .. Array arguments ..
      DOUBLE PRECISION BV1(3,3),RR(3,0:NRD)
!    ..
!    .. Local scalars ..
      DOUBLE PRECISION EPSSHL,R,R1,R2,R3,RMAX,RR2,RS
      INTEGER I,J,K,N1,N2,N3,POS,IPRINT
      INTEGER NINT
      DOUBLE PRECISION DBLE
!..
!.. Local arrays
      DOUBLE PRECISION RABS(NRD),RR1(3,NRD), &
                       V(3),VX(3),VY(3),VZ(3), &
                       VX0(3),VY0(3),VZ0(3)
      INTEGER IND(NRD)
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,MIN,SQRT,NINT
!..
!.. External Subroutines ..
      EXTERNAL DSORT,SCALPR,VADD,VEQ
!..
!.. Data Statements ..
      DATA  EPSSHL /1.0D-5/
!     ..................................................................
WRITE (1337,'(5X,A,/)') '< RRGEN > : generation of real space mesh RR(NR)'

iprint = 0

CALL scalpr(bv1(1,1),bv1(1,1),r1)
CALL scalpr(bv1(1,2),bv1(1,2),r2)
CALL scalpr(bv1(1,3),bv1(1,3),r3)
rmax = 5.d0

r1 = SQRT(r1)
r2 = SQRT(r2)
r3 = SQRT(r3)
r = 1.5D0*rmax + SQRT(r1*r1+r2*r2+r3*r3) + epsshl
rs = r*r
n1 = nint(r/r1)
n2 = nint(r/r2)
IF ( .NOT.lsurf ) n3 = nint(r/r3)

n1 = MIN(12,n1)
n2 = MIN(12,n2)
IF ( .NOT.lsurf ) n3 = MIN(12,n3)

n1 = MAX(2,n1)
n2 = MAX(2,n2)
IF ( .NOT.lsurf ) n3 = MAX(2,n3)

IF ( lsurf ) n3 = 0

WRITE (1337,99001) r
WRITE (1337,99002) rs
IF ( lsurf ) THEN
  WRITE (1337,99003) n1,n2
ELSE
  WRITE (1337,99004) n1,n2,n3
END IF

nr = 0
rr(1,0) = 0.0D0
rr(2,0) = 0.0D0
rr(3,0) = 0.0D0

CALL vmul(bv1(1,1),DBLE(-n1-1),vx0(1))
CALL vmul(bv1(1,2),DBLE(-n2-1),vy0(1))
CALL vmul(bv1(1,3),DBLE(-n3-1),vz0(1))
CALL veq(vx0,vx)
! **********************************************************************
DO i = -n1,n1
  CALL vadd(vx,bv1(1,1),vx)
  CALL veq(vy0,vy)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO j = -n2,n2
    CALL vadd(vy,bv1(1,2),vy)
    CALL veq(vz0,vz)
! ----------------------------------------------------------------------
    DO k = -n3,n3
      CALL vadd(vz,bv1(1,3),vz)
      CALL vadd(vx,vy,v)
      CALL vadd(v,vz,v)
      CALL scalpr(v,v,rr2)
      
      IF ( ((rr2 <= rs) .OR. (ABS(i)+ABS(j)+ABS(k) <= 6))  &
            .AND. (rr2 > epsshl) ) THEN
        nr = nr + 1
        
        IF ( nr > nrd ) THEN
          WRITE (6,*) 'Dimension ERROR. Please, change the ',  &
              'parameter NRD in inc.p to ',nr, nrd
          STOP
        END IF
        
        rr1(1,nr) = v(1)
        rr1(2,nr) = v(2)
        rr1(3,nr) = v(3)
        rabs(nr) = SQRT(rr2)
      END IF
    END DO
! ----------------------------------------------------------------------
  END DO
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END DO
! **********************************************************************

WRITE (1337,99005) nr+1

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
IF ( iprint > 0 ) THEN
  WRITE (1337,99006)
  WRITE (1337,99008) 0,0.0,0.0,0.0,0.0
END IF
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

CALL dsort(rabs,ind,nr,pos)
DO i = 1,nr
  pos = ind(i)
  rr(1,i) = rr1(1,pos)
  rr(2,i) = rr1(2,pos)
  rr(3,i) = rr1(3,pos)
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  IF ( iprint > 0 ) WRITE (1337,99008) i,rr(1,i),rr(2,i), rr(3,i),rabs(pos)
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
END DO

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
IF ( iprint > 0 )  WRITE (1337,99007)
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

99001 FORMAT (10X,'Radius R        : ',f15.6,' (ALAT    units)')
99002 FORMAT (10X,'       R**2     : ',f15.6,' (ALAT**2 units)')
99003 FORMAT (10X,'mesh divisions  : ',5X,2I5)
99004 FORMAT (10X,'mesh divisions  : ',3I5)
99005 FORMAT (10X,'vectors created : ',i15)
99006 FORMAT (/,10X,60('+'),/,18X,  &
    'generated real-space mesh-points (ALAT units)',/, 10X,60('+'),/,13X,  &
    'index      x           y           z          distance  ' ,/,10X,60('-'))
99007 FORMAT (10X,60('+'))
99008 FORMAT (10X,i6,3F12.3,f15.4)
END SUBROUTINE rrgen
