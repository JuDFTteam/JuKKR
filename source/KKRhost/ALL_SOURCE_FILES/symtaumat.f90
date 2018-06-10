SUBROUTINE symtaumat(rotname,rotmat,drot,nsym,isymindex,  &
    symunitary,nqmax,nkmmax,nq,nl,krel,iprint, nsymaxd)
!   ********************************************************************
!   *                                                                  *
!   *  Find the symmetry matrices DROT that act on t, tau, ....        *
!   *  KREL=0: for real spherical harmonics                            *
!   *  KREL=1: for relativistic represntation                          *
!   *                                                                  *
!   *  The NSYM allowed symmetry operations are indicated by ISYMINDEX *
!   *  in the table  ROTMAT. For KREL=1, SYMUNITARY=T/F indicates a    *
!   *  unitary/antiunitary symmetry operation.                         *
!   *                                                                  *
!   *  The routine determines first the Euler angles correponding      *
!   *  to a symmetry operation. Reflections are decomposed into        *
!   *  inversion + rotation for this reason.                           *
!   *                                                                  *
!   ********************************************************************

IMPLICIT NONE

! PARAMETER definitions
DOUBLE COMPLEX CI,C1,C0
PARAMETER (CI=(0.0D0,1.0D0),C1=(1.0D0,0.0D0),C0=(0.0D0,0.0D0))

! Dummy arguments
INTEGER IPRINT,KREL,NKMMAX,NL,NQ,NQMAX,NSYM,NSYMAXD
DOUBLE COMPLEX DROT(NKMMAX,NKMMAX,48)
INTEGER ISYMINDEX(NSYMAXD)
DOUBLE PRECISION ROTMAT(64,3,3)
CHARACTER (len=10) :: ROTNAME(64)
LOGICAL SYMUNITARY(48)

! Local variables
double precision  A,B,CO1,CO2,CO3,DET,FACT(0:100),PI, &
        SI1,SI2,SI3,SK,SYMEULANG(3,48),TET1,TET2,TET3
LOGICAL CHECKRMAT
DOUBLE PRECISION DBLE
double precision  DDET33
DOUBLE COMPLEX DINV(NKMMAX,NKMMAX),DTIM(NKMMAX,NKMMAX), &
           RC(NKMMAX,NKMMAX),W1(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX)
LOGICAL EQUAL
INTEGER I,I1,I2,IND0Q(NQMAX),INVFLAG(48),IQ,IREL,IRELEFF,ISYM, &
        ITOP,J,K,L,LOOP,M,N,NK,NKEFF,NKM,NLM,NOK,NS,RJ,RMJ
INTEGER NINT
DOUBLE PRECISION RMAT(3,3)
double precision  W

equal(a,b) = (DABS(a-b) < 1D-7)

WRITE (1337,99001)

pi = 4D0*ATAN(1D0)

irel = krel*3
nk = (1-krel)*nl + krel*(2*nl-1)
nkm = (1+krel)*nl**2

!-----------------------------------------------------------------------
fact(0) = 1.0D0
DO i = 1,100
  fact(i) = fact(i-1)*DBLE(i)
END DO
!-----------------------------------------------------------------------

ind0q(1) = 0
DO iq = 2,nq
  ind0q(iq) = ind0q(iq-1) + nkm
END DO

! ----------------------------------------------------------------------
!    RC  transforms from  REAL to  COMPLEX (L,M,S) - representation
!                 |LC> = sum[LR] |LR> * RC(LR,LC)
! ----------------------------------------------------------------------
IF ( krel == 0 ) THEN
  nlm = nkm
  
  CALL cinit(nkmmax*nkmmax,rc)
  
  w = 1.0D0/SQRT(2.0D0)
  
  DO l = 0,(nl-1)
    DO m = -l,l
      i = l*(l+1) + m + 1
      j = l*(l+1) - m + 1
      
      IF ( m < 0 ) THEN
        rc(i,i) = -ci*w
        rc(j,i) = w
      END IF
      IF ( m == 0 ) THEN
        rc(i,i) = c1
      END IF
      IF ( m > 0 ) THEN
        rc(i,i) = w*(-1.0D0)**m
        rc(j,i) = ci*w*(-1.0D0)**m
      END IF
    END DO
  END DO
END IF

!=======================================================================
!     The routine determines first the Euler angles correponding
!     to a symmetry operation. Reflections are decomposed into
!     inversion + rotation for this reason.
!=======================================================================

DO isym = 1,nsym
  
  DO i1 = 1,3
    DO i2 = 1,3
      rmat(i1,i2) = rotmat(isymindex(isym),i1,i2)
    END DO
  END DO
  
  det = ddet33(rmat)
  
  invflag(isym) = 0
  IF ( det < 0D0 ) THEN
    CALL dscal(9,-1.0D0,rmat,1)
    invflag(isym) = 1
  END IF
  
!----------------------------------------------------------------------
  co2 = rmat(3,3)
  tet2 = ACOS(co2)
  loop = 0
  50      CONTINUE
  IF ( loop == 1 ) tet2 = -tet2
  si2 = SIN(tet2)
  
  IF ( equal(co2,1.0D0) ) THEN
    tet1 = ACOS(rmat(1,1))
    IF ( .NOT.equal(rmat(1,2),SIN(tet1)) ) THEN
      tet1 = -tet1
      IF ( .NOT.equal(rmat(1,2),SIN(tet1)) ) WRITE (1337,*)  &
          '>>>>>>>>>>>>>>> STRANGE 1'
    END IF
    tet2 = 0D0
    tet3 = 0D0
  ELSE IF ( equal(co2,-1D0) ) THEN
    tet1 = ACOS(-rmat(1,1))
    IF ( .NOT.equal(rmat(1,2),-SIN(tet1)) ) THEN
      tet1 = -tet1
      IF ( .NOT.equal(rmat(1,2),-SIN(tet1)) ) WRITE (1337,*)  &
          '>>>>>>>>>>>>>>> STRANGE 2'
    END IF
    tet2 = pi
    tet3 = 0D0
  ELSE
    tet1 = ACOS(rmat(3,1)/si2)
    IF ( .NOT.equal(rmat(3,2),si2*SIN(tet1)) ) THEN
      tet1 = -tet1
      IF ( .NOT.equal(rmat(3,2),si2*SIN(tet1)) ) WRITE (1337,*)  &
          '>>>>>>>>>>>>>>> STRANGE 3'
    END IF
    
    tet3 = ACOS(-rmat(1,3)/si2)
    IF ( .NOT.equal(rmat(2,3),si2*SIN(tet3)) ) THEN
      tet3 = -tet3
      IF ( .NOT.equal(rmat(2,3),si2*SIN(tet3)) ) WRITE (1337,*)  &
          '>>>>>>>>>>>>>>> STRANGE 4'
    END IF
    
  END IF
  
  co1 = COS(tet1)
  si1 = SIN(tet1)
  co2 = COS(tet2)
  si2 = SIN(tet2)
  co3 = COS(tet3)
  si3 = SIN(tet3)
  
  nok = 0
  DO i1 = 1,3
    DO i2 = 1,3
      IF ( checkrmat(rmat,co1,si1,co2,si2,co3,si3,i1,i2) ) THEN
        nok = nok + 1
      ELSE IF ( loop < 1 ) THEN
        loop = loop + 1
        GO TO 50
      END IF
    END DO
  END DO
  
  symeulang(1,isym) = tet1*(180D0/pi)
  symeulang(2,isym) = tet2*(180D0/pi)
  symeulang(3,isym) = tet3*(180D0/pi)
  
  IF ( nok /= 9 ) WRITE (1337,99009) nok
  WRITE (1337,99008) isym,rotname(isymindex(isym)),invflag(isym),  &
      (symeulang(i,isym),i=1,3),symunitary(isym)
  
END DO
WRITE(1337,'(8X,57(1H-),/)')

!-----------------------------------------------------------------------
!                    initialize all rotation matrices
!-----------------------------------------------------------------------

CALL cinit(nkmmax*nkmmax*nsym,drot)

!-----------------------------------------------------------------------
!                       create rotation matrices
!-----------------------------------------------------------------------

IF ( irel <= 2 ) THEN
  ireleff = 0
  nkeff = nl
ELSE
  ireleff = 3
  nkeff = nk
END IF

DO isym = 1,nsym
  
  CALL calcrotmat(nkeff,ireleff,symeulang(1,isym),  &
      symeulang(2,isym),symeulang(3,isym), drot(1,1,isym),fact,nkmmax)
  
END DO
!-----------------------------------------------------------------------
!                     create matrix for inversion
!-----------------------------------------------------------------------
CALL cinit(nkmmax*nkmmax,dinv)

i = 0
IF ( irel > 2 ) THEN
  ns = 2
ELSE
  ns = 1
END IF
DO l = 0,(nl-1)
  DO m = 1,ns*(2*l+1)
    i = i + 1
    dinv(i,i) = (-1.0D0)**l
  END DO
END DO
itop = i

!-----------------------------------------------------------------------
!                         include inversion
!-----------------------------------------------------------------------
DO isym = 1,nsym
  IF ( invflag(isym) /= 0 ) THEN
    
    CALL zgemm('N','N',nkm,nkm,nkm,c1,drot(1,1,isym),nkmmax,  &
        dinv,nkmmax,c0,w2,nkmmax)
    
    DO j = 1,nkm
      CALL zcopy(nkm,w2(1,j),1,drot(1,j,isym),1)
    END DO
  END IF
END DO

!-----------------------------------------------------------------------
!            add second spin-diagonal block for  IREL=2
!            spin off-diagonal blocks have been initialized before
!-----------------------------------------------------------------------
IF ( irel == 2 ) THEN
  nlm = nkm/2
  IF ( itop /= nlm ) CALL errortrap('SYMTAUMAT',11,1)
  DO isym = 1,nsym
    
    DO j = 1,nlm
      CALL zcopy(nlm,drot(1,j,isym),1,drot(nlm+1,nlm+j,isym),1)
    END DO
  END DO
END IF
!-----------------------------------------------------------------------
!            transform to real spherical representation for  KREL=0
!-----------------------------------------------------------------------
n = nkm
m = nkmmax
IF ( krel == 0 ) THEN
  DO isym = 1,nsym
    CALL zgemm('N','N',n,n,n,c1,rc,m,drot(1,1,isym),m,c0,w1,m)
    CALL zgemm('N','C',n,n,n,c1,w1,m,rc,m,c0,drot(1,1,isym),m)
  END DO
END IF
!-----------------------------------------------------------------------
!                     create matrix for time reversal
!-----------------------------------------------------------------------
IF ( irel > 1 ) THEN
  
  CALL cinit(nkmmax*nkmmax,dtim)
  
  i = 0
  DO k = 1,nk
    l = k/2
    IF ( l*2 == k ) THEN
      sk = -1D0
    ELSE
      sk = +1D0
    END IF
    rj = nint(l + sk*0.5D0)
    DO rmj = -rj, + rj, 1
      i1 = nint(2*l*(rj+0.5D0)+rj+rmj+1)
      i2 = nint(2*l*(rj+0.5D0)+rj-rmj+1)
      dtim(i1,i2) = sk*(-1)**nint(rmj+0.5D0)
    END DO
  END DO
  IF ( iprint > 0 ) THEN
    CALL cmatstr('Inversion     MATRIX',20,dinv,nkm,nkmmax,3,3, 0,1D-8,6)
    CALL cmatstr('Time reversal MATRIX',20,dtim,nkm,nkmmax,3,3, 0,1D-8,6)
  END IF
  
END IF
!=======================================================================
!            set up of transformation matrices completed
!=======================================================================

!=======================================================================
!   include time reversal operation for anti-unitary transformations
!=======================================================================
DO isym = 1,nsym
  IF ( .NOT.symunitary(isym) ) THEN
    IF ( irel == 2 ) CALL errortrap('SYMTAUMAT',14,1)
    
    CALL zgemm('N','N',nkm,nkm,nkm,c1,drot(1,1,isym),nkmmax,  &
        dtim,nkmmax,c0,w2,nkmmax)
    DO j = 1,nkm
      CALL zcopy(nkm,w2(1,j),1,drot(1,j,isym),1)
    END DO
  END IF
END DO

!-----------------------------------------------------------------------
! for testing

!ccc      write (6,*) ' NUMBER OF SYMMETRIES : ', NSYM
!ccc
!ccc      do isym = 1,nsym
!ccc         write(6,*) ' ISYM = ',isym
!ccc         call cmatstr('DROT',4,drot(1,1,isym),nkm,nkmmax,krel*3,krel*3,
!ccc     &        0,1d-12,6)
!ccc         write(6,*)
!ccc      end do

!-----------------------------------------------------------------------

IF ( iprint == 0 ) RETURN

!=======================================================================
!       find the structure of the site-diagonal TAU - matrices  TAUQ
!=======================================================================

CALL taustruct(drot,nsym,symunitary,nkm,nq,nqmax,nkmmax,iprint, irel)

RETURN
99001 FORMAT (5X,'<SYMTAUMAT> : rotation matrices acting on t/G/tau',//,  &
    8X,57(1H-),/,8X,  &
    'ISYM            INV          Euler angles      Unitarity',/, 8X,57(1H-))
99008 FORMAT (8X,i2,3X,a,i3,3F10.5,3X,l1)
99009 FORMAT (50('>'),' trouble in <SYMTAUMAT>',i3,f10.5)
END SUBROUTINE symtaumat
