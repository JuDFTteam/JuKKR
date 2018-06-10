SUBROUTINE bzirr3d(nkp,nkxyz,kpoibz,kp,recbv,bravais,  &
        wtkp,volbz,rsymat,nsymat,isymindex,  &
        symunitary,irr,krel,iprint)
!===========================================================================
!info
!info   find irreducible BZ and create mesh in it.
!info           original version Arthur Ernst, Daresbury, 27/03/98
!                                fcc lattice bug corrected
!                                modified 26/5/99, 2/06/99
! Modified on 20.01.2000 To use the 2d inversion!
! Modified Dec 2001 - Apr 2002 to deal with relativistic case, bugs
!          removed (full BZ integration gives now same results as
!          symetrised calculation)   HE/VP, Munich

!===========================================================================
!  Input:
!         nkxyz : original k-mesh net in the 3 directions of the reciprocal
!                 lattice vectors (not xyz directions).
!         recbv : reciprocal unit cell
!                 (normalized: recbv(i,j)=G(i,j)*alat/(2*pi), i=x,y,z j=1,2,3)
!                 if G is the "normal" rec. lattice vector.
!      bravais  : direct unit cell in a.u.
!      rsymat   : symmetry rotation matrices of real lattice
!      nsymat   : number of rotation matrices of real lattice
!         irr   : if .true. then irreducible BZ else take all BZ

!     symunitary: in case of a FM REL calculation, indicates a unitary/
!                 antiunitary symmetry operation

!  Output:
!         nkp   : number of points in irreducible BZ
!         kp    : k-point mesh
!         wtkp  : weights for k-points
!        volbz  : volume of the BZ
!  Inside:
!         ibk   : Flag showing if the mesh point jx,jy,jz has already been
!                 taken. Could also be a logical variable.
!==========================================================================
      implicit none
      integer MAXK1,MAXK2,MAXK3,NSYMAXD
      PARAMETER (MAXK1=501,MAXK2=501,MAXK3=100,NSYMAXD=48)
      !PARAMETER (MAXK1=350,MAXK2=350,MAXK3=70,NSYMAXD=48)
! i/o
      ! made local ibk array allocatable to not take too much memory
      !integer nkp,nkxyz(3),ibk(0:MAXK1,0:MAXK2,0:MAXK3),kpoibz,nsymat
      integer nkp,nkxyz(3),kpoibz,nsymat
      INTEGER ISYMINDEX(*),nkxyz1(3),krel
      double precision kp(3,*),wtkp(*),volbz,CF(3)
      double precision recbv(3,3),bravais(3,3),RSYMAT(64,3,3)
      logical irr
! local
      integer, allocatable :: ibk(:,:,:)
      REAL*8 BGINV(3,3),BGMAT(3,3),BGP(3,3),BV(3)
      INTEGER NBGP(3,3,48),ISYM
      INTEGER IKTAB(3),IND2(3),IPRINT
      INTEGER I,J,JX,JY,JZ,NSYM,IWS,IWT,IS,JA,JB,JC,IX,IY,IZ,K,NK,N
      INTEGER NDIM
      DOUBLE PRECISION U(3,3,48),GQ(3,3),V1
      LOGICAL LSURF
      LOGICAL SYMUNITARY(NSYMAXD)
      DOUBLE PRECISION MSIGN

      DOUBLE PRECISION DDET33,DDOT
      EXTERNAL DDET33,DDOT

allocate(ibk(0:maxk1,0:maxk2,0:maxk3))


!heck if we are in surface mode

lsurf = .false.
IF (bravais(1,3) == 0.d0.AND.bravais(2,3) == 0.d0.AND.  &
    bravais(3,3) == 0.d0) lsurf = .true.

ndim = 3
IF (lsurf) THEN
  nkxyz(3) = 1
  ndim = 2
END IF

IF (nkxyz(1) > maxk1.OR.nkxyz(2) > maxk2.OR.nkxyz(3) > maxk3) THEN
  WRITE(6,*) 'BZIRR3D : Increase MAXK ',(nkxyz(i),i=1,3), maxk1,maxk2,maxk3
  STOP
END IF
DO i=1,3
  nkxyz1(i) = nkxyz(i)
END DO
!-------------------

! Create small unit cell for the integration in the reciprocal space (gq),

!                                          steps along the basis vectors
DO i = 1, 3
  DO j = 1, 3
    gq(i,j)=recbv(i,j)/DBLE(nkxyz1(j))
  END DO
END DO
!=======================================================================
IF (irr) THEN
  
  nsym = nsymat
  
!-----------------------------------------------------------------------
  
  IF (nsym > 48) STOP 'DIM problem in bzirr3d.'
  
  IF ((lsurf).AND.(nsym > 12)) STOP 'DIM problem in bzirr3d, surf mode.'
  
!----------------------------------------------------------------------
  
  DO n=1,nsym
    isym = isymindex(n)
    msign = 1.d0
    IF ((krel == 1).AND.(.NOT.symunitary(n))) msign = -1.d0
    DO i=1,3
      DO j=1,3
        u(i,j,n) = msign * rsymat(isym,i,j)
      END DO
    END DO
  END DO
!========================================================================
ELSE
!========================================================================
  nsym=1
  u(1,1,1)=1.d0
  u(1,2,1)=0.d0
  u(1,3,1)=0.d0
  u(2,1,1)=0.d0
  u(2,2,1)=1.d0
  u(2,3,1)=0.d0
  u(3,1,1)=0.d0
  u(3,2,1)=0.d0
  u(3,3,1)=1.d0
END IF
!==========================================================================

DO j = 1,3
  DO i = 1,3
    bgmat(i,j) = ddot(3,gq(1,i),1,gq(1,j),1)
  END DO
END DO

CALL rinvgj(bginv,bgmat,3,ndim)

!-----------------------------------------------------------------------
!          rotate the 3 step vectors   GQ  --  it should be possible
!          to express the rotated vectors  B(j) in terms of the old ones
!     B(i) = SUM(j) GQ(j) * n(j,i)  with integer coefficients n(j,i)

!==========================================================================
DO is = 1,nsym
  IF ( iprint > 2 ) WRITE (1337,99004) is
  
  DO i = 1,3
    CALL dgemv('N',3,3,1D0,u(1,1,is),3,gq(1,i),1, 0D0,bgp(1,i),1)
    
    DO j = 1,3
      bv(j) = ddot(3,gq(1,j),1,bgp(1,i),1)
    END DO
    CALL dgemv('N',3,3,1D0,bginv,3,bv,1,0D0,cf,1)
    
    DO j = 1,ndim
      IF ( ABS(nint(cf(j))-cf(j)) > 1D-8 ) WRITE (1337,99006) i,j,cf(j)
      nbgp(j,i,is) = nint(cf(j))
    END DO
    IF ( iprint > 2 ) WRITE (1337,99005) i,  &
        (bgp(j,i),j=1,3), bv, cf, (nbgp(j,i,is),j=1,3)
    
  END DO
  
END DO
!========================================================================

nk=nkxyz1(1)*nkxyz1(2)*nkxyz1(3)

DO i=0,nkxyz1(1)
  DO j=0,nkxyz1(2)
    DO k=0,nkxyz1(3)
      ibk(i,j,k)=0
    END DO
  END DO
END DO

ix=nkxyz1(1)
iy=nkxyz1(2)
iz=nkxyz1(3)
nkp=0
iws=0

DO jx = 0, ix-1
  iktab(1) = jx
  DO jy = 0, iy-1
    iktab(2) = jy
    DO jz = 0, iz-1
      iktab(3) = jz
      
      IF(ibk(jx,jy,jz) == 0) THEN
        nkp=nkp+1
        iwt=0
        
!========================================================================
        DO isym = 1, nsym
          
!               rotate k-vector  NKP  and transform into parallelepiped
!          ROT k = SUM(i) m(i) ROT gq(i)
!                = SUM(i) m(i) SUM(j) n(i,j) gq(j)
!                = SUM(j) [SUM(i) m(i) n(i,j)] gq(j)
          
          DO j = 1,3
            is = 0
            DO i = 1,3
              is = is + iktab(i)*nbgp(j,i,isym)
            END DO
            is = MOD(is,nkxyz1(j))
            IF ( is < 0 ) is = is + nkxyz1(j)
            ind2(j) = is
          END DO
          
          ja=ind2(1)
          jb=ind2(2)
          jc=ind2(3)
          
          IF(ibk(ja,jb,jc) == 0) THEN
            ibk(ja,jb,jc)=nkp
            iwt=iwt+1
          END IF
! Mesh point in the unit cell found.
        END DO
!========================================================================
        
        IF ( nkp <= kpoibz ) THEN
          DO i = 1,3
            kp(i,nkp) = 0D0
            DO j = 1,3
              kp(i,nkp) = kp(i,nkp) + gq(i,j)*DBLE(iktab(j))
            END DO
          END DO
          wtkp(nkp) = DBLE(iwt)/DBLE(nk)
        ELSE
          WRITE (6,*) ' No. of k-points ',nkp, ' > KPOIBZ '
          STOP ' <BZIRR3D> increase parameter KPOIBZ'
        END IF
        
      END IF
      iws=iws+iwt
      
    END DO
  END DO
END DO

volbz = 0.d0

IF (lsurf) THEN
  v1 = DABS(recbv(1,1)*recbv(2,2)-recbv(1,2)*recbv(2,1))
ELSE
  v1 = ddet33(recbv)
END IF

DO i=1,nkp
  wtkp(i) = wtkp(i)*v1/DBLE(nsym)
  volbz = volbz + wtkp(i)*DBLE(nsym)
END DO

deallocate(ibk)

RETURN

99004 FORMAT (5X,'rotated GQ  for IROT=',i3)
99005 FORMAT (5X,i3,3F7.3,2X,3F7.3,2X,3F7.3,2X,3I3)
99006 FORMAT (5X,2I3,3F7.3)
END SUBROUTINE bzirr3d
