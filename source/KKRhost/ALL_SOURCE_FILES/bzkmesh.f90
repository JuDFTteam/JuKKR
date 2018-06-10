SUBROUTINE bzkmesh(nbxin,nbyin,nbzin,maxmesh,lirr,bravais,recbv,  &
        nsymat,rsymat,isymindex,symunitary,  &
        ielast,ez,kmesh,iprint,krel,kpoibz,maxmshd)
use mod_types, only: t_inc
use mod_wunfiles, only: t_params
use mod_rhoqtools, only: rhoq_write_kmesh
IMPLICIT NONE
!..
!.. Scalar Arguments ..
      INTEGER MAXMESH,NBXIN,NBYIN,NBZIN,NSYMAT,IPRINT, &
              KREL,KPOIBZ,IELAST,MAXMSHD
      LOGICAL LIRR
!..
!.. Array Arguments ..
      DOUBLE PRECISION BRAVAIS(3,3),RECBV(3,3)
      DOUBLE PRECISION RSYMAT(64,3,3)
      INTEGER ISYMINDEX(*),KMESH(*)
      DOUBLE COMPLEX EZ(*)
!.. unitary/antiunitary symmetry flag
      LOGICAL SYMUNITARY(*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION VOLBZ
      INTEGER I,ID,KS,L,N,NBX,NBY,NBZ,NOFKS
!..
!.. Local Arrays ..
      DOUBLE PRECISION BZKP(3,KPOIBZ),VOLCUB(KPOIBZ)
      INTEGER NXYZ(3)
!..
!.. External Functions ..
      LOGICAL TEST
      EXTERNAL TEST
!..
!.. External Subroutines ..
      EXTERNAL BZIRR3D
! ---------------------------------------------------------------------

! --> set number of different K-meshes

maxmesh = 1
IF ( test('fix mesh') ) THEN
  DO i = 1,ielast
    kmesh(i) = 1
  END DO
ELSE
  DO i = 1,ielast
    IF (DIMAG(ez(ielast)) /= 0) THEN
      n = INT(1.001D0 + LOG(DIMAG(ez(i))/DIMAG(ez(ielast)))/LOG(2.0D0))
    ELSE
      n = 1
    END IF
    kmesh(i) = n
    maxmesh = MAX(maxmesh,n)
    IF ( kmesh(i) < 1 ) kmesh(i) = 1
  END DO
  kmesh(1) = maxmesh
END IF

IF ( test('fix4mesh') ) THEN
  DO i = 1,ielast
    kmesh(i) = maxmesh
  END DO
END IF

IF (maxmesh > maxmshd) THEN
  WRITE (6,FMT='(5X,A,I2)')  &
      'Dimension ERROR: Please increase MAXMSHD to ',maxmesh
  WRITE (6,FMT='(22X,A,/)') 'in the programs < main0 > and < main1b >'
  STOP '        < BZKMESH >'
END IF
! ---------------------------------------------------------------------
nbx = nbxin
nby = nbyin
nbz = nbzin

IF(test('kptsfile')) OPEN (52,FILE='kpoints',FORM='formatted')
!       OPEN (52,FILE='kpoints',FORM='formatted')


! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,99000)
WRITE (1337,99001) maxmesh,nsymat
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

!save maxmesh and allocate kmesh for later use in t_inc and t_params
t_inc%nkmesh = maxmesh
t_params%kpoibz = kpoibz
t_params%maxmesh = maxmesh
allocate(t_inc%kmesh(maxmesh))
allocate(t_params%bzkp(3,kpoibz,maxmesh), t_params%volcub(kpoibz,maxmesh),  &
    t_params%volbz(maxmesh),t_params%nofks(maxmesh))
! needed for wavefunction saving
allocate(t_inc%kmesh_ie(ielast))
t_inc%kmesh_ie = kmesh(1:ielast)
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
DO l = 1,maxmesh
  IF (l > 1) THEN
    nbx = nbx/1.4
    nby = nby/1.4
    nbz = nbz/1.4
  END IF
  IF (nbx < 1) nbx = 1
  IF (nby < 1) nby = 1
  IF (nbz < 1) nbz = 1
  nxyz(1) = nbx
  nxyz(2) = nby
  nxyz(3) = nbz
  
  CALL bzirr3d(nofks,nxyz,kpoibz,bzkp,recbv,bravais,volcub,volbz,  &
      rsymat,nsymat,isymindex,symunitary, lirr,krel,iprint)
  
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  WRITE(1337,99002) l,nofks,(nxyz(i),i=1,3),volbz
  IF ( l == maxmesh ) WRITE(1337,99003)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  
  IF(test('kptsfile')) WRITE(52,FMT='(I8,F15.10,/,(3F12.8,D20.10))')  &
      nofks,volbz,((bzkp(id,i),id=1,3),volcub(i),i=1,nofks)
  IF( test('rhoqtest') .AND. (l==1) ) THEN
    CALL rhoq_write_kmesh(nofks, nxyz, volbz, bzkp, volcub, recbv, bravais)
  END IF
!       WRITE(52,FMT='(I8,F15.10,/,(3F12.8,D20.10))')
!      +        NOFKS,VOLBZ,((BZKP(ID,I),ID=1,3),VOLCUB(I),I=1,NOFKS)
  t_params%nofks(l) = nofks
  t_params%volbz(l) = volbz
  DO i=1,nofks
    DO id=1,3
      t_params%bzkp(id,i,l) = bzkp(id,i)
    END DO
    t_params%volcub(i,l) = volcub(i)
  END DO
! save nofks for this mesh in t_inc
  t_inc%kmesh(l) = nofks
! ---------------------------------------------------------------------
  
! -->  output of k-mesh
  
  IF (test('k-net   ')) THEN
    DO ks = 1,nofks
      WRITE (1337,FMT=9000) (bzkp(i,ks),i=1,3),volcub(ks)
    END DO
  END IF
! ---------------------------------------------------------------------
END DO
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
IF(test('kptsfile')) CLOSE (52)
!       CLOSE (52)
9000 FORMAT (3F12.5,f15.8)
99000 FORMAT (5X,'< BZKMESH > : creating k-mesh,', ' write to file kpoints',/)
99001 FORMAT (8X,'number of different k-meshes :',i2,/,  &
    8X,'the direct lattice',i3,' symmetries will be used',//,  &
    8X,35(1H-),/,8X,'k-mesh NofKs N kx N ky N kz vol BZ',/, 8X,35(1H-))
99002    FORMAT (8X,2I6,3I5,f8.4)
99003    FORMAT (8X,35(1H-),/)
END SUBROUTINE bzkmesh
