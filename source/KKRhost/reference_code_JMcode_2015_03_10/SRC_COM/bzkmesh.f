      SUBROUTINE BZKMESH(NBXIN,NBYIN,NBZIN,MAXMESH,LIRR,BRAVAIS,RECBV,
     &                   NSYMAT,RSYMAT,ISYMINDEX,SYMUNITARY,
     &                   IELAST,EZ,KMESH,IPRINT,KREL,KPOIBZ,MAXMSHD)
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER MAXMESH,NBXIN,NBYIN,NBZIN,NSYMAT,IPRINT,
     &        KREL,KPOIBZ,IELAST,MAXMSHD
      LOGICAL LIRR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION BRAVAIS(3,3),RECBV(3,3)
      DOUBLE PRECISION RSYMAT(64,3,3)
      INTEGER ISYMINDEX(*),KMESH(*)
      DOUBLE COMPLEX EZ(*)
C     .. unitary/antiunitary symmetry flag
      LOGICAL SYMUNITARY(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION VOLBZ
      INTEGER I,ID,KS,L,N,NBX,NBY,NBZ,NOFKS
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION BZKP(3,KPOIBZ),VOLCUB(KPOIBZ)
      INTEGER NXYZ(3)
C     ..
C     .. External Functions ..
      LOGICAL TEST
      EXTERNAL TEST
C     ..
C     .. External Subroutines ..
      EXTERNAL BZIRR3D
C ---------------------------------------------------------------------
C
C --> set number of different K-meshes 
C
      MAXMESH = 1
      IF ( TEST('fix mesh') ) THEN
         DO I = 1,IELAST
            KMESH(I) = 1
         END DO
      ELSE
         DO I = 1,IELAST
            N = INT(1.001D0 + 
     &              LOG(DIMAG(EZ(I))/DIMAG(EZ(IELAST)))/LOG(2.0D0))
            KMESH(I) = N
            MAXMESH = MAX(MAXMESH,N)
            IF ( KMESH(I).LT.1 ) KMESH(I) = 1
         END DO
         KMESH(1) = MAXMESH
      END IF
C
      IF ( TEST('fix4mesh') ) THEN
         DO I = 1,IELAST
            KMESH(I) = MAXMESH
         END DO
      END IF
C 
      IF (MAXMESH.GT.MAXMSHD) THEN
        WRITE (6,FMT='(5X,A,I2)') 
     &        'Dimension ERROR: Please increase MAXMSHD to ',MAXMESH
        WRITE (6,FMT='(22X,A,/)') 
     &        'in the programs < main0 > and < main1b >'
        STOP '        < BZKMESH >'
      END IF
C ---------------------------------------------------------------------
      NBX = NBXIN
      NBY = NBYIN
      NBZ = NBZIN

      OPEN (52,FILE='kpoints',FORM='formatted')
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,99000)
      WRITE (6,99001) MAXMESH,NSYMAT
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      DO L = 1,MAXMESH
         IF (L.GT.1) THEN
            NBX = NBX/1.4
            NBY = NBY/1.4
            NBZ = NBZ/1.4
         END IF
         IF (NBX.LT.1) NBX = 1
         IF (NBY.LT.1) NBY = 1
         IF (NBZ.LT.1) NBZ = 1
         NXYZ(1) = NBX
         NXYZ(2) = NBY
         NXYZ(3) = NBZ
C     
         CALL BZIRR3D(NOFKS,NXYZ,KPOIBZ,BZKP,RECBV,BRAVAIS,VOLCUB,VOLBZ,
     +                RSYMAT,NSYMAT,ISYMINDEX,SYMUNITARY,
     +                LIRR,KREL,IPRINT)
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
         WRITE(6,99002) L,NOFKS,(NXYZ(I),I=1,3),VOLBZ
         IF ( L.EQ.MAXMESH ) WRITE(6,99003)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
         WRITE (52,FMT='(I8,F15.10,/,(3F12.8,D20.10))') 
     +        NOFKS,VOLBZ,((BZKP(ID,I),ID=1,3),VOLCUB(I),I=1,NOFKS)
C ---------------------------------------------------------------------
C
C -->  output of k-mesh
C
        IF (TEST('k-net   ')) THEN
           DO KS = 1,NOFKS
              WRITE (6,FMT=9000) (BZKP(I,KS),I=1,3),VOLCUB(KS)
           END DO
        END IF
C ---------------------------------------------------------------------
      END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      CLOSE (52)
 9000 FORMAT (3F12.5,F15.8)
99000 FORMAT (5X,'< BZKMESH > : creating k-mesh,',
     &        ' write to file kpoints',/)
99001 FORMAT (8X,'number of different k-meshes :',I2,/,
     &        8X,'the direct lattice',I3,' symmetries will be used',//,
     &     8X,35(1H-),/,8X,'k-mesh NofKs N kx N ky N kz vol BZ',/,
     &     8X,35(1H-))
99002    FORMAT (8X,2I6,3I5,F8.4)
99003    FORMAT (8X,35(1H-),/)
      END