      SUBROUTINE BZKMESH(NBXIN,NBYIN,NBZIN,MAXMESH,LIRR,BRAVAIS,RECBV,
     &                   NSYMAT,RSYMAT,ISYMINDEX,
     &                   IELAST,EZ,KMESH,IPRINT,MAXMSHD,
C                        new after inc.p remove
     &                   IEMXD, KPOIBZ, EKMD, IGUESSD)
C     .. called by ..
C     BZKINT0
C     ..
      IMPLICIT NONE
C     .. Parameters ..
C      INCLUDE 'inc.p'
C     ..

      INTEGER IEMXD
      INTEGER KPOIBZ
      INTEGER EKMD
      INTEGER IGUESSD

C     .. Scalar Arguments ..
      INTEGER MAXMESH,NBXIN,NBYIN,NBZIN,NSYMAT,IPRINT,IELAST,MAXMSHD
      LOGICAL LIRR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION BRAVAIS(3,3),RECBV(3,3)
      DOUBLE PRECISION RSYMAT(64,3,3)
      INTEGER ISYMINDEX(*),KMESH(IEMXD)
      DOUBLE COMPLEX EZ(IEMXD)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION VOLBZ,NEWVOLBZ
      INTEGER          I,ID,KS,L,N,NBX,NBY,NBZ,NOFKS,EKMIN
      LOGICAL          NEWKP
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION BZKP(3,KPOIBZ),VOLCUB(KPOIBZ),
     +                 NEWBZKP(3,KPOIBZ,MAXMSHD),
     +                 NEWVOLCUB(KPOIBZ,MAXMSHD)
      INTEGER          NXYZ(3),NOFKS0(MAXMSHD),NEWNOFKS(MAXMSHD)
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
      IF (TEST('fix mesh')) THEN
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
         KMESH(2) = MAXMESH
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
      WRITE (6,'(79(1H=))')
      WRITE (6,99000)
      WRITE (6,'(79(1H=))')
      WRITE (6,*)
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

         IF (NBX.EQ.0) NBX = 1
         IF (NBY.EQ.0) NBY = 1
         IF (NBZ.EQ.0) NBZ = 1

         NXYZ(1) = NBX
         NXYZ(2) = NBY
         NXYZ(3) = NBZ
C
         CALL BZIRR3D(NOFKS,NXYZ,KPOIBZ,BZKP,RECBV,BRAVAIS,VOLCUB,VOLBZ,
     +                RSYMAT,NSYMAT,ISYMINDEX,
     +                LIRR,IPRINT)
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
C
        NOFKS0(L) = NOFKS
C ---------------------------------------------------------------------
      END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      CLOSE (52)
C
C
C check dimensions of EKMD precond. array in case of IGUESSD=1
C
      IF (IGUESSD.EQ.1) THEN
C
        EKMIN = 0
        DO I = 1,IELAST
            EKMIN = EKMIN + NOFKS0(KMESH(I))
        ENDDO
C
        WRITE (6,'(79(1H=))')
        WRITE (6,'(12X,A)') 
     &       'BZKMESH: checking dimensions of precond. arrays ...'
C
        IF (EKMIN.GT.EKMD) THEN
          WRITE (6,*)
     &    ' ERROR: Dimension EKMD in inc.p too small',
     &    EKMIN, EKMD
          STOP '< BZKMESH >'
        ELSE
          WRITE(6,*) '           EKMIN=',EKMIN,'  EKMD=',EKMD
        ENDIF
C
        INQUIRE(FILE='new.kpoints',EXIST=NEWKP)
C
        IF (NEWKP) THEN
C
          OPEN(53,FILE='new.kpoints',FORM='formatted')
          REWIND(53)
          DO L = 1,MAXMESH
            READ(53,FMT='(I8,f15.10)') NEWNOFKS(L),NEWVOLBZ
            READ(53,FMT=*) (NEWBZKP(ID,1,L),ID=1,3),NEWVOLCUB(1,L)
            DO I=2,NEWNOFKS(L)
              READ(53,FMT=*) (NEWBZKP(ID,I,L),ID=1,3),NEWVOLCUB(I,L)
            END DO
          END DO
          CLOSE(53)
C
          EKMIN = 0
          DO I = 1,IELAST
              EKMIN = EKMIN + NEWNOFKS(KMESH(I))
          ENDDO
C
          IF (EKMIN.GT.EKMD) THEN
            WRITE (6,*)
            WRITE (6,*)
     &      '           WARNING: Dimension EKMD ',
     &      ' too small for use of new.kpoints'
          ENDIF
            WRITE(6,*) '      new: EKMIN=',EKMIN,'  EKMD=',EKMD
C
        ENDIF
C
        WRITE (6,'(79(1H=))')
        WRITE (6,*)
C
      ENDIF
C
C
C
 9000 FORMAT (3F12.5,F15.8)
99000 FORMAT (12X,' BZKMESH : creating k-mesh,',
     &        ' write to file kpoints')
99001 FORMAT (8X,'number of different k-meshes :',I2,/,
     &        8X,'the direct lattice',I3,' symmetries will be used',//,
     &     8X,35(1H-),/,8X,'k-mesh NofKs N kx N ky N kz vol BZ',/,
     &     8X,35(1H-))
99002    FORMAT (8X,2I6,3I5,F8.4)
99003    FORMAT (8X,35(1H-),/)
      END
