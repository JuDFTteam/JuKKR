      subroutine bzirr3d(nkp,nkxyz,kpoibz,kp,recbv,bravais,
     &                   wtkp,volbz,rsymat,nsymat,isymindex,
     &                   irr,iprint)
C===========================================================================
Cinfo 
Cinfo   find irreducible BZ and create mesh in it.
Cinfo           original version Arthur Ernst, Daresbury, 27/03/98
C                                fcc lattice bug corrected
C                                modified 26/5/99, 2/06/99
C Modified on 20.01.2000 To use the 2d inversion!
C Modified Dec 2001 - Apr 2002 to deal with relativistic case, bugs
C          removed (full BZ integration gives now same results as 
C          symetrised calculation)   HE/VP, Munich
C
C===========================================================================
C  Input:
C         nkxyz : original k-mesh net in the 3 directions of the reciprocal
C                 lattice vectors (not xyz directions).
C         recbv : reciprocal unit cell 
C                 (normalized: recbv(i,j)=G(i,j)*alat/(2*pi), i=x,y,z j=1,2,3)
C                 if G is the "normal" rec. lattice vector.
C      bravais  : direct unit cell in a.u.
C      rsymat   : symmetry rotation matrices of real lattice
C      nsymat   : number of rotation matrices of real lattice
C         irr   : if .true. then irreducible BZ else take all BZ
C
C  Output:
C         nkp   : number of points in irreducible BZ
C         kp    : k-point mesh
C         wtkp  : weights for k-points
C        volbz  : volume of the BZ
C  Inside:
C         ibk   : Flag showing if the mesh point jx,jy,jz has already been 
C                 taken. Could also be a logical variable.
C==========================================================================
      implicit none
      integer MAXK1,MAXK2,MAXK3,NSYMAXD
      PARAMETER (MAXK1=350,MAXK2=350,MAXK3=350,NSYMAXD=48)
C i/o
C     integer nkp,nkxyz(3),ibk(0:MAXK1,0:MAXK2,0:MAXK3),kpoibz,nsymat
      integer nkp,nkxyz(3),kpoibz,nsymat
      integer, allocatable :: ibk(:,:,:)
      INTEGER ISYMINDEX(*),nkxyz1(3)
      double precision kp(3,*),wtkp(*),volbz,CF(3)
      double precision recbv(3,3),bravais(3,3),RSYMAT(64,3,3)
      logical irr
C local
      REAL*8 BGINV(3,3),BGMAT(3,3),BGP(3,3),BV(3)
      INTEGER NBGP(3,3,48),ISYM
      INTEGER IKTAB(3),IND2(3),IPRINT
      INTEGER I,J,JX,JY,JZ,NSYM,IWS,IWT,IS,JA,JB,JC,IX,IY,IZ,K,NK,N
      INTEGER I0,N0,INVERS2D(64),NDIM
      DOUBLE PRECISION U(3,3,48),GQ(3,3),V1
      LOGICAL LSURF
      DOUBLE PRECISION MSIGN

      DOUBLE PRECISION DDET33,DDOT
      EXTERNAL DDET33,DDOT

C
Check if we are in surface mode
C
      LSURF = .FALSE.
      IF (BRAVAIS(1,3).EQ.0.D0.AND.BRAVAIS(2,3).EQ.0.D0.AND.
     &     BRAVAIS(3,3).EQ.0.D0) LSURF = .TRUE. 

      NDIM = 3
      IF (LSURF) THEN
         NKXYZ(3) = 1
         NDIM = 2
      END IF

      IF (NKXYZ(1).GT.MAXK1.OR.NKXYZ(2).GT.MAXK2.OR.NKXYZ(3).GT.MAXK3)
     &     THEN
          WRITE(6,*) 'BZIRR3D : Increase MAXK ',(NKXYZ(I),I=1,3),
     &         MAXK1,MAXK2,MAXK3
          STOP
      END IF
      DO I=1,3
          NKXYZ1(I) = NKXYZ(I)
      END DO

c-------------------
c
c Create small unit cell for the integration in the reciprocal space (gq),
c
C                                          steps along the basis vectors
      DO I = 1, 3
         DO J = 1, 3
            GQ(I,J)=RECBV(I,J)/DBLE(NKXYZ1(J))
         ENDDO 
      ENDDO
C=======================================================================
      IF (IRR) THEN

          NSYM = NSYMAT
C
C-----------------------------------------------------------------------
C         
          IF (NSYM.GT.48) STOP 'DIM problem in bzirr3d.'

          IF ((LSURF).AND.(NSYM.GT.12)) STOP
     &         'DIM problem in bzirr3d, surf mode.'
C
C----------------------------------------------------------------------
C
          DO N=1,NSYM
              ISYM = ISYMINDEX(N)
              MSIGN = 1.D0
              DO I=1,3
                  DO J=1,3
                      U(I,J,N) = MSIGN * RSYMAT(ISYM,I,J)
                  END DO
              END DO
          END DO
C========================================================================
      ELSE
C========================================================================
         NSYM=1
         U(1,1,1)=1.D0
         U(1,2,1)=0.D0
         U(1,3,1)=0.D0
         U(2,1,1)=0.D0
         U(2,2,1)=1.D0
         U(2,3,1)=0.D0
         U(3,1,1)=0.D0
         U(3,2,1)=0.D0
         U(3,3,1)=1.D0
      ENDIF
C==========================================================================

      DO J = 1,3
          DO I = 1,3
              BGMAT(I,J) = DDOT(3,GQ(1,I),1,GQ(1,J),1)
          END DO
      END DO
C
      CALL RINVGJ(BGINV,BGMAT,3,NDIM)
C
C-----------------------------------------------------------------------
C          rotate the 3 step vectors   GQ  --  it should be possible
C          to express the rotated vectors  B(j) in terms of the old ones
C     B(i) = SUM(j) GQ(j) * n(j,i)  with integer coefficients n(j,i)
C
C==========================================================================
      DO IS = 1,NSYM
          IF ( IPRINT.GT.2 ) WRITE (*,99004) IS
C     
          DO I = 1,3
              CALL DGEMV('N',3,3,1D0,U(1,1,IS),3,GQ(1,I),1,
     &             0D0,BGP(1,I),1)
C     
              DO J = 1,3
                  BV(J) = DDOT(3,GQ(1,J),1,BGP(1,I),1)
              END DO
              CALL DGEMV('N',3,3,1D0,BGINV,3,BV,1,0D0,CF,1)
C     
              DO J = 1,3
                  IF ( ABS(NINT(CF(J))-CF(J)).GT.1D-8 ) 
     &                 WRITE (*,99006) I,J,CF(J)
                  NBGP(J,I,IS) = NINT(CF(J))
              END DO
              IF ( IPRINT.GT.2 ) WRITE (*,99005) I,(BGP(J,I),J=1,3),BV
     &             ,CF,(NBGP(J,I,IS),J=1,3)
C     
          END DO
C     
      END DO
C========================================================================

      NK=NKXYZ1(1)*NKXYZ1(2)*NKXYZ1(3)

C Fix: allocate memory for ibk array on heap to avoid stack overflow ER
      ALLOCATE(IBK(0:NKXYZ1(1),0:NKXYZ1(2),0:NKXYZ1(3)))

      DO I=0,NKXYZ1(1)
          DO J=0,NKXYZ1(2)
              DO K=0,NKXYZ1(3)
                  IBK(I,J,K)=0
              ENDDO
          ENDDO
      ENDDO

      IX=NKXYZ1(1)
      IY=NKXYZ1(2)
      IZ=NKXYZ1(3)
      NKP=0
      IWS=0

      DO JX = 0, IX-1
          IKTAB(1) = JX
          DO JY = 0, IY-1
              IKTAB(2) = JY
              DO JZ = 0, IZ-1
                  IKTAB(3) = JZ
                  
                  IF(IBK(JX,JY,JZ).EQ.0) THEN
                      NKP=NKP+1
                      IWT=0
                      
C========================================================================
                      DO ISYM = 1, NSYM
C
C               rotate k-vector  NKP  and transform into parallelepiped
C          ROT k = SUM(i) m(i) ROT gq(i)
C                = SUM(i) m(i) SUM(j) n(i,j) gq(j)
C                = SUM(j) [SUM(i) m(i) n(i,j)] gq(j)
C
                          DO J = 1,3
                              IS = 0
                              DO I = 1,3
                                  IS = IS + IKTAB(I)*NBGP(J,I,ISYM)
                              END DO
                              IS = MOD(IS,NKXYZ1(J))
                              IF ( IS.LT.0 ) IS = IS + NKXYZ1(J)
                              IND2(J) = IS
                          END DO
C     
                          JA=IND2(1)
                          JB=IND2(2)
                          JC=IND2(3)

                          IF(IBK(JA,JB,JC).EQ.0) THEN
                              IBK(JA,JB,JC)=NKP
                              IWT=IWT+1
                          ENDIF
c Mesh point in the unit cell found.
                      ENDDO
C========================================================================
                      
                      IF ( NKP.LE.KPOIBZ ) THEN
                          DO I = 1,3
                              KP(I,NKP) = 0D0
                              DO J = 1,3
                                  KP(I,NKP) = KP(I,NKP) + 
     &                                 GQ(I,J)*DBLE(IKTAB(J))
                              END DO
                          END DO
                          WTKP(NKP) = DBLE(IWT)/DBLE(NK)
                      ELSE
                          WRITE (6,*) ' No. of k-points ',NKP,
     &                                ' > KPOIBZ '
                          STOP ' <BZIRR3D> increase parameter KPOIBZ'
                      END IF
                      
                  END IF
                  IWS=IWS+IWT
                  
              ENDDO
          ENDDO
      ENDDO

      VOLBZ = 0.D0

      IF (LSURF) THEN
         V1 = DABS(RECBV(1,1)*RECBV(2,2)-RECBV(1,2)*RECBV(2,1))
      ELSE
         V1 = DDET33(RECBV) 
      ENDIF

      DO I=1,NKP
         WTKP(I) = WTKP(I)*V1/DBLE(NSYM)
         VOLBZ = VOLBZ + WTKP(I)*DBLE(NSYM)
      END DO    
      RETURN

99004 FORMAT (5X,'rotated GQ  for IROT=',I3)
99005 FORMAT (5X,I3,3F7.3,2X,3F7.3,2X,3F7.3,2X,3I3)
99006 FORMAT (5X,2I3,3F7.3)
      END
