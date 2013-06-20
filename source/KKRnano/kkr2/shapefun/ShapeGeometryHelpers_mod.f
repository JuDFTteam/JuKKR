      MODULE ShapeGeometryHelpers_mod
      IMPLICIT NONE

      PUBLIC :: POLCHK
      PUBLIC :: PERP

      CONTAINS
C---------------------------------------------------------------------
C>    THIS SUBROUTINE READS THE COORDINATES OF THE VERTICES OF EACH
C>    (POLYGON)  FACE OF  A CONVEX POLYHEDRON AND  CHECKS  IF THESE
C>    VERTICES ARRANGED  CONSECUTIVELY DEFINE A  POLYGON. THEN  THE
C>    SUBROUTINE  DETERMINES  THE  VERTICES  AND  THE  EDGES OF THE
C>    POLYHEDRON AND CHECKS IF  THE  NUMBER  OF  VERTICES  PLUS THE
C>    NUMBER OF FACES EQUALS THE NUMBER OF EDGES PLUS 2.
C>
C>    The angular sum of the polygons is checked. It has to be (n-2)*PI
C>    n is the number of vertices.
C
C     ----------------------------------------------------------------
C=====================================================================
      SUBROUTINE POLCHK(NFACE,NVERTICES,XVERT,YVERT,ZVERT,TOLVDIST)
      IMPLICIT NONE
C
C     .. PARAMETER STATEMENTS ..
C
C     INTEGER NEDGED
C     PARAMETER (NEDGED=NVRTD+NFACED-2)
c
c     ...Arrays ......
c
c     INTEGER   NVERTICES(NFACED)
c     REAL*8    XVERT(NVERTD,NFACED),YVERT(NVERTD,NFACED),
c    &          ZVERT(NVERTD,NFACED)

      INTEGER   NVERTICES(:)
      REAL*8    XVERT(:,:),YVERT(:,:),
     &          ZVERT(:,:)
c     ...Scalars....
      REAL*8 TOLVDIST
C
C     .. LOCAL SCALARS ..
C
      INTEGER   IVERT,INEW,IVERTP,IVERTM,IVRT,IEDGE,NVRT,NEDGE
      INTEGER   IFACE,NFACE,NVERT
      REAL*8    ARG,A1,A2,DOWN,UP,FISUM,T
      REAL*8    VRTX,VRTY,VRTZ,VRTPX,VRTPY,VRTPZ,VRTMX,VRTMY,VRTMZ
      REAL*8    PI314
C
C     .. LOCAL ARRAYS ..
C
c     REAL*8    V1(3,NEDGED),V2(3,NEDGED),V(3,NVERTD),VRT(3,NVRTD)
      REAL*8, allocatable ::    V1(:,:),V2(:,:),V(:,:),VRT(:,:)
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC ABS,ACOS,SIGN,SQRT
C     ----------------------------------------------------------------
c      READ(7,100) NFACE,LDUM,KDUM,DDUM

      INTEGER NFACED, NEDGED, NVERTD, NVRTD

      NFACED = size(NVERTICES)
      NVERTD = size(XVERT,1)
      NVRTD  = NFACED*NVERTD
      NEDGED = NVRTD+NFACED-2

      allocate(V1(3,NEDGED),V2(3,NEDGED),V(3,NVERTD),VRT(3,NVRTD))

      PI314 = 4.D0*DATAN(1.D0)

      NVRT=0
      NEDGE=0
      DO 10 IFACE=1,NFACE
c      READ(7,101) DUM1,DUM2,DUM3,DUM4,NVERT
      NVERT = NVERTICES(IFACE)
      FISUM=(NVERT-2)*PI314
c     !write(6,*) 'starting ',fisum
      DO 25 IVERT=1,NVERT
c      READ(7,102) (V(I,IVERT),I=1,3)
        V(1,IVERT) = XVERT(IVERT,IFACE)
        V(2,IVERT) = YVERT(IVERT,IFACE)
        V(3,IVERT) = ZVERT(IVERT,IFACE)
   25 CONTINUE
C
C------> T R E A T M E N T   O F   V E R T I C E S
C
      DO 2 IVERT =1,NVERT
      VRTX=V(1,IVERT)
      VRTY=V(2,IVERT)
      VRTZ=V(3,IVERT)
      INEW=1                          ! Save all different vertices
      DO 13 IVRT=1,NVRT
      T=(VRTX-VRT(1,IVRT))**2+(VRTY-VRT(2,IVRT))**2
     & +(VRTZ-VRT(3,IVRT))**2
      IF(T.LT.TOLVDIST) INEW=0
   13 CONTINUE
      IF(INEW.EQ.1)                  THEN
      NVRT=NVRT+1
      IF(NVRT.GT.NVRTD) STOP 'INCREASE NVRTD'
      VRT(1,NVRT)=V(1,IVERT)
      VRT(2,NVRT)=V(2,IVERT)
      VRT(3,NVRT)=V(3,IVERT)
                                     END IF
      IVERTP=IVERT+1
      IF(IVERT.EQ.NVERT) IVERTP=1
      VRTPX=V(1,IVERTP)
      VRTPY=V(2,IVERTP)
      VRTPZ=V(3,IVERTP)
      IVERTM=IVERT-1
      IF(IVERT.EQ.1) IVERTM=NVERT
      VRTMX=V(1,IVERTM)
      VRTMY=V(2,IVERTM)               ! Check if the  consecutive
      VRTMZ=V(3,IVERTM)               ! vertices define a polygon
      A1=SQRT((VRTPX-VRTX)**2+(VRTPY-VRTY)**2+(VRTPZ-VRTZ)**2)
      A2=SQRT((VRTMX-VRTX)**2+(VRTMY-VRTY)**2+(VRTMZ-VRTZ)**2)
      DOWN=A1*A2
      UP=(VRTPX-VRTX)*(VRTMX-VRTX)+(VRTPY-VRTY)*(VRTMY-VRTY)+
     &   (VRTPZ-VRTZ)*(VRTMZ-VRTZ)
      ! write(6,*)
      ! write(6,*) VRTX,VRTY,VRTZ
      ! write(6,*) VRTMX,VRTMY,VRTMZ
      ! write(6,*) VRTPX,VRTPY,VRTPZ
      ! write(6,*) 'fisum ',ivert,a1,a2,up,arg,ACOS(ARG)
      IF(DOWN.GE.TOLVDIST)   THEN
        ARG=UP/DOWN
        IF(ABS(ARG).GE.1.D0) ARG=SIGN(1.D0,ARG)
        FISUM=FISUM-ACOS(ARG)
      ELSE
        WRITE(6,*) 'DOWN',DOWN
        STOP 'IDENTICAL CONSECUTIVE VERTICES'
      ENDIF
C
C------> T R E A T M E N T   O F   E D G E S
C
      INEW=1                          ! Save all different edges
      DO 14 IEDGE=1,NEDGE
      T=(VRTX-V1(1,IEDGE))**2+(VRTY-V1(2,IEDGE))**2
     & +(VRTZ-V1(3,IEDGE))**2
      IF(T.LT.TOLVDIST)         THEN
      T=(VRTPX-V2(1,IEDGE))**2+(VRTPY-V2(2,IEDGE))**2
     & +(VRTPZ-V2(3,IEDGE))**2
      IF(T.LT.TOLVDIST) INEW=0
                             ELSE
      T=(VRTX-V2(1,IEDGE))**2+(VRTY-V2(2,IEDGE))**2
     & +(VRTZ-V2(3,IEDGE))**2
      IF(T.LT.TOLVDIST) THEN
      T=(VRTPX-V1(1,IEDGE))**2+(VRTPY-V1(2,IEDGE))**2
     & +(VRTPZ-V1(3,IEDGE))**2
      IF(T.LT.TOLVDIST) INEW=0
                     END IF
                             END IF
   14 CONTINUE
      IF(INEW.EQ.1)                THEN
      NEDGE=NEDGE+1
      IF(NEDGE.GT.NEDGED) STOP 'INSUFFICIENT NEDGED'
      V1(1,NEDGE)=V(1,IVERT )
      V1(2,NEDGE)=V(2,IVERT )
      V1(3,NEDGE)=V(3,IVERT )
      V2(1,NEDGE)=V(1,IVERTP)
      V2(2,NEDGE)=V(2,IVERTP)
      V2(3,NEDGE)=V(3,IVERTP)
                                   END IF
    2 CONTINUE
      IF(FISUM.GT.1.D-6) THEN
      write(6,*) 'fisum =',fisum
      STOP 'NOT CONSECUTIVE VERTICES OF A POLYGON'
      END IF
   10 CONTINUE
      IF((NVRT+NFACE).NE.(NEDGE+2)) THEN
       WRITE(6,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
       WRITE(6,*) '    SERIOUS WARNING FROM SHAPE      '
       WRITE(6,*) '   >>  STOP ILLEGAL POLYHEDRON      '
       WRITE(6,*) 'NVRT=',NVRT,' ; NFACE=',NFACE,' ; NEDGE=',NEDGE
       STOP
c changed 6.10.2000
c       IF((NVRT+NFACE).NE.(NEDGE+2)) STOP 'ILLEGAL POLYHEDRON'
      END IF
      RETURN
  100 FORMAT(3I5,F10.5)
  101 FORMAT(4F16.8,2I5)
  102 FORMAT(4F16.8)
      END SUBROUTINE

C-----------------------------------------------------------------------
C>    GIVEN  TWO  DISTINCT  POINTS   R1 , R2, THIS  ROUTINE CALCULATES
C>    THE COORDINATES  OF THE FOOT OF  THE  PERPENDICULAR FROM A POINT
C>    R0 TO THE LINE JOINING R1   AND  R2. THE LOGICAL VARIABLE INSIDE
C>    GIVES THE ADDITIONAL INFORMATION WHETHER THE FOOT OF THE PERPEN-
C>    DICULAR LIES WITHIN THE SEGMENT OR NOT.
C-----------------------------------------------------------------------
      SUBROUTINE PERP(R0,R1,R2,RD,INSIDE)
      implicit none
C
C     .. ARRAY ARGUMENTS ..
C
      REAL*8 R0(3),R1(3),R2(3),RD(3)
C
C     .. LOGICAL ARGUMENTS ..
C
      LOGICAL INSIDE
C
C     .. LOCAL SCALARS ..
C
      INTEGER   I
      REAL*8    DX,DY,DZ,S,D,DA,DB,DC,D1,D2,CO
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC DMAX1
C---------------------------------------------------------------------
      DX=R2(1)-R1(1)
      DY=R2(2)-R1(2)
      DZ=R2(3)-R1(3)
      S=R0(1)*DX+R0(2)*DY+R0(3)*DZ
      D=DX*DX+DY*DY+DZ*DZ
      IF(SQRT(D).LT.1.E-6) GO TO 100
      DA=S*DX+DY*(R1(1)*R2(2)-R1(2)*R2(1))+DZ*(R1(1)*R2(3)-R1(3)*R2(1))
      DB=S*DY+DZ*(R1(2)*R2(3)-R1(3)*R2(2))+DX*(R1(2)*R2(1)-R1(1)*R2(2))
      DC=S*DZ+DX*(R1(3)*R2(1)-R1(1)*R2(3))+DY*(R1(3)*R2(2)-R1(2)*R2(3))
      RD(1)=DA/D
      RD(2)=DB/D
      RD(3)=DC/D
      D1=(RD(1)-R1(1))**2+(RD(2)-R1(2))**2+(RD(3)-R1(3))**2
      D2=(RD(1)-R2(1))**2+(RD(2)-R2(2))**2+(RD(3)-R2(3))**2

      CO=D-DMAX1(D1,D2)
      INSIDE=.FALSE.
      IF ( CO.GT.1.E-6)  INSIDE=.TRUE.
      RETURN
  100 WRITE(6,200) (R1(I),I=1,3),(R2(I),I=1,3)
      STOP
  200 FORMAT(///33X,'FROM PERP:   IDENTICAL POINTS'/33X,2('(',3E14.6,')'
     *,3X))
      END SUBROUTINE

      END MODULE ShapeGeometryHelpers_mod
