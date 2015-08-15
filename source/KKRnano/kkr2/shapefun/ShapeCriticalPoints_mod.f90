!> Module to determine radial mesh panel positions.

module ShapeCriticalPoints_mod
  implicit none
  private
  public :: criticalShapePoints
  
  contains

!------------------------------------------------------------------------------
!> Finds critical points for mesh generation and performs geometrical tests.

! *) khcp ? not needed

! TODO: add output parameters to argument list

! TODO: find out whats needed in subsequent part

!    &                 TOLVDIST,   ! Max. tolerance for distance of two vertices
!    &                 TOLEULER,   ! Used in calculation of Euler angles, subr. EULER
!    &                 NMIN,       ! Min. number of points in panel
!    &                 NVERTICES8,XVERT8,YVERT8,ZVERT8,NFACE8,LMAX8,
!    &                 DLT8,KEYPAN8,NM8

!> @param[in]  AFACE coefficients of plane equations of cell faces
!> @param[in]  BFACE
!> @param[in]  CFACE
!> @param[in]  DFACE
!> @param[in]  TOLVDIST numerical tolerance for geometrical constructions
!> @param[in]  TOLEULER numerical tolerance for plane rotations
!> @param[in]  NVERTICES number of vertices of cell
!> @param[in]  XVERT x-coordinates of cell vertices
!> @param[in]  YVERT y-coordinates of cell vertices
!> @param[in]  ZVERT z-coordinates of cell vertices
!> @param[in]  NFACE number of cell faces
!> @param[in]  LMAX calculate Shape-functions up to this cutoff
!> @param[out] NPAN number of panels, depending on critical points on radial mesh
!> @param[out] CRT critical points (where panels have to be placed)
!> Note: the smallest critical point corresponds to the muffin-tin radius

!=====================================================================
subroutine criticalShapePoints(AFACE,BFACE,CFACE,DFACE, &
  TOLVDIST, &
  TOLEULER, &
  NVERTICES,XVERT,YVERT,ZVERT,NFACE,LMAX, &
  NPAN, CRT, &
  NPAND)

  use shape_constants_mod, only: VERBOSITY, CHECK_GEOMETRY, ISUMD, LMAXD1, DP, PI
  use tetrahedra_common, only: NTT, FA, FB, FD, RD, ISIGNU
  use ShapeGeometryHelpers_mod, only: POLCHK

  !integer ::   NVERTICES(NFACED)
  !real(kind=DP) ::    AFACE(NFACED),BFACE(NFACED),CFACE(NFACED),DFACE(NFACED)
  !real(kind=DP) ::    XVERT(NVERTD,NFACED),YVERT(NVERTD,NFACED), &
  !ZVERT(NVERTD,NFACED)

  integer ::   NVERTICES(:)
  real(kind=DP) ::    AFACE(:),BFACE(:),CFACE(:),DFACE(:)
  real(kind=DP) ::    XVERT(:,:),YVERT(:,:), &
  ZVERT(:,:)


  ! output
  integer ::   NPAN
  real(kind=DP) ::    CRT(NPAND)

  integer, intent(in) :: NPAND

  !     .. SCALAR VARIABLES ..

  integer ::     IFACE,IPAN,ISUM,IS0
  integer ::     IV,IVERT,IVTOT,L,LMAX
  integer ::     NFACE,NVERT,NVTOT,IBMAX
  real(kind=DP) ::      SQ3O3,COA
  real(kind=DP) ::      A1,A2,A3,A4
  real(kind=DP) ::      TOLVDIST,TOLEULER
  logical ::    KHCP

  !     .. ARRAY VARIABLES ..
  real(kind=DP), allocatable ::     V(:,:)
  real(kind=DP) ::     Z(3)
  !     .. INTRINSIC FUNCTIONS ..

  intrinsic ABS,ACOS,DMAX1,DMIN1,ATAN,DFLOAT,SQRT
  integer :: NVERTD

  !-----------------------------------------------------------------------

  NVERTD = size(XVERT,1)
  allocate(V(3,NVERTD))

  KHCP=.false.

  !-----------------------------------------
  !  THIS CALL DOES SOME GEOMETRICAL TESTS (N.STEFANOU 98)
  if (CHECK_GEOMETRY .eqv. .true.) then
    call POLCHK(NFACE,NVERTICES,XVERT,YVERT,ZVERT,TOLVDIST)
  end if
  !**********   FOR HCP CASE ONLY    *********
  if (khcp) then
    COA  =SQRT(8.D0/3.D0)
    COA  =2.D0
    SQ3O3=SQRT(3.D0)/3.D0
  end if
  !**********   FOR HCP CASE  ONLY   *********

  IBMAX=(LMAX+1)*(LMAX+1)
  ISUM=0
  do L=0,LMAX
    IS0=(2*L+1)*(2*L+1)
    ISUM=ISUM+IS0
  end do

  if(ISUM > ISUMD .or. LMAX > LMAXD1) goto 200 ! check needed for routine DREAL

  IPAN=0
  IVTOT=0
  !.......................................................................
  !     S T O R A G E            I N    C O M M O N        B L O C K S
  !     C A L C U L A T I O N    O F    R O T A T I O N    M A T R I C E S
  !.......................................................................
  do IFACE=1,NFACE

    A1 = AFACE(IFACE)
    A2 = BFACE(IFACE)
    A3 = CFACE(IFACE)
    A4 = DFACE(IFACE)
    NVERT = NVERTICES(IFACE)

    ! -----------------
    Z(1)=A1/A4
    Z(2)=A2/A4
    Z(3)=A3/A4

    !************    FOR HCP CASE ONLY (TO REMOVE OTHERWISE)   ************
    if (khcp) then
      Z(1)=Z(1)*SQ3O3
      Z(3)=Z(3)*8.D0/COA/3.D0
    end if
    !************    FOR HCP CASE ONLY (TO REMOVE OTHERWISE)   ************

    do IVERT=1,NVERT

      V(1,IVERT) = XVERT(IVERT,IFACE)
      V(2,IVERT) = YVERT(IVERT,IFACE)
      V(3,IVERT) = ZVERT(IVERT,IFACE)

      !************    FOR HCP CASE ONLY (TO REMOVE OTHERWISE)   ************
      if (khcp) then
        V(1,IVERT)=V(1,IVERT)*SQ3O3
        V(3,IVERT)=V(3,IVERT)*COA
      end if
    end do

    call CRIT(IFACE,NVERT,V,Z,IPAN,IVTOT,TOLEULER,TOLVDIST,CRT,NPAND)

    if (VERBOSITY > 0) then
      write(6,105) IFACE,NTT(IFACE)
    end if

  end do ! END of loop over faces

  !.......................................................................
  !     D E F I N I T I O N    O F    T H E    S U I T A B L E    M E S H
  !.......................................................................
  NVTOT=IVTOT
  NPAN=IPAN

  if (VERBOSITY > 1) then
    write(6,102)
    do IV=1,NVTOT
      write(6,101) IV,FA(IV)/PI,FB(IV)/PI,FD(IV)/PI,RD(IV),ISIGNU(IV)
    end do
  end if

  return

200 write(6,107) ISUM,ISUMD,LMAX,LMAXD1
  stop

101 format(I10,4F10.4,I10)
102 format(//15X,'FA/PI',5X,'FB/PI',5X,'FD/PI',6X,'RD',8X,'ISIGNU'/)
105 format(/10X,I3,'-TH PYRAMID SUBDIVIDED IN ',I3,' TETRAHEDRA')
107 format(23X,'FROM MAIN : ISUM=',I7,'  GREATER THAN DIMENSIONED',I7/ &
  23X,'       OR   LMAX=',I7,'  GREATER THAN DIMENSIONED',I7)
110 format(11X,'BUT IS IDENTICAL TO A PREVIOUS ONE.')

end subroutine


!------------------------------------------------------------------------------
!>    Calculates critical points where panels have to be placed.
!>    @brief
!>    Call this routine for each face. Give face index as argument IFACE.
!>
!>    It also performs the subdivision in tetrahedra, which
!>    are then stored in module tetrahedras_common
!>    (needed for shape-function integration)
!>    (should be separate routine)
!>
!>    depends on: EULER, PERP, ROTATE (only needed for this routine!!!)
!>
!>    @param[in]     IFACE face index
!>    @param[in]     NVERT number of vertices
!>    @param[in]     V     coordinates of vertices V(3, NVERTD)
!>    @param[in,out] Z     coefficients determining plane of face Z(3) - changed on output!
!>                         Z(1) * X + Z(2) * Y + Z(3) * Z = 1 (ONE!)
!>    @param[in,out] IPAN  panel counter, pass 0 for 1st call
!>    @param[in,out] IVTOT tetrahedron index: has to be 0 for 1st call
!>    @param[in]     TOLEULER tolerance for Euler angles
!>    @param[in]     TOLVDIST tolerance for distances
!>    @param[in,out] CRT   array of critical points, CRT(NPAND)
!>    @param[in]     NPAND maximal number of panels allowed
subroutine CRIT(IFACE,NVERT,V,Z,IPAN,IVTOT,TOLEULER,TOLVDIST,CRT, &
                NPAND) ! NPAND from inc.geometry
  !-----------------------------------------------------------------------
  !     THIS ROUTINE CALCULATES THE CRITICAL POINTS 'CRT' OF THE SHAPE
  !     FUNCTIONS DUE TO THE FACE: Z(1)*X + Z(2)*Y + Z(3)*Z = 1
  !     THE FACE IS ROTATED THROUGH THE APPROPRIATE EULER ANGLES TO BE
  !     PERPENDICULAR TO THE Z-AXIS. A FURTHER SUBDIVISION OF THE CEN-
  !     TRAL PYRAMID INTO ELEMENTARY TETRAHEDRA IS PERFORMED. THE  NE-
  !     CESSARY QUANTITIES FOR THE CALCULATION ARE STORED IN COMMON.
  !-----------------------------------------------------------------------

  !     .. PARAMETER STATEMENTS ..

  use shape_constants_mod, only: PI, VERBOSITY, DP
  use tetrahedra_common, only: NTT, R0, RD, ISIGNU, FD, FA, FB
  use angles_common, only: ALPHA, BETA, GAMMA
  use ShapeGeometryHelpers_mod, only: PERP

  integer:: NPAND

  !     .. SCALAR ARGUMENTS ..

  integer::IFACE
  integer::NVERT
  integer::IPAN
  integer::IVTOT
  real(kind=DP):: TOLEULER, TOLVDIST

  !     .. ARRAY ARGUMENTS ..

  !real(kind=DP):: V(3,NVERTD)
  real(kind=DP):: V(:,:)
  real(kind=DP)::Z(3)
  real(kind=DP)::CRT(NPAND)

  !     .. LOCAL SCALARS ..

  integer::I
  integer::IX
  integer::ICORN
  integer::IVERT
  integer::INEW
  integer::IP
  integer::IVERTP
  integer::IVERT1
  integer::IBACK
  real(kind=DP)::ARG
  real(kind=DP)::A1
  real(kind=DP)::A2
  real(kind=DP)::A3
  real(kind=DP)::CF1
  real(kind=DP)::CF2
  real(kind=DP)::CF3
  real(kind=DP)::CO
  real(kind=DP)::CRRT
  real(kind=DP)::DD
  real(kind=DP)::DOWN
  real(kind=DP)::D1
  real(kind=DP)::D2
  real(kind=DP)::FF
  real(kind=DP)::F1
  real(kind=DP)::F2
  real(kind=DP)::OMEGA
  real(kind=DP)::RDD
  real(kind=DP)::S
  real(kind=DP)::SF1
  real(kind=DP)::SF2
  real(kind=DP)::SF3
  real(kind=DP)::UP
  real(kind=DP)::XJ
  real(kind=DP)::YJ
  real(kind=DP)::ZMOD2
  real(kind=DP)::ZVMOD

  !     .. LOCAL ARRAYS ..

  !integer::IN(NVERTD)
  !real(kind=DP)::VZ(3,NVERTD)

  integer, allocatable ::IN(:)
  real(kind=DP), allocatable ::VZ(:,:)

  real(kind=DP)::RDV(3)
  real(kind=DP)::ORIGIN(3)

  !     .. LOCAL LOGICAL ..

  logical:: INSIDE

  !     .. INTRINSIC FUNCTIONS ..

  intrinsic ABS,ACOS,DMAX1,DMIN1,ATAN2,SIGN,SQRT

  !     .. DATA STATEMENTS ..

  data ORIGIN/3*0.D0/

  real(kind=DP), parameter :: TOL_SMALL = 1D-6
  real(kind=DP), parameter :: TOL_LARGE = 1D-4

  !-----------------------------------------------------------------------

  allocate(IN(size(V,2)))
  allocate(VZ(3, size(V,2)))

  NTT(IFACE)=0

  if (VERBOSITY > 0) then
    write(6,203) IFACE,(Z(I),I=1,3)
  end if

  ZMOD2=Z(1)*Z(1)+Z(2)*Z(2)+Z(3)*Z(3)

  if(ZMOD2 <= TOL_SMALL) goto 100  ! check if normal vector of plane is valid

  S=2D0*PI

  do I=1,3
    Z(I)=Z(I)/ZMOD2
  end do

  IX=1
  ZVMOD=DSQRT((V(1,1)-Z(1))**2+(V(2,1)-Z(2))**2+(V(3,1)-Z(3))**2)

  if(ZVMOD < TOL_SMALL)  IX=2

  call EULER(Z,V(:,IX),IFACE,TOLEULER)

  if (VERBOSITY > 0) then
    write(6,204) ALPHA(IFACE)/PI,BETA(IFACE)/PI,GAMMA(IFACE)/PI
  end if

  call ROTATE(V,VZ,IFACE,NVERT)

  R0(IFACE)=1D0/DSQRT(ZMOD2)
  ICORN=0

  if (VERBOSITY > 0) then
    write(6,207)
  end if

  do IVERT =1,NVERT
    if(DABS(R0(IFACE)-VZ(3,IVERT)) > TOL_SMALL) goto 101  !check if vertices lie in same plane
    !.......................................................................
    !     D I S T A N C E S   O F   V E R T I C E S   F R O M   C E N T E R
    !.......................................................................
    CRRT=SQRT(VZ(1,IVERT)**2+VZ(2,IVERT)**2+VZ(3,IVERT)**2)

    !     Check if a new panel has been found (at radial point CRRT)
    INEW=1
    do IP=1,IPAN
      if(ABS(CRRT-CRT(IP)) < TOL_SMALL) INEW=0
    end do

    if(INEW == 1) then
      IPAN=IPAN+1

      if(IPAN > NPAND) goto 102

      CRT(IPAN)=CRRT
    end if
    IVERTP=IVERT+1

    if(IVERT == NVERT) IVERTP=1

    !.......................................................................
    !     D I S T A N C E S   O F   E D G E S   F R O M   C E N T E R
    !.......................................................................

    call PERP(ORIGIN,VZ(1,IVERT),VZ(1,IVERTP),RDV,TOLVDIST,INSIDE)

    RDD=SQRT(RDV(1)*RDV(1)+RDV(2)*RDV(2)+RDV(3)*RDV(3)) ! footpoint of line origin-to-edge

    if(INSIDE) then  ! add a new panel only when footpoint is contained on edge
      INEW=1
      do IP=1,IPAN
        if(ABS(RDD-CRT(IP)) < TOL_LARGE) INEW=0
      end do
      if(INEW == 1) then
        IPAN=IPAN+1
        if(IPAN > NPAND) goto 102
        CRT(IPAN)=RDD
      end if
    end if
    A1=SQRT(VZ(1,IVERT )*VZ(1,IVERT )+VZ(2,IVERT )*VZ(2,IVERT ))
    A2=SQRT(VZ(1,IVERTP)*VZ(1,IVERTP)+VZ(2,IVERTP)*VZ(2,IVERTP))
    DOWN=A1*A2
    UP=VZ(1,IVERT)*VZ(1,IVERTP)+VZ(2,IVERT)*VZ(2,IVERTP)

    if(DOWN > TOL_SMALL) then       ! true if not a corner
      ARG=UP/DOWN
      if(ABS(ARG) >= 1D0) ARG=SIGN(1D0,ARG)
      OMEGA=DACOS(ARG)
      S=S-OMEGA

      if(ABS(OMEGA-PI) > TOL_SMALL) then
        !.......................................................................
        !     S U B D I V I S I O N    I N T O    T E T R A H E D R A
        !.......................................................................
        NTT(IFACE)=NTT(IFACE)+1
        IVTOT=IVTOT+1

        if (VERBOSITY > 0) then
          write(6,205) IVTOT,IVERT,(VZ(I,IVERT),I=1,3)
          write(6,206)       IVERTP,(VZ(I,IVERTP),I=1,3)
        end if

        A3=DSQRT(RDV(1)*RDV(1)+RDV(2)*RDV(2))
        RD    (IVTOT)=RDD
        ISIGNU(IVTOT)=1
        CF1=VZ(1,IVERT )/A1
        CF2=VZ(1,IVERTP)/A2
        SF1=VZ(2,IVERT )/A1
        SF2=VZ(2,IVERTP)/A2
        CF3=RDV(1)/A3
        SF3=RDV(2)/A3

        if(ABS(SF1) < TOLEULER .and. ABS(CF1+1D0) < TOLEULER) then
          F1=PI
        else
          F1=2D0*ATAN2(SF1,CF1+1D0)
        end if

        if(ABS(SF2) < TOLEULER .and. ABS(CF2+1D0) < TOLEULER) then
          F2=PI
        else
          F2=2D0*ATAN2(SF2,CF2+1D0)
        end if

        if(ABS(SF3) < TOLEULER .and. ABS(CF3+1D0) < TOLEULER) then
          FD(IVTOT)=PI
        else
          FD(IVTOT)=2D0*ATAN2(SF3,CF3+1D0)
        end if

        FA(IVTOT)=DMIN1(F1,F2)
        FB(IVTOT)=DMAX1(F1,F2)
        if((FB(IVTOT)-FA(IVTOT)) > PI) then
          FF=FA(IVTOT)+2D0*PI
          FA(IVTOT)=FB(IVTOT)
          FB(IVTOT)=FF
        end if
        if((FA(IVTOT)-FD(IVTOT)) > PI)  FD(IVTOT)= 2D0*PI+FD(IVTOT)
        if((FD(IVTOT)-FA(IVTOT)) > PI)  FD(IVTOT)=-2D0*PI+FD(IVTOT)
      end if
    else
      ICORN=1
    end if
  end do   ! end of vertex loop
  !.......................................................................
  !     F O O T   O F   T H E    P E R P END I C U L A R   TO    T H E
  !     F A C E   O U T S I D E   O R  I N S I D E   T H E   P O L Y G O N
  !.......................................................................
  if(S < TOL_SMALL .or. ICORN == 1) then
    INEW=1
    do IP=1,IPAN
      if(ABS(R0(IFACE)-CRT(IP)) < TOL_LARGE) INEW=0
    end do
    if(INEW == 1) then
      IPAN=IPAN+1
      if(IPAN > NPAND) goto 102
      CRT(IPAN)=R0(IFACE)
    end if
  else
    do IVERT1=1,NVERT
      IN(IVERT1)=0
      do IVERT=1,NVERT
        IVERTP=IVERT+1

        if(IVERT == NVERT) IVERTP=1

        if(IVERT == IVERT1 .or. IVERTP == IVERT1) goto 7
        DOWN=VZ(2,IVERT1)*(VZ(1,IVERTP)-VZ(1,IVERT)) &
        -VZ(1,IVERT1)*(VZ(2,IVERTP)-VZ(2,IVERT))

        if(ABS(DOWN) <= TOL_SMALL) goto 7

        UP  =VZ(1,IVERT1)*(VZ(2,IVERT)*(VZ(1,IVERTP)+VZ(1,IVERT)) &
            -VZ(1,IVERT)*(VZ(2,IVERTP)+VZ(2,IVERT)))
        XJ=UP/DOWN
        YJ=XJ*VZ(2,IVERT1)/VZ(1,IVERT1)
        DD=(VZ(1,IVERTP)-VZ(1,IVERT))**2+(VZ(2,IVERTP)-VZ(2,IVERT))**2
        D1=(XJ-VZ(1,IVERT ))**2+(YJ-VZ(2,IVERT ))**2
        D2=(XJ-VZ(1,IVERTP))**2+(YJ-VZ(2,IVERTP))**2

        CO=DD-DMAX1(D1,D2)

        if(CO > TOL_SMALL) then
          IN(IVERT1)=1
          goto 6
        end if

7     continue
      end do
6   continue
    end do
    IBACK=IVTOT-NVERT

    do IVERT=1,NVERT
      IBACK=IBACK+1
      IVERTP=IVERT+1
      if(IVERT == NVERT) IVERTP=1
      if(IN(IVERT) == 0 .and. IN(IVERTP) == 0) ISIGNU(IBACK)=-1
    end do

  end if
  return
100 write(6,200) IFACE,(Z(I),I=1,3)
  stop
101 write(6,201) IFACE,R0(IFACE),((VZ(I,IVERT),I=1,3),IVERT=1,NVERT)
  stop
102 write(6,202) IPAN,NPAND
  stop
200 format(//13X,'FATAL ERROR FROM CRIT: THE',I3,'-TH FACE OF THE POLYHEDRON PASSES THROUGH THE CENTER'  /13X,'(',3e14.7,' )')
201 format(//13X,'FATAL ERROR FROM CRIT: THE VERTICES OF THE',I3,'-TH ROTATED POLYGON DO NOT LIE ON THE PLANE:'  ,E13.6,' *Z = 1'/30(/13X, 3e13.6))
202 format(//13X,'ERROR FROM CRIT: NUMBER OF PANELS=',I5,' GREATER THAN DIMENSIONED='  ,I5)
203 format(//80('*')/3X,'FACE:',I3,' EQUATION:',F10.4,'*X +',F10.4, &
  '*Y +',F10.4,'*Z  =  1')
204 format(3X,'ROTATION ANGLES  :',3(F10.4,4X)/)
205 format(I5,'       VZ(',I2,')  =  (',3F10.4,' )')
206 format(5X,'       VZ(',I2,')  =  (',3F10.4,' )')
207 format(/'TETRAHEDRON',14X,'COORDINATES'/11('*'),14X,11('*')/)
end subroutine CRIT

!-----------------------------------------------------------------------
!>    GIVEN TWO DISTINCT POINTS (Z(1),Z(2),Z(3)) AND (XX(1),XX(2),XX(3))
!>    THIS ROUTINE DEFINES  A LOCAL COORDINATE  SYSTEM WITH THE  Z- AXIS
!>    PASSING THROUGH (Z(1),Z(2),Z(3))  AND THE X- AXIS PARALLEL TO  THE
!>    VECTOR : (XX(1)-Z(1),XX(2)-Z(2),XX(3)-Z(3)).
!>    THE EULER ANGLES ROTATING THIS LOCAL COORDINATE SYSTEM BACK TO THE
!>    ORIGINAL FRAME OF REFERENCE ARE CALCULATED  AND STORED  IN COMMON.
!-----------------------------------------------------------------------
subroutine EULER(Z,XX,IFACE,TOLEULER)

  use shape_constants_mod, only: PI, DP
  use angles_common, only: ALPHA, BETA, GAMMA

  !     .. SCALAR ARGUMENTS ..

  integer ::   IFACE

  !     .. ARRAY ARGUMENTS ..

  real(kind=DP) :: XX(3),Z(3)

  !     .. LOCAL SCALARS ..

  integer ::   I
  real(kind=DP) ::    RX,RZ,S,P,RZP,SA,CA,SG,CG
  real(kind=DP) ::    TOLEULER   ! introduced by Phivos (05.2008) to account for inaccuracies.
  ! Earlier, 1.D-5 was hard-coded at the places in this subr. where TOL is used
  !     DATA TOLEULER /1.D-10/

  !     .. LOCAL ARRAYS ..

  real(kind=DP) ::    X(3),Y(3)

  !     .. INTRINSIC FUNCTIONS ..

  intrinsic SQRT,ACOS,ABS,ATAN2

  real(kind=DP), parameter :: TOLERANCE = 1D-6

  !-----------------------------------------------------------------------
  if(IFACE > size(ALPHA)) goto 20
  do I=1,3
    X(I)=XX(I)-Z(I)
  end do
  RX=DSQRT(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
  RZ=DSQRT(Z(1)*Z(1)+Z(2)*Z(2)+Z(3)*Z(3))

  if (RX < TOLERANCE .or. RZ < TOLERANCE)  goto 30

  S=X(1)*Z(1)+X(2)*Z(2)+X(3)*Z(3)

  if (S > TOLERANCE) goto 30

  P=DSQRT(Z(1)*Z(1)+Z(2)*Z(2))

  do I=1,3
    X(I)=X(I)/RX
    Z(I)=Z(I)/RZ
  end do

  ALPHA(IFACE)=0D0
  BETA(IFACE) =DACOS(Z(3))

  if (P < TOLEULER) goto 10

  RZP=RZ/P
  Y(1)=Z(2)*X(3)-Z(3)*X(2)
  Y(2)=Z(3)*X(1)-Z(1)*X(3)
  Y(3)=Z(1)*X(2)-Z(2)*X(1)
  SA=Y(3)*RZP
  CA=X(3)*RZP

  if(DABS(SA) < TOLEULER .and. DABS(CA+1.D0) < TOLEULER)  then
    ALPHA(IFACE)=PI
  else
    ALPHA(IFACE)=2D0*DATAN2(SA,CA+1D0)
  end if

  SG= Z(2)*RZP
  CG=-Z(1)*RZP

  if(DABS(SG) < TOLEULER .and. DABS(CG+1D0) < TOLEULER)  then
    GAMMA(IFACE)=PI
  else
    GAMMA(IFACE)=2D0*DATAN2(SG,CG+1D0)
  end if

  do I=1,3
    Z(I)=Z(I)*RZ
  end do

  return

10 SG=-Z(3)*X(2)
  CG= Z(3)*X(1)

  if(DABS(SG) < TOLEULER .and. DABS(CG+1D0) < TOLEULER)  then
    GAMMA(IFACE)=PI
  else
    GAMMA(IFACE)=2D0*DATAN2(SG,CG+1D0)
  end if

  do I=1,3
    Z(I)=Z(I)*RZ
  end do

  return

20 write(6,100) IFACE,size(ALPHA)
  stop
30 write(6,101) (X(I),I=1,3),(Z(I),I=1,3)
  stop
100 format(//13X,'NUMBER OF FACES:',I5,' GREATER THAN DIMENSIONED',I5)
101 format(/13X,'FROM EULER,ILLEGAL VECTORS:'/13X,2(' (',3e13.6,' )'))
end subroutine EULER

!-----------------------------------------------------------------------
!> Rotation by Euler angles given in module angles_common.
!>    THIS ROUTINE PERFORMS THE ROTATION OF NVERT VECTORS THROUGH THE
!>    EULER ANGLES: ALPHA(IFACE),BETA(IFACE),GAMMA(IFACE).
!>    V (I,IVERT) : INPUT   VECTORS
!>    VZ(I,IVERT) : ROTATED VECTORS
!-----------------------------------------------------------------------
subroutine ROTATE(V,VZ,IFACE,NVERT)

  use shape_constants_mod, only: PI, DP
  use angles_common, only: ALPHA, BETA, GAMMA


  !     .. SCALAR ARGUMENTS ..

  integer ::   IFACE,NVERT

  !     .. ARRAY ARGUMENTS ..

  !real(kind=DP) :: V(3,NVERTD),VZ(3,NVERTD)
  real(kind=DP) :: V(:,:),VZ(:,:)


  !     .. LOCAL SCALARS ..

  integer ::   I,J,IVERT
  real(kind=DP) ::    SA,SB,SG,CA,CB,CG

  !     .. LOCAL ARRAYS ..

  real(kind=DP) :: A(3,3)

  !-----------------------------------------------------------------------
  CA=DCOS(ALPHA(IFACE))
  SA=DSIN(ALPHA(IFACE))
  CB=DCOS(BETA(IFACE))
  SB=DSIN(BETA(IFACE))
  CG=DCOS(GAMMA(IFACE))
  SG=DSIN(GAMMA(IFACE))
  A(1,1)=CA*CB*CG-SA*SG
  A(2,1)=SA*CB*CG+CA*SG
  A(3,1)=-SB*CG
  A(1,2)=-CA*CB*SG-SA*CG
  A(2,2)=-SA*CB*SG+CA*CG
  A(3,2)=SB*SG
  A(1,3)=CA*SB
  A(2,3)=SA*SB
  A(3,3)=CB
  do IVERT=1,NVERT
    do I=1,3
      VZ(I,IVERT)=0D0
      do J=1,3
        VZ(I,IVERT)=VZ(I,IVERT)+A(I,J)*V(J,IVERT)
      end do
    end do
  end do
  return
end subroutine ROTATE

end module ShapeCriticalPoints_mod
