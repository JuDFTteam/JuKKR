!> Module to calculate shape-functions + mesh given Voronoi-cell data.

module ShapeFunctions_mod
  implicit none
  private
  public :: shapef
  
  
#ifdef DEBUGSHAPEFUNCTIONS
  integer, parameter :: verbosity = 2
#else
  integer, parameter :: verbosity = 0
#endif

  contains

!------------------------------------------------------------------------------
!> Calculates non-zero shape-functions from geometrical (Voronoi-)cell input.
!> @brief
!> Assumes that center of Voronoi cell is at (0,0,0)
!> @todo
!> dynamic array solution - get rid of LMAXD1, MESHND, IPAND
!> the use of modules to store data is a bit ugly - but who cares?
!> @author N. Stefanou
!> @author P. Mavropoulos
!> @author E. Rabel (modernisation, documentation)
!> @author R. Zeller et al.
!> @param[in]  NPOI minimal number of radial mesh points, on outpuit MESHN
!>             contains real, number of mesh points used
!> @param[in]  AFACE coefficients of plane equations of cell faces
!> @param[in]  BFACE
!> @param[in]  CFACE
!> @param[in]  DFACE
!> @param[in]  TOLVDIST numerical tolerance for geometrical constructions
!> @param[in]  TOLEULER numerical tolerance for plane rotations
!> @param[in]  NMIN minimum number of points in critical regions of radial mesh
!> @param[in]  NVERTICES number of vertices of each face
!> @param[in]  XVERT x-coordinates of cell vertices
!> @param[in]  YVERT y-coordinates of cell vertices
!> @param[in]  ZVERT z-coordinates of cell vertices
!> @param[in]  NFACE number of cell faces
!> @param[in]  LMAX calculate Shape-functions up to this cutoff
!> @param[in]  KEYPAN =0 not used
!> @param[in]  DLT step size for Gauss-Legendre integration
!> @param[out] NPAN number of panels, depending on critical points on radial mesh
!> @param[out] NM number of points in each critical region
!> @param[out] XRN radial mesh points
!> @param[out] DRN weights for radial mesh points
!> @param[out] MESHN number of radial mesh points
!> @param[out] THETAS_S on output it contains the non-zero shape-functions
!!                      dimension MESHND x IBMAXD
!> @param[out] LMIFUN_S array with lm-indeices for which the shapefunction is non-zero
!> @param[out] NFUN number of non-zero shape-functions
!> @param[in]  IBMAXD
!> @param[in]  MESHND
!> @param[in]  NPAND

!------------------ OLD HELP TEXT: partially out of date ----------------------
  !----------------------------------------------------------------------
  !                       S H A P E   P R O G R A M
  !        F O R  A R B I T R A R Y  V O R O N O I  P O L Y H E D R A
  !                                                           N.Stefanou
  !----------------------------------------------------------------------
  !     In order to improve the efficency of the code and to
  !     check more, this program was changed in summer 1998 by N.Stefanou
  !     ATTENTION: BUG is removed!
  !                In the old version IPAN=-1 is wrong and has to be
  !                set to IPAN=0.
  !     THIS PROGRAM CALCULATES THE ANGULAR MOMENTUM COMPONENTS OF THE
  !     SHAPE FUNCTION FOR AN ARBITRARY VORONOI POLYHEDRON.
  !     A REAL SPHERICAL HARMONIC BASIS IS USED FOR THE DECOMPOSITION.
  !     ON INPUT WE GIVE :
  !        LMAX          :  MAXIMUM ANGULAR MOMENTUM
  !        DLT           :  DEFINES THE STEP FOR GAUSS-LEGENDRE CALC.
  !  TIME (DEC)  DLT    TOTAL ENERGY (BCC-TEST CdSb in Ge 12 Shells)
  !  1104S      0.002  -.59583341\
  !   315S      0.005  -.59583341 \
  !    51S      0.050  -.59583341  --> NO CHANGES ALSO IN
  !    40S      0.100  -.59583341 /    CHARGES OR FORCES
  !    38S      0.200  -.59583341/
  !    38S      0.300  -.59583354---> CHARGES DIFFER IN 10NTH DIGIT
  !    35S      0.400  -.59583673---> CHARGES DIFFER IN 7NNTH DIGIT
  !                                   FORCES DIFFER IN 5TH DIGIT
  !    --->TO BE AT THE SAVE SIDE USE 0.05 OR 0.1 (SHOULD BE QUITE GOOD)
  !        SIMILAR RESULTS WERE HELD FOR Cu in Fe NN relaxation.
  !        NFACE         :  NUMBER OF FACES OF THE POLYHEDRON
  !        KEYPAN        :  KEY TO DEFINE  THE  RADIAL  MESH.  IF KEYPAN=
  !                         THE DEFAULT  RADIAL  DIVISION  OF PANNELS GIVE
  !                         IN DATA STATEMENT IS USED.OTHERWISE THE  NUMBE
  !                         OF  RADIAL  MESH  POINTS PER PANNEL  (NM(IPAN)
  !                         IS READ IN INPUT
  !                      ** IN THIS VERSION THE MESH IS DETERMINED
  !                         BY SUBROUTINE MESH0.
  !        Z(I)          :  COEFFICIENTS OF THE EQUATION OF A FACE
  !                         Z(1)*X + Z(2)*Y + Z(3)*Z  =  1
  !        NVERT         :  NUMBER OF VERTICES OF A FACE
  !        V(I,IVERT)    :  COORDINATES OF THE VERTICES OF A FACE
  !        NEWSCH(IFACE) :  INTEGER   PARAMETER TO CALCULATE   (=1)
  !                         THE CONTRIBUTION OF THE CORRESPONDING
  !                         PYRAMID TO THE SHAPE FUNCTIONS  .  IF
  !                         NEWSCH.NE.1 THE CONTRIBUTION IS TAKEN
  !                         EQUAL TO THAT OF THE PREVIOUS PYRAMID


  !     IN ORDER TO SAVE MEMORY WE STORE IN LOCAL TEMPORARY FILES IN  UNIT
  !     30+1 , 30+2 , ... , 30+NFACE THE TRANSFORMATION MATRICES ASSOCIATE
  !     WITH THE ROTATION OF EACH PYRAMID. THE TEMPORARY DIRECT ACCESS FIL
  !     IN UNIT 10 CONTAINS THE CALCULATED COMPONENTS OF THE SHAPE FUNCTIO

  !                  ...........I N P U T  C A R D...(Bcc/fcc)

  !                                       if not (Bcc/fcc) change main prg
  ! bcc                          <----- Gives the lattice parameters
  !    16    1   0.05000                lmax,nkey,division
  !                                     LMAX=4*LMAX(KKR),
  !                                     NKEY is not used,
  !                                     DIVISION is DLT
  !   125    0                          number of mesh points,keypan
  !                                     Number of mesh points
  !                                     used for the radial mesh (Depends
  !                                     on the number of pannels).
  !                                     If keypan is 1,
  !                                     then the radial mesh division is
  !                                     taken from the input
  !   -3.30000 -3.30000  -3.30000      relaxation percent

  !    63   32   30    7   21   15   15   15   15   15 \
  !    15   17   15   15   23   15   15    0    0    0  \ This is the
  !     0   17   15   15   23   15   15    0    0    0  / radial mesh info
  !     0    0    0    0    0    0    0    0    0    0 /
  !                  .........................................


  !     THE DEFINITION OF REAL SPHERICAL HARMONICS IS NOT THE STANDARD  ON
  !     REFERED IN THE PAPER:
  !     N.STEFANOU,H.AKAI AND R.ZELLER,COMPUTER PHYS.COMMUN. 60 (1990) 231
  !     IF YOU WANT TO HAVE ANGULAR MOMENTUM COMPONENTS IN THE STANDARD BA
  !     SIS CHANGE THE FOLLOWING STATEMENTS IN THE ROUTINES :
  !     IN CCOEF      ISI=1                        ---->    ISI=1-2*MOD(M,
  !     IN D_REAL     IF(MOD(M+MP),2).EQ.0) D=-D   ---->    DELETE THE LIN


!------------------------------------------------------------------------------
  subroutine shapef(npoi, planes, tolvdist, toleuler, nmin, nvertices, vert, nface, lmax, &
                    dlt, npan, nm, xrn, drn, meshn, &  ! radial mesh ! output parameters
                    thetas_s, lmifun_s, nfun, & ! shape function
                    ibmaxd, meshnd, npand, atom_id)
    
    use ShapeCriticalPoints_mod, only: criticalshapepoints
    use ShapeStandardMesh_mod, only: mesh
    use ShapeIntegration_mod, only: shapeintegration
    use PolygonFaces_mod, only: PolygonFace, destroy
    integer :: ist

    integer, intent(in) :: npoi
    double precision, intent(in) :: tolvdist, toleuler
    integer, intent(in) :: nmin, nface, lmax
    double precision, intent(in) :: dlt

    integer, intent(in) :: ibmaxd
    integer, intent(in) :: meshnd
    integer, intent(in) :: npand

    integer, intent(in) :: nvertices(:) ! (nfaced)
    double precision, intent(in) :: planes(0:,:) ! (0:3,nfaced)
    double precision, intent(in) :: vert(:,:,:) ! (3,nvertd,nfaced)
    ! output radial grid
    double precision, intent(out) :: xrn(meshnd), drn(meshnd)
    integer, intent(out) :: npan, meshn, nm(npand)
    ! output shape functions
    double precision, intent(out) :: thetas_s(meshnd,ibmaxd)
    integer, intent(out) :: lmifun_s(ibmaxd)
    integer, intent(out) :: nfun
    integer, intent(in) :: atom_id

    double precision :: crt(npand) ! critical points
    type(PolygonFace), allocatable :: faces(:)
    
    nm = 0
    xrn = 0.d0
    drn = 0.d0
    crt = 0.d0
    meshn = 0
    npan = 0
    
    thetas_s = 0.d0
    lmifun_s = 0
    nfun = 0

    allocate(faces(size(nvertices)), stat=ist)

    call criticalShapePoints(planes, tolvdist, toleuler, nvertices, vert, nface, faces, npan, crt, atom_id) ! outputs: faces, npan, crt

    ! increase number of mesh points if necessary but use at least 'npoi' points (otherwise mesh0 complains)
    call mesh(crt, npan, nm, xrn, drn, meshn, max(npoi, npan*nmin), 0, nmin, meshnd, npand, verbosity)

    call shapeIntegration(lmax, faces(1:nface), meshn, xrn, dlt, thetas_s, lmifun_s, nfun, meshnd, ibmaxd)

    npan = npan - 1 ! in old code 1 is substracted from npan - why?

    call destroy(faces)
    
  endsubroutine ! shapef

endmodule ! ShapeFunctions_mod

