!*==lattice3d.f    processed by SPAG 6.05Rc at 19:02 on 18 May 2004
!-------------------------------------------------------------------------------
! SUBROUTINE: LATTICE3D
!> @brief Generates the lattice vectors of direct and reciprocal space from
! basic translation vectors for a 3D system
!> @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine LATTICE3D(ALAT,BRAVAIS,RECBV,NGMAX,NRMAX,NSHLG,NSHLR,  &
   NSG,NSR,GN,RM,RMAX,GMAX,IPRINT,NMAXD,ISHLD)
   ! **********************************************************************
   ! *                                                                    *
   ! *  generate lattice vectors of direct and reciprocal space from      *
   ! *  basic translation vectors br                                      *
   ! *                                                                    *
   ! *  alat            : lattice constant                                *
   ! *  br(i,j)         : i=x,y,z j= 1,2,3 bravais vectors                *
   ! *                    *** in a.u. ****                                *
   ! *  rmax            : maximum radius in real space        (input)     *
   ! *  gmax            : maximum radius in reciprocal space  (input)     *
   ! *  ngmax           : Number of reciprocal lattice vectors            *
   ! *  gn(3,nmaxd)     : x,y,z   of reciprocal lattice vectors           *
   ! *  nrmax           : Number of real lattice vectors                  *
   ! *  rm(3,nmaxd)     : x,y,z  of real space vectors                    *
   ! *  nshlg           : shells in reciprocal space                      *
   ! *  nshlr           : shells in real space                            *
   ! *  nsg,nsr         : integer arrays, number of atoms in each shell   *
   ! *                                                                    *
   ! *  Dimension of arrays GN,RM changed from (4,*) to (3,*), the 4th    *
   ! *  one it is used only locally (GNR/RMR)       v.popescu May 2004    *
   ! *                                                                    *
   ! **********************************************************************
   implicit none
   ! ..
   ! .. Input variables
   integer, intent(in) :: NGMAX  !< Number of reciprocal space vectors
   integer, intent(in) :: NRMAX  !< Number of real space vectors rr
   integer, intent(in) :: NSHLG  !< Shells in reciprocal space
   integer, intent(in) :: NSHLR  !< Shells in real space
   integer, intent(in) :: NMAXD  !< Paremeters for the Ewald summations
   integer, intent(in) :: ISHLD  !< Paremeters for the Ewald summations
   integer, intent(in) :: IPRINT
   double precision, intent(in) :: ALAT   !< Lattice constant in a.u.
   double precision, dimension(3,3), intent(in) :: RECBV   !< Reciprocal basis vectors
   DOUBLE PRECISION, dimension(3,3), intent(in) :: BRAVAIS  !< Bravais lattice vectors
   double precision, intent(inout) :: RMAX   !< Ewald summation cutoff parameter for real space summation
   double precision, intent(inout) :: GMAX   !< Ewald summation cutoff parameter for reciprocal space summation
   ! ..
   ! .. Input/Output variables
   integer, dimension(ISHLD), intent(inout) :: NSG
   integer, dimension(ISHLD), intent(inout) :: NSR
   double precision, dimension(3,NMAXD), intent(inout) :: GN   !< x,y,z   of reciprocal lattice vectors
   double precision, dimension(3,NMAXD), intent(inout) :: RM   !< x,y,z  of real space vectors
   ! ..
   ! .. Local scalars ..
   integer :: IDINT
   integer :: I,K,L,M,N,N1,NG,NR,NSH,NSHL,NUMG,NUMGH,NUMR,NUMRH
   double precision :: DBLE
   double precision :: A,ABSGM,ABSRM,AG,AR,B,C,DA,DB,GX,GY,GZ,PI,RX,RY,RZ,VMIN
   ! ..
   ! .. Local arrays ..
   double precision, dimension(NMAXD)     :: GNR
   double precision, dimension(NMAXD)     :: RMR
   double precision, dimension(3)         :: ABSG
   double precision, dimension(3)         :: ABSR
   double precision, dimension(3,3)       :: BG
   double precision, dimension(3,3)       :: BR
   double precision, dimension(4,NMAXD)   :: CJ
   ! ..
   ! .. Intrinsic functions ..
   intrinsic :: ABS,ATAN,DBLE,IDINT,MAX,MOD,SQRT
   ! ..
   ! .. External subroutines ..
   external :: IOINPUT
   !----------------------------------------------------------------------------
   PI = 4.0D0*ATAN(1.0D0)
   !----------------------------------------------------------------------------
   ! OUTPUT
   !----------------------------------------------------------------------------
   write (1337,'(5X,2A,/)') '< LATTICE3D > : ','generating direct/reciprocal lattice vectors'
   !----------------------------------------------------------------------------
   ! OUTPUT
   !----------------------------------------------------------------------------
   RMAX = RMAX*ALAT
   GMAX = GMAX/ALAT
   !----------------------------------------------------------------------------
   ! OUTPUT
   !----------------------------------------------------------------------------
   write (1337,FMT=99001) RMAX,GMAX
   !----------------------------------------------------------------------------
   ! OUTPUT
   !----------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! Basic trans. vectors and basis vectors
   !----------------------------------------------------------------------------
   do I = 1,3
      BR(1,I) = BRAVAIS(1,I)*ALAT
      BR(2,I) = BRAVAIS(2,I)*ALAT
      BR(3,I) = BRAVAIS(3,I)*ALAT
   end do
   !----------------------------------------------------------------------------
   ! Generate primitive vectors BG of reciprocal space
   !----------------------------------------------------------------------------
   do I = 1,3
      BG(1,I) = RECBV(1,I)*2D0*PI/ALAT
      BG(2,I) = RECBV(2,I)*2D0*PI/ALAT
      BG(3,I) = RECBV(3,I)*2D0*PI/ALAT
   end do
   !----------------------------------------------------------------------------
   ! Estimate no. of lattice vectors
   !----------------------------------------------------------------------------
   do I = 1,3
      ABSR(I) = SQRT(BR(1,I)**2+BR(2,I)**2+BR(3,I)**2)
      ABSG(I) = SQRT(BG(1,I)**2+BG(2,I)**2+BG(3,I)**2)
   end do
   !
   ABSRM = MAX(ABSR(1),ABSR(2),ABSR(3))
   ABSGM = MAX(ABSG(1),ABSG(2),ABSG(3))
   ABSRM = 2.0D0*PI/ABSRM
   ABSGM = 2.0D0*PI/ABSGM
   NUMR = 2*(IDINT(RMAX/ABSGM)+1) + 1
   NUMG = 2*(IDINT(GMAX/ABSRM)+1) + 1
   NUMRH = NUMR/2 + 1
   NUMGH = NUMG/2 + 1
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                 generate lattice vectors of real space
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   NR = 0
   !----------------------------------------------------------------------------
   do L = 1,NUMR
      A = DBLE(L-NUMRH)
      do M = 1,NUMR
         B = DBLE(M-NUMRH)
         do N = 1,NUMR
            C = DBLE(N-NUMRH)
            !-------------------------------------------------------------------
            RX = A*BR(1,1) + B*BR(1,2) + C*BR(1,3)
            RY = A*BR(2,1) + B*BR(2,2) + C*BR(2,3)
            RZ = A*BR(3,1) + B*BR(3,2) + C*BR(3,3)
            AR = SQRT(RX*RX+RY*RY+RZ*RZ)
            !-------------------------------------------------------------------
            if ( AR.LE.RMAX ) then
               NR = NR + 1
               if ( NR.gt.NMAXD ) then
                  write (6,*)' ERROR: Dimension NMAXD in the inputcard is too small',NR,NMAXD
                  stop
               end if
               CJ(1,NR) = RX
               CJ(2,NR) = RY
               CJ(3,NR) = RZ
               CJ(4,NR) = AR
            end if
         end do
      end do
   end do
   !----------------------------------------------------------------------------
   NRMAX = NR
   !----------------------------------------------------------------------------
   ! Sort vectors in order of increasing absolute value
   !----------------------------------------------------------------------------
   DA = 1.D-06
   NSH = 0
   NSHL = -1
   !----------------------------------------------------------------------------
   do K = 1,NR
      VMIN = RMAX + 1.0D0
      do N = 1,NR
         if ( CJ(4,N)-VMIN.LT.0D0 ) then
            VMIN = CJ(4,N)
            N1 = N
         end if
      end do
      !
      NSHL = NSHL + 1
      RM(1,K) = CJ(1,N1)
      RM(2,K) = CJ(2,N1)
      RM(3,K) = CJ(3,N1)
      RMR(K)  = CJ(4,N1)
      DB = VMIN
      !-------------------------------------------------------------------------
      if ( DB.GT.DA+1.D-06 ) then
         NSH = NSH + 1
         if ( NSH.GT.ISHLD ) then
            write (6,*) ' ERROR: Dimension ISHLD in the inputcard is too small',   &
            NSH,ISHLD
            stop
         end if
         !
         NSR(NSH) = NSHL
         NSHL = 0
         DA = DB
      end if
      !-------------------------------------------------------------------------
      CJ(4,N1) = RMAX + 1.0D0
   end do
   !----------------------------------------------------------------------------
   NSH = NSH + 1
   NSHL = NSHL + 1
   IF ( NSH.GT.ISHLD ) THEN
      WRITE (6,*) ' ERROR: Dimension ISHLD in the inputcard is too small',NSH,ISHLD
      STOP
   END IF
   !
   NSR(NSH) = NSHL
   NSHLR = NSH
   IF ( NSHLR.LE.1 ) STOP ' ERROR: cut-off radius RMAX too small '
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                 generate lattice vectors of real space
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   NG = 0
   !----------------------------------------------------------------------------
   do L = 1,NUMG
      A = DBLE(L-NUMGH)
      do M = 1,NUMG
         B = DBLE(M-NUMGH)
         do N = 1,NUMG
            C = DBLE(N-NUMGH)
            !-------------------------------------------------------------------
            GX = A*BG(1,1) + B*BG(1,2) + C*BG(1,3)
            GY = A*BG(2,1) + B*BG(2,2) + C*BG(2,3)
            GZ = A*BG(3,1) + B*BG(3,2) + C*BG(3,3)
            AG = SQRT(GX*GX+GY*GY+GZ*GZ)
            !-------------------------------------------------------------------
            if ( AG.LE.GMAX ) then
               NG = NG + 1
               if ( NG.gt.NMAXD ) then
                  write (6,*)' ERROR: Dimension NMAXD in the inputcard is too small',  &
                  NG,NMAXD
                  stop
               end if
               CJ(1,NG) = GX
               CJ(2,NG) = GY
               CJ(3,NG) = GZ
               CJ(4,NG) = AG
            end if
         end do
      end do
   end do
   !----------------------------------------------------------------------------
   NGMAX = NG
   !----------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! Sort vectors in order of increasing abs. value
   !----------------------------------------------------------------------------
   DA = 1.D-06
   NSH = 0
   NSHL = -1
   !----------------------------------------------------------------------------
   do K = 1,NG
      VMIN = GMAX + 1.0D0
      do N = 1,NG
         if ( CJ(4,N).LT.VMIN ) then
            VMIN = CJ(4,N)
            N1 = N
         end if
      end do
      !
      NSHL = NSHL + 1
      GN(1,K) = CJ(1,N1)
      GN(2,K) = CJ(2,N1)
      GN(3,K) = CJ(3,N1)
      GNR(K)  = CJ(4,N1)
      DB = VMIN
      !-------------------------------------------------------------------------
      if ( DB.GT.DA+1.D-07 ) then
         NSH = NSH + 1
         if ( NSH.GT.ISHLD ) then
            write (6,*) ' ERROR: Dimension ISHLD in the inputcard is too small', &
            NSH,ISHLD
            stop
         end if
         !
         NSG(NSH) = NSHL
         NSHL = 0
         DA = DB
      end if
      !-------------------------------------------------------------------------
      CJ(4,N1) = GMAX + 1.0D0
   end do
   !----------------------------------------------------------------------------
   NSH = NSH + 1
   NSHL = NSHL + 1
   if ( NSH.GT.ISHLD ) then
      write (6,*) ' ERROR: Dimension ISHLD in the inputcard is too small',NSH,ISHLD
      stop
   end if
   !
   NSG(NSH) = NSHL
   NSHLG = NSH
   if ( NSHLG.le.1 ) stop ' ERROR: cut-off radius GMAX too small '
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! OUTPUT
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write (1337,FMT=99002)
   write (1337,FMT=99003) 'Direct  lattice',NRMAX,NSHLR,RMR(NRMAX)
   write (1337,FMT=99003) 'Recipr. lattice',NGMAX,NSHLG,GNR(NGMAX)
   write (1337,FMT=99004)
   !
   if ( IPRINT.lt.3 ) return
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! OUTPUT
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !----------------------------------------------------------------------------
   K = 0
   write (1337,FMT=99005) 'real-space'
   do L = 1,NSHLR
      write (1337,99006) L,NSR(L),RMR(K+1),(RM(M,K+1),M=1,3)
      do N = 2,NSR(L)
         write (1337,FMT=99007) (RM(M,K+N),M=1,3)
      end do
      if ( L.ne.NSHLR ) write (1337,99008)
      K = K + NSR(L)
   end do
   write (1337,99009)
   K = 0
   write (1337,FMT=99005) 'reciprocal'
   do L = 1,NSHLG
      write (1337,99006) L,NSG(L),GNR(K+1),(GN(M,K+1),M=1,3)
      do N = 2,NSG(L)
         write (1337,FMT=99007) (GN(M,K+N),M=1,3)
      end do
      if ( L.ne.NSHLG ) write (1337,99008)
      K = K + NSG(L)
   end do
   write (1337,99009)
   !----------------------------------------------------------------------------
   !
   99001 format (10X,'R max =',F9.5,' (a.u.)',/,10X,'G max =',F9.5,' (1/a.u.)',/)
   99002 format (10X,'               vectors  shells  max. R ',/,10X,   &
   '               ------------------------------')
   99003 format (10X,A,I7,2X,I6,2X,F9.5)
   99004 format (10X,'               ------------------------------',/)
   99005 format (10X,55('+'),/,18X,'generated ',A,' lattice vectors',/,10X,   &
   55('+'),/,10X,'shell Nvec    radius          x         y         z',/,10X,55('-'))
   99006 format (10X,I5,I5,F12.6,2X,3F10.5)
   99007 format (34X,3F10.5)
   99008 format (13X,52('-'))
   99009 format (10X,55('+'),/)
end subroutine LATTICE3D
