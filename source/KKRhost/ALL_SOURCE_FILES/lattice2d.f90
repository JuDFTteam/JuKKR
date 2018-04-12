!-------------------------------------------------------------------------------
! SUBROUTINE: lattice2d
!> @brief Generates the lattice vectors of direct and reciprocal space from
!> basic translation vectors for a 2D system
!> @note - V. Popescu May 2004: The routine has been brought to a form which is very similar to
!> LATTICE2D -- from which it has been originally derived. Dimension of arrays GN,RM
!> changed from (4,*) to (2,*), the 4th one it is used only locally (GNR/RMR) -- only GN/RM(2,*) are
!> actually needed in EWALD2D
!> @note - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine LATTICE2D(ALAT,BRAVAIS,RECBV,NGMAX,NRMAX,NSHLG,NSHLR,  &
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
   ! *  gn(2,nmaxd)     : x,y,z   of reciprocal lattice vectors           *
   ! *  nrmax           : Number of real lattice vectors                  *
   ! *  rm(2,nmaxd)     : x,y,z  of real space vectors                    *
   ! *  nshlg           : shells in reciprocal space                      *
   ! *  nshlr           : shells in real space                            *
   ! *  nsg,nsr         : integer arrays, number of atoms in each shell   *
   ! *                                                                    *
   ! *  The routine has been brought to a form which is very similar to   *
   ! *  LATTICE2D -- from which it has been originally derived            *
   ! *  Dimension of arrays GN,RM changed from (4,*) to (2,*), the 4th    *
   ! *  one it is used only locally (GNR/RMR) -- only GN/RM(2,*) are      *
   ! *  actually needed in EWALD2D                  v.popescu May 2004    *
   ! *                                                                    *
   ! **********************************************************************
   implicit none
   ! ..
   ! .. Input variables
   integer, intent(in) :: NMAXD  !< Paremeters for the Ewald summations
   integer, intent(in) :: ISHLD  !< Paremeters for the Ewald summations
   integer, intent(in) :: IPRINT
   double precision, intent(in) :: ALAT   !< Lattice constant in a.u.
   double precision, dimension(3,3), intent(in) :: RECBV   !< Reciprocal basis vectors
   double precision, dimension(3,3), intent(in) :: BRAVAIS  !< Bravais lattice vectors
   ! ..
   ! .. Input/Output variables
   double precision, intent(inout) :: GMAX   !< Ewald summation cutoff parameter for reciprocal space summation
   double precision, intent(inout) :: RMAX   !< Ewald summation cutoff parameter for real space summation
   integer, dimension(ISHLD), intent(inout) :: NSG
   integer, dimension(ISHLD), intent(inout) :: NSR
   double precision, dimension(2,NMAXD), intent(inout) :: GN   !< x,y,z   of reciprocal lattice vectors
   double precision, dimension(2,NMAXD), intent(inout) :: RM   !< x,y,z  of real space vectors
   ! .. Output variables
   integer, intent(out) :: NSHLR  !< Shells in real space
   integer, intent(out) :: NSHLG  !< Shells in reciprocal space
   integer, intent(out) :: NRMAX  !< Number of real space vectors rr
   integer, intent(out) :: NGMAX  !< Number of reciprocal space vectors
   ! ..
   ! .. Local scalars ..
   integer :: IDINT
   integer :: I,K,L,M,N,N1,NG,NR,NSH,NSHL,NUMG,NUMGH,NUMR,NUMRH
   double precision :: DBLE
   double precision :: RX,RY,VMIN
   double precision :: A,ABSGM,ABSRM,AG,AR,B,DA,DB,GX,GY,PI
   ! ..
   ! .. Local arrays ..
   double precision, dimension(NMAXD)  :: GNR
   double precision, dimension(NMAXD)  :: RMR
   double precision, dimension(3)      :: ABSG
   double precision, dimension(3)      :: ABSR
   double precision, dimension(NMAXD)  :: LENGTH
   double precision, dimension(3,3)       :: BG
   double precision, dimension(3,3)       :: BR
   double precision, dimension(4,NMAXD)   :: CJ
   ! ..
   ! .. Intrinsic functions ..
   intrinsic ABS,ATAN,DBLE,IDINT,MAX,MOD,SQRT
   ! ..
   ! .. External subroutines ..
   external IOINPUT
   !----------------------------------------------------------------------------
   PI = 4.0D0*ATAN(1.0D0)
   !----------------------------------------------------------------------------
   ! OUTPUT
   !----------------------------------------------------------------------------
   write (1337,'(5X,2A,/)') '< LATTICE2D > : ','generating direct/reciprocal lattice vectors'
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
   end do
   !----------------------------------------------------------------------------
   ! Generate primitive vectors BG of reciprocal space
   !----------------------------------------------------------------------------
   do I = 1,3
      BG(1,I) = RECBV(1,I)*2D0*PI/ALAT
      BG(2,I) = RECBV(2,I)*2D0*PI/ALAT
   end do
   !----------------------------------------------------------------------------
   ! Estimate no. of lattice vectors
   !----------------------------------------------------------------------------
   do I = 1,3
      ABSR(I) = SQRT(BR(1,I)**2+BR(2,I)**2)
      ABSG(I) = SQRT(BG(1,I)**2+BG(2,I)**2)
   end do
   !
   ABSRM = MAX(ABSR(1),ABSR(2))
   ABSGM = MAX(ABSG(1),ABSG(2))
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
   write(1337,*) 'Real space...'
   NR = 0
   !----------------------------------------------------------------------------
   do L = 1,NUMR
      A = DBLE(L-NUMRH)
     do M = 1,NUMR
         B = DBLE(M-NUMRH)
         !----------------------------------------------------------------------
         RX = A*BR(1,1) + B*BR(1,2)
         RY = A*BR(2,1) + B*BR(2,2)
         AR = SQRT(RX*RX+RY*RY)
         !----------------------------------------------------------------------
         if ( AR.LE.RMAX ) then
            NR = NR + 1
            if ( NR.GT.NMAXD ) then
               write (6,*) &
                  'lattice2d: ERROR: Dimension NMAXD in the inputcard too small',&
                  NR,NMAXD
               stop 'lattice2d'
            end if
            CJ(1,NR) = RX
            CJ(2,NR) = RY
            CJ(3,NR) = 0D0
            CJ(4,NR) = AR
         end if
      end do
   end do
   !----------------------------------------------------------------------------
   NRMAX = NR
   !----------------------------------------------------------------------------
   ! Sort vectors in order of increasing absolute value
   !----------------------------------------------------------------------------
   write(1337,FMT='(A11,I8,A11)') '...sorting ',NRMAX,' vectors...'

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
      RMR(K)  = CJ(4,N1)
      DB = VMIN
      !-------------------------------------------------------------------------
      if ( DB.GT.DA+1.D-06 ) then
         NSH = NSH + 1
         if ( NSH.GT.ISHLD ) then
            write (6,*) ' ERROR: Dimension ISHLD in the inputcard too small', &
            NSH,ISHLD
            stop 'lattice2d'
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
   if ( NSH.GT.ISHLD ) then
      write (6,*) ' ERROR: Dimension ISHLD in the inputcard too small',NSH,ISHLD
      stop 'lattice2d'
   end if
   !
   NSR(NSH) = NSHL
   NSHLR = NSH
   if ( NSHLR.LE.1 ) stop 'lattice2d: ERROR: cut-off radius RMAX too small '
   write(1337,*) '...done.'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                 generate lattice vectors of real space
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   write(1337,*) 'Reciprocal space...'
   NG = 0
   !----------------------------------------------------------------------------
   do L = 1,NUMG
      A = DBLE(L-NUMGH)
      do M = 1,NUMG
         B = DBLE(M-NUMGH)
         !----------------------------------------------------------------------
         GX = A*BG(1,1) + B*BG(1,2)
         GY = A*BG(2,1) + B*BG(2,2)
         AG = SQRT(GX*GX+GY*GY)
         !----------------------------------------------------------------------
         if ( AG.LE.GMAX ) then
            NG = NG + 1
            if ( NG.GT.NMAXD ) then
               write (6,*)' ERROR: Dimension NMAXD in the inputcard too small',&
               NG,NMAXD
               stop 'lattice2d'
            end if
            CJ(1,NG) = GX
            CJ(2,NG) = GY
            CJ(3,NG) = 0D0
            CJ(4,NG) = AG
         end if
      end do
   end do
   !----------------------------------------------------------------------------
   NGMAX = NG
   !----------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! Sort vectors in order of increasing abs. value
   !----------------------------------------------------------------------------
   write(1337,FMT='(A11,I8,A11)') '...sorting ',NGMAX,' vectors...'
   do N = 1,NG
      LENGTH(N) = CJ(4,N)
   enddo
   DA = 1.D-06
   NSH = 0
   NSHL = -1
   !----------------------------------------------------------------------------
   do K = 1,NG
      VMIN = GMAX + 1.0D0
      do N = 1,NG
         if ( LENGTH(N).LT.VMIN ) then ! ( CJ(4,N).LT.VMIN ) THEN
            VMIN = LENGTH(N) ! CJ(4,N)
            N1 = N
         end if
      end do
      !
      NSHL = NSHL + 1
      GN(1,K) = CJ(1,N1)
      GN(2,K) = CJ(2,N1)
      GNR(K)  = LENGTH(N1)  ! CJ(4,N1)
      DB = VMIN
      !-------------------------------------------------------------------------
      if ( DB.GT.DA+1.D-07 ) then
         NSH = NSH + 1      ! Number of shells of different length
         if ( NSH.GT.ISHLD ) then
            write (6,*) ' ERROR: Dimension ISHLD in the inputcard too small', &
            NSH,ISHLD
            stop 'lattice2d'
         end if
         !
         NSG(NSH) = NSHL    ! Number of vectors in shell
         NSHL = 0
         DA = DB
      end if
      !-------------------------------------------------------------------------
      LENGTH(N1) = GMAX + 1.0D0 !  CJ(4,N1) = GMAX + 1.0D0
   end do
   !----------------------------------------------------------------------------
   NSH = NSH + 1
   NSHL = NSHL + 1
   if ( NSH.GT.ISHLD ) then
      write (6,*) ' ERROR: Dimension ISHLD in the inputcard too small',NSH,ISHLD
      stop 'lattice2d'
   end if
   !
   NSG(NSH) = NSHL
   NSHLG = NSH
   if ( NSHLG.LE.1 ) stop 'lattice2dERROR: cut-off radius GMAX too small '

   write(1337,*) '...done.'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! OUTPUT
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write (1337,FMT=99002)
   write (1337,FMT=99003) 'Direct  lattice',NRMAX,NSHLR,RMR(NRMAX)
   write (1337,FMT=99003) 'Recipr. lattice',NGMAX,NSHLG,GNR(NGMAX)
   write (1337,FMT=99004)
   !
   if ( IPRINT.LT.3 ) return
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! OUTPUT
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !----------------------------------------------------------------------------
   K = 0
   write (1337,FMT=99005) 'real-space'
   do L = 1,NSHLR
      write (1337,99006) L,NSR(L),RMR(K+1),(RM(M,K+1),M=1,2)
      do N = 2,NSR(L)
         write (1337,FMT=99007) (RM(M,K+N),M=1,2)
      end do
      if ( L.NE.NSHLR ) write (1337,99008)
      K = K + NSR(L)
   end do
   write (1337,99009)
   K = 0
   write (1337,FMT=99005) 'reciprocal'
   do L = 1,NSHLG
      write (1337,99006) L,NSG(L),GNR(K+1),(GN(M,K+1),M=1,2)
      do N = 2,NSG(L)
         write (1337,FMT=99007) (GN(M,K+N),M=1,2)
      end do
      if ( L.NE.NSHLG ) write (1337,99008)
      K = K + NSG(L)
   end do
   write (1337,99009)
   !----------------------------------------------------------------------------
   !
   99001 format (10X,'R max =',F10.5,' (a.u.)',/,10X,'G max =',F10.5,' (1/a.u.)',/)
   99002 format (10X,'               vectors  shells  max. R ',/,10X,         &
   '               ------------------------------')
   99003 format (10X,A,I7,2X,I6,2X,F9.5)
   99004 format (10X,'               ------------------------------',/)
   99005 format (10X,45('+'),/,13X,'generated ',A,' lattice vectors',/,10X,   &
   45('+'),/,10X,'shell Nvec    radius          x         y',/,10X,45('-'))
   99006 format (10X,I5,I5,F12.6,2X,2F10.5)
   99007 format (34X,2F10.5)
   99008 format (13X,42('-'))
   99009 format (10X,45('+'),/)
end subroutine LATTICE2D
