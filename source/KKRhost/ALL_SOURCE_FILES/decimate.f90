!-------------------------------------------------------------------------------
! SUBROUTINE: DECIMATE
!> @brief Decimation method
!> - Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine DECIMATE(GLLKE,NAEZ,TINVBUP,TINVBDOWN,VACFLAG,FACTL,NLBASIS,NRBASIS,&
   ALM,NDIM,LMMAXD)

   use global_variables

   implicit none

   !     .. Parameters ..
   !
   ! *********************************************************************
   ! * For KREL = 1 (relativistic mode)                                  *
   ! *                                                                   *
   ! *  NPOTD = 2 * NATYPD                                               *
   ! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
   ! *  NSPIND = 1                                                       *
   ! *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
   ! *          function, set up in the spin-independent non-relativstic *
   ! *          (l,m_l)-representation                                   *
   ! *                                                                   *
   ! *********************************************************************
   !
   integer, intent(in) :: ALM       !< NAEZ*LMMAXD
   integer, intent(in) :: NDIM      !< NPRINCD*LMMAXD
   integer, intent(in) :: NAEZ      !< Number of atoms in unit cell
   integer, intent(in) :: LMMAXD    !< (KREL+KORBIT+1)(LMAX+1)^2
   integer, intent(in) :: NLBASIS   !< Number of basis layers of left host (repeated units)
   integer, intent(in) :: NRBASIS   !< Number of basis layers of right host (repeated units)
   logical, dimension(2), intent(in) :: VACFLAG
   double complex, dimension(LMMAXD,LMMAXD), intent(in) :: FACTL
   double complex, dimension(LMMAXD,LMMAXD,*), intent(in) :: TINVBUP
   double complex, dimension(LMMAXD,LMMAXD,*), intent(in) :: TINVBDOWN
   double complex, dimension(ALM,ALM), intent(inout) :: GLLKE
   ! .. Local Scalars
   integer :: LDI1,LDI1T,LDI2,LDI2T,LM1,LM2,NLAYER,ICALL
   integer :: ICHCK,IHOST,II1,II2,IL1,IL2,IP1,IP1T,IP2,IP2T,ITERMAX
   double precision :: ERRMAX
   ! .. Local Arrays
   double complex, dimension(NDIM,NDIM) :: A1,AN,B1,BN,C1,CN,X1,XN
   ! ..
   data ICALL /0/
   ! .. External Subroutines
   external :: BOFM,SURFGF
   ! .. Save statement
   SAVE ICALL,NLAYER,ITERMAX,ERRMAX,ICHCK
   ! .. External Functions
   logical :: OPT
   external :: OPT
   ! .. Intrinsic Functions
   intrinsic :: MOD
   !
   ICALL = ICALL + 1
   !----------------------------------------------------------------------------
   if (ICALL.EQ.1) then
      NLAYER = NAEZ/NPRINCD
      ! Parameters for the "decimation" technique.
      ITERMAX = 300
      ERRMAX = 1.0D-180
      ICHCK = 1
   end if
   !----------------------------------------------------------------------------
   if ( .NOT.VACFLAG(1) ) then
      !-------------------------------------------------------------------------
      ! Get the matrix B1
      !-------------------------------------------------------------------------
      call BOFM(1,1,B1,NDIM,GLLKE,ALM)

      ! Now Subtract t-mat of left host
      do IP1 = 1,NPRINCD
         IHOST = MOD(IP1-1,NLBASIS) + 1
         do LM1 = 1,LMMAXD
            do LM2 = 1,LMMAXD
               IL1 = LMMAXD* (IP1-1) + LM1
               IL2 = LMMAXD* (IP1-1) + LM2
               B1(IL1,IL2) = (B1(IL1,IL2)-TINVBUP(LM1,LM2,IHOST))
            end do
         end do
      end do

      call BOFM(1,2,C1,NDIM,GLLKE,ALM)
      call BOFM(2,1,A1,NDIM,GLLKE,ALM)

      ! It performs the 'space decimation' iterative procedure.
      call SURFGF(NDIM,A1,B1,C1,X1,ITERMAX,ERRMAX,ICHCK,LMMAXD)
      ! Adds to the matrix GLLKE the elements that couples the
      ! interface to the two half-spaces.
      do IP1 = 1,NPRINCD
         do IP2 = 1,NPRINCD
            II1 = IP1
            II2 = IP2
            do LM1 = 1,LMMAXD
               do LM2 = 1,LMMAXD
                  LDI1 = LMMAXD* (IP1-1) + LM1
                  IL1 = LMMAXD* (II1-1) + LM1
                  LDI2 = LMMAXD* (IP2-1) + LM2
                  IL2 = LMMAXD* (II2-1) + LM2
                  GLLKE(IL1,IL2) = GLLKE(IL1,IL2) - X1(LDI1,LDI2)
               end do
            end do
         end do
      end do
   end if

   if ( .NOT.VACFLAG(2) ) then

      !  If 'ONEBULK' is activated then it calculates the xn decimated element
      !  from the x1 element: this is just in the case of equal bulks on the

      if ( .NOT.OPT('ONEBULK ') ) then

         !----------------------------------------------------------------------
         ! Get the matrix BN
         !----------------------------------------------------------------------
         call BOFM(NLAYER,NLAYER,BN,NDIM,GLLKE,ALM)

         ! Now Substract t-mat right host
         ! Notes : the indexing is easier like that
         do IP1 = 1,NPRINCD
            IHOST = NRBASIS - MOD(IP1,NRBASIS)
            IHOST = MOD(IP1-1,NRBASIS) + 1
            do LM1 = 1,LMMAXD
               do LM2 = 1,LMMAXD
                  IL1 = LMMAXD* (IP1-1) + LM1
                  IL2 = LMMAXD* (IP1-1) + LM2
                  BN(IL1,IL2) = (BN(IL1,IL2)-TINVBDOWN(LM1,LM2,IHOST))
               end do
            end do
         end do

         call BOFM(NLAYER,NLAYER-1,AN,NDIM,GLLKE,ALM)
         call BOFM(NLAYER-1,NLAYER,CN,NDIM,GLLKE,ALM)

         ! It performs the 'space decimation' iterative procedure.
         call SURFGF(NDIM,CN,BN,AN,XN,ITERMAX,ERRMAX,ICHCK,LMMAXD)
         !
      else
         !
         do IP1 = 1,NPRINCD
            do IP2 = 1,NPRINCD
               IP1T = (NPRINCD+1) - IP2
               IP2T = (NPRINCD+1) - IP1
               do LM1 = 1,LMMAXD
                  do LM2 = 1,LMMAXD
                     LDI1 = LMMAXD* (IP1-1) + LM1
                     LDI2 = LMMAXD* (IP2-1) + LM2
                     LDI1T = LMMAXD* (IP1T-1) + LM2
                     LDI2T = LMMAXD* (IP2T-1) + LM1
                     XN(LDI1T,LDI2T) = FACTL(LM1,LM2)*X1(LDI1,LDI2)
                  end do
               end do
            end do
         end do
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !             Added on 1.02.2000
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Adds to the matrix GLLKE the elements that couples the
      ! interface to the two half-spaces.
      do IP1 = 1,NPRINCD
         do IP2 = 1,NPRINCD
            II1 = (NLAYER-1)*NPRINCD + IP1
            II2 = (NLAYER-1)*NPRINCD + IP2
            do LM1 = 1,LMMAXD
               do LM2 = 1,LMMAXD
                  LDI1 = LMMAXD* (IP1-1) + LM1
                  IL1 = LMMAXD* (II1-1) + LM1
                  LDI2 = LMMAXD* (IP2-1) + LM2
                  IL2 = LMMAXD* (II2-1) + LM2
                  GLLKE(IL1,IL2) = GLLKE(IL1,IL2) - XN(LDI1,LDI2)
               end do
            end do
         end do
      end do
   end if

   return

end subroutine DECIMATE
