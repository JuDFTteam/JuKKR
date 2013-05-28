module SingleSiteRef_mod

CONTAINS

!**********************************************************************
!> @param E complex energy
!> @param CLEB array of Gaunt coefficients
!> @param ICLEB index array for Gaunt coefficients
!> @param LOFLM array that maps (lm)-index to L-index
!> @param IEND ???
!> @param TREFLL reference T-matrix
!> @param DTREFLL derivative of reference T-matrix
!> @param RATOM real space positions of atoms in ref. cluster
!> @param NATOM number of atoms in reference cluster
!> @param ALAT length of unit vector in Bohr
!> @param GREF0 TODO
!> @param DGDEOUT ??? energy derivative of Green's function
!> @param LLY_G0TR Trace(M^-1 dM/dE) with M = (1 - g0 \Delta t_ref)
!> @param lmaxd angular momentum cutoff (it would be better to rewrite routine to pass lmmaxd)
!> @param naclsd dimension array: maximal number of atoms in reference clusters
!> @param ncleb number of Gaunt coefficients in CLEB
!> @param LLY do Lloyd's formula calculations 1=yes/0=no
subroutine GLL95(E,CLEB,ICLEB,LOFLM,IEND,TREFLL,DTREFLL, &
                 RATOM,NATOM,ALAT,GREF0,DGDEOUT, &
                 LLY_G0TR, &
!                new input parameters after inc.p removal
                 lmaxd, naclsd, ncleb, LLY)

! **********************************************************************
!
!     solution of the DYSON equation for a cluster of potentials
!     (TREFLL) centered at positions RATOM in free space,
!
! ----------------------------------------------------------------------

  implicit none

  integer, intent(in) :: lmaxd
  integer, intent(in) :: naclsd
  integer, intent(in) :: ncleb
  integer, intent(in) :: LLY

  !
  !     INTEGER LMGF0D,NGD
  !     PARAMETER (LMGF0D= (LMAXD+1)**2,NGD=LMGF0D*NACLSD)
  !     INTEGER LLYNGD
  !     PARAMETER (LLYNGD=LLY*(LMGF0D*NACLSD-1)+1)

  double complex :: CONE
  double complex :: CZERO
  parameter (CONE= (1.D0,0.D0),CZERO= (0.D0,0.D0))
  !     ..
  !     .. Scalar Arguments ..
  double complex :: E
  double complex :: LLY_G0TR
  double precision :: ALAT

  integer :: NATOM
  integer :: INFO
  !     ..
  !     .. Array Arguments ..
  integer :: IPVT(NACLSD*(LMAXD+1)**2)

  !     DOUBLE COMPLEX GREF0(NGD,LMGF0D)
  double complex :: GREF0(NACLSD*(LMAXD+1)**2,(LMAXD+1)**2)

  double complex :: TREFLL((LMAXD+1)**2,(LMAXD+1)**2, naclsd)
  double complex :: DTREFLL((LMAXD+1)**2,(LMAXD+1)**2, naclsd)

  double complex :: DGDEOUT(LLY*(NACLSD*(LMAXD+1)**2-1)+1, &
                            (LMAXD+1)**2)

  double precision :: CLEB(ncleb)
  integer :: ICLEB(NCLEB,3)
  integer :: LOFLM(:)
  integer :: IEND

  double precision :: RATOM(3,*)  ! first dim: 3

  !     ..
  !     .. Local Scalars ..
  integer :: N1
  integer :: N2
  integer :: NDIM
  integer :: site_lm_index2

! ---------------------------------------------------------------------
!     The following arrays can be very large (> 60 MB in typical cases)
!     therefore they are allocated on the heap
! ---------------------------------------------------------------------

  double complex, allocatable, dimension(:,:) :: GREF
  double complex, allocatable, dimension(:,:) :: GTREF

  double complex, allocatable, dimension(:,:) ::  DGTDE
  double complex, allocatable, dimension(:,:) :: DGTDE0
  double complex, allocatable, dimension(:,:) :: DGDE

  !     ..
  !     .. External Subroutines ..
  external GFREE,GREFSY,ZCOPY,ZGEMM

  !     ..
  !     .. Intrinsic Functions ..
  intrinsic ABS,DBLE
  !     ..

  integer :: memory_stat
  logical :: memory_fail

  integer :: LMGF0D
  integer :: NGD
  integer :: LMMAXD

  LMMAXD = (LMAXD+1)**2
  LMGF0D = LMMAXD
  NGD = LMMAXD*NACLSD

  !     Allocate arrays
  memory_stat = 0
  memory_fail = .false.

  allocate(GREF(NGD,NGD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  allocate(GTREF(NGD,LMMAXD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  allocate(DGTDE(LLY*(LMMAXD*NACLSD-1)+1,LMMAXD), stat = memory_stat)  !FIX: workaround, need to pass it to GREFSY
  if (memory_stat /= 0) memory_fail = .true.

  if (LLY == 1) then
    !allocate(DGTDE(NGD,LMMAXD), stat = memory_stat) ! Wrong: need to allocate in any case because
    !if (memory_stat /= 0) memory_fail = .true.      ! it has to be passed to GREFSY

    allocate(DGTDE0(NGD,NGD), stat = memory_stat)
    if (memory_stat /= 0) memory_fail = .true.

    allocate(DGDE(NGD,NGD), stat = memory_stat)
    if (memory_stat /= 0) memory_fail = .true.
  end if

  if (memory_fail .eqv. .true.) then
    write(*,*) "GLL95: FATAL Error, failure to allocate memory."
    write(*,*) "       Probably out of memory."
    stop
  end if

  NDIM = LMMAXD*NATOM     ! Warning: NDIM can be smaller than NGD=LMMAXD*NACLSD
  call calcFreeGreens(GREF, E, LMMAXD, NATOM, RATOM, ALAT, CLEB, ICLEB, ncleb, IEND, LOFLM, .false.)

  if (LLY==1) then
    call calcFreeGreens(DGDE, E, LMMAXD, NATOM, RATOM, ALAT, CLEB, ICLEB, ncleb, IEND, LOFLM, .true.)
  endif

! construct right hand side of linear equation system for GREFSY
! the first LMMAXD columns of GREF are copied into GREF0
! GREF0 then contains g0^{(1)N'}_{LL'}, the free space structural
! Green's function for the central cluster atom (E.R.)
! --------------------------------------------------------------
   call ZCOPY(NGD*LMMAXD,GREF,1,GREF0,1)
! --------------------------------------------------------------

   if (LLY==1) then

     do N2 = 1,NATOM
       site_lm_index2 = (N2-1)*LMMAXD + 1

       ! -dG_0/dE * \Delta t_ref    -- stored in GTREF
       call ZGEMM('N','N',NDIM,LMMAXD,LMMAXD,-CONE,DGDE(1,site_lm_index2),NGD, &
                  TREFLL(1,1,N2), LMMAXD, &
                  CZERO,GTREF,NGD)

       !   - G_0 * d(\Delta t_ref)/dE + GTREF  -- stored again in GTREF
       ! = -dG_0/dE * \Delta t_ref - G_0 * d(\Delta t_ref)/dE

       call ZGEMM('N','N',NDIM,LMMAXD,LMMAXD,-CONE,GREF(1,site_lm_index2),NGD, &
                  DTREFLL(1,1,N2), LMMAXD, &
                  CONE,GTREF,NGD)

       ! copy GTREF to DGTDE0 - GTREF is reused
       call ZCOPY(NGD*LMMAXD,GTREF,1,DGTDE0(1,site_lm_index2),1)
     end do
   end if  ! (LLY==1)

   do N2 = 1,NATOM
     site_lm_index2 = (N2-1)*LMMAXD + 1

     ! -G_ref * \Delta t_ref  -- stored in GTREF
     call ZGEMM('N','N',NDIM,LMMAXD,LMMAXD,-CONE,GREF(1,site_lm_index2),NGD, &
     TREFLL(1,1,N2), LMMAXD, &
     CZERO,GTREF,NGD)

     call ZCOPY(NGD*LMMAXD,GTREF,1,GREF(1,site_lm_index2),1)
   end do

   if (LLY==1) then
     do N2 = 1,LMMAXD
       do N1 = 1,NGD
         DGTDE(N1,N2) = DGTDE0(N1,N2)
       enddo
     enddo
   endif


   ! Solve Dyson-Equation for reference system
   ! Solves (1 - g0 \Delta t) G_ref = g0 for G_ref.
   LLY_G0TR = CZERO
   call GREFSY(GREF,GREF0,IPVT,NDIM,DGTDE, &
               LLY_G0TR, &
               NACLSD, LMMAXD, LLY)

   if (LLY==1) then

     call ZGEMM('N','N',NDIM,LMMAXD,NDIM,-CONE,DGTDE0,NGD, &
                GREF0,NGD,CONE,DGDE,NGD)

     call ZGETRS('N',NDIM,LMMAXD,GREF,NGD,IPVT,DGDE,NGD,INFO)

     do N2 = 1,LMMAXD
       do N1 = 1,NGD
         DGDEOUT(N1,N2) = DGDE(N1,N2)
       enddo
     enddo

   endif

   !     Deallocate arrays

   deallocate(GREF)
   deallocate(GTREF)

   if (LLY == 1) then
     deallocate(DGTDE)
     deallocate(DGTDE0)
     deallocate(DGDE)
   end if

 end subroutine GLL95

 !> Calculates the Free-Space Green-Function or the Derivative of the Free-Space Green-Function
 !> @param[out] greenFree   the free space Green-Function
 !> @param      NATOM       number of atoms in reference cluster
 !> @param[in]  derivative  .false. = calculate Free-Space-Greens Function
 !>                         .true.  = calculate Derivative of Free-Space-Greens Function
 subroutine calcFreeGreens(greenFree, energy, lmmaxd, NATOM, RATOM, ALAT, CLEB, ICLEB, ncleb, IEND, LOFLM, derivative)
   use kkr_helpers_mod
   implicit none
   integer, intent(in) :: ncleb
   double precision :: ALAT
   double precision :: CLEB(:)
   double complex :: energy
   double complex :: greenFree(:,:)
   integer :: ICLEB(NCLEB,3)
   integer :: lmmaxd
   integer :: LOFLM(:)
   integer :: NATOM
   double precision :: RATOM(3,*)
   logical, intent(in) :: derivative

   ! local
   double complex, parameter :: CZERO = (0.0d0, 0.0d0)
   double precision :: RDIFF(3)
   integer :: ind
   integer :: LM1
   integer :: LM2
   integer :: IEND
   double complex :: GLL(lmmaxd,lmmaxd) ! automatic array
   integer :: N1
   integer :: N2
   integer :: site_lm_index1
   integer :: site_lm_index2
   integer :: lmaxd

   lmaxd = lmmaxToLmax(lmmaxd)

   !
   ! ---> construct free Green's function
   ! The free space structural Green's function g0 is a matrix of dimension LMMAXD x LMMAXD
   ! (for a certain pair of reference cluster atoms N and N')
   ! It is calculated in routine GFREE and stored in GLL
   !
   ! Then the matrix GREF^{NN'}_{LL'} is constructed
   ! The blocks N /= N' contain g0,
   ! the other N=N' blocks are set to zero (Green's function is not defined for r=r')
   ! (E.R.)
   do N1 = 1,NATOM
     do N2 = 1,NATOM
       do ind = 1,3
         !            RDIFF(I) = (RATOM(I,N1) - RATOM(I,N2))*ALAT
         !           changed P.Z. 4.7.97
         RDIFF(ind) = - (RATOM(ind,N1)-RATOM(ind,N2))*ALAT
       end do

       if (N1/=N2) then

         if (derivative) then
           call DGFREE(RDIFF,Energy,GLL,CLEB,ICLEB,LOFLM,IEND, lmaxd, ncleb)
         else
           call GFREE(RDIFF,Energy,GLL,CLEB,ICLEB,LOFLM,IEND, lmaxd, ncleb)
         end if

         do LM2 = 1,lmmaxd
           site_lm_index2 = (N2-1)*lmmaxd + LM2
           do LM1 = 1,lmmaxd
             site_lm_index1 = (N1-1)*lmmaxd + LM1
             greenFree(site_lm_index1,site_lm_index2) = GLL(LM1,LM2)
           end do
         end do
       else
         do LM2 = 1,lmmaxd
           site_lm_index2 = (N2-1)*lmmaxd + LM2
           do LM1 = 1,lmmaxd
             site_lm_index1 = (N1-1)*lmmaxd + LM1
             greenFree(site_lm_index1,site_lm_index2) = CZERO
           end do
         end do

       end if  ! (N1 /= N2)
     end do
    end do
 end subroutine


 end module
