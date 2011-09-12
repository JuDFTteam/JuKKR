module inputcard_reader

  implicit none

  type InputcardParams
    double precision :: ALAT
    double precision :: ABASIS
    double precision :: BBASIS
    double precision :: CBASIS

    integer :: NCLS
    double precision :: E1
    double precision :: E2
    double precision :: TK
    integer :: NPOL
    integer :: NPNT1
    integer :: NPNT2
    integer :: NPNT3

    integer :: NSTEPS ! number of self consistent field steps
    integer :: IMIX
    double precision :: MIXING
    double precision :: QBOUND
    double precision :: FCM
    integer :: ITDBRY

    integer :: NAEZ   ! Number of atoms in unit cell
    integer :: IRM

    integer :: NREF
    integer :: ICST
    integer :: IFILE
    integer :: IPE
    integer :: IPF
    integer :: IPFE
    integer :: KHFELD
    integer :: KPRE
    integer :: KTE    ! Calculate total energy (0/1)
    integer :: KVMAD
    integer :: KVREL
    integer :: KXC
    integer :: LMAX
    integer :: LMPOT
    integer :: LPOT 
    integer :: NSPIN
    integer :: ISHIFT
    integer :: INTERVX ! Division of the Brillouin zone in x direction
    integer :: INTERVY
    integer :: INTERVZ
    double precision :: HFIELD
    double precision :: VCONST

    character(len=40) :: I13
    character(len=40) :: I19
    double precision :: RCUTZ
    double precision :: RCUTXY
    double precision :: RCUTJIJ
    logical :: JIJ
    double precision :: RCUTTRC
    logical :: LDAU

    integer :: KFORCE  ! calculate forces (0/1)
    integer :: IGUESS
    integer :: BCP
    double precision :: QMRBOUND
    logical :: LCARTESIAN
    double precision :: RMAX
    double precision :: GMAX
  end type

   type InputcardArrays
     ! location of basis atoms
     double precision, dimension(:,:), pointer :: RBASIS
     integer, dimension(:), pointer :: CLS
     integer, dimension(:), pointer :: IRNS
     integer, dimension(:), pointer :: NTCELL
     double precision, dimension(:), pointer :: Z
     integer, dimension(:), pointer :: REFPOT
     ! initial spin polarisation
     integer, dimension(:), pointer :: INIPOL
     double precision, dimension(:), pointer :: RMTREF

     double precision, dimension(2) :: VBC
   end type

  CONTAINS

  subroutine createInputcardArrays(input_arrays, NAEZD, NREFD)
    implicit none
    type (InputcardArrays), intent(inout) :: input_arrays
    integer, intent(in) :: NAEZD, NREFD

    allocate(input_arrays%RBASIS(3,NAEZD))
    allocate(input_arrays%CLS(NAEZD))
    allocate(input_arrays%IRNS(NAEZD))
    allocate(input_arrays%NTCELL(NAEZD))
    allocate(input_arrays%Z(NAEZD))
    allocate(input_arrays%REFPOT(NAEZD))
    allocate(input_arrays%INIPOL(NAEZD))
    allocate(input_arrays%RMTREF(NREFD))

  end subroutine createInputcardArrays

  subroutine destroyInputcardArrays(input_arrays)
    implicit none
    type (InputcardArrays), intent(inout) :: input_arrays

    deallocate(input_arrays%RBASIS)
    deallocate(input_arrays%CLS)
    deallocate(input_arrays%IRNS)
    deallocate(input_arrays%NTCELL)
    deallocate(input_arrays%Z)
    deallocate(input_arrays%REFPOT)
    deallocate(input_arrays%INIPOL)
    deallocate(input_arrays%RMTREF)

  end subroutine destroyInputcardArrays

  subroutine readInput(input_params, input_arrays)

  use inc_p_replace
  implicit none

  type (InputcardParams), intent(inout) :: input_params
  type (InputcardArrays), intent(inout) :: input_arrays

  call inc_p_replace_init()

  call RINPUT99( &
  & input_params%ALAT, &
  & input_arrays%RBASIS, &
  & input_params%ABASIS, &
  & input_params%BBASIS, &
  & input_params%CBASIS, &
  & input_arrays%CLS, &
  & input_params%NCLS, &
  & input_params%E1, &
  & input_params%E2, &
  & input_params%TK, &
  & input_params%NPOL, &
  & input_params%NPNT1, &
  & input_params%NPNT2, &
  & input_params%NPNT3, &
  & input_params%NSTEPS, &
  & input_params%IMIX, &
  & input_params%MIXING, &
  & input_params%QBOUND, &
  & input_params%FCM, &
  & input_params%ITDBRY, &
  & input_arrays%IRNS, &
  & input_arrays%NTCELL, &
  & input_params%NAEZ, &
  & input_params%IRM, &
  & input_arrays%Z, &
  & input_params%NREF, &
  & input_params%ICST, &
  & input_params%IFILE, &
  & input_params%IPE, &
  & input_params%IPF, &
  & input_params%IPFE, &
  & input_params%KHFELD, &
  & input_params%KPRE, &
  & input_params%KTE, &
  & input_params%KVMAD, &
  & input_params%KVREL, &
  & input_params%KXC, &
  & input_params%LMAX, &
  & input_params%LMPOT, &
  & input_params%LPOT , &
  & input_params%NSPIN, &
  & input_arrays%REFPOT, &
  & input_params%ISHIFT, &
  & input_params%INTERVX, &
  & input_params%INTERVY, &
  & input_params%INTERVZ, &
  & input_params%HFIELD, &
  & input_arrays%VBC, &
  & input_params%VCONST, &
  & input_arrays%INIPOL, &
  & input_params%I13, &
  & input_params%I19, &
  & input_params%RCUTZ, &
  & input_params%RCUTXY, &
  & input_params%RCUTJIJ, &
  & input_params%JIJ, &
  & input_params%RCUTTRC, &
  & input_params%LDAU, &
  & input_arrays%RMTREF, &
  & input_params%KFORCE, &
  & input_params%IGUESS, &
  & input_params%BCP, &
  & input_params%QMRBOUND, input_params%LCARTESIAN, &
  & input_params%RMAX, input_params%GMAX, &
  & LMAXD, IRNSD, TRC, LPOTD, NSPIND, &
  & IRMD, NAEZD)

  end subroutine readinput

end module inputcard_reader

! program test_it
! 
!   use inputcard_reader
! 
!   implicit none
! 
!   integer :: NAEZD
!   integer :: NREFD
! 
!   type (InputcardParams) :: input_params
!   type (InputcardArrays) :: input_arrays
! 
!   NAEZD = 4
!   NREFD = 1
! 
!   call createInputcardArrays(input_arrays, NAEZD, NREFD)
!   call readInput(input_params, input_arrays)
! 
!   write(*,*) "BASIS ====="
!   write(*,*) input_arrays%RBASIS
!   write(*,*) "Nuclear charge ======"
!   write(*,*) input_arrays%Z
!   write(*,*) "IRNS ======"
!   write(*,*) input_arrays%IRNS
!   write(*,*) "INIPOL ======"
!   write(*,*) input_arrays%INIPOL
! 
!   call destroyInputcardArrays(input_arrays)
! 
! end program
