!> multiple scattering k-loop and symmetrisation
module kloopz1_mod

CONTAINS

    subroutine KLOOPZ1_new(GMATN,ALAT, &
    NOFKS, VOLBZ, BZKP, VOLCUB, &
    RR, GINP_LOCAL, &
    NSYMAT,DSYMLL, &
    TMATLL, atom_indices, &
    QMRBOUND, lmmaxd,  &
    nrd, trunc2atom_index, communicator, &
    iguess_data, cluster_info)

! **********************************************************************

! Only part of arrays for corresponding spin direction is passed
! (GMATN, TSST_LOCAL, DTDE_LOCAL, LLY_GRDT, TR_ALPH, GMATXIJ)
!
! NOFKS .. number of k-points, integer
! VOLBZ .. Brillouin zone volume, double
! BZKP ... k-points of used k-mesh ... dimension (3, KPOIBZ)
! VOLCUB . array of Brillouin zone integration weights for each k-point ... dimension (KPOIBZ)

! GINP_LOCAL ... reference Green's function
! DGINP ...      derivative of reference Green's function
! TSST_LOCAL ..  t-matrix

    use kkrmat_new_mod
    use TEST_lcutoff_mod
    use InitialGuess_mod
    use ClusterInfo_mod

    use jij_calc_mod, only: global_jij_data, symjij

    implicit none

    integer, intent(in) :: lmmaxd
    integer, intent(in) :: nrd
    !> mapping trunc. index -> atom index
    integer, intent(in) :: trunc2atom_index(:)
    integer, intent(in) :: communicator
    type(InitialGuess), intent(inout) :: iguess_data
    type(ClusterInfo), target :: cluster_info
    !     .. Parameters ..

    integer, parameter :: NSYMAXD = 48
    double complex, parameter :: CONE = ( 1.0D0,0.0D0)
    double complex, parameter :: CZERO = ( 0.0D0,0.0D0)
    !     ..
    !     .. Scalar Arguments ..
    !     ..
    double precision:: ALAT
    double precision::VOLBZ
    double precision::QMRBOUND
    integer::NOFKS
    integer::NSYMAT
    !     .. Array Arguments ..

    double complex :: DSYMLL(LMMAXD,LMMAXD,NSYMAXD)

    double complex :: GMATN(:,:,:)

    double complex :: GINP_LOCAL(:, :, :, :)
    double complex, intent(inout), dimension(:,:,:) :: TMATLL

    double precision::RR(3,0:NRD)
    double precision::BZKP(:,:)
    double precision::VOLCUB(:) ! dim kpoibz

    !     .. Local Scalars ..
    !     ..
    double complex :: TAUVBZ
    double precision:: RFCTOR
    integer::LM1
    integer::LM2

    integer::INFO    ! for LAPACK calls
    integer::symmetry_index

    integer :: ispin
    !     ..
    !     .. Local Arrays ..
    !     ..
    double complex :: GLL  (LMMAXD, LMMAXD)
    double complex, allocatable :: GS(:,:,:,:)

    !     effective (site-dependent) Delta_t^(-1) matrix
    double complex, allocatable ::  MSSQ (:, :, :)
    double complex :: work_array(LMMAXD, LMMAXD)    ! work array for LAPACK ZGETRI
    double complex ::       TPG (LMMAXD, LMMAXD)
    double complex ::        XC (LMMAXD, LMMAXD)      ! to store temporary matrix-matrix mult. result

    integer::         IPVT(LMMAXD)                 ! work array for LAPACK

    !-----------------------------------------------------------------------

    integer, intent(in) :: atom_indices(:)
    integer :: num_local_atoms
    integer :: ilocal

    integer :: memory_stat
    logical :: memory_fail

    num_local_atoms = size(atom_indices)

! -------------------------------------------------------------------
! Allocate Arrays
! -------------------------------------------------------------------
    memory_stat = 0
    memory_fail = .false.

    allocate(GS(LMMAXD, LMMAXD, NSYMAXD, num_local_atoms))
    if (memory_stat /= 0) memory_fail = .true.

    allocate(MSSQ(LMMAXD, LMMAXD, num_local_atoms))
    if (memory_stat /= 0) memory_fail = .true.

    if (memory_fail .eqv. .true.) then
      write(*,*) "KLOOPZ1: FATAL Error, failure to allocate memory."
      write(*,*) "       Probably out of memory."
      stop
    end if

!     RFCTOR=A/(2*PI) conversion factor to p.u.
    RFCTOR = ALAT/(8.D0*ATAN(1.0D0))           ! = ALAT/(2*PI)

    do ilocal = 1, num_local_atoms

      do LM2 = 1,LMMAXD
          do LM1 = 1,LMMAXD
              MSSQ(LM1,LM2,ilocal) =  TMATLL(LM1,LM2, atom_indices(ilocal))
          end do
      end do

  ! ---> inversion

  !     The (local) Delta_t matrix is inverted and stored in MSSQ

      call ZGETRF(LMMAXD,LMMAXD,MSSQ(:,:,ilocal),LMMAXD,IPVT,INFO)
      call ZGETRI(LMMAXD,MSSQ(:,:,ilocal),LMMAXD,IPVT,work_array, &
      LMMAXD*LMMAXD,INFO)

    end do !ilocal

!=======================================================================
!     Note: the actual k-loop is in kkrmat01 (it is not parallelized)
!     The integration over k is also performed in kkrmat01

    TAUVBZ = 1.D0/VOLBZ
    ! cutoffmode:
    ! 0 no cutoff, 1 T-matrix cutoff, 2 full matrix cutoff, &
    ! 3 T-matrix cutoff with new solver, 4 T-matrix cutoff with direct solver
    if (cutoffmode > 2) then
      call KKRMAT01_new(BZKP,NOFKS,GS,VOLCUB,TMATLL, &
      ALAT, NSYMAT, RR, GINP_LOCAL, atom_indices, &
      QMRBOUND, lmmaxd, trunc2atom_index, communicator, &
      iguess_data, cluster_info)
    else
      write(*,*) "cutoffmode<3 not supported."
      STOP
    endif
!-------------------------------------------------------- SYMMETRISE GLL


!      kkrmat01 returns GS (local) which contains NSYMAT -copies- (!)
!      (see 3rd index) of the scattering path operator
!      (already integrated over the irreducible wedge in k-space)
!      scattering path operator: ((Delta_T)^-1 - G_ref)^-1

!      All the symmetry operations are applied on GS and summed over
!     - the result is stored in GLL
!      Note: the symmetry operations apply on the (LL')-space

!------------------------------------------------------------------------------
    do ilocal = 1, num_local_atoms
!------------------------------------------------------------------------------

        do symmetry_index = 1,NSYMAT

        ! --->    GLL = sum(i=1,iumax)(tauvbz * DLL(i) * GS * DLL(i)^H)

            if ( symmetry_index == 1 ) then

            ! --->    ull(1) is equal to unity matrix

                call ZCOPY(LMMAXD*LMMAXD,GS(1,1,1,ilocal),1,GLL,1)
                call ZSCAL(LMMAXD*LMMAXD,TAUVBZ,GLL,1)

            else

            ! --->      tpg = tauvbz * DLL * GS

                call ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,TAUVBZ, &
                DSYMLL(1,1,symmetry_index),LMMAXD,GS(1,1,symmetry_index,ilocal),LMMAXD, &
                CZERO,TPG,LMMAXD)

            ! --->    GLL = GLL + TPG * DLL(i)^H
            !                           C  ! dsymll might be complex in REL case

                call ZGEMM('N','C',LMMAXD,LMMAXD,LMMAXD,CONE,TPG,LMMAXD, &
                DSYMLL(1,1,symmetry_index),LMMAXD,CONE,GLL,LMMAXD)
            end if

        end do
    !-------------------------------------------------------- IU = 1,NSYMAT


    ! --->  XC = Delta_t^(-1) * GLL

        call ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,MSSQ(:,:,ilocal), &
        LMMAXD,GLL,LMMAXD,CZERO,XC,LMMAXD)

    !       GLL is overwritten with the following expression:
    !       (for the in configuration space diagonal structural Green's
    !       function of the REAL system - already integrated over k)

    ! --->  GLL = - Delta_t^(-1) - Delta_t^(-1) * GLL * Delta_t^(-1)

    !       copy overwrite GLL with content of MSSQ
        call ZCOPY(LMMAXD**2,MSSQ(:,:,ilocal),1,GLL,1)


    !       GLL = (-1) * XC                     *    MSSQ     + (-1) * GLL
    !                    |                            |                 |
    !            Delta_t^-1 * scat. path op.      Delta_t^-1    Delta_t^-1

        call ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,-CONE,XC,LMMAXD, &
        MSSQ(:,:,ilocal),LMMAXD,-CONE,GLL,LMMAXD)


    ! --->  GMATN = GMATLL = GLL/RFCTOR...............rescaled and copied into output array

        do LM1 = 1,LMMAXD
            do LM2 = 1,LMMAXD
                GMATN(LM2,LM1,ilocal) = GLL(LM2,LM1)/RFCTOR
            end do
        end do

!------------------------------------------------------------------------------
    end do ! ilocal
!------------------------------------------------------------------------------

    if (global_jij_data%do_jij_calculation) then

      ispin = global_jij_data%active_spin

      call SYMJIJ( &
      ALAT,TAUVBZ, &
      NSYMAT,DSYMLL, &
      global_jij_data%NXIJ,global_jij_data%IXCP, &
      TMATLL,MSSQ, &
      global_jij_data%GSXIJ, &
      global_jij_data%GMATXIJ(:,:,:,ispin), &  ! Result
      cluster_info%naez_trc, lmmaxd, global_jij_data%nxijd)
    end if

! -------------------------------------------------------------------
! Deallocate Arrays
! -------------------------------------------------------------------

    deallocate(GS)
    deallocate(MSSQ)

  end subroutine KLOOPZ1_new

end module
