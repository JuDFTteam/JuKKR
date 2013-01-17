module kloopz1_mod

CONTAINS

    subroutine KLOOPZ1_new(GMATN,ALAT,ITER, &
    NAEZ,NOFKS,VOLBZ,BZKP,VOLCUB,CLS, &
    NACLS,RR,EZOA,ATOM,GINP_LOCAL,DGINP, &
    NSYMAT,DSYMLL, &
    TMATLL,DTDE_LOCAL, &
    NUMN0,INDN0,atom_indices, &
    PRSC,EKM,NOITER, &
    QMRBOUND,IGUESS,BCP, &
    NXIJ,XCCPL,IXCP,ZKRXIJ, &            ! > input
    LLY_GRDT,TR_ALPH,GMATXIJ, &          ! < output
    communicator, comm_size, &           ! > input
    ! new parameters after inc.p removal
    lmmaxd, naclsd, nclsd, xdim, ydim, zdim, natbld, LLY, &
    nxijd, nguessd, kpoibz, nrd, ekmd)

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
    use kkrmat_mod !TODO: remove
    use TEST_lcutoff_mod !TODO: remove
    implicit none
    include 'mpif.h'

    integer, intent(in) :: communicator
    integer, intent(in) :: comm_size

    integer, intent(in) :: lmmaxd
    integer, intent(in) :: naclsd  ! max. number of atoms in reference cluster
    integer, intent(in) :: nclsd   ! number of reference clusters
    integer, intent(in) :: xdim
    integer, intent(in) :: ydim
    integer, intent(in) :: zdim
    integer, intent(in) :: natbld  ! number of atoms in preconditioning blocks
    integer, intent(in) :: LLY
    integer, intent(in) :: nxijd   ! max. number of atoms in cluster for exchange coupling-calculation
    integer, intent(in) :: nguessd
    integer, intent(in) :: kpoibz
    integer, intent(in) :: nrd
    integer, intent(in) :: ekmd

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
    integer::ITER
    integer::NOFKS
    integer::NXIJ
    integer::NAEZ
    integer::IGUESS
    integer::BCP
    integer::EKM
    integer::NOITER
    integer::NSYMAT
    double complex :: LLY_GRDT
    double complex :: TR_ALPH
    double complex :: BZTR2
    !     ..
    !     .. Array Arguments ..
    !     ..

    !----- Initial Guess arrays-----------------------------------------------
    complex::          PRSC(NGUESSD*LMMAXD,EKMD)
    !-------------------------------------------------------------------------

    double complex :: DSYMLL(LMMAXD,LMMAXD,NSYMAXD)

    double complex :: GMATN(:,:,:)

    double complex :: DGINP(LMMAXD,LMMAXD, NACLSD, NCLSD)
    double complex :: GINP_LOCAL(LMMAXD, LMMAXD, NACLSD, NCLSD)
    double complex :: GMATXIJ   (LMMAXD, LMMAXD, NXIJD)
    double complex, intent(inout), dimension(lmmaxd,lmmaxd,naez) :: TMATLL
    double complex :: DTDE_LOCAL(LMMAXD, LMMAXD)

    double precision::RR(3,0:NRD)
    double precision::ZKRXIJ(48,3,NXIJD)
    double precision::BZKP(3,KPOIBZ)
    double precision::VOLCUB(KPOIBZ)
    integer:: ATOM(NACLSD,naez)
    integer:: CLS(naez)
    integer:: EZOA(NACLSD,naez)
    integer:: NACLS(naez)
    !     ..
    integer::NUMN0(NAEZ)
    integer::INDN0(NAEZ,NACLSD)
    integer::IXCP(NXIJD)
    !     ..
    !     .. Local Scalars ..
    !     ..
    double complex :: TAUVBZ
    double precision:: RFCTOR
    integer::LM1
    integer::LM2

    integer::INFO    ! for LAPACK calls
    integer::IERR
    integer::symmetry_index
    !     ..
    !     .. Local Arrays ..
    !     ..

    double complex :: GLL  (LMMAXD, LMMAXD)
    double complex, allocatable :: GS(:,:,:,:)

    !double complex :: GSXIJ(LMMAXD,LMMAXD,NSYMAXD,NXIJD)
    double complex, allocatable, dimension(:,:,:,:) :: GSXIJ

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

    logical::XCCPL

    external ZCOPY,ZAXPY,ZGETRF,ZGETRS,ZGETRI,ZGEMM,ZSCAL
!     ..
!     .. Intrinsic Functions ..
    intrinsic ATAN
!     ..
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

    allocate(GSXIJ(LMMAXD, LMMAXD, NSYMAXD, NXIJD))
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
    ! 0 no cutoff, 1 T-matrix cutoff, 2 full matrix cutoff, 3 T-matrix cutoff with new solver, 4 T-matrix cutoff with direct solver
    if (cutoffmode > 2) then
      call KKRMAT01_new(BZKP,NOFKS,GS,VOLCUB,TMATLL,MSSQ(:,:,1), &
      ITER, &
      ALAT,NSYMAT,NAEZ,CLS,NACLS,RR,EZOA,ATOM, &
      GINP_LOCAL,DGINP, &
      NUMN0,INDN0, atom_indices, &
      PRSC, &
      EKM,NOITER, &
      QMRBOUND,IGUESS,BCP, &
      DTDE_LOCAL, &
      GSXIJ, &
      NXIJ,XCCPL,IXCP,ZKRXIJ, &
      BZTR2, &
      communicator, comm_size, &
      lmmaxd, naclsd, nclsd, xdim, ydim, zdim, natbld, LLY, &
      nxijd, nguessd, kpoibz, nrd, ekmd)
    else
      call KKRMAT01(BZKP,NOFKS,GS(:,:,:,1),VOLCUB,TMATLL,MSSQ(:,:,1), &
      ITER, &
      ALAT,NSYMAT,NAEZ,CLS,NACLS,RR,EZOA,ATOM, &
      GINP_LOCAL,DGINP, &
      NUMN0,INDN0,atom_indices(1), &
      PRSC, &
      EKM,NOITER, &
      QMRBOUND,IGUESS,BCP, &
      DTDE_LOCAL, &
      GSXIJ, &
      NXIJ,XCCPL,IXCP,ZKRXIJ, &
      BZTR2, &
      communicator, comm_size, &
      lmmaxd, naclsd, nclsd, xdim, ydim, zdim, natbld, LLY, &
      nxijd, nguessd, kpoibz, nrd, ekmd)
    endif

    ! TODO: move out
    !=========== Lloyd's Formula =====================================

    if(LLY == 1)  then
      BZTR2 = BZTR2*NSYMAT/VOLBZ + TR_ALPH

      call MPI_ALLREDUCE(BZTR2,LLY_GRDT,1,MPI_DOUBLE_COMPLEX,MPI_SUM, &
                         communicator,IERR)

    endif
    !========== END Lloyd's Formula ==================================


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

!================================
    if (XCCPL) then
!================================
    ! todo: loop over local atoms
    if (size(atom_indices) /= 1) then
      write(*,*) "XCCPL and num_local_atoms>1 not possible."
      STOP
    endif
        call SYMJIJ( &
        ALAT,TAUVBZ, &
        NSYMAT,DSYMLL, &
        NXIJ,IXCP, &
        TMATLL,MSSQ(:,:,1), &
        GSXIJ, &
        GMATXIJ, &
        naez, lmmaxd, nxijd)

!================================
    endif
!================================

! -------------------------------------------------------------------
! Deallocate Arrays
! -------------------------------------------------------------------

    deallocate(GSXIJ)
    deallocate(GS)
    deallocate(MSSQ)

  end subroutine KLOOPZ1_new

end module
