#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

! JUST FOR TESTING purposes
! replace by proper implementation
module TEST_lcutoff_mod
  implicit none
  private

  public :: initLcutoffNew

  integer,                protected, public :: lm_low
  double precision,       protected, public :: cutoff_radius
  integer,                protected, public :: lm_low2
  double precision,       protected, public :: cutoff_radius2
  integer, allocatable,   protected, public :: lmarray(:)
  integer,                protected, public :: cutoffmode
  logical,                protected, public :: DEBUG_dump_matrix = .false.
  integer,                protected, public :: num_untruncated
  integer,                protected, public :: num_truncated
  integer,                protected, public :: num_truncated2
  logical,                protected, public :: real_space_cutoff

  contains

  !----------------------------------------------------------------------------
  subroutine initLcutoffNew(trunc_zone, atom_ids, arrays)

    use Main2Arrays_mod, only: Main2Arrays
    use lcutoff_mod, only: calcCutoffarray
    use TruncationZone_mod, only: TruncationZone, createTruncationZone

    type(TruncationZone), intent(inout) :: trunc_zone
    type(Main2Arrays), intent(in) :: arrays
    integer, intent(in) :: atom_ids(:)

    integer :: lmmaxd, atomindex, ilocal, ii, ind, ios, num_local_atoms
    integer, allocatable :: lmarray_temp(:), lmarray_full(:)

    lmmaxd = arrays%lmmaxd

    allocate(lmarray_full(size(arrays%rbasis,2)))
    allocate(lmarray_temp(size(arrays%rbasis,2)))

    real_space_cutoff = .false.
    open(91, file='lcutoff', form='formatted', action='read', status='old', iostat=ios)
    if (ios == 0) then ! opening was successful
      read(91,*) cutoff_radius
      read(91,*) lm_low
      read(91,*) cutoffmode

      lm_low2 = -1
      cutoff_radius2 = 0.d0
      if (cutoffmode > 4) then
        cutoffmode = cutoffmode - 2
        ! cutoff-mode 5: iterative solver with 2 cutoffs,
        !             6: direct solver with 2 cutoffs
        read(91,*) cutoff_radius2
        read(91,*) lm_low2
        real_space_cutoff = .true.
      endif
      close(91)
    else
      write(6,*) 'No file "lcutoff" found, use defaults.' ! todo: convert to warning
      cutoff_radius =  9.d9 ! effectively infinity
      cutoff_radius2 = 9.d9
      lm_low  = 999 
      lm_low2 = 999
      cutoffmode = 4 ! 3:iterative solver, 4:full solver
    endif

    lmarray_full      = lmmaxd
    num_local_atoms = size(atom_ids)

    do ilocal = 1, num_local_atoms
      atomindex = atom_ids(ilocal)

      lmarray_temp = lmmaxd

      ! inner truncation zone
      call calcCutoffarray(lmarray_temp, arrays%rbasis, arrays%rbasis(:,atomindex), &
                           arrays%bravais, cutoff_radius, lm_low, .false.)

      ! outer truncation zone
      if (real_space_cutoff) then
        call calcCutoffarray(lmarray_temp, arrays%rbasis, arrays%rbasis(:,atomindex), &
                             arrays%bravais, cutoff_radius2, lm_low2, .false.)
      endif

      if (ilocal == 1) then
        lmarray_full = lmarray_temp
      else
        do ii = 1, size(lmarray_full)
          ! merge truncation zones of local atoms
          lmarray_full(ii) = max(lmarray_full(ii), lmarray_temp(ii))
        enddo
      endif
    enddo

    ! TODO: a bit confusing, is never deallocated
    call createTruncationZone(trunc_zone, lmarray_full, arrays)

    num_untruncated = 0
    num_truncated = 0
    num_truncated2 = 0
    do ii = 1, size(lmarray_full)
      if (lmarray_full(ii) == lmmaxd) then
        num_untruncated = num_untruncated + 1
      else if (lmarray_full(ii) == lm_low) then
        num_truncated = num_truncated + 1
      else if (lmarray_full(ii) == lm_low2) then
        num_truncated2 = num_truncated2 + 1
      endif
    enddo

    ind = 0
    do ii = 1, size(lmarray_full)
      if (lmarray_full(ii) > 0) then
        ind = ind + 1
      endif
    enddo

    allocate(lmarray(ind)) ! never deallocated - who cares

    ind = 0
    do ii = 1, size(lmarray_full)
      if (lmarray_full(ii) > 0) then
        ind = ind + 1
        lmarray(ind) = lmarray_full(ii)
        CHECKASSERT(trunc_zone%index_map(ii) == ind) ! TODO
      endif
    enddo

    deallocate(lmarray_temp, lmarray_full)
  endsubroutine

  !----------------------------------------------------------------------------
  subroutine cropGLLH(GLLH, lmmaxd, naclsd, naezd, lmarray, numn0, indn0)
    double complex, intent(inout) :: GLLH(LMMAXD,NACLSD*LMMAXD,NAEZD)
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: naclsd
    integer, intent(in) :: naezd
    integer, intent(in) :: lmarray(naezd)
    integer, intent(in) :: numn0(naezd)
    integer, intent(in) :: indn0(naezd,naclsd)

    integer ii, jj, lmmax1, lmmax2, lm1, lm2, clustersitelm
    double complex, parameter :: CZERO = (0.d0, 0.d0)

    do ii = 1, naezd
      do jj = 1, numn0(ii)

        lmmax1 = lmarray(ii)
        lmmax2 = lmarray(indn0(ii,jj))

        do lm2 = 1, lmmaxd
          do lm1 = 1, lmmaxd

            clustersitelm = (jj-1)*lmmaxd + lm2

            if (lm1 > lmmax1 .or. lm2 > lmmax2)  GLLH(lm1,clustersitelm,ii) = CZERO

          enddo ! lm1
        enddo ! lm2
      enddo ! jj
    enddo ! ii
    
  endsubroutine


!------------------------------------------------------------------------------
!> Generates matrix (\Delta T G_ref - 1) BUT PERFORMS L-CUTOFF BY TRUNCATING
!> T-MATRIX
!> on input: GLLH contains G_ref, on output: GLLH contains coefficient matrix
subroutine generateCoeffMatrixCROPPED(gllh, numn0, indn0, tmatll, naez, lmmaxd, naclsd, lmarray)
  integer, intent(in) :: lmmaxd
  integer, intent(in) :: naclsd
  integer, intent(in) :: naez
  double complex, intent(inout) :: gllh(lmmaxd,naclsd*lmmaxd,naez)
  integer, intent(in) :: indn0(naez,naclsd)
  integer, intent(in) :: numn0(naez)
  doublecomplex, intent(in) :: tmatll(lmmaxd,lmmaxd,naez)
  integer, intent(in) :: lmarray(:)

  double complex, parameter :: cone  = ( 1.d0,0.d0)
  double complex, parameter :: czero = ( 0.d0,0.d0)
  double complex :: tgh(lmmaxd)
  integer :: il1b
  integer :: il2b
  integer :: lm1
  integer :: lm2
  integer :: lm3
  integer :: site_index
  integer :: site_lm_index
  integer :: cluster_site_index
  integer :: cluster_site_lm_index
  integer :: lmmax1, lmmax2, lmmax3

  ! -------------- Calculation of (Delta_t * G_ref - 1) ---------------
  !
  !
  ! NUMN0(site_index) is the number of atoms in the reference cluster
  ! of atom/site 'site_index' (inequivalent atoms only!)
  ! INDN0 stores the index of the atom in the basis corresponding to
  ! the reference cluster atom (inequivalent atoms only!)
  ! -------------------------------------------------------------------

  !$omp parallel do private(site_index, site_lm_index, cluster_site_index, &
  !$omp                     cluster_site_lm_index, il1b, il2b, &
  !$omp                     lm1, lm2, lm3, tgh)
  do site_index = 1, naez
    il1b = lmmaxd*(site_index-1)
    do cluster_site_index = 1, numn0(site_index)

      lmmax2 = lmarray(indn0(site_index,cluster_site_index))
      lmmax1 = lmarray(site_index)
      lmmax3 = lmmax1

      do lm2 = 1, lmmaxd
        cluster_site_lm_index = lmmaxd*(cluster_site_index-1)+lm2
        il2b = lmmaxd*(indn0(site_index,cluster_site_index)-1)+lm2

        tgh = czero
        do lm1 = 1, lmmax1
          do lm3 = 1, lmmax3
            tgh(lm1) = tgh(lm1) + tmatll(lm1,lm3,site_index)*gllh(lm3,cluster_site_lm_index,site_index)
          enddo ! lm3
        enddo ! lm1

        do lm1 = 1, lmmaxd
          site_lm_index = il1b+lm1
          gllh(lm1,cluster_site_lm_index,site_index) = tgh(lm1)

          if (site_lm_index == il2b) then
            ! substract 1 only at the 'diagonal'
            gllh(lm1,cluster_site_lm_index,site_index) = &
            gllh(lm1,cluster_site_lm_index,site_index) - cone
          endif

          if (lm1 > lmmax1 .or. lm2 > lmmax2) gllh(lm1,cluster_site_lm_index,site_index) = czero

        enddo ! lm1
      enddo ! lm2

    enddo ! cluster_site_index
  enddo ! site_index
  !$omp endparallel do
endsubroutine

endmodule TEST_lcutoff_mod
