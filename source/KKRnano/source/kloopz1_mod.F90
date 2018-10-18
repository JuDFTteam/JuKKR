!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> multiple scattering k-loop and symmetrisation
module kloopz1_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: kloopz1

  double complex, parameter :: CONE=(1.d0,0.d0), ZERO=(0.d0,0.d0)
  
  contains

  subroutine kloopz1(GmatN, solv, op, precond, alat, NofKs, volBZ, Bzkp, k_point_weights, rr, Ginp_local, &
                     dsymLL, tmatLL, global_atom_id, communicator, xTable, iguess_data, ienergy, ispin, &
                     tr_alph, Lly_grdt, global_atom_idx_Lly, Lly, & ! LLY 
                     korbit, & ! NOCO
                     solver_type, kpoint_timer, kernel_timer)

! only part of arrays for corresponding spin direction is passed
! (GmatN, tsst_local, dtde_local, Lly_grdt, tr_alph, gmatxij)
!
! NofKs .. number of k-points, integer
! volBZ .. brillouin zone volume, double
! Bzkp ... k-points of used k-mesh ... dimension (3, kpoibz)
! k_point_weights . array of brillouin zone integration weights for each k-point ... dimension (kpoibz)

! Ginp_local ... reference green's function
! dginp ...      derivative of reference green's function
! tsst_local ..  t-matrix

    use KKRmat_mod, only: MultipleScattering
    use InitialGuess_mod, only: InitialGuess
    use IterativeSolver_mod, only: IterativeSolver
    use BCPOperator_mod, only: BCPOperator
    use KKROperator_mod, only: KKROperator
    use ExchangeTable_mod, only: ExchangeTable 
    use Constants_mod, only: pi
    use TimerMpi_mod, only: TimerMpi
    use jij_calc_mod, only: global_jij_data, symjij
    
    type(IterativeSolver), intent(inout) :: solv
    type(KKROperator), intent(inout) :: op
    type(BCPOperator), intent(inout) :: precond
    integer, intent(in) :: global_atom_id(:) !> mapping trunc. index -> atom index
    integer, intent(in) :: communicator
    type(ExchangeTable), intent(in) :: xTable
    type(InitialGuess), intent(inout) :: iguess_data
    integer, intent(in) :: ienergy, ispin
    double precision, intent(in) :: alat
    double complex, intent(in) :: dsymLL(:,:,:) !< dim(lmmaxd,lmmaxd,nsymat)
    double complex, intent(out) :: GmatN(:,:,:) !< dim(lmsd,lmsd,num_local_atoms) ! result
    double complex, intent(inout) :: Ginp_local(:,:,0:,:,:) !< dim(lmmaxd,lmmaxd,0:Lly,naclsd,num_local_atoms) reference green function 
    double complex, intent(in) :: tmatLL(:,:,:,0:) !< t-matrices (lmsd,lmsd,num_trunc_atoms,0:Lly)
    double precision, intent(in) :: rr(:,0:) !< lattice vectors(1:3,0:nrd)
    integer, intent(in) :: NofKs
    double precision, intent(in) :: volBZ
    double precision, intent(in) :: Bzkp(:,:) ! dim(3,kpoibz)
    double precision, intent(in) :: k_point_weights(:) ! dim kpoibz

    integer, intent(in) :: korbit ! NOCO
    
    ! LLY
    double complex, intent(in)    :: tr_alph(:)
    double complex, intent(out)   :: Lly_grdt
    integer       , intent(in)    :: global_atom_idx_Lly
    integer       , intent(in)    :: Lly
    
    integer, intent(in) :: solver_type
    type(TimerMpi), intent(inout) :: kpoint_timer, kernel_timer

    external :: zgetri, zgetrf, zgemm ! LAPACK routines
 
    ! locals
    double precision :: mrfctori, tauvBZ
    integer :: N, isym, num_local_atoms, naclsd, num_trunc_atoms, ila, nsymat, ist ! status for LAPACK calls
    double complex, allocatable :: GS(:,:,:)
    integer, allocatable :: ipvt(:), info(:,:) ! work array for LAPACK
    double complex, allocatable :: temp(:), gll(:,:), tpg(:,:), xc(:,:), mssq(:,:,:) ! effective (site-dependent) delta_t^(-1) matrix

    N = size(tmatLL, 2) ! ToDo: perpare for non-collinear
    assert( N == size(tmatLL, 1) )
    nsymat = size(dsymLL, 3)
    num_local_atoms = size(op%atom_indices)
    num_trunc_atoms = size(tmatLL, 3)
    naclsd = size(Ginp_local, 4)
    
    assert( all(shape(dsymLL) == [N,N,nsymat]) )
    assert( all(shape(GmatN) == [N,N,num_local_atoms]) )
    assert( all(shape(Ginp_local) == [N/(korbit+1),N/(korbit+1),1+Lly,naclsd,num_local_atoms]) )
    assert( all(shape(tmatLL) == [N,N,num_trunc_atoms,1+Lly]) )

    allocate(GS(N,N,num_local_atoms), mssq(N,N,num_local_atoms), info(2,num_local_atoms), stat=ist)
    if (ist /= 0) stop "KLOOPZ1: FATAL Error, failure to allocate memory, probably out of memory."
    
!$omp parallel private(ipvt, temp)
    allocate(ipvt(N), temp(N*N), stat=ist)
!$omp do private(ila)
    do ila = 1, num_local_atoms
      mssq(:,:,ila) = tmatLL(1:N,1:N,op%atom_indices(ila),0)
      ! inversion:
      !     The (local) Delta_t matrix is inverted and stored in mssq
      call zgetrf(N,N,mssq(:,:,ila),N,ipvt,info(1,ila)) ! LU-factorize
      call zgetri(N,mssq(:,:,ila),N,ipvt,temp,N*N,info(2,ila)) ! compute inverse
    enddo ! ila
!$omp end do
    deallocate(ipvt, temp, stat=ist)
!$omp end parallel

    if (any(info /= 0)) write(*, '(a,999("  f",i0," i",i0))') "zgetrf or zgetri returned an error:", info ! warn
    deallocate(info, stat=ist)
    
!     rfctor=A/(2*PI) conversion factor to p.u.
    mrfctori = -(2.d0*pi)/alat ! = inverse of -alat/(2*PI)

!=======================================================================
!     Note: the actual k-loop is in MultipleScattering (it is not parallelized)
!     The integration over k is also performed in MultipleScattering
    
    
    
    ! solver_type=3 T-matrix cutoff with new solver
    ! solver_type=4 T-matrix cutoff with direct solver
    call MultipleScattering(solv, op, precond, Bzkp, NofKs, k_point_weights, GS, tmatLL, alat, nsymat, rr, &
                      Ginp_local, global_atom_id, communicator, xTable, iguess_data, ienergy, ispin, &
                      mssq, tr_alph, Lly_grdt, volBZ, global_atom_idx_Lly, Lly, & !LLY
                      solver_type, kpoint_timer, kernel_timer)
                      
!-------------------------------------------------------- SYMMETRISE gll

!      MultipleScattering returns GS (local) which contains the scattering path operator
!      (already integrated over the irreducible wedge in k-space)
!      scattering path operator: ((Delta_T)^-1 - G_ref)^-1

!      All the symmetry operations are applied on GS and summed over
!     - the result is stored in gll
!      Note: the symmetry operations apply on the (LL')-space
    tauvBZ = 1.d0/volBZ

!$omp parallel private(gll, tpg, xc)
    allocate(gll(N,N), tpg(N,N), xc(N,N), stat=ist)
    if (ist /= 0) stop "KLOOPZ1: FATAL Error2, failure to allocate memory, probably out of memory."
!$omp do private(ila, isym)
    do ila = 1, num_local_atoms

      gll(:,:) = GS(:N,:N,ila) ! 1st symmetry matrix is the unity operation
      do isym = 2, nsymat

    !     gll = sum(i=1,iumax)(tauvBZ * DLL(i) * GS * DLL(i)^H)
    !     ull(1) is equal to unity matrix

        !       tpg = tauvBZ * DLL * GS

        call zgemm('n','n',N,N,N,CONE,dsymLL(:,:,isym),size(dsymLL, 1),GS(:,:,ila),N,ZERO,tpg(:,:),N)

        !     gll = gll + tpg * DLL(i)^H
        !                           C  ! dsymLL might be complex in REL case

        call zgemm('n','c',N,N,N,CONE,tpg(:,:),N,dsymLL(:,:,isym),size(dsymLL, 1),CONE,gll(:,:),N)

      enddo ! isym
      gll(:,:) = gll(:,:)*tauvBZ


!     xc = Delta_t^(-1) * gll
      call zgemm('n','n',N,N,N,CONE,mssq(:,:,ila),N,gll(:,:),N,ZERO,xc(:,:),N)

  !       gll is overwritten with the following expression:
  !       (for the in configuration space diagonal structural Green's
  !       function of the REAL system - already integrated over k)

  !   gll = - Delta_t^(-1) - Delta_t^(-1) * gll * Delta_t^(-1)

  !       copy overwrite gll with content of mssq
      gll(:,:) = mssq(:,:,ila) ! call ZCOPY(N*N,mssq(:,1),1,gll,1)


  !       gll = (-1) * xc                     *    mssq     + (-1) * gll
  !                    |                            |                 |
  !            Delta_t^-1 * scat. path op.      Delta_t^-1    Delta_t^-1

      call zgemm('n','n',N,N,N,CONE,xc,N,mssq(:,:,ila),N,CONE,gll,N)

!       gll(:,:) = xc(:,:) + gll(:,:)

  !   GmatN = GMATLL = -gll/rfctor...............rescaled and copied into output array

      GmatN(:N,:N,ila) = gll(:,:)*mrfctori

    enddo ! ila
!$omp end do
    deallocate(gll, tpg, xc, stat=ist)
!$omp end parallel

    if (global_jij_data%do_jij_calculation) then
      call SYMJIJ(alat, tauvBZ, nsymat, dsymLL, global_jij_data%NXIJ, global_jij_data%IXCP, &
                  tmatLL, mssq, global_jij_data%GSXIJ, global_jij_data%GMATXIJ(:,:,:,global_jij_data%active_spin), &  ! Result
                  num_trunc_atoms, N, global_jij_data%nxijd)
    endif ! jij
    
    deallocate(GS, mssq, stat=ist)
    
  endsubroutine ! kloopz1

endmodule ! kloopz1_mod
