!> multiple scattering k-loop and symmetrisation
module kloopz1_mod
  implicit none
  private
  public :: kloopz1_new

  double complex, parameter :: CONE=(1.d0,0.d0), ZERO=(0.d0,0.d0)
  
  contains

  subroutine kloopz1_new(Gmatn, solv, op, precond, alat, NofKs, volBZ, Bzkp, k_point_weights, rr, Ginp_local, &
                         nsymat, dsymll, tmatLL, lmmaxd, trunc2atom_index, communicator, iguess_data, &
                         DGinp_local, dtde, tr_alph, lly_grdt, &
                         global_atom_idx_lly, lly) ! LLY 

! only part of arrays for corresponding spin direction is passed
! (Gmatn, tsst_local, dtde_local, lly_grdt, tr_alph, gmatxij)
!
! NofKs .. number of k-points, integer
! volBZ .. brillouin zone volume, double
! Bzkp ... k-points of used k-mesh ... dimension (3, kpoibz)
! k_point_weights . array of brillouin zone integration weights for each k-point ... dimension (kpoibz)

! Ginp_local ... reference green's function
! dginp ...      derivative of reference green's function
! tsst_local ..  t-matrix

    use kkrmat_new_mod, only: kkrmat01_new
    
    use TEST_lcutoff_mod, only: cutoffmode
    use InitialGuess_mod, only: InitialGuess
    use IterativeSolver_mod, only: IterativeSolver
    use BCPOperator_mod, only: BCPOperator
    use KKROperator_mod, only: KKROperator
    use jij_calc_mod, only: global_jij_data, symjij
    integer, parameter :: nsymaxd = 48

    type(IterativeSolver), intent(inout) :: solv
    type(KKROperator), intent(inout) :: op
    type(BCPOperator), intent(inout) :: precond

    integer, intent(in) :: lmmaxd
#define N lmmaxd
    !> mapping trunc. index -> atom index
    integer, intent(in) :: trunc2atom_index(:)
    integer, intent(in) :: communicator
    type(InitialGuess), intent(inout) :: iguess_data
    double precision, intent(in) :: alat
    integer, intent(in) :: nsymat
    double complex, intent(in) :: dsymll(N,N,nsymat) !<
    double complex, intent(out) :: Gmatn(:,:,:) !< (N,N,num_local_atoms)
    double complex, intent(inout) :: Ginp_local(:,:,:,:) !< reference green function
    double complex, intent(in) :: tmatLL(:,:,:) !< t-matrices (lmmaxd,lmmaxd,naez)
    double precision, intent(in) :: rr(:,0:) !< lattice vectors(1:3,0:nrd)
    integer, intent(in) :: NofKs
    double precision, intent(in) :: volBZ
    double precision, intent(in) :: Bzkp(:,:) ! dim (3,kpoibz)
    double precision, intent(in) :: k_point_weights(:) ! dim kpoibz
 
    ! LLY
    double complex, intent(inout)   :: DGinp_local(:,:,:,:) !< energy derivative of reference green function
    double complex, intent(in)   :: tr_alph(:)
    double complex, intent(in)   :: dtde(:,:,:) 
    double complex, intent(out)  :: lly_grdt
    integer       , intent(in)   :: global_atom_idx_lly
    integer       , intent(in)   :: lly

    external :: zgetri, zgetrf, zgemm ! LAPACK routines
 
    ! locals
    double complex :: tauvBZ ! could be real but it enters a zBLAS call
    double precision :: rfctori, tauvBZr
    integer :: ist, info ! status for LAPACK calls
    integer :: ispin, isym
    integer :: num_local_atoms, ila
    double complex, allocatable :: mssq(:,:,:), GS(:,:,:) ! effective (site-dependent) delta_t^(-1) matrix
    integer, allocatable :: ipvt(:) ! work array for LAPACK
    double complex, allocatable :: temp(:,:) ! work array for LAPACK
    double complex, allocatable :: gll(:,:), tpg(:,:), xc(:,:)

#define cluster_info op%cluster_info

    num_local_atoms = size(op%atom_indices)


!     rfctor=A/(2*PI) conversion factor to p.u.
!    mrfctori = -(8.d0*atan(1.0d0))/alat ! = inverse of -alat/(2*PI)
     rfctori = (8.d0*atan(1.0d0))/alat ! = inverse of -alat/(2*PI)

    allocate(temp(N,N), ipvt(N), mssq(N,N,num_local_atoms), stat=ist)
    do ila = 1, num_local_atoms

      mssq(:,:,ila) = tmatLL(1:N,1:N,op%atom_indices(ila))

      ! inversion:
      !     The (local) Delta_t matrix is inverted and stored in mssq
      call zgetrf(N,N,mssq(:,:,ila),N,ipvt,info)
      call zgetri(N,mssq(:,:,ila),N,ipvt,temp,N*N,info)

    enddo ! ila
    deallocate(temp, ipvt, stat=ist) ! ignore status

!=======================================================================
!     Note: the actual k-loop is in kkrmat01 (it is not parallelized)
!     The integration over k is also performed in kkrmat01

    ! cutoffmode:
    ! 0 no cutoff
    ! 1 T-matrix cutoff, 2 full matrix cutoff (not supported anymore)
    if (cutoffmode == 1 .or. cutoffmode == 2) stop "0 < cutoffmode < 3 not supported."
    
    allocate(GS(N,N,num_local_atoms), stat=ist)
    if (ist /= 0) stop "KLOOPZ1: FATAL Error, failure to allocate memory, probably out of memory."
    
    
    ! 3 T-matrix cutoff with new solver
    ! 4 T-matrix cutoff with direct solver
    call kkrmat01_new(solv, op, precond, Bzkp, NofKs, k_point_weights, GS, tmatLL, alat, nsymat, rr, Ginp_local, lmmaxd, trunc2atom_index, communicator, iguess_data, &
                      mssq, DGinp_local, dtde, tr_alph, lly_grdt, k_point_weights, volBZ, global_atom_idx_lly, lly) !LLY
!-------------------------------------------------------- SYMMETRISE gll

!      kkrmat01 returns GS (local) which contains the scattering path operator
!      (already integrated over the irreducible wedge in k-space)
!      scattering path operator: ((Delta_T)^-1 - G_ref)^-1

!      All the symmetry operations are applied on GS and summed over
!     - the result is stored in gll
!      Note: the symmetry operations apply on the (LL')-space
    tauvBZ  = CONE/volBZ
    tauvBZr = 1.d0/volBZ
    
    allocate(gll(N,N), tpg(N,N), xc(N,N), stat=ist)
    if (ist /= 0) stop "KLOOPZ1: FATAL Error2, failure to allocate memory, probably out of memory."

    do ila = 1, num_local_atoms

      do isym = 1, nsymat

    !     gll = sum(i=1,iumax)(tauvBZ * DLL(i) * GS * DLL(i)^H)

        if (isym == 1) then

        !     ull(1) is equal to unity matrix

          gll(:,:) = tauvBZr*GS(1:N,1:N,ila)

        else

        !       tpg = tauvBZ * DLL * GS

          call zgemm('n','n',N,N,N,tauvBZ,dsymll(:,:,isym),N,GS(:,:,ila),N,ZERO,tpg(:,:),N)

        !     gll = gll + tpg * DLL(i)^H
        !                           C  ! dsymll might be complex in REL case

          call zgemm('n','c',N,N,N,CONE,tpg(:,:),N,dsymll(:,:,isym),N,CONE,gll(:,:),N)
        endif ! isym == 1

      enddo ! isym


!     xc = Delta_t^(-1) * gll

     call zgemm('n','n',N,N,N,CONE,mssq(:,:,ila),N,gll(:,:),N,ZERO,xc(:,:),N)

  !       gll is overwritten with the following expression:
  !       (for the in configuration space diagonal structural Green's
  !       function of the REAL system - already integrated over k)

  !   gll = - Delta_t^(-1) - Delta_t^(-1) * gll * Delta_t^(-1)

  !       copy overwrite gll with content of mssq
       gll(:,:) = mssq(:,:,ila) ! call ZCOPY(N**2,mssq(:,:,ila),1,gll,1)


  !       gll = (-1) * xc                     *    mssq     + (-1) * gll
  !                    |                            |                 |
  !            Delta_t^-1 * scat. path op.      Delta_t^-1    Delta_t^-1


       call zgemm('n','n',N,N,N,-CONE,xc,N, mssq(:,:,ila),N,-CONE,gll,N)

!       gll(:,:) = - xc(:,:) - gll(:,:)

  !   Gmatn = GMATLL = gll/rfctor...............rescaled and copied into output array

       Gmatn(1:N,1:N,ila) = gll(1:N,1:N)*rfctori
!      Gmatn(1:N,1:N,ila) = (mssq(:,:,ila) + xc(:,:))*mrfctori

    enddo ! ila

    if (global_jij_data%do_jij_calculation) then

      ispin = global_jij_data%active_spin

      call SYMJIJ(alat, tauvBZ, nsymat, dsymll, global_jij_data%NXIJ, global_jij_data%IXCP, &
                  tmatLL, mssq, global_jij_data%GSXIJ, global_jij_data%GMATXIJ(:,:,:,ispin), &  ! Result
                  cluster_info%naez_trc, lmmaxd, global_jij_data%nxijd)
    endif ! jij

    deallocate(GS, mssq, gll, tpg, xc, stat=ist)

#undef cluster_info
#undef ms
  endsubroutine ! kloopz1_new

endmodule ! kloopz1_mod
