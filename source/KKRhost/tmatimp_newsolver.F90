!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculate and write down impurity t-matrix and delta matrix first calculate t-matrix for the host corresponding to imp. cluster
!> Author: N. H. Long
!> Calculate and write down impurity t-matrix and delta matrix
!> first calculate t-matrix for the host corresponding to imp. cluster
!------------------------------------------------------------------------------------
!> @note - Adapted to new routines (mainly changed interfaces) to work in KKRcode
!> also added MPI parallelization. Philipp Rüssmann, Juelich, 09.2017
!> @endnote
!------------------------------------------------------------------------------------
module mod_tmatimp_newsolver

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculate and write down impurity t-matrix and delta matrix first calculate t-matrix for the host corresponding to imp. cluster
  !> Author: N. H. Long
  !> Category: single-site, KKRhost 
  !> Deprecated: False 
  !> Calculate and write down impurity t-matrix and delta matrix
  !> first calculate t-matrix for the host corresponding to imp. cluster
  !-------------------------------------------------------------------------------
  !> @note - Adapted to new routines (mainly changed interfaces) to work in KKRcode
  !> also added MPI parallelization. Philipp Rüssmann, Juelich, 09.2017
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine tmatimp_newsolver(irm,ksra,lmax,iend,irid,lpot,natyp,ncleb,ipand,irnsd,&
    nfund,ihost,ntotd,nspin,lmpot,ncheb,lmmax0d,korbit,nspotd,ielast,irmind,npan_eq, &
    npan_log,natomimp,r_log,vins,vm2z,ipan,irmin,hostimp,ipanimp,irwsimp,atomimp,   &
    irminimp,icleb,ircut,ircutimp,zat,zimp,rmesh,cleb,rimp,rclsimp,eryd,vm2zimp,    &
    vinsimp,dtmtrx,lmmaxd)

#ifdef CPP_MPI
    use :: mpi
    use :: mod_mympi, only: myrank, master, nranks, distribute_linear_on_tasks
#else
    use :: mod_mympi, only: myrank, master, nranks
#endif
    use :: mod_types, only: t_inc, t_imp
    use :: mod_runoptions, only: disable_tmat_sratrick, write_green_imp, write_pkkr_operators
    use :: mod_create_newmesh, only: create_newmesh 
    use :: mod_version_info, only: 
    use :: mod_wunfiles, only: t_params
    use :: mod_constants, only: czero, pi, cone, cvlight
    use :: mod_profiling, only: memocc
    use :: mod_datatypes, only: dp
    use :: mod_calcsph, only: calcsph
    use :: mod_interpolate_poten, only: interpolate_poten
    use :: mod_intcheb_cell, only: intcheb_cell
    use :: mod_rllsll, only: rllsll
    use :: mod_rllsllsourceterms, only: rllsllsourceterms
    use :: mod_spinorbit_ham, only: spinorbit_ham
    use :: mod_rotatespinframe, only: rotatematrix
    use :: mod_vllmat, only: vllmat
    use :: mod_vllmatsra, only: vllmatsra

    implicit none

    ! .. Input variables
    integer, intent (in) :: irm    !! Maximum number of radial points
    integer, intent (in) :: ksra
    integer, intent (in) :: lmax   !! Maximum l component in wave function expansion
    integer, intent (in) :: iend
    integer, intent (in) :: irid
    integer, intent (in) :: lpot   !! Maximum l component in potential expansion
    integer, intent (in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: ncleb  !! Number of Clebsch-Gordon coefficients
    integer, intent (in) :: ipand  !! Number of panels in non-spherical part
    integer, intent (in) :: irnsd
    integer, intent (in) :: nfund  !! Shape functions parameters in non-spherical part
    integer, intent (in) :: ntotd
    integer, intent (in) :: ihost
    integer, intent (in) :: nspin  !! Counter for spin directions
    integer, intent (in) :: lmpot  !! (LPOT+1)**2
    integer, intent (in) :: ncheb  !! Number of Chebychev pannels for the new solver
    integer, intent (in) :: korbit !! Spin-orbit/non-spin-orbit (1/0) added to the Schroedinger or SRA equations. Works with FP. KREL and KORBIT cannot be both non-zero.
    integer, intent (in) :: lmmax0d !! (LMAX+1)^2
    integer, intent (in) :: lmmaxd !! (KREL+KORBIT+1)*(LMAX+1)^2
    integer, intent (in) :: nspotd !! Number of potentials for storing non-sph. potentials
    integer, intent (in) :: ielast
    integer, intent (in) :: irmind !! IRM-IRNSD
    integer, intent (in) :: npan_eq !! Number of intervals from [R_LOG] to muffin-tin radius Used in conjunction with runopt NEWSOSOL
    integer, intent (in) :: npan_log !! Number of intervals from nucleus to [R_LOG] Used in conjunction with runopt NEWSOSOL
    integer, intent (in) :: natomimp !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
    real (kind=dp), intent (in) :: r_log !! Radius up to which log-rule is used for interval width. Used in conjunction with runopt NEWSOSOL
    integer, dimension (natyp), intent (in) :: ipan !! Number of panels in non-MT-region
    integer, dimension (natyp), intent (in) :: irmin !! Max R for spherical treatment
    integer, dimension (ihost), intent (in) :: hostimp
    integer, dimension (natomimp), intent (in) :: ipanimp
    integer, dimension (natomimp), intent (in) :: irwsimp
    integer, dimension (natomimp), intent (in) :: atomimp
    integer, dimension (natomimp), intent (in) :: irminimp
    integer, dimension (ncleb, 4), intent (in) :: icleb !! Pointer array
    integer, dimension (0:ipand, natyp), intent (in) :: ircut !! R points of panel borders
    integer, dimension (0:ipand, natomimp), intent (in) :: ircutimp
    real (kind=dp), dimension (natyp), intent (in) :: zat !! Nuclear charge
    real (kind=dp), dimension (natomimp), intent (in) :: zimp
    real (kind=dp), dimension (irm, natyp), intent (in) :: rmesh !! Radial mesh ( in units a Bohr)
    real (kind=dp), dimension (ncleb, 2), intent (in) :: cleb !! GAUNT coefficients (GAUNT)
    real (kind=dp), dimension (irm, natomimp), intent (in) :: rimp
    real (kind=dp), dimension (3, natomimp), intent (in) :: rclsimp
    complex (kind=dp), intent (in) :: eryd
    ! .. In/Out variables
    real (kind=dp), dimension (irm, nspin*natyp), intent (inout) :: vm2z
    real (kind=dp), dimension (irmind:irm, lmpot, nspotd*natyp), intent (inout) :: vins
    real (kind=dp), dimension (irm, nspin*natomimp), intent (inout) :: vm2zimp
    real (kind=dp), dimension (irmind:irm, lmpot, nspotd*natomimp), intent (inout) :: vinsimp
    complex (kind=dp), dimension ((korbit+1)*lmmax0d*natomimp, (korbit+1)*lmmax0d*natomimp), intent (inout) :: dtmtrx
    ! .. Local variables
    integer :: ipot
    integer :: i1, ir, nsra, use_sratrick, nvec, lm1, lm2, i2, il1, il2, irmdnewd
    integer :: i_stat, i_all
#ifdef CPP_MPI
    integer :: ierr
#endif
    real (kind=dp) :: theta, phi
    complex (kind=dp) :: gmatprefactor
    real (kind=dp), dimension (natomimp) :: phiimp
    real (kind=dp), dimension (natyp) :: phihost
    real (kind=dp), dimension (natomimp) :: thetaimp
    real (kind=dp), dimension (natyp) :: thetahost
    complex (kind=dp), dimension (2*(lmax+1)) :: tmatsph
    complex (kind=dp), dimension (2*(lmax+1)) :: dummy_alpha
    complex (kind=dp), dimension ((korbit+1)*lmmax0d, (korbit+1)*lmmax0d, ihost) :: tmatll
    complex (kind=dp), dimension ((korbit+1)*lmmax0d, (korbit+1)*lmmax0d) :: dummy_alphaget
    ! .. Allocatable variables
    integer, allocatable :: npan_tot(:)
    integer, allocatable :: npan_log_at(:)
    integer, allocatable :: npan_eq_at(:)
    integer, allocatable :: npan_inst(:)
    integer, dimension (:), allocatable :: irmdnew
    integer, dimension (:), allocatable :: jlk_index
    integer, dimension (:, :), allocatable :: ipan_intervall
    real (kind=dp), dimension (:, :), allocatable :: rnew, rpan_intervall
    real (kind=dp), dimension (:, :, :), allocatable :: vinsnew
    complex (kind=dp), dimension (:), allocatable :: deltatmp
    complex (kind=dp), dimension (:, :), allocatable :: hlk, jlk, hlk2, jlk2
    complex (kind=dp), dimension (:, :), allocatable :: radialhost, radialimp
    complex (kind=dp), dimension (:, :), allocatable :: vllimp, deltav, deltaimp
    complex (kind=dp), dimension (:, :, :), allocatable :: vnsimp
    complex (kind=dp), dimension (:, :, :), allocatable :: rll, sll
    complex (kind=dp), dimension (:, :, :), allocatable :: deltabg, deltasm
    complex (kind=dp), dimension (:, :, :), allocatable :: tmatllimp, deltamtr
    complex (kind=dp), dimension (:, :, :), allocatable :: vnspll0, vnspll, vnspll1
    complex (kind=dp), dimension (:, :, :, :), allocatable :: rllhost
    complex (kind=dp), dimension (:, :, :, :), allocatable :: vnshost
    ! .. Parallelization variables
    integer :: i1_start, i1_end, i1_start_imp, i1_end_imp
#ifdef CPP_MPI
    integer, dimension (0:nranks-1) :: ntot_pt, ioff_pt
    complex (kind=dp), dimension (:, :, :), allocatable :: temp
    complex (kind=dp), dimension (:, :, :, :), allocatable :: temp2 ! needed for MPI communication
#endif

    if (myrank==master) write (6, *) 'in tmatimp'
    if (ksra>=1) then
      nsra = 2
    else
      nsra = 1
    end if

    ! these are used by both host and impuity and therefore need the maximal
    ! aray size
    allocate (npan_tot(max(natyp, natomimp)), stat=i_stat)
    call memocc(i_stat, product(shape(npan_tot))*kind(npan_tot), 'npan_tot', 'tmatimp_newsolver')
    allocate (npan_eq_at(max(natyp, natomimp)), stat=i_stat)
    call memocc(i_stat, product(shape(npan_eq_at))*kind(npan_eq_at), 'npan_eq', 'tmatimp_newsolver')
    allocate (npan_log_at(max(natyp, natomimp)), stat=i_stat)
    call memocc(i_stat, product(shape(npan_log_at))*kind(npan_log_at), 'npan_log_at', 'tmatimp_newsolver')
    allocate (npan_inst(max(natyp, natomimp)), stat=i_stat)
    call memocc(i_stat, product(shape(npan_inst))*kind(npan_inst), 'npan_inst', 'tmatimp_newsolver')

    allocate (jlk_index(2*lmmaxd), stat=i_stat)
    call memocc(i_stat, product(shape(jlk_index))*kind(jlk_index), 'JLK_INDEX', 'tmatimp_newsolver')
    allocate (deltav(lmmaxd,lmmaxd), stat=i_stat)
    call memocc(i_stat, product(shape(deltav))*kind(deltav), 'DELTAV', 'tmatimp_newsolver')
    allocate (deltamtr(lmmaxd,lmmaxd,natomimp), stat=i_stat)
    call memocc(i_stat, product(shape(deltamtr))*kind(deltamtr), 'DELTAMTR', 'tmatimp_newsolver')
    allocate (tmatllimp(lmmaxd,lmmaxd,natomimp), stat=i_stat)
    call memocc(i_stat, product(shape(tmatllimp))*kind(tmatllimp), 'TMATLLIMP', 'tmatimp_newsolver')
    allocate (rllhost(nsra*lmmaxd,lmmaxd,ihost,ntotd*(ncheb+1)), stat=i_stat)
    call memocc(i_stat, product(shape(rllhost))*kind(rllhost), 'RLLHOST', 'tmatimp_newsolver')
    allocate (vnshost(nsra*lmmaxd,nsra*lmmaxd,ihost,ntotd*(ncheb+1)), stat=i_stat)
    call memocc(i_stat, product(shape(vnshost))*kind(vnshost), 'VNSHOST', 'tmatimp_newsolver')
    tmatll = czero
    rllhost = czero
    vnshost = czero
    tmatllimp = czero
    deltamtr = czero
    deltav = czero
    tmatsph = czero

    if (myrank==master) then
      write(*,'("**************************************************************************************************")')
      write(*,'("***  WARNING: tmatimp_newsolver still uses the old rllsll!                                     ***")')
      write(*,'("**************************************************************************************************")')
      ! read angles from nonco_ange files
      open (unit=12, file='nonco_angle.dat', form='FORMATTED')
      do i1 = 1, natyp
        read (12, *) thetahost(i1), phihost(i1)
        thetahost(i1) = thetahost(i1)/360.0_dp*2.0_dp*pi
        phihost(i1) = phihost(i1)/360.0_dp*2.0_dp*pi
      end do
      close (12)
      open (unit=13, file='nonco_angle_imp.dat', form='FORMATTED')
      do i1 = 1, natomimp
        read (13, *) thetaimp(i1), phiimp(i1)
        thetaimp(i1) = thetaimp(i1)/360.0_dp*2.0_dp*pi
        phiimp(i1) = phiimp(i1)/360.0_dp*2.0_dp*pi
      end do
      close (13)
    end if

#ifdef CPP_MPI
    ! broadcast read-in values to all ranks (MPI_COMM_WORLD since
    ! atom dimension is solely used without energy parallelization)
    call mpi_bcast(thetahost, natyp, mpi_double_precision, master, mpi_comm_world, ierr)
    if (ierr/=0) stop 'Error MPI_Bcast THETAhost in tmatimp'
    call mpi_bcast(phihost, natyp, mpi_double_precision, master, mpi_comm_world, ierr)
    if (ierr/=0) stop 'Error MPI_Bcast PHIhost in tmatimp'
    call mpi_bcast(thetaimp, natomimp, mpi_double_precision, master, mpi_comm_world, ierr)
    if (ierr/=0) stop 'Error MPI_Bcast THETAimp in tmatimp'
    call mpi_bcast(phiimp, natomimp, mpi_double_precision, master, mpi_comm_world, ierr)
    if (ierr/=0) stop 'Error MPI_Bcast PHIimp in tmatimp'

    ! set start/end of parallel atom loops
    if (t_inc%i_write>0) write (1337, *) 'Parallelization host atoms:'
    call distribute_linear_on_tasks(nranks,myrank,master,ihost,ntot_pt,ioff_pt,.true.)
    i1_start = ioff_pt(myrank) + 1
    i1_end = ioff_pt(myrank) + ntot_pt(myrank)
    if (t_inc%i_write>0) write (1337, *) 'Parallelization imp. atoms:'
    call distribute_linear_on_tasks(nranks,myrank,master,natomimp,ntot_pt,ioff_pt,.true.)
    i1_start_imp = ioff_pt(myrank) + 1
    i1_end_imp = ioff_pt(myrank) + ntot_pt(myrank)
#else
    i1_start = 1
    i1_end = ihost
    i1_start_imp = 1
    i1_end_imp = natomimp
#endif

    if (write_green_imp) then

      !------------------------------------------------------------------------------
      ! START calculate tmat and radial wavefunctions of host atoms
      !------------------------------------------------------------------------------

      ! create new mesh before loop starts
      ! data for the new mesh
      allocate (irmdnew(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(irmdnew))*kind(irmdnew), 'irmdnew', 'tmatimp_newsolver')
      irmdnewd = 0
      do i1 = 1, natyp
        npan_inst(i1) = ipan(i1) - 1
        npan_tot(i1) = npan_log + npan_eq + npan_inst(i1)
        if (npan_tot(i1)*(ncheb+1)>irmdnewd) then
          irmdnewd = npan_tot(i1)*(ncheb+1)
        end if
        irmdnew(i1) = npan_tot(i1)*(ncheb+1)
      end do
      ! new mesh
      allocate (rnew(irmdnewd,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(rnew))*kind(rnew), 'RNEW', 'tmatimp_newsolver')
      allocate (rpan_intervall(0:ntotd,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(rpan_intervall))*kind(rpan_intervall), 'RPAN_INTERVALL', 'tmatimp_newsolver')
      allocate (ipan_intervall(0:ntotd,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(ipan_intervall))*kind(ipan_intervall), 'IPAN_INTERVALL', 'tmatimp_newsolver')
      allocate (vinsnew(irmdnewd,lmpot,nspotd*natyp), stat=i_stat) ! NSPIND*max(NATYP,NATOMIMP)))
      call memocc(i_stat, product(shape(vinsnew))*kind(vinsnew), 'VINSNEW', 'tmatimp_newsolver')

      call create_newmesh(natyp,irm,ipand,irid,ntotd,nfund,ncheb,irmdnewd,nspin,    &
        rmesh(:,:),irmin(:),ipan(:),ircut(0:ipand,:),r_log,npan_log,npan_eq,        &
        npan_log_at(:),npan_eq_at(:),npan_tot(:),rnew(:,:),                         &
        rpan_intervall(0:ntotd,:),ipan_intervall(0:ntotd,:),1)

      ! in second step interpolate potential (gain atom by atom with NATYPD==1)
      call interpolate_poten(lpot,irm,irnsd,natyp,ipand,lmpot,nspotd*natyp,ntotd,irmdnewd,&
        nspin,rmesh(:,:),irmin(:),t_params%irws(:),ircut(0:ipand,:),                &
        vins(irmind:irm,1:lmpot,:),vm2z(:,:),npan_log_at(:),npan_eq_at(:),          &
        npan_tot(:),rnew(:,:),ipan_intervall(0:ntotd,:),vinsnew)

      ! calculate tmat and radial wavefunctions of host atoms
      ! parallelized with MPI over atoms
      do i2 = i1_start, i1_end
        i1 = hostimp(i2)

        theta = thetahost(i1)
        phi = phihost(i1)
        ipot = nspin*(i1-1) + 1
        write (6, *) 'HOST', i2, i1, irmdnew(i1)

        ! set up the non-spherical ll' matrix for potential VLL'
        if (nsra==2) then
          use_sratrick = 1
          if (disable_tmat_sratrick) use_sratrick = 0
        else if (nsra==1) then
          use_sratrick = 0
        end if

        allocate (vnspll0(lmmaxd,lmmaxd,irmdnew(i1)), stat=i_stat)
        call memocc(i_stat, product(shape(vnspll0))*kind(vnspll0), 'VNSPLL0', 'tmatimp_newsolver')
        vnspll0 = czero
        ! output potential onto which SOC is added
        allocate (vnspll1(lmmaxd,lmmaxd,irmdnew(i1)), stat=i_stat)
        call memocc(i_stat, product(shape(vnspll1))*kind(vnspll1), 'VNSPLL1', 'tmatimp_newsolver')
        vnspll1 = czero

        call vllmat(1,irmdnew(i1),irmdnew(i1),lmmax0d,lmmaxd,vnspll0,               &
          vinsnew(1:irmdnew(i1),1:lmpot,ipot:ipot+nspin-1),lmpot,cleb,icleb,iend,   &
          nspin,zat(i1),rnew(1:irmdnew(i1),i1),use_sratrick,ncleb)
        ! contruct the spin-orbit coupling hamiltonian and add to potential
        call spinorbit_ham(lmax,lmmax0d,vinsnew(1:irmdnew(i1),1:lmpot,ipot:ipot+nspin-1),&
          rnew(1:irmdnew(i1),i1),eryd,zat(i1),cvlight,t_params%socscale(i1),nspin,  &
          lmpot,theta,phi,ipan_intervall(0:ntotd,i1),rpan_intervall(0:ntotd,i1),    &
          npan_tot(i1),ncheb,irmdnew(i1),irmdnew(i1),vnspll0,vnspll1,'1')
        ! extend matrix for the SRA treatment
        if (nsra==2) then
          allocate (vnspll(2*lmmaxd,2*lmmaxd,irmdnew(i1)), stat=i_stat)
          call memocc(i_stat, product(shape(vnspll))*kind(vnspll), 'VNSPLL', 'tmatimp_newsolver')
          if (use_sratrick==0) then
            call vllmatsra(vnspll1,vnspll,rnew(1:irmdnew(i1),i1),lmmaxd,           &
              irmdnew(i1),irmdnew(i1),eryd,lmax,0,'Ref=0')
          else if (use_sratrick==1) then
            call vllmatsra(vnspll1,vnspll,rnew(1:irmdnew(i1),i1),lmmaxd,           &
              irmdnew(i1),irmdnew(i1),eryd,lmax,0,'Ref=Vsph')
          end if
        else
          allocate (vnspll(lmmaxd,lmmaxd,irmdnew(i1)), stat=i_stat)
          call memocc(i_stat, product(shape(vnspll))*kind(vnspll), 'VNSPLL', 'tmatimp_newsolver')
          vnspll(:, :, :) = vnspll1(:, :, :)
        end if

        ! calculate the source terms in the Lippmann-Schwinger equation
        ! these are spherical hankel and bessel functions
        allocate (hlk(1:4*(lmax+1),irmdnew(i1)), stat=i_stat)
        call memocc(i_stat, product(shape(hlk))*kind(hlk), 'HLK', 'tmatimp_newsolver')
        allocate (jlk(1:4*(lmax+1),irmdnew(i1)), stat=i_stat)
        call memocc(i_stat, product(shape(jlk))*kind(jlk), 'JLK', 'tmatimp_newsolver')
        allocate (hlk2(1:4*(lmax+1),irmdnew(i1)), stat=i_stat)
        call memocc(i_stat, product(shape(hlk2))*kind(hlk2), 'HLK2', 'tmatimp_newsolver')
        allocate (jlk2(1:4*(lmax+1),irmdnew(i1)), stat=i_stat)
        call memocc(i_stat, product(shape(jlk2))*kind(jlk2), 'JLK2', 'tmatimp_newsolver')
        hlk = czero
        jlk = czero
        hlk2 = czero
        jlk2 = czero
        gmatprefactor = czero
        call rllsllsourceterms(nsra,nvec,eryd,rnew(1:irmdnew(i1),i1),irmdnew(i1),   &
          irmdnew(i1),lmax,lmmaxd,1,jlk_index,hlk,jlk,hlk2,jlk2,gmatprefactor)

        ! using spherical potential as reference
        if (use_sratrick==1) then
          call calcsph(nsra,irmdnew(i1),irmdnew(i1),lmax,nspin,zat(i1),eryd,lmpot,  &
            lmmaxd,rnew(1:irmdnew(i1),i1),                                         &
            vinsnew(1:irmdnew(i1),1:lmpot,ipot:ipot+nspin-1),ncheb,npan_tot(i1),    &
            rpan_intervall(0:ntotd,i1),jlk_index,hlk,jlk,hlk2,jlk2,gmatprefactor,   &
            tmatsph,dummy_alpha,use_sratrick,.true.)
        end if

        ! calculate the tmat and wavefunctions
        allocate (rll(nvec*lmmaxd,lmmaxd,irmdnewd), stat=i_stat)
        call memocc(i_stat, product(shape(rll))*kind(rll), 'RLL', 'tmatimp_newsolver')
        allocate (sll(nvec*lmmaxd,lmmaxd,irmdnewd), stat=i_stat)
        call memocc(i_stat, product(shape(sll))*kind(sll), 'SLL', 'tmatimp_newsolver')
        rll = czero
        sll = czero

        ! right solutions
        call rllsll(rpan_intervall(0:ntotd,i1),rnew(1:irmdnew(i1),i1),vnspll,rll,   &
          sll,tmatll(:,:,i2),ncheb,npan_tot(i1),lmmaxd,nvec*lmmaxd,4*(lmax+1),    &
          irmdnew(i1),nsra,jlk_index,hlk,jlk,hlk2,jlk2,gmatprefactor,'1','1','0',   &
          use_sratrick,dummy_alphaget)
        if (nsra==2) then
          do ir = 1, irmdnew(i1)
            do lm1 = 1, lmmaxd
              do lm2 = 1, lmmaxd
                rll(lm1+lmmaxd, lm2, ir) = rll(lm1+lmmaxd, lm2, ir)/cvlight
                sll(lm1+lmmaxd, lm2, ir) = sll(lm1+lmmaxd, lm2, ir)/cvlight
              end do
            end do
          end do
        end if
        ! save radial wavefunction for a host
        do ir = 1, irmdnew(i1)
          do lm1 = 1, nvec*lmmaxd
            do lm2 = 1, lmmaxd
              rllhost(lm1, lm2, i2, ir) = rll(lm1, lm2, ir)
            end do
          end do
        end do

        ! add spherical contribution of tmatrix
        if (use_sratrick==1) then
          do lm1 = 1, (korbit+1)*lmmax0d
            tmatll(lm1, lm1, i2) = tmatll(lm1, lm1, i2) + tmatsph(jlk_index(lm1))
          end do
        end if

        ! rotate tmatrix and radial wavefunction to global frame
        call rotatematrix(tmatll(1,1,i2), theta, phi, lmmax0d, 0)

        ! create SRA potential for host
        ! set up the non-spherical ll' matrix for potential VLL'
        vnspll0 = czero
        vnspll1 = czero
        call vllmat(1,irmdnew(i1),irmdnew(i1),lmmax0d,lmmaxd,vnspll0,               &
          vinsnew(1:irmdnew(i1),1:lmpot,ipot:ipot+nspin-1),lmpot,cleb,icleb,iend,   &
          nspin,zat(i1),rnew(1:irmdnew(i1),i1),0,ncleb)

        ! contruct the spin-orbit coupling hamiltonian and add to potential
        call spinorbit_ham(lmax,lmmax0d,vinsnew(1:irmdnew(i1),1:lmpot,ipot:ipot+nspin-1),&
          rnew(1:irmdnew(i1),i1),eryd, zat(i1),cvlight,t_params%socscale(i1),nspin, &
          lmpot,theta,phi,ipan_intervall(0:ntotd,i1),rpan_intervall(0:ntotd,i1),    &
          npan_tot(i1),ncheb,irmdnew(i1),irmdnew(i1),vnspll0,vnspll1,'1')

        ! save potential for a host
        do ir = 1, irmdnew(i1)
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              vnshost(lm1, lm2, i2, ir) = vnspll1(lm1, lm2, ir)
              if (nsra==2) then
                vnshost(lm1+lmmaxd, lm2+lmmaxd, i2, ir) = vnspll1(lm1, lm2, ir)
              end if
            end do
          end do
        end do

        i_all = -product(shape(vnspll0))*kind(vnspll0)
        deallocate (vnspll0, stat=i_stat)
        call memocc(i_stat, i_all, 'VNSPLL0', 'tmatimp_newsolver')
        i_all = -product(shape(vnspll1))*kind(vnspll1)
        deallocate (vnspll1, stat=i_stat)
        call memocc(i_stat, i_all, 'VNSPLL1', 'tmatimp_newsolver')
        i_all = -product(shape(vnspll))*kind(vnspll)
        deallocate (vnspll, stat=i_stat)
        call memocc(i_stat, i_all, 'VNSPLL', 'tmatimp_newsolver')
        i_all = -product(shape(hlk))*kind(hlk)
        deallocate (hlk, stat=i_stat)
        call memocc(i_stat, i_all, 'HLK', 'tmatimp_newsolver')
        i_all = -product(shape(jlk))*kind(jlk)
        deallocate (jlk, stat=i_stat)
        call memocc(i_stat, i_all, 'JLK', 'tmatimp_newsolver')
        i_all = -product(shape(hlk2))*kind(hlk2)
        deallocate (hlk2, stat=i_stat)
        call memocc(i_stat, i_all, 'HLK2', 'tmatimp_newsolver')
        i_all = -product(shape(jlk2))*kind(jlk2)
        deallocate (jlk2, stat=i_stat)
        call memocc(i_stat, i_all, 'JLK2', 'tmatimp_newsolver')
        i_all = -product(shape(sll))*kind(sll)
        deallocate (sll, stat=i_stat)
        call memocc(i_stat, i_all, 'SLL', 'tmatimp_newsolver')
        i_all = -product(shape(rll))*kind(rll)
        deallocate (rll, stat=i_stat)
        call memocc(i_stat, i_all, 'RLL', 'tmatimp_newsolver')
      end do                       ! I2

      i_all = -product(shape(rnew))*kind(rnew)
      deallocate (rnew, stat=i_stat)
      call memocc(i_stat, i_all, 'RNEW', 'tmatimp_newsolver')
      i_all = -product(shape(vinsnew))*kind(vinsnew)
      deallocate (vinsnew, stat=i_stat)
      call memocc(i_stat, i_all, 'VINSNEW', 'tmatimp_newsolver')
      i_all = -product(shape(rpan_intervall))*kind(rpan_intervall)
      deallocate (rpan_intervall, stat=i_stat)
      call memocc(i_stat, i_all, 'RPAN_INTERVALL', 'tmatimp_newsolver')
      i_all = -product(shape(ipan_intervall))*kind(ipan_intervall)
      deallocate (ipan_intervall, stat=i_stat)
      call memocc(i_stat, i_all, 'IPAN_INTERVALL', 'tmatimp_newsolver')

      !------------------------------------------------------------------------------

#ifdef CPP_MPI
      ! collect results and write out only on master
      ! communicate VNSHOST, RLLHOST, TMATLL
      i1 = nsra*lmmaxd*lmmaxd*ihost*ntotd*(ncheb+1)
      ! Allocation of temp2 for RLLHOST
      allocate (temp2(nsra*lmmaxd,lmmaxd,ihost,ntotd*(ncheb+1)), stat=i_stat)
      call memocc(i_stat, product(shape(temp2))*kind(temp2), 'temp2', 'tmatimp_newsolver')
      temp2 = czero

      call mpi_allreduce(rllhost,temp2,i1,mpi_double_complex,mpi_sum,mpi_comm_world,ierr)
      if (ierr/=0) stop 'Error in MPI_Allreduce for RLLHOST in tmatimp'
      rllhost = temp2
      ! Deallocation of temp2 for RLLHOST
      i_all = -product(shape(temp2))*kind(temp2)
      deallocate (temp2, stat=i_stat)
      call memocc(i_stat, i_all, 'temp2', 'tmatimp_newsolver')

      i1 = nsra*lmmaxd*nsra*lmmaxd*ihost*ntotd*(ncheb+1)
      ! Allocation of temp2 for VNSHOST
      allocate (temp2(nsra*lmmaxd,nsra*lmmaxd,ihost,ntotd*(ncheb+1)), stat=i_stat)
      call memocc(i_stat, product(shape(temp2))*kind(temp2), 'temp2', 'tmatimp_newsolver')
      temp2 = czero

      call mpi_allreduce(vnshost,temp2,i1,mpi_double_complex,mpi_sum,mpi_comm_world,ierr)
      if (ierr/=0) stop 'Error in MPI_Allreduce for VNSHOST in tmatimp'
      vnshost = temp2
      ! Deallocation of temp2 for RLLHOST
      i_all = -product(shape(temp2))*kind(temp2)
      deallocate (temp2, stat=i_stat)
      call memocc(i_stat, i_all, 'temp2', 'tmatimp_newsolver')

      i1 = (korbit+1)*lmmax0d*(korbit+1)*lmmax0d*ihost
      ! Allocation of temp for TMATLL
      allocate (temp(lmmaxd,lmmaxd,natomimp), stat=i_stat)
      call memocc(i_stat, product(shape(temp))*kind(temp), 'temp', 'tmatimp_newsolver')
      temp = czero
      call mpi_allreduce(tmatll,temp,i1,mpi_double_complex,mpi_sum,mpi_comm_world,ierr)
      if (ierr/=0) stop 'Error in MPI_Allreduce for TMATLL in tmatimp'
      tmatll = temp
      ! Deallocation of temp for TMATLL
      i_all = -product(shape(temp))*kind(temp)
      deallocate (temp, stat=i_stat)
      call memocc(i_stat, i_all, 'temp', 'tmatimp_newsolver')
#endif

      ! write out DTMTRX file containgin Delta_T and Delta-matrices
      if (myrank==master) then
        if (ielast==1) then
          open (unit=20, file='DTMTRX', form='FORMATTED')
          write (20, '(I5)') natomimp
          do i1 = 1, natomimp
            write (20, '(3e17.9)')(rclsimp(i2,i1), i2=1, 3)
          end do
        end if
      end if

      ! cleanup allocation
      i_all = -product(shape(irmdnew))*kind(irmdnew)
      deallocate (irmdnew, stat=i_stat)
      call memocc(i_stat, i_all, 'irmdnew', 'tmatimp_newsolver')

      !------------------------------------------------------------------------------
      ! END  calculate tmat and radial wavefunctions of host atoms
      !------------------------------------------------------------------------------

    else if (myrank==master) then

      write (*, *) 'skipping host atom loop in tmatimp_newsolver'

    end if                         ! (write_green_imp)

    !--------------------------------------------------------------------------------
    ! START calculate tmat and radial wavefunctions of impurity atoms
    !--------------------------------------------------------------------------------

    ! create new mesh before loop starts
    ! data for the new mesh
    allocate (irmdnew(natomimp), stat=i_stat)
    call memocc(i_stat, product(shape(irmdnew))*kind(irmdnew), 'irmdnew', 'tmatimp_newsolver')
    irmdnewd = 0
    do i1 = 1, natomimp
      npan_inst(i1) = ipanimp(i1) - 1
      npan_tot(i1) = npan_log + npan_eq + npan_inst(i1)
      if (npan_tot(i1)*(ncheb+1)>irmdnewd) then
        irmdnewd = npan_tot(i1)*(ncheb+1)
      end if
      irmdnew(i1) = npan_tot(i1)*(ncheb+1)
    end do
    ! new mesh
    allocate (rnew(irmdnewd,natomimp), stat=i_stat)
    call memocc(i_stat, product(shape(rnew))*kind(rnew), 'RNEW', 'tmatimp_newsolver')
    allocate (rpan_intervall(0:ntotd,natomimp), stat=i_stat)
    call memocc(i_stat, product(shape(rpan_intervall))*kind(rpan_intervall), 'RPAN_INTERVALL', 'tmatimp_newsolver')
    allocate (ipan_intervall(0:ntotd,natomimp), stat=i_stat)
    call memocc(i_stat, product(shape(ipan_intervall))*kind(ipan_intervall), 'IPAN_INTERVALL', 'tmatimp_newsolver')
    allocate (vinsnew(irmdnewd,lmpot,nspotd*natomimp), stat=i_stat)
    call memocc(i_stat, product(shape(vinsnew))*kind(vinsnew), 'VINSNEW', 'tmatimp_newsolver')

    ! initialize with zeros
    tmatllimp = czero
    tmatsph = czero

    call create_newmesh(natomimp,irm,ipand,irid,ntotd,nfund,ncheb,irmdnewd,nspin,   &
      rimp(:,1:natomimp),irminimp(1:natomimp),ipanimp(1:natomimp),                  &
      ircutimp(0:ipand,1:natomimp),r_log,npan_log,npan_eq,npan_log_at(1:natomimp),  &
      npan_eq_at(1:natomimp),npan_tot(1:natomimp),rnew(1:irmdnewd,1:natomimp),      &
      rpan_intervall(0:ntotd,1:natomimp),ipan_intervall(0:ntotd,1:natomimp),1)

    ! In second step interpolate potential
    call interpolate_poten(lpot,irm,irnsd,natomimp,ipand,lmpot,nspotd*natomimp,ntotd,        &
      irmdnewd,nspin,rimp(:,1:natomimp),irminimp(1:natomimp),irwsimp(1:natomimp),   &
      ircutimp(0:ipand,1:natomimp),vinsimp(irmind:irm,1:lmpot,1:nspin*natomimp),          &
      vm2zimp(1:irm,1:nspin*natomimp),npan_log_at(1:natomimp),npan_eq_at(1:natomimp),     &
      npan_tot(1:natomimp),rnew(1:irmdnewd,1:natomimp),                             &
      ipan_intervall(0:ntotd,1:natomimp),vinsnew)

    ! now start loop over atoms
    do i1 = i1_start_imp, i1_end_imp
      theta = thetaimp(i1)
      phi = phiimp(i1)
      ipot = nspin*(i1-1) + 1
      write (6, *) 'IMP', i1

      allocate (vnsimp(nsra*lmmaxd,nsra*lmmaxd,irmdnew(i1)), stat=i_stat)
      call memocc(i_stat, product(shape(vnsimp))*kind(vnsimp), 'VNSIMP', 'tmatimp_newsolver')
      vnsimp = czero
      ! set up the non-spherical ll' matrix for potential VLL'
      if (nsra==2) then
        use_sratrick = 1
      else if (nsra==1) then
        use_sratrick = 0
      end if
      allocate (vnspll0(lmmaxd,lmmaxd,irmdnew(i1)), stat=i_stat)
      call memocc(i_stat, product(shape(vnspll0))*kind(vnspll0), 'VNSPLL0', 'tmatimp_newsolver')
      vnspll0 = czero
      ! output potential onto which SOC is added
      allocate (vnspll1(lmmaxd,lmmaxd,irmdnew(i1)), stat=i_stat)
      call memocc(i_stat, product(shape(vnspll1))*kind(vnspll1), 'VNSPLL1', 'tmatimp_newsolver')
      vnspll1 = czero

      call vllmat(1,irmdnew(i1),irmdnew(i1),lmmax0d,lmmaxd,vnspll0,                 &
        vinsnew(1:irmdnew(i1),1:lmpot,ipot:ipot+nspin-1),lmpot,cleb,icleb,iend,     &
        nspin,zimp(i1),rnew(1:irmdnew(i1),i1),use_sratrick,ncleb)

      ! Contruct the spin-orbit coupling hamiltonian and add to potential
      call spinorbit_ham(lmax,lmmax0d,vinsnew(1:irmdnew(i1),1:lmpot,ipot:ipot+nspin-1),&
        rnew(1:irmdnew(i1),i1),eryd,zimp(i1),cvlight,t_imp%socscale(i1),nspin,   &
        lmpot,theta,phi,ipan_intervall(0:ntotd,i1),rpan_intervall(0:ntotd,i1),      &
        npan_tot(i1),ncheb,irmdnew(i1),irmdnew(i1),vnspll0,vnspll1,'1')

      ! extend matrix for the SRA treatment
      if (nsra==2) then
        allocate (vnspll(2*lmmaxd,2*lmmaxd,irmdnew(i1)), stat=i_stat)
        call memocc(i_stat, product(shape(vnspll))*kind(vnspll), 'VNSPLL', 'tmatimp_newsolver')
        if (use_sratrick==0) then
          call vllmatsra(vnspll1,vnspll,rnew(1:irmdnew(i1),i1),lmmaxd,irmdnew(i1), &
            irmdnew(i1),eryd,lmax,0,'Ref=0')
        else if (use_sratrick==1) then
          call vllmatsra(vnspll1,vnspll,rnew(1:irmdnew(i1),i1),lmmaxd,irmdnew(i1), &
            irmdnew(i1),eryd,lmax,0,'Ref=Vsph')
        end if
      else
        allocate (vnspll(lmmaxd,lmmaxd,irmdnew(i1)), stat=i_stat)
        call memocc(i_stat, product(shape(vnspll))*kind(vnspll), 'VNSPLL', 'tmatimp_newsolver')
        vnspll(:, :, :) = vnspll1(:, :, :)
      end if

      ! calculate the source terms in the Lippmann-Schwinger equation
      ! these are spherical hankel and bessel functions
      allocate (hlk(1:4*(lmax+1),irmdnew(i1)), stat=i_stat)
      call memocc(i_stat, product(shape(hlk))*kind(hlk), 'HLK', 'tmatimp_newsolver')
      allocate (jlk(1:4*(lmax+1),irmdnew(i1)), stat=i_stat)
      call memocc(i_stat, product(shape(jlk))*kind(jlk), 'JLK', 'tmatimp_newsolver')
      allocate (hlk2(1:4*(lmax+1),irmdnew(i1)), stat=i_stat)
      call memocc(i_stat, product(shape(hlk2))*kind(hlk2), 'HLK2', 'tmatimp_newsolver')
      allocate (jlk2(1:4*(lmax+1),irmdnew(i1)), stat=i_stat)
      call memocc(i_stat, product(shape(jlk2))*kind(jlk2), 'JLK2', 'tmatimp_newsolver')
      hlk = czero
      jlk = czero
      hlk2 = czero
      jlk2 = czero
      gmatprefactor = czero
      call rllsllsourceterms(nsra,nvec,eryd,rnew(1:irmdnew(i1),i1),irmdnew(i1),     &
        irmdnew(i1),lmax,lmmaxd,1,jlk_index,hlk, jlk,hlk2,jlk2,gmatprefactor)
      ! using spherical potential as reference
      if (use_sratrick==1) then
        call calcsph(nsra,irmdnew(i1),irmdnew(i1),lmax,nspin,zimp(i1),eryd,lmpot,   &
          lmmaxd,rnew(1:irmdnew(i1),i1),                                           &
          vinsnew(1:irmdnew(i1),1:lmpot,ipot:ipot+nspin-1),ncheb,npan_tot(i1),      &
          rpan_intervall(0:ntotd,i1),jlk_index,hlk,jlk,hlk2,jlk2,gmatprefactor,     &
          tmatsph,dummy_alpha,use_sratrick,.true.)
      end if

      ! calculate the tmat and wavefunctions
      allocate (rll(nvec*lmmaxd,lmmaxd,irmdnewd), stat=i_stat)
      call memocc(i_stat, product(shape(rll))*kind(rll), 'RLL', 'tmatimp_newsolver')
      allocate (sll(nvec*lmmaxd,lmmaxd,irmdnewd), stat=i_stat)
      call memocc(i_stat, product(shape(sll))*kind(sll), 'SLL', 'tmatimp_newsolver')
      rll = czero
      sll = czero

      ! Right solutions
      call rllsll(rpan_intervall(0:ntotd,i1),rnew(1:irmdnew(i1),i1),vnspll,rll,sll, &
        tmatllimp(:,:,i1),ncheb,npan_tot(i1),lmmaxd,nvec*lmmaxd,4*(lmax+1),       &
        irmdnew(i1),nsra,jlk_index,hlk,jlk,hlk2,jlk2,gmatprefactor,'1','1','0',     &
        use_sratrick,dummy_alphaget)
      if (nsra==2) then
        do ir = 1, irmdnew(i1)
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              rll(lm1+lmmaxd, lm2, ir) = rll(lm1+lmmaxd, lm2, ir)/cvlight
              sll(lm1+lmmaxd, lm2, ir) = sll(lm1+lmmaxd, lm2, ir)/cvlight
            end do
          end do
        end do
      end if

      ! for OPERATOR option save impurity wavefuncitons
      if (write_pkkr_operators) then
        t_imp%rllimp(:, :, :, i1) = rll(:, :, :)
      end if

      ! add spherical contribution of tmatrix
      if (use_sratrick==1) then
        do lm1 = 1, (korbit+1)*lmmax0d
          tmatllimp(lm1, lm1, i1) = tmatllimp(lm1, lm1, i1) + tmatsph(jlk_index(lm1))
        end do
      end if

      ! rotate tmatrix and radial wavefunction to global frame
      call rotatematrix(tmatllimp(:,:,i1), theta, phi, lmmax0d, 0)

      ! create SRA potential for impurity
      ! set up the non-spherical ll' matrix for potential VLL'
      vnspll0 = czero
      call vllmat(1,irmdnew(i1),irmdnew(i1),lmmax0d,lmmaxd,vnspll0,                 &
        vinsnew(1:irmdnew(i1),1:lmpot,ipot:ipot+nspin-1),lmpot,cleb,icleb,iend,     &
        nspin,zimp(i1),rnew(1:irmdnew(i1),i1),0,ncleb)
      ! +             ZIMP(I1),RNEW(:,I1),USE_SRATRICK)

      ! contruct the spin-orbit coupling hamiltonian and add to potential
      call spinorbit_ham(lmax,lmmax0d,vinsnew(1:irmdnew(i1),1:lmpot,ipot:ipot+nspin-1),&
        rnew(1:irmdnew(i1),i1),eryd,zimp(i1),cvlight,t_imp%socscale(i1),nspin,   &
        lmpot,theta,phi,ipan_intervall(0:ntotd,i1),rpan_intervall(0:ntotd,i1),      &
        npan_tot(i1),ncheb,irmdnew(i1),irmdnew(i1),vnspll0,vnspll1,'1')
      do ir = 1, irmdnew(i1)
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            vnsimp(lm1, lm2, ir) = vnspll1(lm1, lm2, ir)
            if (nsra==2) then
              vnsimp(lm1+lmmaxd, lm2+lmmaxd, ir) = vnspll1(lm1, lm2, ir)
            end if
          end do
        end do
      end do

      ! calculate delta_t_imp matrix written in TMATLLIMP
      do i2 = 1, ihost
        if (atomimp(i1)==hostimp(i2)) then
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              tmatllimp(lm1, lm2, i1) = tmatllimp(lm1, lm2, i1) - tmatll(lm1, lm2, i2)
            end do
          end do
          do lm1 = 1, nsra*lmmaxd
            do lm2 = 1, nsra*lmmaxd
              do ir = 1, irmdnew(i1)
                vnsimp(lm1, lm2, ir) = vnsimp(lm1, lm2, ir) - vnshost(lm1, lm2, i2, ir)
              end do
            end do
          end do
        end if
      end do

      ! calculate delta matrix \delta=int{R_imp*\deltaV*R_host}
      if (ielast==1) then
        allocate (deltabg(lmmaxd,lmmaxd,irmdnew(i1)), stat=i_stat)
        call memocc(i_stat, product(shape(deltabg))*kind(deltabg), 'DELTABG', 'tmatimp_newsolver')
        allocate (deltasm(lmmaxd,lmmaxd,irmdnew(i1)), stat=i_stat)
        call memocc(i_stat, product(shape(deltasm))*kind(deltasm), 'DELTASM', 'tmatimp_newsolver')
        deltabg = czero
        deltasm = czero
        allocate (deltatmp(irmdnew(i1)), stat=i_stat)
        call memocc(i_stat, product(shape(deltatmp))*kind(deltatmp), 'DELTATMP', 'tmatimp_newsolver')
        allocate (radialhost(lmmaxd,lmmaxd), stat=i_stat)
        call memocc(i_stat, product(shape(radialhost))*kind(radialhost), 'RADIALHOST', 'tmatimp_newsolver')
        allocate (radialimp(lmmaxd,lmmaxd), stat=i_stat)
        call memocc(i_stat, product(shape(radialimp))*kind(radialimp), 'RADIALIMP', 'tmatimp_newsolver')
        allocate (vllimp(lmmaxd,lmmaxd), stat=i_stat)
        call memocc(i_stat, product(shape(vllimp))*kind(vllimp), 'VLLIMP', 'tmatimp_newsolver')
        deltatmp = czero

        ! big component for SRA stored in DELTABG
        do ir = 1, irmdnew(i1)
          radialhost = czero
          radialimp = czero
          vllimp = czero
          deltav = czero
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              do i2 = 1, ihost
                if (atomimp(i1)==hostimp(i2)) then
                  radialhost(lm1, lm2) = rllhost(lm1, lm2, i2, ir)
                end if
              end do               ! I2
              radialimp(lm1, lm2) = rll(lm1, lm2, ir)
              vllimp(lm1, lm2) = vnsimp(lm1, lm2, ir)
            end do                 ! LM2
          end do                   ! LM1

          call zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cone,vllimp,lmmaxd,radialimp, &
            lmmaxd,czero,deltav,lmmaxd)
          call zgemm('C','N',lmmaxd,lmmaxd,lmmaxd,cone,radialhost,lmmaxd,deltav,&
            lmmaxd,czero,deltabg(:,:,ir),lmmaxd)

          ! small component for SRA stored in DELTASM
          if (nsra==2) then
            radialhost = czero
            radialimp = czero
            vllimp = czero
            deltav = czero
            do lm1 = 1, lmmaxd
              do lm2 = 1, lmmaxd
                do i2 = 1, ihost
                  if (atomimp(i1)==hostimp(i2)) then
                    radialhost(lm1, lm2) = rllhost(lm1+lmmaxd, lm2, i2, ir)
                  end if
                end do
                radialimp(lm1, lm2) = rll(lm1+lmmaxd, lm2, ir)
                vllimp(lm1, lm2) = vnsimp(lm1+lmmaxd, lm2+lmmaxd, ir)
              end do
            end do
            call zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cone,vllimp,lmmaxd,         &
              radialimp,lmmaxd,czero,deltav,lmmaxd)
            call zgemm('C','N',lmmaxd,lmmaxd,lmmaxd,cone,radialhost,lmmaxd,     &
              deltav,lmmaxd,czero,deltasm(:,:,ir),lmmaxd)

            ! sum up big and small component stored in DELTABG
            do lm1 = 1, lmmaxd
              do lm2 = 1, lmmaxd
                deltabg(lm1, lm2, ir) = deltabg(lm1, lm2, ir) + deltasm(lm1, lm2, ir)
              end do
            end do

          end if                   ! NSRA
        end do                     ! IR

        ! integrate
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            do ir = 1, irmdnew(i1)
              deltatmp(ir) = deltabg(lm1, lm2, ir)
            end do
            call intcheb_cell(deltatmp,deltamtr(lm1,lm2,i1),                        &
              rpan_intervall(0:ntotd,i1),ipan_intervall(0:ntotd,i1),npan_tot(i1),   &
              ncheb,irmdnew(i1))
          end do
        end do

        i_all = -product(shape(deltatmp))*kind(deltatmp)
        deallocate (deltatmp, stat=i_stat)
        call memocc(i_stat, i_all, 'DELTATMP', 'tmatimp_newsolver')
        i_all = -product(shape(radialhost))*kind(radialhost)
        deallocate (radialhost, stat=i_stat)
        call memocc(i_stat, i_all, 'RADIALHOST', 'tmatimp_newsolver')
        i_all = -product(shape(radialimp))*kind(radialimp)
        deallocate (radialimp, stat=i_stat)
        call memocc(i_stat, i_all, 'RADIALIMP', 'tmatimp_newsolver')
        i_all = -product(shape(vllimp))*kind(vllimp)
        deallocate (vllimp, stat=i_stat)
        call memocc(i_stat, i_all, 'VLLIMP', 'tmatimp_newsolver')
        i_all = -product(shape(deltabg))*kind(deltabg)
        deallocate (deltabg, stat=i_stat)
        call memocc(i_stat, i_all, 'DELTABG', 'tmatimp_newsolver')
        i_all = -product(shape(deltasm))*kind(deltasm)
        deallocate (deltasm, stat=i_stat)
        call memocc(i_stat, i_all, 'DELTASM', 'tmatimp_newsolver')

      end if                       ! IELAST.EQ.1

      i_all = -product(shape(vnspll0))*kind(vnspll0)
      deallocate (vnspll0, stat=i_stat)
      call memocc(i_stat, i_all, 'VNSPLL0', 'tmatimp_newsolver')
      i_all = -product(shape(vnspll1))*kind(vnspll1)
      deallocate (vnspll1, stat=i_stat)
      call memocc(i_stat, i_all, 'VNSPLL1', 'tmatimp_newsolver')
      i_all = -product(shape(vnspll))*kind(vnspll)
      deallocate (vnspll, stat=i_stat)
      call memocc(i_stat, i_all, 'VNSPLL', 'tmatimp_newsolver')
      i_all = -product(shape(hlk))*kind(hlk)
      deallocate (hlk, stat=i_stat)
      call memocc(i_stat, i_all, 'HLK', 'tmatimp_newsolver')
      i_all = -product(shape(jlk))*kind(jlk)
      deallocate (jlk, stat=i_stat)
      call memocc(i_stat, i_all, 'JLK', 'tmatimp_newsolver')
      i_all = -product(shape(hlk2))*kind(hlk2)
      deallocate (hlk2, stat=i_stat)
      call memocc(i_stat, i_all, 'HLK2', 'tmatimp_newsolver')
      i_all = -product(shape(jlk2))*kind(jlk2)
      deallocate (jlk2, stat=i_stat)
      call memocc(i_stat, i_all, 'JLK2', 'tmatimp_newsolver')
      i_all = -product(shape(rll))*kind(rll)
      deallocate (rll, stat=i_stat)
      call memocc(i_stat, i_all, 'RLL', 'tmatimp_newsolver')
      i_all = -product(shape(sll))*kind(sll)
      deallocate (sll, stat=i_stat)
      call memocc(i_stat, i_all, 'SLL', 'tmatimp_newsolver')
      i_all = -product(shape(vnsimp))*kind(vnsimp)
      deallocate (vnsimp, stat=i_stat)
      call memocc(i_stat, i_all, 'VNSIMP', 'tmatimp_newsolver')

    end do                         ! I1 impurity

    i_all = -product(shape(rnew))*kind(rnew)
    deallocate (rnew, stat=i_stat)
    call memocc(i_stat, i_all, 'RNEW', 'tmatimp_newsolver')
    i_all = -product(shape(vinsnew))*kind(vinsnew)
    deallocate (vinsnew, stat=i_stat)
    call memocc(i_stat, i_all, 'VINSNEW', 'tmatimp_newsolver')
    i_all = -product(shape(rpan_intervall))*kind(rpan_intervall)
    deallocate (rpan_intervall, stat=i_stat)
    call memocc(i_stat, i_all, 'RPAN_INTERVALL', 'tmatimp_newsolver')
    i_all = -product(shape(ipan_intervall))*kind(ipan_intervall)
    deallocate (ipan_intervall, stat=i_stat)
    call memocc(i_stat, i_all, 'IPAN_INTERVALL', 'tmatimp_newsolver')

    !--------------------------------------------------------------------------------
    ! END calculate tmat and radial wavefunctions of impurity atoms
    !--------------------------------------------------------------------------------

    ! final writeout only on master
#ifdef CPP_MPI
    ! collect results and write out only on master
    ! collect TMATLLIMP, DELTAMTR
    ! communicate TMATLLIMP, DELTAMTR
    i1 = lmmaxd*lmmaxd*natomimp
    allocate (temp(lmmaxd,lmmaxd,natomimp), stat=i_stat)
    call memocc(i_stat, product(shape(temp))*kind(temp), 'temp', 'tmatimp_newsolver')
    temp = czero
    call mpi_allreduce(tmatllimp, temp, i1, mpi_double_complex, mpi_sum, mpi_comm_world, ierr)
    if (ierr/=0) stop 'Error in MPI_Allreduce for TMATLLIMP in tmatimp'
    tmatllimp = temp
    i_all = -product(shape(temp))*kind(temp)
    deallocate (temp, stat=i_stat)
    call memocc(i_stat, i_all, 'temp', 'tmatimp_newsolver')

    i1 = lmmaxd*lmmaxd*natomimp
    allocate (temp(lmmaxd,lmmaxd,natomimp), stat=i_stat)
    call memocc(i_stat, product(shape(temp))*kind(temp), 'temp', 'tmatimp_newsolver')
    temp = czero
    call mpi_allreduce(deltamtr, temp, i1, mpi_double_complex, mpi_sum, mpi_comm_world, ierr)
    if (ierr/=0) stop 'Error in MPI_Allreduce for DELTAMTR in tmatimp'
    deltamtr = temp
    i_all = -product(shape(temp))*kind(temp)
    deallocate (temp, stat=i_stat)
    call memocc(i_stat, i_all, 'temp', 'tmatimp_newsolver')
#endif

    ! collect results and writeout only for GREENIMP option
    if (write_green_imp .and. myrank==master) then

      do i1 = 1, natomimp
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            il1 = lmmaxd*(i1-1) + lm1
            il2 = lmmaxd*(i1-1) + lm2
            dtmtrx(il1, il2) = tmatllimp(lm1, lm2, i1)
          end do
        end do
      end do

      ! write down to the file DTMTRX
      if (ielast==1) then
        allocate (deltaimp((korbit+1)*lmmax0d*natomimp,(korbit+1)*lmmax0d*natomimp), stat=i_stat)
        call memocc(i_stat, product(shape(deltaimp))*kind(deltaimp), 'DELTAIMP', 'tmatimp_newsolver')
        deltaimp = czero
        do i1 = 1, natomimp
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              il1 = lmmaxd*(i1-1) + lm1
              il2 = lmmaxd*(i1-1) + lm2
              deltaimp(il1, il2) = deltamtr(lm1, lm2, i1)
            end do
          end do
        end do
        do lm1 = 1, lmmaxd*natomimp
          do lm2 = 1, lmmaxd*natomimp
            write (20, '((2I5),(4e17.9))') lm2, lm1, dtmtrx(lm2, lm1), deltaimp(lm2, lm1)
          end do
        end do
        i_all = -product(shape(deltaimp))*kind(deltaimp)
        deallocate (deltaimp, stat=i_stat)
        call memocc(i_stat, i_all, 'DELTAIMP', 'tmatimp_newsolver')
        if (myrank==master) write (6, *) 'end of delta t'
      end if                       ! IELAST.EQ.1

      close (20)                   ! output file DTMTRX

    end if                         ! myrank==master

    i_all = -product(shape(vnshost))*kind(vnshost)
    deallocate (vnshost, stat=i_stat)
    call memocc(i_stat, i_all, 'VNSHOST', 'tmatimp_newsolver')
    i_all = -product(shape(rllhost))*kind(rllhost)
    deallocate (rllhost, stat=i_stat)
    call memocc(i_stat, i_all, 'RLLHOST', 'tmatimp_newsolver')
    i_all = -product(shape(tmatllimp))*kind(tmatllimp)
    deallocate (tmatllimp, stat=i_stat)
    call memocc(i_stat, i_all, 'TMATLLIMP', 'tmatimp_newsolver')
    i_all = -product(shape(deltamtr))*kind(deltamtr)
    deallocate (deltamtr, stat=i_stat)
    call memocc(i_stat, i_all, 'DELTAMTR', 'tmatimp_newsolver')
    i_all = -product(shape(deltav))*kind(deltav)
    deallocate (deltav, stat=i_stat)
    call memocc(i_stat, i_all, 'DELTAV', 'tmatimp_newsolver')
    i_all = -product(shape(irmdnew))*kind(irmdnew)
    deallocate (irmdnew, stat=i_stat)
    call memocc(i_stat, i_all, 'irmdnew', 'tmatimp_newsolver')
  end subroutine tmatimp_newsolver

end module mod_tmatimp_newsolver
