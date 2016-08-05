!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module mod_calconfs

  use mod_ioformat, only: MODE_INT, MODE_VIS
  implicit none

  private
  public :: calc_on_fsurf_inputcard, calc_on_fsurf, get_nsqa, calc_spinvalue_state, calculate_spinmixing_int, calculate_spinmixing_vis, calculate_response_functions_CRTA_int, calculate_response_functions_CRTA_vis, order_lines, calculate_dos_int

  integer, parameter :: SPCZMAX=1, SPCXY0=2, ROT_NO=0, ROT_FULLBZ=1, ROT_SPEC=2

  type :: cfg_TYPE

    integer :: N1 = 16
    integer :: lspin=-1
    integer :: lfvel=-1
    integer :: lrashba=-1
    integer :: lspinperatom=-1
    integer :: ltorqperatom=-1
    integer :: ltorq=-1
    integer :: lspinflux=-1
    integer :: lalpha=-1
    integer :: nsqa=-1
    integer :: mode=-1
    integer :: rotatemode=-1
    integer :: ilayer=-1
    double precision :: dk_fv = -1d0
    logical :: saveeigv=.false.
    logical :: simpson=.false.

    integer :: N2 = 3
    integer,          allocatable :: ispincomb(:)
    double precision, allocatable :: nvect(:,:)

  end type cfg_TYPE

  logical, save :: cfg_read=.false.

  type(cfg_type), save :: cfg




contains



  subroutine calc_on_fsurf_inputcard(inc,lattice,cluster,tgmatrx,nkpts,kpoints)

    use type_inc,       only: inc_type
    use type_data,      only: lattice_type, cluster_type, tgmatrx_type
    use mod_symmetries, only: symmetries_type, set_symmetries, rotate_kpoints, expand_visarrays, expand_areas
    use mod_read,       only: read_kpointsfile_vis, read_kpointsfile_int
    use mod_parutils,   only: distribute_linear_on_tasks
    use mod_mympi,      only: myrank, nranks, master
    use mod_ioformat,   only: filemode_vis, filemode_int, filemode_rot, filename_vtkonfs, ext_vtkxml, fmt_fn_sub_ext
    use mod_iohelp,     only: getBZvolume
    use mod_vtkxml,     only: write_pointdata_rot
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    type(tgmatrx_type), intent(in) :: tgmatrx

    integer,                       intent(inout) :: nkpts
    double precision, allocatable, intent(inout) :: kpoints(:,:)

    integer                       :: nkpts1, nkpts_all1, nkpts_all
    integer,          allocatable :: kpt2irr1(:), irr2kpt1(:), kpt2irr(:), irr2kpt(:), kpt2irr_ord(:), band_indices(:)
    double precision, allocatable :: kpoints1(:,:)


    type(symmetries_type) :: symmetries
    double precision, allocatable :: fermivel(:,:), spinval(:,:,:), spinmix(:), torqval(:,:,:),&
                                   & torqval_atom(:,:,:,:), spinvec(:,:,:,:), spinvec_atom(:,:,:,:),&
                                   & spinflux_atom(:,:,:,:), alphaval(:,:,:)
    double precision, allocatable :: areas1(:), areas(:)

    integer :: nsym, ii, iatom, ierr, ikp, isqa
    double precision :: integ, dos, vabs, BZVol
    integer, allocatable :: isym(:)
    character(len=256) :: filename

    integer                       :: nscalar, nvector, iscalar, ivector
    double precision, allocatable :: scalardata(:,:), &
                                   & vectordata(:,:,:)
    character(len=256) :: filemode
    character(len=256), allocatable :: scalarstring(:), vectorstring(:)

    call set_symmetries(inc, lattice, symmetries)
    BZVol = getBZvolume(lattice%recbv)

    if(.not.cfg_read) call read_cfg()

    if(cfg%lspin/=1 .and. cfg%lfvel/=1) return

    !read in the (irreducible) k-points and apply symmetry operations to expand to full BZ
    select case (cfg%mode)

    case (MODE_VIS)
      filemode = filemode_vis
      if(cfg%rotatemode==ROT_FULLBZ)then
        call read_kpointsfile_vis(nkpts_all1, nkpts1, kpoints1, nsym, isym, kpt2irr1, irr2kpt1)
        call rotate_kpoints(symmetries%rotmat, nkpts1, kpoints1, nsym, isym, nkpts, kpoints)
        call expand_visarrays(nsym, nkpts_all1, nkpts1, kpt2irr1, irr2kpt1, kpt2irr, irr2kpt)
        nkpts_all = nsym*nkpts_all1
        deallocate(kpt2irr1, irr2kpt1,kpoints1)
      elseif(cfg%rotatemode==ROT_NO)then
        call read_kpointsfile_vis(nkpts_all, nkpts, kpoints, nsym, isym, kpt2irr, irr2kpt)
      else
        stop 'cfg%rotatemode not known!'
      end if

    case (MODE_INT)
      filemode = filemode_int
      if(cfg%rotatemode==ROT_FULLBZ)then
        call read_kpointsfile_int(nkpts1, kpoints1, areas1, nsym, isym)
        call rotate_kpoints(symmetries%rotmat, nkpts1, kpoints1, nsym, isym, nkpts, kpoints)
        call expand_areas(nsym,nkpts1,areas1,areas)
        deallocate(areas1,kpoints1)
      elseif(cfg%rotatemode==ROT_NO)then
        call read_kpointsfile_int(nkpts, kpoints, areas, nsym, isym)
      else
        stop 'cfg%rotatemode not known!'
      end if

    case default; stop 'case not known in select case (cfg%mode)'
    end select

    if(cfg%rotatemode==ROT_FULLBZ)then
      !redefine symmetries (set to unit transformation only)
      nsym = 1
      deallocate(isym)
      allocate(isym(nsym))
      isym = (/ 1 /)
    end if

    !calculate the properties
    nvector=0
    nscalar=0

    if(cfg%lfvel==1)then
      nvector=nvector+1
      allocate(fermivel(3,nkpts), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating fermivel'
    endif!cfg%lfvel==1

    if(cfg%lalpha==1)then
      nvector=nvector+1
      allocate(alphaval(3,inc%ndegen,nkpts), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating alphaval'
    endif!cfg%lalpha==1

    if(cfg%lspin==1)then
      nscalar=nscalar+cfg%nsqa
      nvector=nvector+cfg%nsqa
      allocate( spinval(inc%ndegen,cfg%nsqa,nkpts),   &
              & spinvec(3,inc%ndegen,cfg%nsqa,nkpts), &
              & STAT=ierr                             )
      if(ierr/=0) stop 'Problem allocating spinval and spinvec'
    endif!cfg%lspin==1

    if(cfg%lspinperatom==1)then
      nvector = nvector+inc%natypd
      allocate( spinvec_atom(3,inc%natypd,inc%ndegen,nkpts), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating spinvec_atom etc.'
    endif!cfg%lspinperatom==1

    if(cfg%ltorq==1)then
      nvector = nvector+1
      allocate(torqval(3,inc%ndegen,nkpts), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating torqval'
    endif!cfg%ltorq==1

    if(cfg%ltorqperatom==1)then
      nvector = nvector+inc%natypd
      allocate( torqval_atom(3,inc%natypd,inc%ndegen,nkpts), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating torqval_atom'
    endif!cfg%ltorqperatom==1

    if(cfg%lspinflux==1)then
      nvector = nvector+inc%natypd
      allocate( spinflux_atom(3,inc%natypd,inc%ndegen,nkpts), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating spinflux_atom'
    endif!cfg%lspinflux==1

    call calc_on_fsurf( inc,lattice,cluster,tgmatrx,nkpts,kpoints,                       &
                      & fermivelocity=fermivel,spinvalue=spinval,save_eigv=cfg%saveeigv, &
                      & spinvec=spinvec,spinvec_atom=spinvec_atom,torqvalue=torqval,     &
                      & torqvalue_atom=torqval_atom,spinflux_atom=spinflux_atom,         &
                      & alphavalue=alphaval                                              )

    allocate( scalarstring(nscalar), vectorstring(nvector), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating scalarstring etc.'


    !save the calculated properties to a file
    if(cfg%lfvel==1   .and. myrank==master) call save_fermivelocity(trim(filemode), nkpts, fermivel, nsym, isym)
    if(cfg%lalpha==1  .and. myrank==master) call save_alpha(trim(filemode), nkpts, inc%ndegen, alphaval, nsym, isym)
    if(cfg%lfvel==1   .and. myrank==master .and. cfg%mode==MODE_INT) call calculate_and_save_weights(nkpts,nsym,isym,areas,fermivel)
    if(cfg%lspin==1   .and. myrank==master)then
       call save_spinvalue(trim(filemode), nkpts, cfg%nsqa, inc%ndegen, cfg%ispincomb, cfg%nvect, spinval, nsym, isym)
       call save_spinvec  (trim(filemode), nkpts, cfg%nsqa, inc%ndegen, cfg%ispincomb, cfg%nvect, spinvec, nsym, isym)
    end if!cfg%lspin==1 .and. myrank==master
    if(cfg%lspinperatom==1 .and. myrank==master) call save_spinvec_atom(trim(filemode), nkpts, inc%ndegen, inc%natypd, spinvec_atom, nsym, isym)
    if(cfg%ltorq==1        .and. myrank==master) call save_torqvalue(trim(filemode), nkpts, inc%ndegen, torqval, nsym, isym)
    if(cfg%ltorqperatom==1 .and. myrank==master) call save_torqvalue_atom(trim(filemode), nkpts, inc%ndegen, inc%natypd, torqval_atom, nsym, isym)
    if(cfg%lspinflux==1    .and. myrank==master) call save_spinflux_atom(trim(filemode), nkpts, inc%ndegen, inc%natypd, spinflux_atom, nsym, isym)


    !calculate averages
    if(myrank==master .and. cfg%lfvel==1 .and. cfg%lspin==1) then

      allocate(spinmix(cfg%nsqa), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating spinmix'

      if(cfg%mode==MODE_INT) call calculate_spinmixing_int(nsym,inc%ndegen,cfg%nsqa,nkpts,areas,fermivel,spinval,spinmix,dos,.true.,BZVol)
      if(cfg%mode==MODE_VIS) call calculate_spinmixing_vis(nsym,inc%ndegen,cfg%nsqa,nkpts,nkpts_all,kpt2irr,kpoints,fermivel,spinval,spinmix,dos,.true.,BZVol,inc%nBZdim)

    elseif(myrank==master .and. cfg%lfvel==1 .and. cfg%lspin==0) then

      if(cfg%mode==MODE_INT) call calculate_dos_int(nsym,isym,symmetries%rotmat,lattice%alat,BZVol,nkpts,areas,fermivel,inc%nBZdim)
      if(cfg%mode==MODE_VIS) call calculate_dos_vis(nsym,nkpts,nkpts_all,kpt2irr,kpoints,fermivel,dos,.true.,BZVol)

    end if

!    if(myrank==master .and. cfg%lfvel==1 .and. cfg%ltorq==1) then
    if(myrank==master.and.cfg%lfvel==1) then
      if(cfg%mode==MODE_INT) call calculate_response_functions_CRTA_int(nsym,isym,symmetries%rotmat,lattice%alat,BZVol,inc%ndegen,inc%natypd,nkpts,areas,fermivel,torqval,torqval_atom,spinvec_atom,spinflux_atom)
!BZcomment: order changed to avoid simpson=T and 3D-case
      if(cfg%mode==MODE_VIS)then
        if(cfg%simpson .and. inc%nBZdim==2)then
          write(*,*) 'Simpson rule is implemented only for torkance and damping. Only these quantities will be computed.'
          call calculate_torkance_CRTA_vis_simpson2D(nsym,isym,symmetries%rotmat,lattice%alat,BZVol,inc%ndegen,nkpts,nkpts_all,kpt2irr,irr2kpt,kpoints,fermivel,torqval)
        else!cfg%simpson==.true.
          if(cfg%simpson) write(*,*) 'Simpson rule is not implemented for inc%nBZdim =/= 2. Taking standard (=linear) integration.'
          call calculate_response_functions_CRTA_vis(inc%nBZdim,nsym,isym,symmetries%rotmat,lattice%alat,BZVol,inc%ndegen,inc%natypd,nkpts,nkpts_all,kpt2irr,irr2kpt,kpoints,fermivel,torqval,torqval_atom,spinvec_atom,spinflux_atom)
        end if!cfg%simpson==.true.
      end if!cfg%mode==MODE_VIS
    end if!myrank==master

    !output 3D visualization data
    if(myrank==master .and. cfg%mode==MODE_VIS .and. inc%nBZdim==3) then

      !===================================================================
      !allocate arrays
      if(nvector>0) then
        allocate(vectordata(3,nkpts,nvector), STAT=ierr)
        if(ierr/=0) stop 'Problem allocating vectordata'
      end if

      if(nscalar>0) then
        allocate(scalardata(nkpts,nscalar), STAT=ierr)
        if(ierr/=0) stop 'Problem allocating scalardata'
      end if
      !===================================================================

      iscalar=0
      ivector=0

      !===================================================================
      if(cfg%lfvel==1) then
        ivector = ivector+1
        vectordata(:,:,ivector) = fermivel
        vectorstring(ivector) = 'fvelocity'
      end if!cfg%lfvel
      !===================================================================

      !===================================================================
      if(cfg%lalpha==1) then
        ivector = ivector+1
        vectordata(:,:,ivector) = alphaval(:,1,:)
        vectorstring(ivector) = 'alpha'
      end if!cfg%lalpha
      !===================================================================


      !===================================================================
      if(cfg%lspin==1)then


        !*****************************!
        !*** save scalar spinvalue ***!
        if(inc%ndegen==2) then
          do isqa=1,cfg%nsqa
            iscalar = iscalar+1
            scalardata(:,iscalar) = (1d0-spinval(1,isqa,:))/2
            write(scalarstring(iscalar),'(A,I0)') 'Eyaf_',isqa
          end do!isqa
        else!inc%ndegen==2
          do isqa=1,cfg%nsqa
            iscalar = iscalar+1
            scalardata(:,iscalar) = spinval(1,isqa,:)
            write(scalarstring(iscalar),'(A,I0)') 'Spin_',isqa
          end do!isqa
        end if!inc%ndegen==2
        !*** save scalar spinvalue ***!
        !*****************************!


        !************************!
        !*** save spin vector ***!
        do isqa=1,cfg%nsqa
          ivector = ivector+1
          vectordata(:,:,ivector) = spinvec(:,1,isqa,:)
          write(vectorstring(ivector),'(A,I0)') 'Spinv_',isqa
        end do!isqa
        !*** save spin vector ***!
        !************************!

      end if!cfg%lspin
      !===================================================================



      !===================================================================
      if(cfg%lspinperatom==1)then
        do iatom=1,inc%natyp
          ivector = ivector+1
          vectordata(:,:,ivector) = spinvec_atom(:,iatom,1,:)
          write(vectorstring(ivector),'(A,I0)') 'Spinvec_atom_', iatom
        end do!iatom
      end if!cfg%lspinperatom==1
      !===================================================================



      !===================================================================
      if(cfg%ltorq==1)then
        ivector = ivector+1
        vectordata(:,:,ivector) = torqval(:,1,:)
        vectorstring(ivector) = 'torq'
      end if!cfg%ltorq==1
      !===================================================================


      !===================================================================
      if(cfg%ltorqperatom==1)then
        do iatom=1,inc%natyp
          ivector = ivector+1
          vectordata(:,:,ivector) = torqval_atom(:,iatom,1,:)
          write(vectorstring(ivector),'(A,I0)') 'torq_atom_', iatom
        end do!iatom
      end if!cfg%ltorqperatom==1
      !===================================================================

      !===================================================================
      if(cfg%lspinflux==1)then
        do iatom=1,inc%natyp
          ivector = ivector+1
          vectordata(:,:,ivector) = spinflux_atom(:,iatom,1,:)
          write(vectorstring(ivector),'(A,I0)') 'spinflux_atom_', iatom
        end do!iatom
      end if!cfg%lspinflux==1
      !===================================================================

      write(filename,fmt_fn_sub_ext) filename_vtkonfs, filemode_rot, ext_vtkxml
      call write_pointdata_rot( trim(filename),nkpts,kpoints,   &
                              & nscalar,scalardata,scalarstring,&
                              & nvector,vectordata,vectorstring,&
                              & nsym,symmetries%rotmat,isym,    &
                              & nkpts_all, kpt2irr              )

    end if!myrank==master .and. cfg%mode==MODE_VIS

  end subroutine calc_on_fsurf_inputcard





  subroutine calc_on_fsurf(inc,lattice,cluster,tgmatrx,nkpts,kpoints,fermivelocity,spinvalue,save_eigv,spinvec,spinvec_atom,torqvalue,torqvalue_atom,spinflux_atom,alphavalue)

    use type_inc,       only: inc_type
    use type_data,      only: lattice_type, cluster_type, tgmatrx_type
    use mod_symmetries, only: symmetries_type
    use mod_parutils,   only: distribute_linear_on_tasks
    use mod_mympi,      only: myrank, nranks, master
    use mod_kkrmat,     only: compute_kkr_eigenvectors2
    use mod_mathtools,  only: findminindex, bubblesort
    use mod_ioformat,   only: filename_eigvect, filemode_int, filemode_vis
#ifdef CPP_MPI
    use mod_iohelp,     only: open_mpifile_setview, close_mpifile
    use mpi
#endif
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    type(tgmatrx_type), intent(in) :: tgmatrx

    integer,          intent(in) :: nkpts
    double precision, intent(in) :: kpoints(3,nkpts)
    double precision, allocatable, intent(inout) :: fermivelocity(:,:),&!fermivelocity(3,nkpts),&
                                   & spinvalue(:,:,:),       &!spinvalue(inc%ndegen,cfg%nsqa,nkpts),&
                                   & spinvec(:,:,:,:),       &!spinvec(3,inc%ndegen,cfg%nsqa,nkpts),&
                                   & spinvec_atom(:,:,:,:),  &!spinvec_atom(3,inc%natypd,inc%ndegen,nkpts),  &
                                   & torqvalue(:,:,:),       &!torqvalue(3,inc%ndegen,nkpts),                &
                                   & torqvalue_atom(:,:,:,:),&!torqvalue_atom(3,inc%natypd,inc%ndegen,nkpts),&
                                   & spinflux_atom(:,:,:,:), &!spinflux_atom(3,inc%natypd,inc%ndegen,nkpts), &
                                   & alphavalue(:,:,:)        !alphavalue(3,inc%ndegen,nkpts)
    logical,          intent(in), optional  :: save_eigv
!   double complex,   intent(out), optional :: eigenvectors(inc%lmmaxso,inc%natypd,inc%ndegen,cfg%nsqa,nkpts)

    double precision, allocatable :: fermivel_tmp(:,:),&
                                   & spinval_tmp(:,:,:),&
                                   & spinvec_tmp(:,:,:,:),&
                                   & spinvec_atom_tmp(:,:,:,:),&
                                   & torqval_tmp(:,:,:),       &
                                   & torqval_atom_tmp(:,:,:,:),&
                                   & spinflux_atom_tmp(:,:,:,:),&
                                   & alphaval_tmp(:,:,:)

    integer :: nkpt, ioff, ntot_pT(0:nranks-1), ioff_pT(0:nranks-1)
    double precision :: kpoint(3)

    integer          :: lm_fs, lm_fs2
    double precision, allocatable :: spinvec_tmp1(:,:,:),      &!spinvec_tmp1(3,inc%ndegen,cfg%nsqa),        &
                                   & spinvec_atom_tmp1(:,:,:), &!spinvec_atom_tmp1(3,inc%natypd,inc%ndegen), &
                                   & torqval_tmp1(:,:),        &!torqval_tmp1(3,inc%ndegen),                 &
                                   & torqval_atom_tmp1(:,:,:), &!torqval_atom_tmp1(3,inc%natypd,inc%ndegen), &
                                   & spinflux_atom_tmp1(:,:,:),&!spinflux_atom_tmp1(3,inc%natypd,inc%ndegen),&
                                   & alphaval_tmp1(:,:)         !alphaval_tmp1(3,inc%ndegen) 

    double complex   :: eigwEF(inc%almso),           &
                      & LVeigEF(inc%almso,inc%almso),&
                      & RVeigEF(inc%almso,inc%almso),&
                      & rveig(inc%almso,inc%ndegen)

    integer          :: printstep, ikp, ierr, i3dim(3), itmp, irecv(0:nranks-1), idispl(0:nranks-1)
    integer          :: indsort(inc%almso)
    double precision :: deigsort(inc%almso)
    character(len=256) :: filename

    logical :: l_save_eigv
#ifdef CPP_MPI
    integer :: fh_eigv, irveigel, dimens(4)
    double complex, allocatable :: rveigrot(:,:,:,:)
    integer :: mpistatus(MPI_STATUS_SIZE)
#endif

    !Checks
    if(.not.cfg_read)                                         stop 'cfg not read when entering calc_on_fsurf'
    if(allocated(torqvalue).and.cfg%nsqa>1)                     stop 'cfg%nsqa must be 1 if torquevalue is present in calc_on_fsurf'
    if(allocated(torqvalue).and..not.inc%ltorq)                 stop 'TBKKR_torq-file not present'
    if(allocated(spinvalue).and..not.inc%lrhod)                 stop 'TBKKR_rhod-file not present'
    if(allocated(torqvalue).and.(.not.allocated(spinvalue)))      stop 'For torqvalue also spinvalue must be present in calc_on_fsurf'
    if(allocated(spinvec).and.(.not.allocated(spinvalue)))        stop 'For spinvec also spinvalue must be present in calc_on_fsurf'
    if(allocated(spinvec_atom).and.(.not.allocated(spinvalue)))   stop 'For spinvec_atom also spinvalue must be present in calc_on_fsurf'
    if(allocated(torqvalue_atom).and.(.not.allocated(torqvalue))) stop 'For torqvalue_atom also torqvalue must be present in calc_on_fsurf'
    if(allocated(spinflux_atom).and..not.inc%lspinflux)         stop 'TBkkr_spinflux-file not present'
    if(cfg%lrashba==1 .and. cfg%nsqa>1)                       stop 'cfg%nsqa must be 1 if cfg%lrashba=1 in calc_on_fsurf'
    if(allocated(spinvec_atom) .and. cfg%nsqa>1)                stop 'cfg%nsqa must be 1 for spinvec_atom in calc_on_fsurf'
    if(allocated(alphavalue).and..not.inc%lalpha)               stop 'TBkkr_alpha-file not present'

    !Parallelize
    call distribute_linear_on_tasks(nranks,myrank,master,nkpts,ntot_pT,ioff_pT,.true.)
    nkpt = ntot_pT(myrank)
    ioff = ioff_pT(myrank)

    !flag whether eigenvalues shall be saved to file
    l_save_eigv = .false.
#ifdef CPP_MPI
    if(present(save_eigv))then
      if(save_eigv) l_save_eigv = .true.
    else 
      l_save_eigv = cfg%saveeigv
    end if
#endif

    !allocate task-local arrays
    if(allocated(fermivelocity))then
      allocate( fermivel_tmp(3,nkpt), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating fermivel_tmp'
      fermivelocity = 0d0
      fermivel_tmp  = 0d0
    end if

    if(allocated(spinvalue))then
      allocate( spinval_tmp(inc%ndegen,cfg%nsqa,nkpt), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating spinval_tmp'
      spinvalue   = 0d0
      spinval_tmp = 0d0
    end if

    if(allocated(spinvec))then
      allocate( spinvec_tmp(3,inc%ndegen,cfg%nsqa,nkpt), &
              & spinvec_tmp1(3,inc%ndegen,cfg%nsqa), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating spinvec_tmp'
      spinvec      = 0d0
      spinvec_tmp  = 0d0
      spinvec_tmp1 = 0d0
    end if

    if(allocated(spinvec_atom))then
      allocate( spinvec_atom_tmp(3,inc%natypd,inc%ndegen,nkpt), &
              & spinvec_atom_tmp1(3,inc%natypd,inc%ndegen), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating spinvec_atom_tmp'
      spinvec_atom      = 0d0
      spinvec_atom_tmp  = 0d0
      spinvec_atom_tmp1 = 0d0
    end if

    if(allocated(torqvalue))then
      allocate( torqval_tmp(3,inc%ndegen,nkpt), &
              & torqval_tmp1(3,inc%ndegen),STAT=ierr )
      if(ierr/=0) stop 'Problem allocating torqval_tmp'
      torqvalue    = 0d0
      torqval_tmp  = 0d0
      torqval_tmp1 = 0d0
    end if

    if(allocated(torqvalue_atom))then
      allocate( torqval_atom_tmp(3,inc%natypd,inc%ndegen,nkpt), &
              & torqval_atom_tmp1(3,inc%natypd,inc%ndegen), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating torqval_atom_tmp'
      torqvalue_atom    = 0d0
      torqval_atom_tmp  = 0d0
      torqval_atom_tmp1 = 0d0
    end if

    if(allocated(spinflux_atom))then
      allocate( spinflux_atom_tmp(3,inc%natypd,inc%ndegen,nkpt), &
              & spinflux_atom_tmp1(3,inc%natypd,inc%ndegen), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating spinflux_atom_tmp'
      spinflux_atom      = 0d0
      spinflux_atom_tmp  = 0d0
      spinflux_atom_tmp1 = 0d0
    end if

    if(allocated(alphavalue))then
      allocate( alphaval_tmp(3,inc%ndegen,nkpt), &
              & alphaval_tmp1(3,inc%ndegen), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating alphaval_tmp'
      alphavalue    = 0d0
      alphaval_tmp  = 0d0
      alphaval_tmp1 = 0d0
    end if


#ifdef CPP_MPI
    !open the mpi-io file to store the eigenvector
    if(l_save_eigv)then
      allocate( rveigrot(inc%lmmaxso,inc%natypd,inc%ndegen,cfg%nsqa), STAT=ierr ) 
      if(ierr/=0) stop 'Problem allocating rveigrot'
      dimens = (/ inc%lmmaxso,inc%natypd,inc%ndegen,cfg%nsqa /)
      irveigel = product(dimens)
      if(cfg%mode==MODE_INT) write(filename,'(A,A)') filename_eigvect, filemode_int
      if(cfg%mode==MODE_VIS) write(filename,'(A,A)') filename_eigvect, filemode_vis
      call open_mpifile_setview( trim(filename), 'write', 4, dimens, ntot_pT, MPI_DOUBLE_COMPLEX, &
                               & myrank, nranks, MPI_COMM_WORLD, fh_eigv                          )
    end if!l_save_eigv
#endif

    !******************************************!
    !*************** B E G I N ****************!
    !******************************************!
    !*** (parallelized) loop over FS points ***!
    !******************************************!
    printstep = nkpt/50
    if(printstep==0) printstep=1
    !print header of statusbar
    if(myrank==master) write(*,'("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
    if(myrank==master) write(*,FMT=190) !beginning of statusbar
    do ikp=1,nkpt
      !update statusbar
      if(mod(ikp,printstep)==0 .and. myrank==master) write(*,FMT=200)

      kpoint = kpoints(:,ikp+ioff)

      !===============================================================!
      !=== find the eigenvector and eigenvalue at the fermi energy ===!
      !===============================================================!
      call compute_kkr_eigenvectors2( inc, lattice, cluster, tgmatrx%ginp(:,:,:,1),&
                                    & tgmatrx%tmat(:,:,1), kpoint,                 &
                                    & eigwEF, LVeigEF, RVeigEF                     )
      deigsort = abs(eigwEF)
      call findminindex(inc%almso, deigsort, lm_fs)
!     write(*,'(A,2ES25.16)') 'minimal eigenvector=', eigwEF(lm_fs)



      !====================================!
      !=== Calculate the fermi-velocity ===!
      !====================================!
      if(allocated(fermivelocity))then

        call calc_fermivel( inc,lattice,cluster,tgmatrx,kpoint,   &
                              & eigwEF(lm_fs), LVeigEF(:,lm_fs),      &
                              & RVeigEF(:,lm_fs), fermivel_tmp(:,ikp) )

      end if!present(fermivelocity)



      !============================================!
      !===   Calculate the expectation values   ===!
      !============================================!

      if(allocated(spinvalue).or.allocated(spinvec).or.allocated(spinvec_atom).or.&
        &allocated(torqvalue).or.allocated(torqvalue_atom).or.allocated(alphavalue))then
        !copy the band
        rveig(:,1) = RVeigEF(:,lm_fs)

        !find degenerate band for degnerate cases
        if(inc%ndegen==2)then
          deigsort = abs(eigwEF-eigwEF(lm_fs))
          call bubblesort(inc%almso, deigsort, indsort)
          lm_fs2 = indsort(2)
          rveig(:,2) = RVeigEF(:,lm_fs2)
        end if

#ifdef CPP_MPI
        if(l_save_eigv)then
          call calc_spinvalue_state_generalized( inc, tgmatrx, rveig, &
                                     & spinval_tmp(:,:,ikp),          &
                                     & spinvec_tmp1,                  &
                                     & spinvec_atom_tmp1,             &
                                     & torqval_tmp1,                  &
                                     & torqval_atom_tmp1,             &
                                     & spinflux_atom_tmp1,            &
                                     & alphaval_tmp1,                 &
                                     & eigvect_rot=rveigrot           )
          call MPI_File_write( fh_eigv, rveigrot, irveigel,        &
                             & MPI_DOUBLE_COMPLEX, mpistatus, ierr )

        else!l_save_eigv
#endif
          call calc_spinvalue_state_generalized( inc, tgmatrx, rveig,  &
                                     & spinval_tmp(:,:,ikp),             &
                                     & spinvec=spinvec_tmp1,             &
                                     & spinvec_atom=spinvec_atom_tmp1,   &
                                     & torq_value=torqval_tmp1,          &
                                     & torq_value_atom=torqval_atom_tmp1,&
                                     & spinflux_atom=spinflux_atom_tmp1, &
                                     & alpha_value=alphaval_tmp1         )
#ifdef CPP_MPI
        end if!l_save_eigv
#endif
        if(allocated(spinvec))        spinvec_tmp(:,:,:,ikp) = spinvec_tmp1
        if(allocated(spinvec_atom))   spinvec_atom_tmp(:,:,:,ikp) = spinvec_atom_tmp1
        if(allocated(torqvalue))      torqval_tmp(:,:,ikp) = torqval_tmp1
        if(allocated(torqvalue_atom)) torqval_atom_tmp(:,:,:,ikp) = torqval_atom_tmp1
        if(allocated(spinflux_atom))  spinflux_atom_tmp(:,:,:,ikp) = spinflux_atom_tmp1
        if(allocated(alphavalue))     alphaval_tmp(:,:,ikp) = alphaval_tmp1 

      end if!allocated(spinvalue).or.allocated(spinvec).or.allocated(spinvec_atom).or.
            !allocated(torqvalue).or.allocated(torqvalue_atom).or.allocated(alphavalue)
 
    end do!ikp
    if(myrank==master) write(*,*) ''
!   write(*,*) 'after loop, myrank=', myrank
    !******************************************!
    !*** (parallelized) loop over FS points ***!
    !******************************************!
    !***************** E N D ******************!
    !******************************************!


    !close the eigenvector-file
#ifdef CPP_MPI
    if(l_save_eigv) call close_mpifile(fh_eigv)
#endif


    !gather the fermi-velocity from the processes
    if(allocated(fermivelocity))then
#ifdef CPP_MPI
      irecv  = 3*ntot_pT
      idispl = 3*ioff_pT
      call MPI_Allgatherv( fermivel_tmp, 3*nkpt, MPI_DOUBLE_PRECISION, &
                         & fermivelocity, irecv, idispl,               &
                         & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr  )
      if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv fermivelocity'
#else
      fermivelocity = fermivel_tmp
#endif
    end if!present(fermivelocity)


    !gather the spinvalue from the processes
    if(allocated(spinvalue))then
#ifdef CPP_MPI
      itmp   = inc%ndegen*cfg%nsqa
      irecv  = itmp*ntot_pT
      idispl = itmp*ioff_pT
      call MPI_Allgatherv( spinval_tmp, itmp*nkpt, MPI_DOUBLE_PRECISION, &
                         & spinvalue, irecv, idispl,                     &
                         & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr    )
      if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv spinvalue'
#else
      spinvalue= spinval_tmp
#endif
    end if!present(spinvalue)


    !gather the spinvec from the processes
    if(allocated(spinvec))then
#ifdef CPP_MPI
      itmp   = 3*inc%ndegen*cfg%nsqa
      irecv  = itmp*ntot_pT
      idispl = itmp*ioff_pT
      call MPI_Allgatherv( spinvec_tmp, itmp*nkpt, MPI_DOUBLE_PRECISION, &
                         & spinvec, irecv, idispl,                       &
                         & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr    )
      if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv spinvec'
#else
      spinvec= spinvec_tmp
#endif
    end if!present(spinvec)


    !gather the spinvec_atom from the processes
    if(allocated(spinvec_atom))then
#ifdef CPP_MPI
      itmp   = 3*inc%ndegen*inc%natypd
      irecv  = itmp*ntot_pT
      idispl = itmp*ioff_pT
      call MPI_Allgatherv( spinvec_atom_tmp, itmp*nkpt, MPI_DOUBLE_PRECISION, &
                         & spinvec_atom, irecv, idispl,                       &
                         & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr         )
      if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv spinvec_atom'
#else
      spinvec_atom= spinvec_atom_tmp
#endif
    end if!present(spinvec_atom)


    !gather the torquevalue from the processes
    if(allocated(torqvalue))then
#ifdef CPP_MPI
      itmp   = 3*inc%ndegen
      irecv  = itmp*ntot_pT
      idispl = itmp*ioff_pT
      call MPI_Allgatherv( torqval_tmp, itmp*nkpt, MPI_DOUBLE_PRECISION, &
                         & torqvalue, irecv, idispl,                     &
                         & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr    )
      if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv torqvalue'
#else
      torqvalue= torqval_tmp
#endif
    end if!present(torqvalue)

    !gather the torquevalue_atom from the processes
    if(allocated(torqvalue_atom))then
#ifdef CPP_MPI
      itmp   = 3*inc%ndegen*inc%natypd
      irecv  = itmp*ntot_pT
      idispl = itmp*ioff_pT
      call MPI_Allgatherv( torqval_atom_tmp, itmp*nkpt, MPI_DOUBLE_PRECISION, &
                         & torqvalue_atom, irecv, idispl,                     &
                         & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr    )
      if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv torqvalue_atom'
#else
      torqvalue_atom= torqval_atom_tmp
#endif
    end if!present(torqvalue_atom)

    !gather the spinflux_atom from the processes
    if(allocated(spinflux_atom))then
#ifdef CPP_MPI
      itmp   = 3*inc%ndegen*inc%natypd
      irecv  = itmp*ntot_pT
      idispl = itmp*ioff_pT
      call MPI_Allgatherv( spinflux_atom_tmp, itmp*nkpt, MPI_DOUBLE_PRECISION, &
                         & spinflux_atom, irecv, idispl,                     &
                         & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr    )
      if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv spinflux_atom'
#else
      spinflux_atom= spinflux_atom_tmp
#endif
    end if!present(spinflux_atom)

    !gather the alphavalue from the processes
    if(allocated(alphavalue))then
#ifdef CPP_MPI
      itmp   = 3*inc%ndegen
      irecv  = itmp*ntot_pT
      idispl = itmp*ioff_pT
      call MPI_Allgatherv( alphaval_tmp, itmp*nkpt, MPI_DOUBLE_PRECISION, &
                         & alphavalue, irecv, idispl,                     &
                         & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr     )
      if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv alphavalue'
#else
      alphavalue= alphaval_tmp
#endif
    end if!present(alphavalue)


190 FORMAT('                 |'$)
200 FORMAT('|'$)
  end subroutine calc_on_fsurf

  subroutine calc_spinvalue_state(inc,tgmatrx,rveig_in,spin_value,eigvect_rot)

    use type_inc,       only: inc_type
    use type_data,      only: tgmatrx_type
    use mod_eigvects,   only: rewrite_eigv_atom, orthgonalize_wavefunc, normalize_wavefunc, calc_proj_wavefunc, calc_norm_wavefunc
    use mod_spintools,  only: spin_expectationvalue, spin_crossterm, Spinrot_AlphaBeta, Spinrot_AlphaBeta_Rashba, rotate_wavefunction
    use mod_mympi, only: myrank
    implicit none
    type(inc_type),     intent(in)  :: inc
    type(tgmatrx_type), intent(in)  :: tgmatrx
    double complex,     intent(in)  :: rveig_in(inc%almso,inc%ndegen)
    double precision,   intent(out) :: spin_value(inc%ndegen,cfg%nsqa)
    double complex,     intent(out), optional :: eigvect_rot(inc%lmmaxso,inc%natypd,inc%ndegen,cfg%nsqa)

    double complex :: rveig_atom_unrot(inc%lmmaxso,inc%natypd,inc%ndegen), &
                    & rveig_atom_rot(inc%lmmaxso,inc%natypd,inc%ndegen)
    double complex :: Spin_ini(3,2), Spin_ini_atom(3,inc%natypd,2), Scross(3), Scross_atom(3,inc%natypd), Spin_final(3,inc%ndegen), Spin_final_atom(3,inc%natypd,inc%ndegen)
    double precision :: alpha, beta, Spin_estimated
    integer :: isqa, ient

    if(.not.cfg_read) call read_cfg()

    call rewrite_eigv_atom(inc,inc%ndegen,rveig_in,rveig_atom_unrot)

    select case (inc%ndegen)
    case( 1 ); !do nothing
    case( 2 ); call orthgonalize_wavefunc(inc, tgmatrx%rhod(:,:,:,1), rveig_atom_unrot)
    case default; stop 'inc%ndegen /= (1,2) not implemented'
    end select
    call normalize_wavefunc(inc, inc%ndegen, tgmatrx%rhod(:,:,:,1), rveig_atom_unrot)


    if(inc%ndegen==1)then
      !************************************************************************
      !* in this case, the wave-functions are already uniquely defined        *
      !* (except an overall phase, which does not influence the final result) *
      !************************************************************************

      !calculate spin-expectation value (vector-form)
      spin_value = 0d0
      call spin_expectationvalue(inc, inc%ndegen, tgmatrx%rhod, rveig_atom_unrot, Spin_final, Spin_atom=Spin_final_atom)

      !make a projection onto the SQA
      do isqa=1,cfg%nsqa
        spin_value(1,isqa) = sum(cfg%nvect(:,isqa)*dble(Spin_final(:,1)))

        !save the eigenvector into the array
        if(present(eigvect_rot)) then
          eigvect_rot(:,:,:,isqa) = rveig_atom_unrot
        end if!present(eigvect_rot)

      end do!isqa

    elseif(inc%ndegen==2)then
      !******************************************************************************
      !* In this case, all states are two-fold degenerate and any linear            *
      !* combination is also an eigenstat of the Hamiltonian. Therefore, we have    *
      !* to find the (correct) linear combination first to find the wave functions. *
      !******************************************************************************
      
      !calculate spin-expectation values of the (initial) wavefunctions
      call spin_expectationvalue(inc, inc%ndegen, tgmatrx%rhod, rveig_atom_unrot, Spin_ini, Spin_ini_atom)
      call spin_crossterm(inc, tgmatrx%rhod, rveig_atom_unrot, Scross, Scross_atom)

      do isqa=1,cfg%nsqa

        ! find parameter alpha and beta to perform the linear combination
        if(cfg%lrashba==1)then
          call Spinrot_AlphaBeta_Rashba(Spin_ini_atom(:,cfg%ilayer,:), Scross_atom(:,cfg%ilayer), alpha, beta)
        else!cfg%lrashba
          call Spinrot_AlphaBeta( cfg%ispincomb(isqa), cfg%nvect(:,isqa), Spin_ini, &
                                & Scross, alpha, beta, Spin_estimated               )
        end if!cfg%lrashba

        ! perform the linear combination of the degenerate wavefunctions
        call rotate_wavefunction(inc%lmmaxso, inc%natypd, alpha, beta, rveig_atom_unrot, rveig_atom_rot)

        ! calculate the true (=final) spin-expectation value
        call spin_expectationvalue(inc, inc%ndegen, tgmatrx%rhod, rveig_atom_rot, Spin_final, Spin_final_atom)

        ! project onto direction n and save to output-array
        do ient=1,2
          spin_value(ient,isqa) = sum(cfg%nvect(:,isqa)*dble(Spin_final(:,ient)))
        end do!ient

        !save the eigenvector into the array
        if(present(eigvect_rot)) then
          eigvect_rot(:,:,:,isqa) = rveig_atom_rot
        end if!present(eigvect_rot)

      end do!isqa

    end if!inc%ndegen==2

  end subroutine calc_spinvalue_state



  subroutine calc_spinvalue_state_generalized(inc,tgmatrx,rveig_in,spin_value,spinvec,spinvec_atom,torq_value,torq_value_atom,spinflux_atom,alpha_value,eigvect_rot)

    use type_inc,       only: inc_type
    use type_data,      only: tgmatrx_type
    use mod_eigvects,   only: rewrite_eigv_atom, orthgonalize_wavefunc, normalize_wavefunc, calc_proj_wavefunc, calc_norm_wavefunc
    use mod_spintools,  only: spin_expectationvalue, spin_crossterm, Spinrot_AlphaBeta, Spinrot_AlphaBeta_Rashba, rotate_wavefunction,&
                            & torq_expectationvalue, spinflux_expectationvalue, alpha_expectationvalue
    use mod_mympi, only: myrank
    implicit none
    type(inc_type),     intent(in)  :: inc
    type(tgmatrx_type), intent(in)  :: tgmatrx
    double complex,     intent(in)  :: rveig_in(inc%almso,inc%ndegen)
    double precision,   intent(out) :: spin_value(inc%ndegen,cfg%nsqa)
    double complex,     intent(inout), optional    :: eigvect_rot(:,:,:,:)  !eigvect_rot(inc%lmmaxso,inc%natypd,inc%ndegen,cfg%nsqa)
    double precision,   allocatable, intent(inout) :: spinvec(:,:,:)        !spinvec(3,inc%ndegen,cfg%nsqa)
    double precision,   allocatable, intent(inout) :: spinvec_atom(:,:,:)   !spinvec_atom(3,inc%natypd,inc%ndegen)
    double precision,   allocatable, intent(inout) :: torq_value(:,:)       !torq_value(3,inc%ndegen)
    double precision,   allocatable, intent(inout) :: torq_value_atom(:,:,:)!torq_value_atom(3,inc%natypd,inc%ndegen)
    double precision,   allocatable, intent(inout) :: spinflux_atom(:,:,:)  !spinflux_atom(3,inc%natypd,inc%ndegen)
    double precision,   allocatable, intent(inout) :: alpha_value(:,:)      !alpha_value(3,inc%ndegen)

    double complex :: rveig_atom_unrot(inc%lmmaxso,inc%natypd,inc%ndegen), &
                    & rveig_atom_rot(inc%lmmaxso,inc%natypd,inc%ndegen)
    double complex :: Spin_ini(3,2), Spin_ini_atom(3,inc%natypd,2), Scross(3), Scross_atom(3,inc%natypd), Spin_final(3,inc%ndegen), Spin_final_atom(3,inc%natypd,inc%ndegen), Torq_final(3,inc%ndegen), Torq_atom(3,inc%natypd,inc%ndegen), Spinfluxes_atom(3,inc%natypd,inc%ndegen), Alpha_Final(3,inc%ndegen)
    double precision :: alpha, beta, Spin_estimated
    integer :: isqa, ient

    if(.not.cfg_read) call read_cfg()

    call rewrite_eigv_atom(inc,inc%ndegen,rveig_in,rveig_atom_unrot)

    select case (inc%ndegen)
    case( 1 ); !do nothing
    case( 2 ); call orthgonalize_wavefunc(inc, tgmatrx%rhod(:,:,:,1), rveig_atom_unrot)
    case default; stop 'inc%ndegen /= (1,2) not implemented'
    end select
    call normalize_wavefunc(inc, inc%ndegen, tgmatrx%rhod(:,:,:,1), rveig_atom_unrot)


    if(inc%ndegen==1)then
      !************************************************************************
      !* in this case, the wave-functions are already uniquely defined        *
      !* (except an overall phase, which does not influence the final result) *
      !************************************************************************

      !calculate spin-expectation value (vector-form)
      spin_value = 0d0
      call spin_expectationvalue(inc, inc%ndegen, tgmatrx%rhod, rveig_atom_unrot, Spin_final, Spin_atom=Spin_final_atom)


      !make a projection onto the SQA
      do isqa=1,cfg%nsqa
        spin_value(1,isqa) = sum(cfg%nvect(:,isqa)*dble(Spin_final(:,1)))

        !save the eigenvector into the array
        if(present(eigvect_rot)) then
          eigvect_rot(:,:,:,isqa) = rveig_atom_unrot
        end if!present(eigvect_rot)

        if(allocated(spinvec)) spinvec(:,:,isqa)=dble(Spin_final)
        if(allocated(spinvec_atom) .and. isqa==1) spinvec_atom(:,:,:)=dble(Spin_final_atom)
      end do!isqa

      !calculate torque expectation value
      if(allocated(torq_value)) then
        call torq_expectationvalue(inc, inc%ndegen, tgmatrx%torq, rveig_atom_unrot, Torq_final, Torq_atom=Torq_atom)
        torq_value(:,:)=dble(Torq_final(:,:))
        if(allocated(torq_value_atom)) torq_value_atom(:,:,:)=dble(Torq_atom(:,:,:))
      end if!present(torq_value)

      !calculate spinflux expectation value
      if(allocated(spinflux_atom)) then
        call spinflux_expectationvalue(inc, inc%ndegen, tgmatrx%spinflux, rveig_atom_unrot, Spinflux_atom=Spinfluxes_atom)
        spinflux_atom(:,:,:)=dble(Spinfluxes_atom(:,:,:))
      end if!present(spinflux_atom)

      !calculate alpha expectation value
      if(allocated(alpha_value)) then
        call alpha_expectationvalue(inc, inc%ndegen, tgmatrx%alpha, rveig_atom_unrot, Alpha_final)
        alpha_value(:,:)=dble(Alpha_final(:,:))
      end if!present(alpha_value)

    elseif(inc%ndegen==2)then
      !******************************************************************************
      !* In this case, all states are two-fold degenerate and any linear            *
      !* combination is also an eigenstat of the Hamiltonian. Therefore, we have    *
      !* to find the (correct) linear combination first to find the wave functions. *
      !******************************************************************************
      
      !calculate spin-expectation values of the (initial) wavefunctions
      call spin_expectationvalue(inc, inc%ndegen, tgmatrx%rhod, rveig_atom_unrot, Spin_ini, Spin_ini_atom)
      call spin_crossterm(inc, tgmatrx%rhod, rveig_atom_unrot, Scross, Scross_atom)

      do isqa=1,cfg%nsqa

        ! find parameter alpha and beta to perform the linear combination
        if(cfg%lrashba==1)then
          call Spinrot_AlphaBeta_Rashba(Spin_ini_atom(:,cfg%ilayer,:), Scross_atom(:,cfg%ilayer), alpha, beta)
        else!cfg%lrashba
          call Spinrot_AlphaBeta( cfg%ispincomb(isqa), cfg%nvect(:,isqa), Spin_ini, &
                                & Scross, alpha, beta, Spin_estimated               )
        end if!cfg%lrashba

        ! perform the linear combination of the degenerate wavefunctions
        call rotate_wavefunction(inc%lmmaxso, inc%natypd, alpha, beta, rveig_atom_unrot, rveig_atom_rot)

        ! calculate the true (=final) spin-expectation value
        call spin_expectationvalue(inc, inc%ndegen, tgmatrx%rhod, rveig_atom_rot, Spin_final, Spin_final_atom)

        ! project onto direction n and save to output-array
        do ient=1,2
          spin_value(ient,isqa) = sum(cfg%nvect(:,isqa)*dble(Spin_final(:,ient)))
        end do!ient

        !save the result to the putput-arrays
        if(allocated(spinvec)) spinvec(:,:,isqa)=dble(Spin_final)
        if(allocated(spinvec_atom) .and. isqa==1) spinvec_atom(:,:,:)=dble(Spin_final_atom)

        !save the eigenvector into the array
        if(present(eigvect_rot)) then
          eigvect_rot(:,:,:,isqa) = rveig_atom_rot
        end if!present(eigvect_rot)

        ! calculate the torque expectation value
        ! ATTENTION!!!!! This is not properly implemented. If it made sense to calculate the torque
        !                also when there is the two-fold conjugation degeneracy (AFMs???), this implementation
        !                would only store the torque of the last SQA. Therefore, in "read_cfg", we check 
        !                that cfg%nsqa=1 if ltorq=1.
        if(allocated(torq_value)) then
          call torq_expectationvalue(inc, inc%ndegen, tgmatrx%torq, rveig_atom_unrot, Torq_final, Torq_atom)
          torq_value(:,:)=dble(Torq_final(:,:))
          if(allocated(torq_value_atom)) torq_value_atom(:,:,:)=dble(Torq_atom(:,:,:))
        end if!present(torq_value)

        if(allocated(torq_value_atom)) then
          call spinflux_expectationvalue(inc, inc%ndegen, tgmatrx%spinflux, rveig_atom_unrot, Spinfluxes_atom)
          spinflux_atom(:,:,:)=dble(Spinfluxes_atom(:,:,:))
        end if!present(torq_value_atom)

        if(allocated(alpha_value)) then
          call alpha_expectationvalue(inc, inc%ndegen, tgmatrx%spinflux, rveig_atom_unrot, Alpha_final)
          alpha_value(:,:)=dble(Alpha_final(:,:))
        end if!present(alpha_value)

      end do!isqa

    end if!inc%ndegen==2

  end subroutine calc_spinvalue_state_generalized




  subroutine calc_fermivel_new(inc,lattice,cluster,tgmatrx,kpoint,eigw_in,LVin,RVin,fermi_velocity)

    use type_inc,       only: inc_type
    use type_data,      only: lattice_type, cluster_type, tgmatrx_type
    use mod_mympi,      only: myrank, nranks, master
    use mod_kkrmat,     only: compute_kkr_eigenvectors2, compute_kkr_eigenvectors2_dk
    use mod_eigvects,   only: compare_eigv
    use mod_mathtools,  only: bubblesort
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    type(tgmatrx_type), intent(in) :: tgmatrx

    double precision, intent(in) :: kpoint(3)
    double complex,   intent(in) :: eigw_in, LVin(inc%almso), RVin(inc%almso)
    
    double precision, intent(out) :: fermi_velocity(3)

    double precision :: Delta_k(3), kpoint_dk(3), fac_fv
    double complex   :: eigwPT(inc%almso), LVeigPT(inc%almso,inc%almso), RVeigPT(inc%almso,inc%almso)

    double complex   :: Delta_lambda_kxyz(3), Delta_lambda_E(2)

    integer :: ikxyz, ipm, LMout

    logical, save :: first=.true.

!   integer          :: sortind(inc%almso)
!   double precision :: dtmp(inc%almso)

!   integer, parameter :: dimenk=3

    !==================================================!
    !+++         FERMI-VELOCITY-CALULATION          +++!
    !==================================================!
    !  a) perturb the eigenvalue w.r.t. the k-point    !
    !  b) perturb the eigenvalue w.r.t. the energy     !
    !  c) determine the fermi-velocity by              !
    !                     \nabla_k lambda              !
    !        v_F = --------------------------------    !
    !                (\delta lambda) / (\delta E)      !
    !==================================================!


    !++++++++++++++++!
    !+++ BEGIN a) +++!
    !++++++++++++++++!

    call compute_kkr_eigenvectors2_dk( inc, lattice, cluster, tgmatrx%ginp(:,:,:,1),&
                                     & tgmatrx%tmat(:,:,1), kpoint_dk, LVin, RVin,  &
                                     & Delta_lambda_kxyz                            )

    write(6000,'(6ES25.16)') Delta_lambda_kxyz

    !++++++++++++++++!
    !+++ BEGIN b) +++!
    !++++++++++++++++!
    do ipm=1,2

      !calculate the eigenvalue at the perturbed k-point
      call compute_kkr_eigenvectors2( inc, lattice, cluster, tgmatrx%ginp(:,:,:,1+ipm),&
                                    & tgmatrx%tmat(:,:,1+ipm), kpoint,                 &
                                    & eigwPT, LVeigPT, RVeigPT                         )

!     dtmp = abs(eigwPT)
!     call bubblesort(inc%almso, dtmp, sortind)
!     write(*,'(4X,"sorted abs(eigw)=", 6ES18.9)') abs(eigwPT(sortind(1:6)))

      call compare_eigv(inc, LVin, RVeigPT, LMout)

!     write(*,'(2X,"pert e= (",2ES18.9,"), lmo= ",I0)') eigwPT(LMout), LMout

      Delta_lambda_E(ipm) = eigwPT(LMout) - eigw_in
      
    end do!ipm
 

    !++++++++++++++++!
    !+++ BEGIN c) +++!
    !++++++++++++++++!
    fac_fv = 2*dble(tgmatrx%energies(2) - tgmatrx%energies(1))
    fermi_velocity(:) = -dble(   Delta_lambda_kxyz(:) / (Delta_lambda_E(1)-Delta_lambda_E(2)) ) *fac_fv
    write(7000,'(6ES25.16)') Delta_lambda_kxyz(:) / (Delta_lambda_E(1)-Delta_lambda_E(2))*fac_fv

  end subroutine calc_fermivel_new





  subroutine calc_fermivel(inc,lattice,cluster,tgmatrx,kpoint,eigw_in,LVin,RVin,fermi_velocity)

    use type_inc,       only: inc_type
    use type_data,      only: lattice_type, cluster_type, tgmatrx_type
    use mod_mympi,      only: myrank, nranks, master
    use mod_kkrmat,     only: compute_kkr_eigenvectors2
    use mod_eigvects,   only: compare_eigv
    use mod_mathtools,  only: bubblesort
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    type(tgmatrx_type), intent(in) :: tgmatrx

    double precision, intent(in) :: kpoint(3)
    double complex,   intent(in) :: eigw_in, LVin(inc%almso), RVin(inc%almso)
    
    double precision, intent(out) :: fermi_velocity(3)

    double precision :: Delta_k(3), kpoint_dk(3), fac_fv
    double complex   :: eigwPT(inc%almso), LVeigPT(inc%almso,inc%almso), RVeigPT(inc%almso,inc%almso)

    double complex   :: Delta_lambda_kxyz(3,2), Delta_lambda_E(2)

    integer :: ikxyz, ipm, LMout
!   integer          :: sortind(inc%almso)
!   double precision :: dtmp(inc%almso)

!   integer, parameter :: dimenk=3

    !==================================================!
    !+++         FERMI-VELOCITY-CALULATION          +++!
    !==================================================!
    !  a) perturb the eigenvalue w.r.t. the k-point    !
    !  b) perturb the eigenvalue w.r.t. the energy     !
    !  c) determine the fermi-velocity by              !
    !                     \nabla_k lambda              !
    !        v_F = --------------------------------    !
    !                (\delta lambda) / (\delta E)      !
    !==================================================!


    !++++++++++++++++!
    !+++ BEGIN a) +++!
    !++++++++++++++++!
    Delta_lambda_kxyz = (0d0,0d0)
    do ikxyz=1,inc%nBZdim!dimenk
      do ipm=1,2
        !determine the perturbed k-point
        Delta_k = 0d0
        Delta_k(ikxyz) = dble(3-2*ipm)*cfg%dk_fv ! = +-Delta_k_pert depending on ipm=i_plus_minus
        kpoint_dk = kpoint + Delta_k

        !calculate the eigenvalue at the perturbed k-point
        call compute_kkr_eigenvectors2( inc, lattice, cluster, tgmatrx%ginp(:,:,:,1),&
                                      & tgmatrx%tmat(:,:,1), kpoint_dk,              &
                                      & eigwPT, LVeigPT, RVeigPT                     )

        call compare_eigv(inc, LVin, RVeigPT, LMout)

!       write(*,'(2X,"pert k= (",2ES18.9,"), lmo= ",I0)') eigwPT(LMout), LMout

        Delta_lambda_kxyz(ikxyz,ipm) = eigwPT(LMout) - eigw_in

      end do!ipm
    end do!ikxyz


!   write(6001,'(6ES25.16)') Delta_lambda_kxyz(:,1)-Delta_lambda_kxyz(:,2)/(2*cfg%dk_fv)

    !++++++++++++++++!
    !+++ BEGIN b) +++!
    !++++++++++++++++!
    do ipm=1,2

      !calculate the eigenvalue at the perturbed k-point
      call compute_kkr_eigenvectors2( inc, lattice, cluster, tgmatrx%ginp(:,:,:,1+ipm),&
                                    & tgmatrx%tmat(:,:,1+ipm), kpoint,                 &
                                    & eigwPT, LVeigPT, RVeigPT                         )

!     dtmp = abs(eigwPT)
!     call bubblesort(inc%almso, dtmp, sortind)
!     write(*,'(4X,"sorted abs(eigw)=", 6ES18.9)') abs(eigwPT(sortind(1:6)))

      call compare_eigv(inc, LVin, RVeigPT, LMout)

!     write(*,'(2X,"pert e= (",2ES18.9,"), lmo= ",I0)') eigwPT(LMout), LMout

      Delta_lambda_E(ipm) = eigwPT(LMout) - eigw_in
      
    end do!ipm
 

    !++++++++++++++++!
    !+++ BEGIN c) +++!
    !++++++++++++++++!
    fac_fv = dble(tgmatrx%energies(2) - tgmatrx%energies(1))/cfg%dk_fv
    fermi_velocity(:) = -dble( (Delta_lambda_kxyz(:,1) - Delta_lambda_kxyz(:,2))/(Delta_lambda_E(1)-Delta_lambda_E(2)))*fac_fv
!   write(7001,'(6ES25.16)') -(Delta_lambda_kxyz(:,1) - Delta_lambda_kxyz(:,2))/(Delta_lambda_E(1)-Delta_lambda_E(2))*fac_fv

  end subroutine calc_fermivel





  subroutine read_cfg(force_spinread)

    use mod_ioinput, only: IoInput
    use mod_mympi,   only: myrank, nranks, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    logical, intent(in), optional :: force_spinread

    integer :: nspincomb, angrep, idir, inum, ierr
    double precision :: dtmpin(3), theta, phi
    character(len=80) :: uio
    double precision, allocatable :: vectintmp(:,:)
    logical :: force_spin

    if(myrank==master) then
      force_spin = .false.
      if(present(force_spinread)) force_spin=force_spinread

      call IoInput('LFVEL     ',uio,1,7,ierr)
      read(unit=uio,fmt=*) cfg%lfvel
      if(cfg%lfvel==1)then
        call IoInput('DKFVEL    ',uio,1,7,ierr)
        read(unit=uio,fmt=*) cfg%dk_fv
      end if!lfvel==1



      call IoInput('LSPIN     ',uio,1,7,ierr)
      read(unit=uio,fmt=*) cfg%lspin
      if(cfg%lspin==1 .or. force_spin) then

        call IoInput('LSPINATOM ',uio,1,7,ierr)
        read(unit=uio,fmt=*) cfg%lspinperatom

        call IoInput('LTORQ     ',uio,1,7,ierr)
        read(unit=uio,fmt=*) cfg%ltorq

        call IoInput('LTORQATOM ',uio,1,7,ierr)
        read(unit=uio,fmt=*) cfg%ltorqperatom

        call IoInput('LSPINFLUX ',uio,1,7,ierr)
        if(ierr==0) then
          read(unit=uio,fmt=*) cfg%lspinflux
        else ! ensures compatibility of old format input files
          write(*,*) "Warning : LSPINFLUX set to 0 by default !"
          cfg%lspinflux=0
        end if!ierr==0

        call IoInput('LALPHA    ',uio,1,7,ierr)
        if(ierr==0) then
          read(unit=uio,fmt=*) cfg%lalpha
        else ! ensures compatibility of old format input files
          write(*,*) "Warning : LALPHA set to 0 by default !"
          cfg%lalpha=0
        end if!ierr==0

        call IoInput('SIMPSON   ',uio,1,7,ierr)
        read(unit=uio,fmt=*) cfg%simpson

        call IoInput('LSAVEEIGV ',uio,1,7,ierr)
        read(unit=uio,fmt=*) inum
        if(inum==1)then
           cfg%saveeigv=.true.
        else
           cfg%saveeigv=.false.
        end if

        call IoInput('LRASHBA   ',uio,1,7,ierr)
        read(unit=uio,fmt=*) cfg%lrashba
        if(cfg%lrashba==1)then
          call IoInput('ILAYER    ',uio,1,7,ierr)
          read(unit=uio,fmt=*) cfg%ilayer
        end if!lrashba==1

        call IoInput('NSQA      ',uio,1,7,ierr)
        read(unit=uio,fmt=*) cfg%nsqa
        write(6,*)"cfg%nsqa=",cfg%nsqa
        allocate(vectintmp(3,cfg%nsqa), STAT=ierr)
        if(ierr/=0) stop 'Problem allocating vectintmp'

        !check consistency of NSQA and LRASHBA
        if(cfg%lrashba==1 .and. cfg%nsqa>1)then
          write(*,'(A)') 'WARNING! setting NSQA=1 for option LRASHBA'
          cfg%nsqa=1
        end if

        !check consistency of NSQA and LSPINATOM
        if(cfg%lspinperatom==1 .and. cfg%nsqa>1)then
          write(*,'(A)') 'WARNING! setting NSQA=1 for option LSPINATOM'
          cfg%nsqa=1
        end if

        call IoInput('SPINCOMB  ',uio,1,7,ierr)
        read(unit=uio,fmt=*) nspincomb

        call IoInput('ANGREP    ',uio,1,7,ierr)
        read(unit=uio,fmt=*) angrep

        !read in vectors
        do idir=1,cfg%nsqa
          if(angrep==0) then
            call IoInput('ANISOVECS ',uio,1+idir,7,ierr)
            read (unit=uio,fmt=*) inum, dtmpin(1:2)
          elseif(angrep==1) then
            call IoInput('ANISOVECS ',uio,1+idir,7,ierr)
            read (unit=uio,fmt=*) inum, dtmpin(1:3)
          end if
          call convert_nvect(angrep,dtmpin,vectintmp(:,idir))
        end do!idir

        !check consistency of NSQA and LTORQ
        if(cfg%ltorq==1 .and. cfg%nsqa>1)then
          write(*,'(A)') 'WARNING! setting NSQA=1 for option LTORQ'
          cfg%nsqa=1
          if(nspincomb==0)then
            write(*,'(A)') 'WARNING! option SPINCOMB=0 not allaowed for LTORQ. Assunming SPINCOMB=1'
            nspincomb=1
          end if!nspincomb==0
        end if!cfg%ltorq==1 .and. cfg%nsqa>1

        !combine the spin-sqa-data
        if( nspincomb==0 )then
          !case that both conditions shall be treated

          allocate(cfg%ispincomb(2*cfg%nsqa), cfg%nvect(3,2*cfg%nsqa), STAT=ierr)
          if(ierr/=0) stop 'Problem allocating cfg%ispincomb etc. on master'

          do idir=1,cfg%nsqa
            cfg%nvect(:,2*idir-1) = vectintmp(:,idir)
            cfg%nvect(:,2*idir  ) = vectintmp(:,idir)
            cfg%ispincomb(2*idir-1) = 1
            cfg%ispincomb(2*idir  ) = 2
          end do!idir

          cfg%nsqa = 2*cfg%nsqa

        else!nspincomb
          !case that only one condition shall be treated

          allocate(cfg%ispincomb(cfg%nsqa), cfg%nvect(3,cfg%nsqa), STAT=ierr)
          if(ierr/=0) stop 'Problem allocating cfg%ispincomb etc. on master'

          cfg%nvect = vectintmp
          cfg%ispincomb = nspincomb

        end if!nspincomb

      end if!lspin==1


      if(cfg%lfvel==1 .or. cfg%lspin==1)then
        call IoInput('ONFSMODE  ',uio,1,7,ierr)
        read(unit=uio,fmt=*) cfg%mode
        call IoInput('ROTMODE   ',uio,1,7,ierr)
        read(unit=uio,fmt=*) cfg%rotatemode
      end if

    end if!myrank==master

#ifdef CPP_MPI
    call myBcast_cfg(cfg)
#endif

    cfg_read = .true.

  end subroutine read_cfg





  subroutine convert_nvect(angrep,dtmpin,nvect)
    use mod_mathtools, only: pi
    implicit none
    integer, intent(in) :: angrep
    double precision, intent(in) :: dtmpin(3)
    double precision, intent(out) :: nvect(3)

    double precision :: theta, phi, norm

    if(angrep==0) then
      theta = dtmpin(1)
      phi   = dtmpin(2)
      nvect(1) = dsin(theta*pi)*dcos(phi*pi)
      nvect(2) = dsin(theta*pi)*dsin(phi*pi)
      nvect(3) = dcos(theta*pi)
    elseif(angrep==1) then
      norm = sqrt(sum(dtmpin**2))
      nvect = dtmpin/norm
    end if
    
  end subroutine convert_nvect





#ifdef CPP_MPI
  subroutine myBcast_cfg(cfg)

    use mpi
    use mod_mympi,   only: myrank, nranks, master
    implicit none

    type(cfg_type), intent(inout)  :: cfg

    integer :: blocklen1(cfg%N1), etype1(cfg%N1), myMPItype1
    integer :: blocklen2(cfg%N2), etype2(cfg%N2), myMPItype2
    integer :: ierr
    integer(kind=MPI_ADDRESS_KIND) :: disp1(cfg%N1), disp2(cfg%N2), base

    call MPI_Get_address(cfg%N1,           disp1(1), ierr)
    call MPI_Get_address(cfg%lspin,        disp1(2), ierr)
    call MPI_Get_address(cfg%lfvel,        disp1(3), ierr)
    call MPI_Get_address(cfg%lrashba,      disp1(4), ierr)
    call MPI_Get_address(cfg%lspinperatom, disp1(5), ierr)
    call MPI_Get_address(cfg%ltorqperatom, disp1(6), ierr)
    call MPI_Get_address(cfg%ltorq,        disp1(7), ierr)
    call MPI_Get_address(cfg%lspinflux,    disp1(8), ierr)
    call MPI_Get_address(cfg%lalpha,       disp1(9), ierr)
    call MPI_Get_address(cfg%nsqa,         disp1(10), ierr)
    call MPI_Get_address(cfg%mode,         disp1(11), ierr)
    call MPI_Get_address(cfg%rotatemode,   disp1(12), ierr)
    call MPI_Get_address(cfg%ilayer,       disp1(13), ierr)
    call MPI_Get_address(cfg%dk_fv,        disp1(14), ierr)
    call MPI_Get_address(cfg%saveeigv,     disp1(15), ierr)
    call MPI_Get_address(cfg%simpson,      disp1(16), ierr)

    base  = disp1(1)
    disp1 = disp1 - base

    blocklen1(1:16)=1

    etype1(1:13) = MPI_INTEGER
    etype1(14)   = MPI_DOUBLE_PRECISION
    etype1(15:16)= MPI_LOGICAL

    call MPI_Type_create_struct(cfg%N1, blocklen1, disp1, etype1, myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_cfg_1'

    call MPI_Type_commit(myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error commiting create_mpimask_cfg_1'

    call MPI_Bcast(cfg%N1, 1, myMPItype1, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting cfg_1'

    call MPI_Type_free(myMPItype1, ierr)

    if(cfg%lspin==1)then

      if(myrank/=master) then
        allocate( cfg%ispincomb(cfg%nsqa), cfg%nvect(3,cfg%nsqa), STAT=ierr )
        if(ierr/=0) stop 'Problem allocating cfg_arrays on slaves'
      end if

      call MPI_Get_address(cfg%N2,        disp2(1), ierr)
      call MPI_Get_address(cfg%ispincomb, disp2(2), ierr)
      call MPI_Get_address(cfg%nvect,     disp2(3), ierr)

      base  = disp2(1)
      disp2 = disp2 - base

      blocklen2(1)=1
      blocklen2(2)=size(cfg%ispincomb)
      blocklen2(3)=size(cfg%nvect)

      etype2(1:2) = MPI_INTEGER
      etype2(3)   = MPI_DOUBLE_PRECISION

      call MPI_Type_create_struct(cfg%N2, blocklen2, disp2, etype2, myMPItype2, ierr)
      if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_cfg_2'

      call MPI_Type_commit(myMPItype2, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error commiting create_mpimask_cfg_2'

      call MPI_Bcast(cfg%N2, 1, myMPItype2, master, MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error brodcasting cfg_2'

      call MPI_Type_free(myMPItype2, ierr)

    end if!cfg%lspin==1

  end subroutine myBcast_cfg
#endif





  subroutine save_fermivelocity(filemode, nkpts, fermivel, nsym, isym)

    use mod_ioformat,   only: fmt_fn_sub_ext, ext_formatted, filename_fvel
    implicit none

    character(len=*), intent(in) :: filemode
    integer,          intent(in) :: nkpts, nsym, isym(nsym)
    double precision, intent(in) :: fermivel(3,nkpts)

    integer :: ii
    character(len=256) :: filename

    write(filename,fmt_fn_sub_ext) filename_fvel, filemode, ext_formatted
    open(unit=326528, file=trim(filename), form='formatted', action='write')
    write(326528,'(2I8)') nkpts, nsym
    write(326528,'(12I8)') isym
    write(326528,'(3ES25.16)') fermivel
    close(326528)

  end subroutine save_fermivelocity





  subroutine save_spinvalue(filemode, nkpts, nsqa, ndegen, ispincomb, nvect, spinval, nsym, isym)

    use mod_ioformat,   only: fmt_fn_sub_ext, ext_formatted, filename_spin
    implicit none

    character(len=*), intent(in) :: filemode
    integer,          intent(in) :: nkpts, nsqa, ndegen, nsym, isym(nsym), ispincomb(nsqa)
    double precision, intent(in) :: nvect(3,nsqa), spinval(ndegen,nsqa,nkpts)


    integer :: ii
    character(len=256) :: filename

    write(filename,fmt_fn_sub_ext) filename_spin, filemode, ext_formatted
    open(unit=326529, file=trim(filename), form='formatted', action='write')
    write(326529,'(2I8)') nkpts, nsym, nsqa, ndegen
    write(326529,'(12I8)') isym
    write(326529,'(10I8)') ispincomb
    write(326529,'(3ES25.16)') nvect
    write(326529,'(ES25.16)') spinval
    close(326529)

  end subroutine save_spinvalue



  subroutine save_spinvec(filemode, nkpts, nsqa, ndegen, ispincomb, nvect, spinvec, nsym, isym)

    use mod_ioformat,   only: fmt_fn_sub_ext, ext_formatted, filename_spinvec
    implicit none

    character(len=*), intent(in) :: filemode
    integer,          intent(in) :: nkpts, nsqa, ndegen, nsym, isym(nsym), ispincomb(nsqa)
    double precision, intent(in) :: nvect(3,nsqa), spinvec(3,ndegen,nsqa,nkpts)

    integer :: ii
    character(len=256) :: filename

    write(filename,fmt_fn_sub_ext) filename_spinvec, filemode, ext_formatted
    open(unit=326529, file=trim(filename), form='formatted', action='write')
    write(326529,'(2I8)') nkpts, nsym, nsqa, ndegen
    write(326529,'(12I8)') isym
    write(326529,'(10I8)') ispincomb
    write(326529,'(3ES25.16)') nvect
    write(326529,'(3ES25.16)') spinvec
    close(326529)

  end subroutine save_spinvec



  subroutine save_spinvec_atom(filemode, nkpts, ndegen, natyp, spinvec_atom, nsym, isym)

    use mod_ioformat,   only: fmt_fn_atom_sub_ext, ext_formatted, filename_spinvec
    implicit none

    character(len=*), intent(in) :: filemode
    integer,          intent(in) :: nkpts, ndegen, nsym, isym(nsym), natyp
    double precision, intent(in) :: spinvec_atom(3,natyp,ndegen,nkpts)

    integer :: ii
    character(len=256) :: filename

    do ii=1,natyp
      write(filename,fmt_fn_atom_sub_ext) filename_spinvec, ii, filemode, ext_formatted
      open(unit=326529, file=trim(filename), form='formatted', action='write')
      write(326529,'(3I8)') nkpts, nsym, ndegen
      write(326529,'(12I8)') isym
      write(326529,'(3ES25.16)') spinvec_atom(:,ii,:,:)
      close(326529)
   end do!ii=1,natyp

  end subroutine save_spinvec_atom



  subroutine save_torqvalue(filemode, nkpts, ndegen, torqval, nsym, isym)

    use mod_ioformat,   only: fmt_fn_sub_ext, ext_formatted, filename_torq
    implicit none

    character(len=*), intent(in) :: filemode
    integer,          intent(in) :: nkpts, ndegen, nsym, isym(nsym)
    double precision, intent(in) :: torqval(3,ndegen,nkpts)

    integer :: ii
    character(len=256) :: filename

    write(filename,fmt_fn_sub_ext) filename_torq, filemode, ext_formatted
    open(unit=326529, file=trim(filename), form='formatted', action='write')
    write(326529,'(3I8)') nkpts, nsym, ndegen
    write(326529,'(12I8)') isym
    write(326529,'(3ES25.16)') torqval
    close(326529)

  end subroutine save_torqvalue


  subroutine save_torqvalue_atom(filemode, nkpts, ndegen, natyp, torqval_atom, nsym, isym)

    use mod_ioformat,   only: fmt_fn_atom_sub_ext, ext_formatted, filename_torq
    implicit none

    character(len=*), intent(in) :: filemode
    integer,          intent(in) :: nkpts, ndegen, nsym, isym(nsym), natyp
    double precision, intent(in) :: torqval_atom(3,natyp,ndegen,nkpts)

    integer :: ii
    character(len=256) :: filename

    do ii=1,natyp
      write(filename,fmt_fn_atom_sub_ext) filename_torq, ii, filemode, ext_formatted
      open(unit=326529, file=trim(filename), form='formatted', action='write')
      write(326529,'(3I8)') nkpts, nsym, ndegen
      write(326529,'(12I8)') isym
      write(326529,'(3ES25.16)') torqval_atom(:,ii,:,:)
      close(326529)
    end do!ii=1,natyp

  end subroutine save_torqvalue_atom

  subroutine save_spinflux_atom(filemode, nkpts, ndegen, natyp, spinflux_atom, nsym, isym)

    use mod_ioformat,   only: fmt_fn_atom_sub_ext, ext_formatted, filename_spinflux
    implicit none

    character(len=*), intent(in) :: filemode
    integer,          intent(in) :: nkpts, ndegen, nsym, isym(nsym), natyp
    double precision, intent(in) :: spinflux_atom(3,natyp,ndegen,nkpts)

    integer :: ii
    character(len=256) :: filename

    do ii=1,natyp
      write(filename,fmt_fn_atom_sub_ext) filename_spinflux, ii, filemode, ext_formatted
      open(unit=326529, file=trim(filename), form='formatted', action='write')
      write(326529,'(3I8)') nkpts, nsym, ndegen
      write(326529,'(12I8)') isym
      write(326529,'(3ES25.16)') spinflux_atom(:,ii,:,:)
      close(326529)
    end do!ii=1,natyp

  end subroutine save_spinflux_atom

  subroutine save_alpha(filemode, nkpts, ndegen, alpha, nsym, isym)

    use mod_ioformat,   only: fmt_fn_sub_ext, ext_formatted, filename_alpha
    implicit none

    character(len=*), intent(in) :: filemode
    integer,          intent(in) :: nkpts, ndegen, nsym, isym(nsym)
    double precision, intent(in) :: alpha(3,ndegen,nkpts)

    integer :: ii
    character(len=256) :: filename

    write(filename,fmt_fn_sub_ext) filename_alpha, filemode, ext_formatted
    open(unit=326529, file=trim(filename), form='formatted', action='write')
    write(326529,'(A)') '# To adapt this file to the fermivelocity.vis.txt format:'
    write(326529,'(A)') '# if ndegen==1 : just remove the 1 on the second line'
    write(326529,'(A)') '# if ndegen==2 : remove the 2 on the second line'
    write(326529,'(A)') '#                and remove second half of the data'
    write(326529,'(2I8)') nkpts, nsym, ndegen
    write(326529,'(12I8)') isym
    do ii=1,ndegen
      write(326529,'(3ES25.16)') alpha(:,ii,:)
    end do
    close(326529)

  end subroutine save_alpha

  subroutine calculate_response_functions_CRTA_int(nsym,isym,rotmat,alat,BZVol,ndeg,natyp,nkpts,areas,fermivel,torqval,torqval_atom,spinvec_atom,spinflux)
    use mpi
    use mod_mympi,      only: myrank, master
    use mod_mathtools,  only: tpi
    use mod_symmetries, only: rotate_kpoints, expand_areas
    implicit none

    integer,          intent(in)  :: nsym, ndeg, natyp, nkpts, isym(nsym)
    double precision, intent(in)  :: alat, BZVol
    double precision, intent(in)  :: rotmat(64,3,3), areas(nkpts), fermivel(3,nkpts)
    double precision, allocatable, intent(in)  :: torqval(:,:,:)      !torqval(3,ndeg,nkpts)
    double precision, allocatable, intent(in)  :: torqval_atom(:,:,:,:) !torqval_atom(3,natypd,ndeg,nkpts)
    double precision, allocatable, intent(in)  :: spinvec_atom(:,:,:,:) !spinvec_atom(3,natypd,ndeg,nkpts)
    double precision, allocatable, intent(in)  :: spinflux(:,:,:,:)   !spinflux(3,natypd,ndegen,nkpts)

    double precision, allocatable :: fermivel_r(:,:), areas_r(:), arrtmp1(:,:), arrtmp2(:,:)
    double precision, allocatable :: torqval_r(:,:,:), torqval_atom_r(:,:,:,:), spinvec_atom_r(:,:,:,:), spinflux_r(:,:,:,:)

    integer :: ierr, ikp, i1, i2, ideg, itmp, nkpts_tmp, nkpts_r
    double precision :: torkance(3,3), torkance_atom(3,3,natyp), damping(3,3), spinacc_atom(3,3,natyp), spinflux_resp(3,3,natyp)
    double precision :: integ_torkance(3,3), integ_torkance_atom(3,3,natyp), integ_damping(3,3), integ_spinacc_atom(3,3,natyp), integ_spinflux(3,3,natyp)
    double precision :: vabs
    double precision, parameter :: RyToinvfs = 20.67068667282055d0
    integer, parameter :: iounit=13517

    if(nsym>1)then
      !apply symmetries to fermivel and areas
      call rotate_kpoints(rotmat, nkpts, fermivel, nsym, isym, nkpts_r, fermivel_r)
      call expand_areas(nsym,nkpts,areas,areas_r)

      if(allocated(torqval))then
        !next apply symmetries to torqval; as a first step, flatten the array
        itmp = size(torqval)/3
        allocate(arrtmp1(3,itmp), STAT=ierr)
        if(ierr/=0) stop 'problem alloaction arrtmp1'
        arrtmp1 = reshape(torqval, (/3, itmp/))

        !now apply the symmetries
        call rotate_kpoints(rotmat, nkpts*ndeg, arrtmp1, nsym, isym, nkpts_tmp, arrtmp2)
        if(nkpts_r*ndeg /= nkpts_tmp) stop 'nkpts_r*ndeg /= nkpts_tmp'

        !transform back to old shape
        allocate(torqval_r(3,ndeg,nkpts_r), STAT=ierr)
        if(ierr/=0) stop 'problem allocating torqval_r'
        torqval_r = reshape(arrtmp2,(/3,ndeg,nkpts_r/))

        deallocate(arrtmp1, arrtmp2)
      end if!allocated(torqval)

      if(allocated(torqval_atom))then
        !next apply symmetries to torqval_atom; as a first step, flatten the array
        itmp = size(torqval_atom)/3
        allocate(arrtmp1(3,itmp), STAT=ierr)
        if(ierr/=0) stop 'problem alloaction arrtmp1'
        arrtmp1 = reshape(torqval_atom, (/3, itmp/))

        !now apply the symmetries
        call rotate_kpoints(rotmat, nkpts*ndeg*natyp, arrtmp1, nsym, isym, nkpts_tmp, arrtmp2)
        if(nkpts_r*ndeg*natyp /= nkpts_tmp) stop 'nkpts_r*ndeg /= nkpts_tmp'

        !transform back to old shape
        allocate(torqval_atom_r(3,natyp,ndeg,nkpts_r), STAT=ierr)
        if(ierr/=0) stop 'problem allocating torqval_atom_r'
        torqval_atom_r = reshape(arrtmp2,(/3,natyp,ndeg,nkpts_r/))

        deallocate(arrtmp1, arrtmp2)
      end if!allocated(torqval_atom)

      if(allocated(spinvec_atom))then
        !next apply symmetries to spinvec_atom; as a first step, flatten the array
        itmp = size(spinvec_atom)/3
        allocate(arrtmp1(3,itmp), STAT=ierr)
        if(ierr/=0) stop 'problem alloaction arrtmp1'
        arrtmp1 = reshape(spinvec_atom, (/3, itmp/))

        !now apply the symmetries
        call rotate_kpoints(rotmat, nkpts*ndeg*natyp, arrtmp1, nsym, isym, nkpts_tmp, arrtmp2)
        if(nkpts_r*ndeg*natyp /= nkpts_tmp) stop 'nkpts_r*ndeg /= nkpts_tmp'

        !transform back to old shape
        allocate(spinvec_atom_r(3,natyp,ndeg,nkpts_r), STAT=ierr)
        if(ierr/=0) stop 'problem allocating spinvec_atom_r'
        spinvec_atom_r = reshape(arrtmp2,(/3,natyp,ndeg,nkpts_r/))

        deallocate(arrtmp1, arrtmp2)
      end if!allocated(spinvec_atom)

      if(allocated(spinflux))then
        !next apply symmetries to spinflux; as a first step, flatten the array
        itmp = size(spinflux)/3
        allocate(arrtmp1(3,itmp), STAT=ierr)
        if(ierr/=0) stop 'problem alloaction arrtmp1'
        arrtmp1 = reshape(spinflux, (/3, itmp/))

        !now apply the symmetries
        call rotate_kpoints(rotmat, nkpts*ndeg*natyp, arrtmp1, nsym, isym, nkpts_tmp, arrtmp2)
        if(nkpts_r*ndeg*natyp /= nkpts_tmp) stop 'nkpts_r*ndeg*natyp /= nkpts_tmp'

        !transform back to old shape
        allocate(spinflux_r(3,natyp,ndeg,nkpts_r), STAT=ierr)
        if(ierr/=0) stop 'problem allocating spinflux_r'
        spinflux_r = reshape(arrtmp2,(/3,natyp,ndeg,nkpts_r/))

        deallocate(arrtmp1, arrtmp2)
      end if!allocated(spinflux)

    else
      nkpts_r = nkpts

      allocate(fermivel_r(3,nkpts_r), areas_r(nkpts_r))
      fermivel_r = fermivel
      areas_r = areas

      if(allocated(torqval)) then
        allocate(torqval_r(3,ndeg,nkpts_r))
        torqval_r  = torqval
      end if!(allocated(torqval))

      if(allocated(torqval_atom)) then
        allocate(torqval_atom_r(3,natyp,ndeg,nkpts_r))
        torqval_atom_r  = torqval_atom
      end if!(allocated(torqval))

      if(allocated(spinvec_atom)) then
        allocate(spinvec_atom_r(3,natyp,ndeg,nkpts_r))
        spinvec_atom_r  = spinvec_atom
      end if!(allocated(spinvec))

      if(allocated(spinflux)) then
        allocate(spinflux_r(3,natyp,ndeg,nkpts_r))
        spinflux_r = spinflux
      end if!(allocated(spinflux))
    end if!nsym>1

    !compute response functions
    integ_torkance      = 0d0
    integ_torkance_atom = 0d0
    integ_damping       = 0d0
    integ_spinacc_atom  = 0d0
    integ_spinflux      = 0d0

    do ikp=1,nkpts_r
     vabs = sqrt(sum(fermivel_r(:,ikp)**2))
     do ideg=1,ndeg
      do i1=1,3
       do i2=1,3
        if(allocated(torqval))      &
        &  integ_torkance(i2,i1)        = integ_torkance(i2,i1)        + areas_r(ikp)*torqval_r(i1,ideg,ikp)*fermivel_r(i2,ikp)/(vabs)
        if(allocated(torqval_atom)) &
        &  integ_torkance_atom(i2,i1,:) = integ_torkance_atom(i2,i1,:) + areas_r(ikp)*torqval_atom_r(i1,:,ideg,ikp)*fermivel_r(i2,ikp)/(vabs)
        if(allocated(torqval))      &
        &  integ_damping(i2,i1)         = integ_damping(i2,i1)         + areas_r(ikp)*torqval_r(i1,ideg,ikp)*torqval_r(i2,ideg,ikp)/(vabs)
        if(allocated(spinvec_atom)) &
        &  integ_spinacc_atom(i2,i1,:)  = integ_spinacc_atom(i2,i1,:)  + areas_r(ikp)*spinvec_atom_r(i1,:,ideg,ikp)*fermivel_r(i2,ikp)/(vabs)
        if(allocated(spinflux))     &
        &  integ_spinflux(i2,i1,:)      = integ_spinflux(i2,i1,:)      + areas_r(ikp)*spinflux_r(i1,:,ideg,ikp)*fermivel_r(i2,ikp)/(vabs)
       end do!i2
      end do!i1
     end do!ideg
    end do!ikp

    if(allocated(torqval)) then
      torkance = integ_torkance/BZVol*alat/tpi

      if(myrank==master)then
        if(nsym/=1) write(*,*) 'ATTENTION: nsym =/= 1, make sure the following output is correct'
        write(*,'(A)')
        write(*,'(A)') "Torkance / relaxation time [in constant relaxation time approximation]:"
        write(*,'(A)') "  in units of e*abohr Rydberg:"
        write(*,'(3(ES25.16))') torkance
        write(*,'(A)')   
        write(*,'(A)') "  in units of e*abohr / fs:     "
        write(*,'(3(ES25.16))') torkance*RyToinvfs
        write(*,'(A)')
        write(*,'(A)') "Torkance at room temperature [in constant relaxation time approximation]:"
        write(*,'(A)') "  in units of e*abohr using hbar/(2*tau) = 25 meV:"
        write(*,'(3(ES25.16))') torkance*13.60569253/(2*25*0.001)
        write(*,'(A)')
      end if!myrank==master
    end if!allocated(torqval)

    if(allocated(torqval_atom)) then
      torkance_atom = integ_torkance_atom/BZVol*alat/tpi

      if(myrank==master)then
        open(unit=iounit,file="torkance.CRTA.int.dat",form='formatted',action='write')
        write(iounit,'(1I8)') natyp
        write(iounit,'(9(ES25.16))') torkance_atom*13.60569253/(2*25*0.001)
        close(iounit)
      end if!myrank==master
    end if!allocated(torqval_atom)

    if(allocated(torqval)) then
      damping = integ_damping/BZVol

      if(myrank==master)then
        if(nsym/=1) write(*,*) 'ATTENTION: nsym =/= 1, make sure the following output is correct'
        write(*,'(A)') "Gilbert damping in constant relaxation time approximation:"
        write(*,'(A)') "In units of  Rydberg:"
        write(*,'(3(ES25.16))') damping
        write(*,'(A)') "In units of eV:"
        write(*,'(3(ES25.16))') damping*13.60569253
        write(*,'(A)') "For room temperature (hbar/(2*tau) = 25 meV):"
        write(*,'(3(ES25.16))') damping*13.60569253/(2*25*0.001)
      end if!myrank==master
    end if!allocated(torqval)

    if(allocated(spinvec_atom)) then
      spinacc_atom = -integ_spinacc_atom/BZVol*alat/tpi ! units of (e a_0 mu_B)/(Rd)

      if(myrank==master)then
        open(unit=iounit,file="spinacc_resp.CRTA.int.dat",form='formatted',action='write')
        write(iounit,'(1I8)') natyp
        write(iounit,'(9(ES25.16))') spinacc_atom*13.60569253/(2*25*0.001)
        close(iounit)
      end if!myrank==master
    end if!allocated(spinvec_atom)

    if(allocated(spinflux)) then
      spinflux_resp = integ_spinflux/BZVol*alat/tpi

      if(myrank==master)then
        open(unit=iounit,file="spinflux_resp.CRTA.int.dat",form='formatted',action='write')
        write(iounit,'(1I8)') natyp
        write(iounit,'(9(ES25.16))') spinflux_resp*13.60569253/(2*25*0.001)
        close(iounit)
      end if!myrank==master
    end if!allocated(spinflux)

  end subroutine calculate_response_functions_CRTA_int


  subroutine calculate_response_functions_CRTA_vis(nBZdim,nsym,isym,rotmat,alat,BZVol,ndeg,natyp,nkpts,nkpts_all,kpt2irr,irr2kpt,kpoints,fermivel,torqval,torqval_atom,spinvec_atom,spinflux)
    use mpi
    use mod_mympi,      only: myrank, master
    use mod_mathtools,  only: tpi, crossprod
    use mod_symmetries, only: rotate_kpoints, expand_visarrays
    use mod_mathtools,  only: simple_integration_general
    implicit none

    integer,          intent(in)  :: nBZdim, nsym, ndeg, natyp, nkpts, nkpts_all, kpt2irr(nkpts_all), irr2kpt(nkpts), isym(nsym)
    double precision, intent(in)  :: alat, BZVol, kpoints(3,nkpts)
    double precision, intent(in)  :: rotmat(64,3,3), fermivel(3,nkpts)
    double precision, allocatable, intent(in)  :: torqval(:,:,:)      !torqval(3,ndeg,nkpts)
    double precision, allocatable, intent(in)  :: torqval_atom(:,:,:,:) !torqval_atom(3,natypd,ndeg,nkpts)
    double precision, allocatable, intent(in)  :: spinvec_atom(:,:,:,:) !spinvec_atom(3,natypd,ndeg,nkpts)
    double precision, allocatable, intent(in)  :: spinflux(:,:,:,:)   !spinflux(3,natypd,ndegen,nkpts)

    integer,          allocatable :: irr2kpt_r(:), kpt2irr_r(:)
    double precision, allocatable :: kpoints_r(:,:), fermivel_r(:,:), arrtmp1(:,:), arrtmp2(:,:)
    double precision, allocatable :: torqval_r(:,:,:), torqval_atom_r(:,:,:,:), spinvec_atom_r(:,:,:,:), spinflux_r(:,:,:,:)

    integer :: ierr, ikp, i1, i2, i3, iat, ideg, itmp, nkpts_tmp, nkpts_r, nkpts_all_r, iline
    integer :: pointspick(nBZdim)
    double precision :: kpoints_corners(3,nBZdim), fermi_velocity_corners(3,nBZdim), torq_corners(3,ndeg,nBZdim), torq_corners_atom(3,natyp,ndeg,nBZdim), spinvec_corners_atom(3,natyp,ndeg,nBZdim), spinflux_corners(3,natyp,ndeg,nBZdim)
    double precision :: k21(3), k31(3), kcross(3), area, integrand(nBZdim), dinteg, d_dos, dos
    double precision :: torkance(3,3), torkance_atom(3,3,natyp), damping(3,3), spinacc_atom(3,3,natyp), spinflux_resp(3,3,natyp)
    double precision :: integ_torkance(3,3), integ_torkance_atom(3,3,natyp), integ_damping(3,3), integ_spinacc_atom(3,3,natyp), integ_spinflux(3,3,natyp)
    double precision :: vabs
    double precision, parameter :: RyToinvfs = 20.67068667282055d0
    integer, parameter :: iounit=13518

    if(nsym>1)then
      !apply symmetries to kpoints and visarray
      call rotate_kpoints(rotmat, nkpts, kpoints, nsym, isym, nkpts_r, kpoints_r)
      call expand_visarrays(nsym, nkpts_all, nkpts, kpt2irr, irr2kpt, kpt2irr_r, irr2kpt_r)

      !apply symmetries to fermivel and areas
      call rotate_kpoints(rotmat, nkpts, fermivel, nsym, isym, nkpts_r, fermivel_r)

      if(allocated(torqval))then
        !next apply symmetries to torqval; as a first step, flatten the array
        itmp = size(torqval)/3
        allocate(arrtmp1(3,itmp), STAT=ierr)
        if(ierr/=0) stop 'problem alloaction arrtmp1'
        arrtmp1 = reshape(torqval, (/3, itmp/))

        !now apply the symmetries
        call rotate_kpoints(rotmat, nkpts*ndeg, arrtmp1, nsym, isym, nkpts_tmp, arrtmp2)
        if(nkpts_r*ndeg /= nkpts_tmp) stop 'nkpts_r*ndeg /= nkpts_tmp'

        !transform back to old shape
        allocate(torqval_r(3,ndeg,nkpts_r), STAT=ierr)
        if(ierr/=0) stop 'problem allocating torqval_r'
        torqval_r = reshape(arrtmp2,(/3,ndeg,nkpts_r/))

        deallocate(arrtmp1, arrtmp2)
      end if!allocated(torqval)

      if(allocated(torqval_atom))then
        !next apply symmetries to torqval_atom; as a first step, flatten the array
        itmp = size(torqval_atom)/3
        allocate(arrtmp1(3,itmp), STAT=ierr)
        if(ierr/=0) stop 'problem alloaction arrtmp1'
        arrtmp1 = reshape(torqval_atom, (/3, itmp/))

        !now apply the symmetries
        call rotate_kpoints(rotmat, nkpts*ndeg*natyp, arrtmp1, nsym, isym, nkpts_tmp, arrtmp2)
        if(nkpts_r*ndeg*natyp /= nkpts_tmp) stop 'nkpts_r*ndeg /= nkpts_tmp'

        !transform back to old shape
        allocate(torqval_atom_r(3,natyp,ndeg,nkpts_r), STAT=ierr)
        if(ierr/=0) stop 'problem allocating torqval_atom_r'
        torqval_atom_r = reshape(arrtmp2,(/3,natyp,ndeg,nkpts_r/))

        deallocate(arrtmp1, arrtmp2)
      end if!allocated(torqval_atom)

      if(allocated(spinvec_atom))then
        !next apply symmetries to spinvec_atom; as a first step, flatten the array
        itmp = size(spinvec_atom)/3
        allocate(arrtmp1(3,itmp), STAT=ierr)
        if(ierr/=0) stop 'problem alloaction arrtmp1'
        arrtmp1 = reshape(spinvec_atom, (/3, itmp/))

        !now apply the symmetries
        call rotate_kpoints(rotmat, nkpts*ndeg*natyp, arrtmp1, nsym, isym, nkpts_tmp, arrtmp2)
        if(nkpts_r*ndeg*natyp /= nkpts_tmp) stop 'nkpts_r*ndeg /= nkpts_tmp'

        !transform back to old shape
        allocate(spinvec_atom_r(3,natyp,ndeg,nkpts_r), STAT=ierr)
        if(ierr/=0) stop 'problem allocating spinvec_atom_r'
        spinvec_atom_r = reshape(arrtmp2,(/3,natyp,ndeg,nkpts_r/))

        deallocate(arrtmp1, arrtmp2)
      end if!allocated(spinvec_atom)

      if(allocated(spinflux))then
        !next apply symmetries to spinflux; as a first step, flatten the array
        itmp = size(spinflux)/3
        allocate(arrtmp1(3,itmp), STAT=ierr)
        if(ierr/=0) stop 'problem alloaction arrtmp1'
        arrtmp1 = reshape(spinflux, (/3, itmp/))

        !now apply the symmetries
        call rotate_kpoints(rotmat, nkpts*ndeg*natyp, arrtmp1, nsym, isym, nkpts_tmp, arrtmp2)
        if(nkpts_r*ndeg*natyp /= nkpts_tmp) stop 'nkpts_r*ndeg*natyp /= nkpts_tmp'

        !transform back to old shape
        allocate(spinflux_r(3,natyp,ndeg,nkpts_r), STAT=ierr)
        if(ierr/=0) stop 'problem allocating spinflux_r'
        spinflux_r = reshape(arrtmp2,(/3,natyp,ndeg,nkpts_r/))

        deallocate(arrtmp1, arrtmp2)
      end if!allocated(spinflux)

      nkpts_all_r = nkpts_all*nsym

    else
      nkpts_r = nkpts
      allocate(fermivel_r(3,nkpts_r), kpoints_r(3,nkpts), kpt2irr_r(nkpts_all))
      kpoints_r  = kpoints
      kpt2irr_r  = kpt2irr
      nkpts_all_r = nkpts_all*nsym
      fermivel_r = fermivel

      if(allocated(torqval)) then
        allocate(torqval_r(3,ndeg,nkpts_r))
        torqval_r  = torqval
      end if!(allocated(torqval))

      if(allocated(torqval_atom)) then
        allocate(torqval_atom_r(3,natyp,ndeg,nkpts_r))
        torqval_atom_r  = torqval_atom
      end if!(allocated(torqval))

      if(allocated(spinvec_atom)) then
        allocate(spinvec_atom_r(3,natyp,ndeg,nkpts_r))
        spinvec_atom_r  = spinvec_atom
      end if!(allocated(spinvec))

      if(allocated(spinflux)) then
        allocate(spinflux_r(3,natyp,ndeg,nkpts_r))
        spinflux_r = spinflux
      end if!(allocated(spinflux))
    end if

    !compute response functions
    integ_torkance      = 0d0
    integ_torkance_atom = 0d0
    integ_damping       = 0d0
    integ_spinacc_atom  = 0d0
    integ_spinflux      = 0d0

    dos    = 0d0
    do iline=1,nkpts_all_r/nBZdim

      do i3=1,nBZdim
        pointspick(i3) = kpt2irr_r((iline-1)*nBZdim+i3)
      end do!i3

      !rewrite corner point data
      kpoints_corners(:,:)        = kpoints_r(:,pointspick(:))
      fermi_velocity_corners(:,:) = fermivel_r(:,pointspick(:))
      if(allocated(torqval)) torq_corners(:,:,:) = torqval_r(:,:,pointspick(:))
      if(allocated(torqval_atom)) torq_corners_atom(:,:,:,:) = torqval_atom_r(:,:,:,pointspick(:))
      if(allocated(spinvec_atom)) spinvec_corners_atom(:,:,:,:) = spinvec_atom_r(:,:,:,pointspick(:))
      if(allocated(spinflux)) spinflux_corners(:,:,:,:) = spinflux_r(:,:,:,pointspick(:))

      if(nBZdim==3)then
        !compute area
        k21 = kpoints_corners(:,2) - kpoints_corners(:,1)
        k31 = kpoints_corners(:,3) - kpoints_corners(:,1)
        call crossprod(k21, k31, kcross)
        area = 0.5d0*sqrt(sum(kcross**2))
      elseif(nBZdim==2)then
        area = sqrt(sum((kpoints_corners(:,2) - kpoints_corners(:,1))**2))
      else!nBZdim
        stop 'nBZdim is neither 2 nor 3'
      end if!nBZdim

      do ideg=1,ndeg
       do i1=1,3
        do i2=1,3
         if(allocated(torqval))then
           integrand(:) = torq_corners(i1,ideg,:)*fermi_velocity_corners(i2,:)
           call simple_integration_general(nBZdim, area, fermi_velocity_corners, integrand, dinteg, d_dos)
           integ_torkance(i2,i1) = integ_torkance(i2,i1) + dinteg
         end if!allocated(torqval)
         if(allocated(torqval_atom))then
           do iat=1,natyp
             integrand(:) = torq_corners_atom(i1,iat,ideg,:)*fermi_velocity_corners(i2,:)
             call simple_integration_general(nBZdim, area, fermi_velocity_corners, integrand, dinteg, d_dos)
             integ_torkance_atom(i2,i1,iat) = integ_torkance_atom(i2,i1,iat) + dinteg
           end do!iat
         end if!allocated(torqval_atom)
         if(allocated(torqval))then
           integrand(:) = torq_corners(i1,ideg,:)*torq_corners(i2,ideg,:)
           call simple_integration_general(nBZdim, area, fermi_velocity_corners, integrand, dinteg, d_dos)
           integ_damping(i2,i1) = integ_damping(i2,i1) + dinteg
         end if!allocated(torqval)
         if(allocated(spinvec_atom))then
           do iat=1,natyp
             integrand(:) = spinvec_corners_atom(i1,iat,ideg,:)*fermi_velocity_corners(i2,:)
             call simple_integration_general(nBZdim, area, fermi_velocity_corners, integrand, dinteg, d_dos)
             integ_spinacc_atom(i2,i1,iat) = integ_spinacc_atom(i2,i1,iat) + dinteg
           end do!iat
         end if!allocated(spinvec_atom)
         if(allocated(spinflux))then
           do iat=1,natyp
             integrand(:) = spinflux_corners(i1,iat,ideg,:)*fermi_velocity_corners(i2,:)
             call simple_integration_general(nBZdim, area, fermi_velocity_corners, integrand, dinteg, d_dos)
             integ_spinflux(i2,i1,iat) = integ_spinflux(i2,i1,iat) + dinteg
           end do!iat
         end if!allocated(spinflux)
        end do!i2
       end do!i1
      end do!ideg
      dos = dos+d_dos
    end do!iline
    if(allocated(torqval)) then
      torkance = integ_torkance/BZVol*alat/tpi

      if(myrank==master)then
        if(nsym/=1) write(*,*) 'ATTENTION: nsym =/= 1, make sure the following output is correct'
        write(*,'(A)')
        write(*,'(A)') "Torkance / relaxation time [in constant relaxation time approximation]:"
        write(*,'(A)') "  in units of e*abohr Rydberg:"
        write(*,'(3(ES25.16))') torkance
        write(*,'(A)')
        write(*,'(A)') "In units of e*abohr / fs:"
        write(*,'(3(ES25.16))') torkance*RyToinvfs
        write(*,'(A)')
        write(*,'(A)') "Torkance at room temperature [in constant relaxation time approximation]:"
        write(*,'(A)') "  in units of e*abohr using hbar/(2*tau) = 25 meV:"
        write(*,'(3(ES25.16))') torkance*13.60569253/(2*25*0.001)
      end if!myrank==master
    end if!allocated(torqval)
    
    if(allocated(torqval_atom)) then
      torkance_atom = integ_torkance_atom/BZVol*alat/tpi

      if(myrank==master)then
        open(unit=iounit,file="torkance.CRTA.vis.dat",form='formatted',action='write')
        write(iounit,'(1I8)') natyp
        write(iounit,'(9(ES25.16))') torkance_atom*13.60569253/(2*25*0.001)
        close(iounit)
      end if!myrank==master
    end if!allocated(torqval_atom)

    if(allocated(torqval)) then
      damping = integ_damping/BZVol

      if(myrank==master)then
        write(*,'(A)') "Gilbert damping in constant relaxation time approximation:"
        write(*,'(A)') "In units of  Rydberg:"
        write(*,'(3(ES25.16))') damping
        write(*,'(A)') "In units of eV:"
        write(*,'(3(ES25.16))') damping*13.60569253
        write(*,'(A)') "For room temperature (hbar/(2*tau) = 25 meV):"
        write(*,'(3(ES25.16))') damping*13.60569253/(2*25*0.001)
      end if!myrank==master
    end if!allocated(torqval)

    if(allocated(spinvec_atom)) then
      spinacc_atom = -integ_spinacc_atom*(2)**0.5/BZVol*alat/tpi ! sqrt(2) for mu_B

      if(myrank==master)then
        open(unit=iounit,file="spinacc_resp.CRTA.vis.dat",form='formatted',action='write')
        write(iounit,'(1I8)') natyp
        write(iounit,'(9(ES25.16))') spinacc_atom*13.60569253/(2*25*0.001)
        close(iounit)
      end if!myrank==master
    end if!allocated(spinvec_atom)

    if(allocated(spinflux)) then
      spinflux_resp = integ_spinflux/BZVol*alat/tpi

      if(myrank==master)then
        open(unit=iounit,file="spinflux_resp.CRTA.vis.dat",form='formatted',action='write')
        write(iounit,'(1I8)') natyp
        write(iounit,'(9(ES25.16))') spinflux_resp*13.60569253/(2*25*0.001)
        close(iounit)
      end if!myrank==master
    end if!allocated(spinflux)

  end subroutine calculate_response_functions_CRTA_vis

  subroutine calculate_torkance_CRTA_vis_simpson2D(nsym,isym,rotmat,alat,BZVol,ndeg,nkpts,nkpts_all,kpt2irr,irr2kpt,kpoints,fermivel,torqval)
    use mpi
    use mod_mympi,      only: myrank, master
    use mod_mathtools,  only: tpi, crossprod
    use mod_symmetries, only: rotate_kpoints, expand_visarrays
    use mod_mathtools,  only: simple_integration_general, simpson2D_integration
    implicit none

    integer,          intent(in)  :: nsym, ndeg, nkpts, nkpts_all, kpt2irr(nkpts_all), irr2kpt(nkpts), isym(nsym)
    double precision, intent(in)  :: alat, BZVol, kpoints(3,nkpts)
    double precision, intent(in)  :: rotmat(64,3,3), fermivel(3,nkpts), torqval(3,ndeg,nkpts)

    integer,          allocatable :: irr2kpt_r(:), kpt2irr_r(:),  kpt2irr_r_ord(:), band_indices_r(:), nkpts_band(:)
    double precision, allocatable :: kpoints_r(:,:), fermivel_r(:,:), torqval_r(:,:,:), arrtmp1(:,:), arrtmp2(:,:)

    integer :: ierr, ikp, i1, i2, i3, ideg, itmp, nkpts_tmp, nkpts_r, nkpts_all_r, nparts, ipart, iband, nbands, offset
    integer :: pointspick(3)
    logical :: fourmultiple
    double precision :: kpoints_corners(3,3), fermi_velocity_corners(3,3), torq_corners(3,ndeg,3), k21(3), k31(3), kcross(3)
    double precision :: d(2), integrand(3), integrand_damping(3), dinteg, dinteg_damping, d_dos, dos
    double precision :: torkance(3,3), damping(3,3)
    double precision :: integ(3,3), integ_damping(3,3), vabs
    double precision, parameter :: RyToinvfs = 20.67068667282055d0

    nkpts_all_r = nkpts_all*nsym

    if(nsym>1)then
      !apply symmetries to kpoints and visarray
      call rotate_kpoints(rotmat, nkpts, kpoints, nsym, isym, nkpts_r, kpoints_r)
      call expand_visarrays(nsym, nkpts_all, nkpts, kpt2irr, irr2kpt, kpt2irr_r, irr2kpt_r)

      allocate(kpt2irr_r_ord(nkpts_all_r),band_indices_r(nkpts_all_r), STAT=ierr)
      if(ierr/=0) stop 'problem allocating kpt2irr_r_ord and band_indices_r'

      !order kpoints according to lines
      call order_lines(kpt2irr_r, nkpts_all_r, kpt2irr_r_ord, band_indices_r, nkpts_band, nbands)
      write(*,*) 'ordered k-points according to bands'
      write(*,*) 'number of bands that were found : ', nbands

      !apply symmetries to fermivel and areas
      call rotate_kpoints(rotmat, nkpts, fermivel, nsym, isym, nkpts_r, fermivel_r)

      !next apply symmetries to torqval; as a first step, flatten the array
      itmp = size(torqval)/3
      allocate(arrtmp1(3,itmp), STAT=ierr)
      if(ierr/=0) stop 'problem alloaction arrtmp1'
      arrtmp1 = reshape(torqval, (/3, itmp/))

      !now apply the symmetries
      call rotate_kpoints(rotmat, nkpts*ndeg, arrtmp1, nsym, isym, nkpts_tmp, arrtmp2)
      if(nkpts_r*ndeg /= nkpts_tmp) stop 'nkpts_r*ndeg /= nkpts_tmp'

      !transform back to old shape
      allocate(torqval_r(3,ndeg,nkpts_r), STAT=ierr)
      if(ierr/=0) stop 'problem allocating torqval_r'
      torqval_r = reshape(arrtmp2,(/3,ndeg,nkpts_r/))

      deallocate(arrtmp1, arrtmp2)

      !nkpts_all_r = nkpts_all*nsym

    else
      write(*,*) 'it seems you want to use Simpson rule integration with nsym=1'
      write(*,*) 'this is not implemented yet (closed bands can not be found)'
      nkpts_r = nkpts
      allocate(fermivel_r(3,nkpts_r), torqval_r(3,ndeg,nkpts_r), kpoints_r(3,nkpts), kpt2irr_r(nkpts_all))
      kpoints_r  = kpoints
      kpt2irr_r  = kpt2irr
      fermivel_r = fermivel
      torqval_r  = torqval
      !nkpts_all_r = nkpts_all*nsym
    end if

    !torkance and damping
    integ  = 0d0
    integ_damping = 0d0
    dos    = 0d0


    do iband=1, nbands     
      if(MOD(nkpts_band(iband), 4).ne.0) fourmultiple=.true.
      if(MOD(nkpts_band(iband), 4) == 0) fourmultiple=.false.
      !find how many parts of 3 or 2 kpts are in the band
      if(fourmultiple) nparts=nkpts_band(iband)/4
      if(.not.fourmultiple) nparts=(nkpts_band(iband)+2)/4

      do ipart=1,nparts
        if(iband>1)offset=sum(nkpts_band(1:iband-1))
        if(iband==1)offset=0
        pointspick(1) = kpt2irr_r_ord(offset+(ipart-1)*4+1)
        pointspick(2) = kpt2irr_r_ord(offset+(ipart-1)*4+2)
        if(.not.fourmultiple .and. ipart==nparts) then 
          !if no more kpt to save, save the same twice (but won't be used)
          pointspick(3) = kpt2irr_r_ord(offset+(ipart-1)*4+2)
        else
          pointspick(3) = kpt2irr_r_ord(offset+(ipart-1)*4+4)
        end if

        !rewrite corner point data
        kpoints_corners(:,:)        = kpoints_r(:,pointspick(:))
        fermi_velocity_corners(:,:) = fermivel_r(:,pointspick(:))
        torq_corners(:,:,:)         = torqval_r(:,:,pointspick(:))

        d(1) = sqrt(sum((kpoints_corners(:,2) - kpoints_corners(:,1))**2))
        d(2) = sqrt(sum((kpoints_corners(:,3) - kpoints_corners(:,2))**2))

        do ideg=1,ndeg
         do i1=1,3
          do i2=1,3
           integrand(:) = torq_corners(i1,ideg,:)*fermi_velocity_corners(i2,:)
           integrand_damping(:) = torq_corners(i1,ideg,:)*torq_corners(i2,ideg,:)
         
           if (.not.fourmultiple .and. ipart==nparts) then
             call simple_integration_general(2, d(1), fermi_velocity_corners(:,1:2), integrand(1:2), dinteg, d_dos)
             call simple_integration_general(2, d(1), fermi_velocity_corners(:,1:2), integrand_damping(1:2), dinteg_damping, d_dos)
           else
             call simpson2D_integration(d,fermi_velocity_corners,integrand, dinteg, d_dos)
             call simpson2D_integration(d,fermi_velocity_corners,integrand_damping, dinteg_damping, d_dos)
           end if

           integ(i2,i1) = integ(i2,i1) + dinteg
           integ_damping(i2,i1) = integ_damping(i2,i1) + dinteg_damping
          end do!i2
         end do!i1
        end do!ideg
        dos = dos+d_dos
      end do!ipart
    end do!iband

    deallocate(kpt2irr_r_ord,band_indices_r)

    torkance = integ/BZVol*alat/tpi

    if(myrank==master)then
      if(nsym/=1) write(*,*) 'ATTENTION: nsym =/= 1, make sure the following output is correct'
!      write(*,'(A)') "Torkance per relaxation time in constant relaxation time approximation:"
!      write(*,'(A)') "In units of e*abohr Rydberg:"
      write(*,'(A)')
      write(*,'(A)') "Torkance / relaxation time [in constant relaxation time approximation]:"
      write(*,'(A)') "  in units of e*abohr Rydberg:"
      write(*,'(3(ES25.16))') torkance
      write(*,'(A)')
      write(*,'(A)') "  in units of e*abohr / fs:     "
!      write(*,'(A)') "In units of e*abohr / fs:"
      write(*,'(3(ES25.16))') torkance*RyToinvfs
!      write(*,'(A)') "In units of e*abohr for room temperature (hbar/(2*tau) = 25 meV) :"
      write(*,'(A)')
      write(*,'(A)') "Torkance at room temperature [in constant relaxation time approximation]:"
      write(*,'(A)') "  in units of e*abohr using hbar/(2*tau) = 25 meV:"
      write(*,'(3(ES25.16))') torkance*13.60569253/(2*25*0.001)
    end if!myrank==master

    damping = integ_damping/BZVol

    if(myrank==master)then
      if(nsym/=1) write(*,*) 'ATTENTION: nsym =/= 1, make sure the following output is correct'
      write(*,'(A)') "Gilbert damping in constant relaxation time approximation:"
      write(*,'(A)') "In units of  Rydberg:"
      write(*,'(3(ES25.16))') damping
      write(*,'(A)') "In units of eV:"
      write(*,'(3(ES25.16))') damping*13.60569253
      write(*,'(A)') "For room temperature (hbar/(2*tau) = 25 meV):"
      write(*,'(3(ES25.16))') damping*13.60569253/(2*25*0.001)
    end if!myrank==master

  end subroutine calculate_torkance_CRTA_vis_simpson2D


  subroutine order_lines(kpt2irr, nkpts_all, kpt2irr_ord, band_indices, nkpts_band, nbands)
    integer,         intent(in)       :: kpt2irr(nkpts_all)
    integer,         intent(in)       :: nkpts_all
    integer,         intent(out)      :: band_indices(nkpts_all), kpt2irr_ord(nkpts_all), nbands
    integer,         allocatable      :: nkpts_band(:)
    logical                           :: beginning_of_band
    integer                           :: i1, i2, i3, i_ord, i_band, ierr, i_ord_new

    ! This subroutine order k-pts according to the band to which they belong
    ! Output  : kpt2irr_ord contains the ordered k-pts
    !           band_indices contains the band indices
    !           nkpts_band(nbands) contain the nb of kpts in each band
    ! Comment : At the moment works only if no band makes a loop in the IBZ

    i_ord =1
    i_band=1    
    !loop over the k-points to find the beginning of bands
    do i1=1,nkpts_all
      beginning_of_band=.true.

      ! check if i1 is the beginning of a band
      do i2=1,nkpts_all
        if(kpt2irr(i2)==kpt2irr(i1) .and. i1/=i2)then
          beginning_of_band=.false.
          exit
        end if!kpt2irr(i2)==kpt2irr(i1)
      end do!i2=i1+1,nkpts_all

      ! check if i1 was not already treated (it could be the end of a band)
      if (i_ord>1)then
        do i3=1, i_ord-1
          if(kpt2irr_ord(i3)==kpt2irr(i1))then
            beginning_of_band=.false.
            exit
          end if!kpt2irr_ord(i3)==kpt2irr(i1)
        end do!i3=1, i_ord-1
      end if!i_ord>1

      !save all kpts from the band
      if(beginning_of_band)then
        call traceback_band(i_ord,i_band,i1,kpt2irr,nkpts_all,kpt2irr_ord,band_indices,i_ord_new)
        i_ord=i_ord_new
        i_band=i_band+1
      end if!end_of_band=.true.

    end do!i1=1,nkpts_all    

    nbands=i_band-1
    allocate(nkpts_band(nbands), STAT=ierr)
    if(ierr/=0) stop 'problem allocating nkpts_band'

    ! count number of kpts for each band
    nkpts_band = 0
    do i1=1,nkpts_all
      do i_band=1, nbands
        if(band_indices(i1)==i_band) nkpts_band(i_band) = nkpts_band(i_band)+1
      end do!i_band=1, nbands
    end do! i1=1,nkpts_all

    if (sum(nkpts_band(:))/=nkpts_all) stop 'Some kpts were not found in order_lines(). &
                                 & Probably there is band that makes a loop in the IBZ. '
    
  end subroutine

  recursive subroutine traceback_band(i_ord, i_band, i, kpt2irr, nkpts_all, kpt2irr_ord, band_indices, i_ord_new)
    integer,         intent(in)   :: kpt2irr(nkpts_all)
    integer,         intent(in)   :: nkpts_all
    integer,         intent(in)   :: i_ord, i_band, i
    integer,         intent(out)  :: band_indices(nkpts_all), kpt2irr_ord(nkpts_all)
    integer                       :: i_ord_new
    integer                       :: i_sameline, j          
    logical                       :: end_of_band

    ! This routine calls itself recursively until the end of the band is reached 

    if(MOD(i, 2).ne.0) i_sameline = i+1 ! i is odd
    if(MOD(i, 2) == 0) i_sameline = i-1 ! i is even
    
    ! save the 2 kpts of the current line and the band indices 
    kpt2irr_ord(i_ord)=kpt2irr(i)
    kpt2irr_ord(i_ord+1)=kpt2irr(i_sameline) 
    band_indices(i_ord)=i_band
    band_indices(i_ord+1)=i_band

    end_of_band=.true.
    ! check if i_sameline is the end of a band
    do j=1,nkpts_all
      if(kpt2irr(j)==kpt2irr(i_sameline) .and. j.ne.i_sameline)then
        end_of_band=.false.
        exit
      end if!kpt2irr(j)==kpt2irr(i_sameline) .and. j.ne.i_sameline
    end do!j=1,nkpts_all

    if(.not.end_of_band)then
      call traceback_band(i_ord+2, i_band, j, kpt2irr, nkpts_all, kpt2irr_ord, band_indices, i_ord_new) 
    else
      i_ord_new=i_ord+2
    end if!end_of_band==.false. 
 
  end subroutine



  subroutine calculate_spinmixing_int(nsym, ndeg,nsqa,nkpts,areas,fermivel,spinval,spinmix,dos,printout,BZVol)
    use mpi
    use mod_mympi, only: myrank, master
    implicit none

    integer,          intent(in)  :: nsym,ndeg, nsqa, nkpts
    double precision, intent(in)  :: areas(nkpts), fermivel(3,nkpts), spinval(ndeg,nsqa,nkpts)
    double precision, intent(out) :: spinmix(nsqa), dos
    logical,          intent(in)  :: printout
    double precision, intent(in), optional :: BZVol

    integer :: isqa, ikp
    double precision :: integ(nsqa), vabs

    integ  = 0d0
    dos  = 0d0
    do ikp=1,nkpts
      vabs = sqrt(sum(fermivel(:,ikp)**2))
      dos  = dos +areas(ikp)/vabs
      integ(:) = integ(:)+areas(ikp)*(1d0-abs(spinval(1,:,ikp)))/(2*vabs)
    end do!ikp

    spinmix = integ/dos

    if(printout .and. myrank==master)then
      write(*,'(A,ES25.16,A)') "DOS(E_F)=", nsym*dos, " [unnormalized]"
      if(present(BZVol)) write(*,'(A,ES25.16,A)') "DOS(E_F)=", nsym*dos/BZVol, " [normalized]"
      do isqa=1,nsqa
        write(*,'(A,I4,3(A,ES25.16))') "Elliott-Yafet parameter: isqa= ", isqa,", b^2= ", nsym*integ(isqa)," [unnormalized], ", spinmix(isqa), " [normalized by DOS]"
      end do!isqa
    end if!printout

  end subroutine calculate_spinmixing_int





  subroutine calculate_dos_int(nsym,isym,rotmat,alat,BZVol,nkpts,areas,fermivel,BZdim)
    use mpi
    use mod_mympi,      only: myrank, master
    use mod_mathtools,  only: tpi
    use mod_symmetries, only: rotate_kpoints, expand_areas
    implicit none

    integer,          intent(in)  :: nsym,nkpts, isym(nsym), BZdim
    double precision, intent(in)  :: rotmat(64,3,3), alat, BZVol, areas(nkpts), fermivel(3,nkpts)

    integer :: ikp, i1, i2, nkpts_r
    double precision :: dos, vabs, conductivity(3,3), dtmp

    double precision, allocatable :: fermivel_r(:,:), areas_r(:)

    double precision, parameter :: e2byhbar = 2.434134807664281d-4, abohr = 0.52917721092d-10, RyToinvfs = 20.67068667282055d0

    dos  = 0d0
    do ikp=1,nkpts
      vabs = sqrt(sum(fermivel(:,ikp)**2))
      dos  = dos +areas(ikp)/vabs
    end do!ikp

    if(myrank==master)then
      write(*,'(A,ES25.16,A)') "DOS(E_F)=", nsym*dos, " [unnormalized]"
      write(*,'(A,ES25.16,A)') "DOS(E_F)=", nsym*dos/BZVol, " [normalized]"
    end if!myrank==master

    if(myrank==master) write(*,'(A,I0,A,48I4)') 'nsym=', nsym, ', isym=', isym

    if(nsym>1)then
      call rotate_kpoints(rotmat, nkpts, fermivel, nsym, isym, nkpts_r, fermivel_r)
      call expand_areas(nsym,nkpts,areas,areas_r)
    else
      nkpts_r = nkpts
      allocate(fermivel_r(3,nkpts_r), areas_r(nkpts_r))
      fermivel_r = fermivel
      areas_r    = areas
    end if

    conductivity=0d0
    do ikp=1,nkpts_r
      vabs = sqrt(sum(fermivel_r(:,ikp)**2))
      dtmp = areas_r(ikp)/vabs
      do i1=1,3
        do i2=1,3
          conductivity(i2,i1) = conductivity(i2,i1) + dtmp*fermivel_r(i2,ikp)*fermivel_r(i1,ikp)
        end do!i2
      end do!i1
    end do!ikp

    if(BZdim==3)then    
      conductivity = conductivity/(alat*tpi**2)

      if(myrank==master)then
        write(*,'(A)') "Conductivity / relaxation time  [in constant relaxation time approximation]:"
        write(*,'(A)') "  in Rydberg:"
        write(*,'(3ES25.16)') conductivity*2 !factor 2 because e^2 = 2 in Rydberg units
        write(*,'(A)') "  in siemens/(meter * femtosecond):"
        write(*,'(3ES25.16)') conductivity*e2byhbar/abohr*RyToinvfs
        write(*,'(A)') "  in siemens/meter for room temperature (hbar/(2*tau) = 25 meV) :"
        write(*,'(3ES25.16)') conductivity*e2byhbar/abohr*13.60569253/(2*25*0.001)
      end if!myrank==master

    elseif (BZdim==2)then
      conductivity = conductivity*e2byhbar/(tpi**2)/alat/abohr

      if(myrank==master)then
        write(*,'(A)')
        write(*,'(A)') "2D MODE : The following values have to be divided by the thickness of the film"
        write(*,'(A)') "          in unit of the lattice constant (ALAT) to obtain the conductivity"
        write(*,'(A)') "          in the corresponding units."
        write(*,'(A)')
        write(*,'(A)') "Conductivity / relaxation time [in constant relaxation time approximation]:"
        write(*,'(A)') "  in (siemens/m)*Rydberg:"
        write(*,'(3ES25.16)') conductivity
        write(*,'(A)')
        write(*,'(A)') "  in (siemens/m)/(femtosecond):"
        write(*,'(3ES25.16)') conductivity*RyToinvfs
        write(*,'(A)')
        write(*,'(A)') "Conductivity at room temperature [in constant relaxation time approximation]:"
        write(*,'(A)') "  in (siemens/m) using hbar/(2*tau) = 25 meV :"
        write(*,'(3ES25.16)') conductivity*13.60569253/(2*25*0.001)
        write(*,'(A)')
      end if!myrank==master
    else
      stop 'BZdim is neiter 2 nor 3 !'
    endif

  end subroutine calculate_dos_int





  subroutine calculate_spinmixing_vis(nsym, ndeg,nsqa,nkpts,nkpts_all,kpt2irr,kpoints,fermivel,spinval,spinmix,dos,printout,BZVol,nBZdim)
    use mpi
    use mod_mympi, only: myrank, master
    use mod_mathtools, only: crossprod, simple_integration_general
    implicit none

    integer,          intent(in)  :: nsym,ndeg, nsqa, nkpts, nkpts_all, kpt2irr(nkpts_all), nBZdim
    double precision, intent(in)  :: kpoints(3,nkpts), fermivel(3,nkpts), spinval(ndeg,nsqa,nkpts)
    double precision, intent(out) :: spinmix(nsqa), dos
    logical,          intent(in)  :: printout
    double precision, intent(in)  :: BZVol

    integer :: isqa, iii, itri, pointspick(nBZdim)
    double precision :: integ(nsqa), vabs, kpoints_triangle(3,nBZdim), fermi_velocity_triangle(3,nBZdim), eyaf_triangle(nBZdim), k21(3), k31(3), kcross(3), area, d_integ, d_dos

    integ  = 0d0
    dos  = 0d0
    do itri=1,nkpts_all/nBZdim

      do iii=1,nBZdim
        pointspick(iii) = kpt2irr((itri-1)*nBZdim+iii)
      end do!iii

      !rewrite corner point data
      kpoints_triangle(:,:)        = kpoints(:,pointspick(:))
      fermi_velocity_triangle(:,:) = fermivel(:,pointspick(:))

      !compute area
      if(nBZdim==3)then
        k21 = kpoints_triangle(:,2) - kpoints_triangle(:,1)
        k31 = kpoints_triangle(:,3) - kpoints_triangle(:,1)
        call crossprod(k21, k31, kcross)
        area = 0.5d0*sqrt(sum(kcross**2))
      else if(nBZdim==2)then
        area = sqrt(sum((kpoints_triangle(:,2) - kpoints_triangle(:,1))**2))
      else
        stop 'nBZdim must be 2 or 3 in calculate_spinmixing_vis'
      end if

      do isqa=1,nsqa
        eyaf_triangle(:) = (1d0-abs(spinval(1,isqa,pointspick(:))))/2
        call simple_integration_general(nBZdim, area, fermi_velocity_triangle, eyaf_triangle, d_integ, d_dos)
        integ(isqa) = integ(isqa)+d_integ
      end do!isqa

      dos = dos+d_dos
    end do!ikp

    spinmix = integ/dos

    if(printout .and. myrank==master)then
      write(*,'(A,ES25.16,A)') "DOS(E_F)=", nsym*dos, " [unnormalized]"
      write(*,'(A,ES25.16,A)') "DOS(E_F)=", nsym*dos/BZVol, " [normalized]"
      do isqa=1,nsqa
        write(*,'(A,I4,3(A,ES25.16))') "Elliott-Yafet parameter: isqa= ", isqa,", b^2= ", nsym*integ(isqa)," [unnormalized], ", spinmix(isqa), " [normalized by DOS]"
      end do!isqa
    end if!printout

  end subroutine calculate_spinmixing_vis




  subroutine calculate_dos_vis(nsym,nkpts,nkpts_all,kpt2irr,kpoints,fermivel,dos,printout,BZVol)
    use mod_mympi, only: myrank, master
    use mod_mathtools, only: crossprod, simple_integration
    use mpi
    implicit none

    integer,          intent(in)  :: nsym, nkpts, nkpts_all, kpt2irr(nkpts_all)
    double precision, intent(in)  :: kpoints(3,nkpts), fermivel(3,nkpts)
    double precision, intent(out) :: dos
    logical,          intent(in)  :: printout
    double precision, intent(in), optional :: BZVol

    integer :: itri, pointspick(3)
    double precision :: integ, vabs, kpoints_triangle(3,3), fermi_velocity_triangle(3,3), dtmp3(3), k21(3), k31(3), kcross(3), area, dtmp, d_dos

    dtmp3 = 0d0
    dos  = 0d0

    do itri=1,nkpts_all/3

      pointspick(1) = kpt2irr((itri-1)*3+1)
      pointspick(2) = kpt2irr((itri-1)*3+2)
      pointspick(3) = kpt2irr((itri-1)*3+3)

      !rewrite corner point data
      kpoints_triangle(:,:)        = kpoints(:,pointspick(:))
      fermi_velocity_triangle(:,:) = fermivel(:,pointspick(:))

      !compute area
      k21 = kpoints_triangle(:,2) - kpoints_triangle(:,1)
      k31 = kpoints_triangle(:,3) - kpoints_triangle(:,1)
      call crossprod(k21, k31, kcross)
      area = 0.5d0*sqrt(sum(kcross**2))

      call simple_integration(area, fermi_velocity_triangle, dtmp3, dtmp, d_dos)
      dos = dos+d_dos

    end do!ikp

    if(printout .and. myrank==master)then
      write(*,'(A,ES25.16,A)') "DOS(E_F)=", nsym*dos, " [unnormalized]"
      if(present(BZVol)) write(*,'(A,ES25.16,A)') "DOS(E_F)=", nsym*dos/BZVol, " [normalized]"
    end if!printout

  end subroutine calculate_dos_vis



  integer function get_nsqa()
    implicit none
    if(.not.cfg_read .or. cfg%nsqa==-1 ) then
      call read_cfg(force_spinread=.true.)
    end if
    get_nsqa =cfg%nsqa
  end function



  subroutine calculate_and_save_weights(nkpts,nsym,isym,areas,fermivel)
    use mod_ioformat, only: fmt_fn_ext, filename_weights, ext_formatted
    implicit none
    integer, intent(in) :: nkpts, nsym, isym(nsym)
    double precision, intent(in) :: areas(nkpts), fermivel(3,nkpts)

    integer :: ikp
    double precision :: weights(nkpts)
    character(len=256) :: filename
    integer, parameter :: iounit=13516

    do ikp=1,nkpts
      weights(ikp) = areas(ikp)/sqrt(sum(fermivel(:,ikp)**2))
    end do

    write(filename,fmt_fn_ext) filename_weights, ext_formatted
    open(unit=iounit,file=trim(filename),form='formatted',action='write')
    write(iounit,'(2I8)') nkpts, nsym
    write(iounit,'(12I8)') isym
    write(iounit,'(10ES25.16)') weights
    close(iounit)

  end subroutine calculate_and_save_weights


end module mod_calconfs
