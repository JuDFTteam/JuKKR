!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module mod_fermisurf_3D

  implicit none

  private
  public :: find_fermisurface_3D

  integer, save :: tetcorners(4,6)=-1, tetdiags(2,4)=-1, cubedges(2,19)=-1, tetedges(6,6)=-1


  integer, parameter :: nkpmax = 512 ! shall be >> 32 to contain enough triangles per cube and band

contains

  subroutine find_fermisurface_3D( inc, lattice, cluster, tgmatrx, symmetries, nCub3, nFSiter, nROOTiter, nstepsconnect, &
                                 & nCut_iter, roottype, rooteps, lrefine, nrefinenew, nkpts_int, kpoints_int, areas_int  )

    use type_inc,       only: inc_type
    use type_data,      only: lattice_type, cluster_type, tgmatrx_type
    use mod_symmetries, only: symmetries_type, get_IBZwedge_faces
    use mod_mympi,      only: myrank, nranks, master
    use mod_ioformat,   only: filename_cubesinfo, ext_formatted, fmt_fn_ext
    use mod_iohelp,     only: file_present
    use mod_fermisurf_basic, only: get_cubesinfo_filename, read_cubesfile, save_cubesfile, mark_cubes_FScross, find_kpoints_irredset, save_kpointsfile_vis
#ifdef CPP_TIMING
    use mod_timing, only: timing_start, timing_stop
#endif
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    type(tgmatrx_type), intent(inout) :: tgmatrx
    type(symmetries_type), intent(in) :: symmetries
    integer,               intent(inout) :: nCub3(3)
    integer,               intent(in) :: nFSiter, nROOTiter, lrefine, nrefinenew, nCut_iter(:), roottype(:), nstepsconnect(:)
    double precision,      intent(in) :: rooteps(:)

    integer, intent(out) :: nkpts_int
    double precision, allocatable, intent(out) :: kpoints_int(:,:), areas_int(:)

    integer               :: nfaces
    double precision      :: bounds(3,2)
    double precision, allocatable :: nvec(:,:), dscal(:)

    integer :: nmarked, ntotal
    integer, allocatable :: imarked(:)

    integer :: nkpts_vis, nkpts_irr
    integer,          allocatable :: kpt2irr(:), irr2kpt(:), vis2int(:)
    double precision, allocatable :: kpoints_vis(:,:)

    integer :: icube, iter, iter2, ierr, ii, nCub3test(3), iterstart
    double precision :: dtmp
    character(len=256) :: filename, filetest
    logical :: l_cubesfile, l_exist

#ifdef CPP_TIMING
    call timing_start('  FS 3D - get cubes')
#endif

    call get_IBZwedge_faces(lattice%recbv, symmetries%nsym_used, symmetries%rotmat, symmetries%isym_used, nfaces, nvec, dscal, bounds)

    !initialize tetraeder indices arrays
    call init_cube2tetralines()


    !=====================================!
    !=== generate or read a cubesfile ====!
    !=====================================!
    write(filename,fmt_fn_ext) filename_cubesinfo, ext_formatted
    if(.not.file_present(trim(filename)).and.lrefine==1) stop 'cannot refine dataset: no cubesinfo present'
    if(file_present(trim(filename)))then

      call read_cubesfile(nCub3test, nmarked, imarked)
      if(any(nCub3 /= nCub3test).and.(myrank==master)) write(*,*) 'Warning!!! nCub3 inconsistent with cubesfile. Taking the values from cubesfile'
      nCub3 = nCub3test

      if(lrefine==1)then
        deallocate(imarked)
        call read_cubesrefine(nCub3test, nmarked, imarked)
        if(any(nCub3/=nCub3test)) stop 'nCub3 not consistent between cubesinfo and cubesrefine'
        if(myrank==master)then
          write(filename,'("cubes_refine_step=",I0,"_",A".vtp")') 0, 'read'
          call cubes2VTK(filename, nCub3, nmarked, imarked, bounds)
        end if!myrank==master
      end if!lrefine==1
!     write(*,*) 'readin of cubesfile done on rank', myrank

    else!cubesfile_present


      !search for existing cubesfiles that match the grid of an iteration within the iterative refinements
      iterstart=1
      !try to take the densest file, therefore count from top to bottom:
      iter_loop: do iter=nFSiter,1,-1

        !find the density of the cubes of iteration #iter
        nCub3test = nCub3
        do iter2=2,iter
          nCub3test = nCub3test*nCut_iter(iter2-1)
        end do!iter2

        !search for the existing cubesfile
        filetest = get_cubesinfo_filename(nCub3test,.true.)
        l_exist = file_present(trim(filetest))

        
        if(l_exist)then
          iterstart=iter+1
          call read_cubesfile(nCub3, nmarked, imarked, nCub3test)
          call cut_and_update_cubes(nCut_iter(iter), nCub3, nmarked, imarked)
          ntotal = nmarked
          if(myrank==master) then
            write(filename,'("cubes_iter=",I0,"_step=",I0,"_",(A),".vtp")') iter, 3, 'cut'
            call cubes2VTK(filename, nCub3, nmarked, imarked, bounds)
          end if!myrank==master

          !cubesfile found, so exit loop
          exit iter_loop
        end if!l_exist
      end do iter_loop


      if(iterstart==1)then
        !initialize cubes indices
        nmarked = product(nCub3)
        ntotal  = product(nCub3)
        allocate(imarked(nmarked), STAT=ierr)
        if(ierr/=0) stop 'Problem allocating imarked in fermisurface'
        do icube=1,nmarked
          imarked(icube) = icube
        end do!icube
      end if!iterstart==1

#ifdef CPP_TIMING
    call timing_stop('  FS 3D - get cubes')
#endif

      !======
      != scan the cubes for intersections
      !======
      do iter=iterstart,nFSiter

#ifdef CPP_TIMING
        call timing_start('  FS 3D - FS iter - mark cubes in IBZ')
#endif
        !=== mark the cubes that lie (at least with one corner) within the first BZ ===!
        call mark_cubes_in_IBZ(nCub3, nfaces, nvec, dscal, bounds, nmarked, imarked)
        if(myrank==master) then
          write(filename,'("cubes_iter=",I0,"_step=",I0,"_",(A),".vtp")') iter, 1, 'IBZ'
          call cubes2VTK(filename, nCub3, nmarked, imarked, bounds)
          write(*,'("***** Iteration ",I0," *****")') iter
          write(*,'(2X,2(A,I0),(A,F5.1,A))') 'Cubes found within the IBZ: ', nmarked, ' of ', ntotal, ' (= ', real(nmarked*100)/ntotal,' %)'
        end if!myrank==master
#ifdef CPP_TIMING
        call timing_stop('  FS 3D - FS iter - mark cubes in IBZ')
#endif

#ifdef CPP_TIMING
        call timing_start('  FS 3D - FS iter - find FS crossings')
#endif
        !========= mark the cubes that cross the Fermi surface =========!
        !=== (only searching across the four diagonals of the cubes) ===!
        ntotal = nmarked
        call mark_cubes_FScross( inc, lattice, cluster, tgmatrx, nCub3,          &
                               & nstepsconnect(iter), 8, 4, tetdiags, bounds,    &
                               & roottype(iter), rooteps(iter), nmarked, imarked )
        if(myrank==master) then
          write(filename,'("cubes_iter=",I0,"_step=",I0,"_",(A),".vtp")') iter, 2, 'mark'
          call cubes2VTK(filename, nCub3, nmarked, imarked, bounds)
          write(*,'(2X,2(A,I0),(A,F5.1,A))') 'Cubes found intersecting with FS: ', nmarked, ' of ', ntotal, ' (= ', real(nmarked*100)/ntotal,' %)'
          call save_cubesfile(nCub3, nmarked, imarked, lintermediate=.true.)
        end if!myrank==master
#ifdef CPP_TIMING
        call timing_stop('  FS 3D - FS iter - find FS crossings')
#endif

#ifdef CPP_TIMING
        call timing_start('  FS 3D - FS iter - divide cubes')
#endif
        !=== cut the remaining cubes into smaller pieces and update the indices ===!
!       write(1000+myrank,*) iter, nCut_iter(iter)
        call cut_and_update_cubes(nCut_iter(iter), nCub3, nmarked, imarked)
        ntotal = nmarked
        if(myrank==master) then
          write(filename,'("cubes_iter=",I0,"_step=",I0,"_",(A),".vtp")') iter, 3, 'cut'
          call cubes2VTK(filename, nCub3, nmarked, imarked, bounds)
        end if!myrank==master
#ifdef CPP_TIMING
        call timing_stop('  FS 3D - FS iter - divide cubes')
#endif

      end do!iter

      !=== mark the cubes that lie (maybe only partly) within the first BZ ===!
      iter = nFSiter+1
      call mark_cubes_in_IBZ(nCub3, nfaces, nvec, dscal, bounds, nmarked, imarked)
      if(myrank==master) then
        write(filename,'("cubes_iter=",I0,"_step=",I0,"_",(A),".vtp")') iter, 1, 'IBZ'
        call cubes2VTK(filename, nCub3, nmarked, imarked, bounds)
        call save_cubesfile(nCub3, nmarked, imarked)
      end if!myrank==master

    end if!cubesfile_present
    !=====================================!
    !=== generate or read a cubesfile ====!
    !=====================================!



    iter = nFSiter+1


    !=======================================!
    !=== refine the grid on a cubesfile ====!
    !=======================================!
    if(lrefine==1)then
      if(myrank==master) write(*,*) 'In REFINE mode...'

!      write(*,'("-->NREFINE= ",I0)') nrefinenew
      call cut_and_update_cubes(nrefinenew, nCub3, nmarked, imarked)
      if(myrank==master)then
        write(filename,'("cubes_refine_step=",I0,"_",A".vtp")') 1, 'cut'
        call cubes2VTK(filename, nCub3, nmarked, imarked, bounds)
      end if!myrank==master

      ntotal = nmarked
      call mark_cubes_in_IBZ(nCub3, nfaces, nvec, dscal, bounds, nmarked, imarked)
      if(myrank==master)then
        write(filename,'("cubes_refine_step=",I0,"_",A,".vtp")') 2, 'IBZ'
        call cubes2VTK(filename, nCub3, nmarked, imarked, bounds)
        write(*,'(2X,2(A,I0),(A,F5.1,A))') 'Cubes found within the IBZ: ', nmarked, ' of ', ntotal, ' (= ', real(nmarked*100)/ntotal,' %)'
      end if!myrank==master

      ntotal = nmarked
      call mark_cubes_FScross( inc, lattice, cluster, tgmatrx, nCub3,                &
                             & nstepsconnect(nFSiter), 8, 4, tetdiags, bounds,       &
                             & roottype(nFSiter), rooteps(nFSiter), nmarked, imarked )
      if(myrank==master) then
        write(filename,'("cubes_refine_step=",I0,"_",A,".vtp")') 3, 'mark'
        call cubes2VTK(filename, nCub3, nmarked, imarked, bounds)
        write(*,'(2X,2(A,I0),(A,F5.1,A))') 'Cubes found intersecting with FS: ', nmarked, ' of ', ntotal, ' (= ', real(nmarked*100)/ntotal,' %)'
      end if!myrank==master

    end if!lrefine==1
    !=======================================!
    !=== refine the grid on a cubesfile ====!
    !=======================================!

!   write(3360+myrank,*) iter, nstepsconnect(iter), nROOTiter
!   write(3360+myrank,*) roottype(iter), rooteps(iter)
!   write(3360+myrank,*) nfaces
!   write(3360+myrank,*) nvec
!   write(3360+myrank,*) dscal
!   dtmp=0d0
!   do ii=1,10000000
!     dtmp=dtmp+datan(ii*0.1d0)/ii
!   end do!ii

#ifdef CPP_TIMING
    call timing_start('  FS 3D - find FS intersections')
#endif

!   write(*,*) 'before sub find_intesection_triangles, myrank=',myrank
    call find_intesection_triangles( inc, lattice, cluster, tgmatrx, symmetries,  &
                                   & nCub3, bounds, nmarked, imarked,             &
                                   & nstepsconnect(iter), nROOTiter,              &
                                   & roottype(iter), rooteps(iter), nfaces,       &
                                   & nvec, dscal, nkpts_vis, nkpts_int,           &
                                   & kpoints_vis, kpoints_int, areas_int, vis2int )

    call find_kpoints_irredset( bounds, nkpts_vis, kpoints_vis, nkpts_irr, kpt2irr, irr2kpt)

#ifdef CPP_TIMING
    call timing_stop('  FS 3D - find FS intersections')
#endif

    !save the visualization k-points to a file
    if(myrank==master) call save_kpointsfile_vis(nkpts_vis, nkpts_irr, kpoints_vis, symmetries%nsym_used, symmetries%isym_used, kpt2irr, irr2kpt, vis2int)


  end subroutine find_fermisurface_3D





  subroutine find_intesection_triangles( inc, lattice, cluster, tgmatrx, symmetries, &
                                       & nCub3, bounds, nmarked, imarked,            &
                                       & nsteps, niter, roottype, rooteps,           &
                                       & nfaces, nvec, dscal, npoints_vis_tot,       &
                                       & npoints_int_tot, kpoints_vis_all,           &
                                       & kpoints_int_all, areas_int_all, vis2int_all )

    use type_inc,       only: inc_type
    use type_data,      only: lattice_type, cluster_type, tgmatrx_type
    use mod_symmetries, only: symmetries_type
    use mod_parutils,   only: distribute_linear_on_tasks
    use mod_mympi,      only: myrank, nranks, master
    use mod_mathtools,  only: bubblesort, findminindex
    use mod_vtkxml,     only: write_pointdata_rot
    use mod_ioformat,   only: fmt_fn_ext, fmt_fn_sub_ext, ext_vtkxml, filemode_rot, filename_vtktest, fmt_fn_rank_ext, filename_outinfo, ext_formatted
    use mod_fermisurf_basic, only: roots_along_edge, compare_two_eigv_in_substeps, ROOT_IMAG, ROOT_REAL, generate_cubevertices
#ifdef CPP_MPI
    use mpi
#endif
#ifdef CPP_TIMING
    use mod_timing,     only: timing_start, timing_stop
#endif
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    type(tgmatrx_type), intent(in) :: tgmatrx
    type(symmetries_type), intent(in) :: symmetries
    integer,            intent(in) :: nCub3(3), nmarked, imarked(nmarked),&
                                    & nsteps, roottype, niter, nfaces
    double precision,   intent(in) :: bounds(3,2), rooteps, nvec(3,nfaces), dscal(nfaces)

    integer,            intent(out) :: npoints_vis_tot, npoints_int_tot
    integer,          allocatable, intent(out) :: vis2int_all(:)
    double precision, allocatable, intent(out) :: kpoints_vis_all(:,:),&
                                                & kpoints_int_all(:,:),&
                                                & areas_int_all(:)

    integer :: irank, npoints_vis,&
             & npoints_vis_tmpa(0:nranks-1),&
             & npoints_int_tmpa(0:nranks-1),&
             & ioffs_tmpa(0:nranks-1),&
             & ioffs_save(0:nranks-1),&
             & isendarr(0:nranks-1),&
             & ioffsarr(0:nranks-1)

    double precision, allocatable :: kpoints_vis(:,:),&
                                   & kpoints_int(:,:),&
                                   & areas_int(:),    &
                                   & scalardata(:,:), &
                                   & vectordata(:,:,:)

    character(len=256) :: scalarstring(3), vectorstring(1)

    integer,          allocatable :: vis2int(:)
    double complex,   allocatable :: eigw_vis(:), eigw_vis_all(:)

    integer :: curid, nbands, nroots(19), lmid(inc%nrootmax,19), lmroot(inc%nrootmax,19)

    double precision :: kverts(3,8),&
                      & kends(3,2),&
                      & kmidcube(3),&
                      & kroot_in_cube(3,inc%nrootmax,19),&
                      & ktriangle(3,3),&
                      & kroot_lm_tet(3,nkpmax)   ! nkpmax >> (3 intersections per triangle) * (2 triangles per tetrahedron) * (6 tetrahedra) x
                                                 !           (multiplication-factor for splitting the triangle if BZ-face is crossed)
    double complex   :: eigw_lm_tet(nkpmax),&
                      & eigw_triangle(3)

    double complex, allocatable :: LVroot(:,:,:,:), &
                                 & RVroot(:,:,:,:), &
                                 & eigwroot(:,:,:)

    double precision :: ksub(3,nsteps+1)
    double complex   :: eigwends(2,inc%nrootmax,19)

    integer :: iproblems, icubeproblems(2,10*nmarked/nranks), iproblems_tmpa(0:nranks-1), iproblems_tot
    integer, allocatable :: icubeproblems_all(:,:)

    integer :: ioff, nkpt, ntot_pT(0:nranks-1), ioff_pT(0:nranks-1), lb, ub
    integer :: ierr, ii, ivc, icub, itmp, itet, iedge, itetedge, itetroot,&
             & icount, iedge1, iedge2, iroot1, iroot2, lm0, lm1, lm2,     &
             & nfound_band, nfound_tet, lookup_tet(2,4), iverts_picked(4),&
             & sorted(4), sorted_tmp(4), i4edge(4), i4tet(4), isave,      &
             & kcounter_lm, kcounter_cub, printstep, itmparr(1)
    double precision :: dtmp_ri(4), dist(nkpmax), weight_lm
    double complex   :: eigw_picked(4)
    logical :: match
    character(len=256) :: filename
    integer, parameter :: iofile=69821

!   write(*,*) 'in sub find_intesection_triangles, myrank=',myrank

    !Parallelize
    call distribute_linear_on_tasks(nranks,myrank,master,nmarked,ntot_pT,ioff_pT,.true.)
    nkpt = ntot_pT(myrank)
    ioff = ioff_pT(myrank)

    !initialize
    icubeproblems = 0
    iproblems = 0
    npoints_vis = 0

    !allocate arrays
    allocate( LVroot(inc%almso,inc%almso,inc%nrootmax,19),&
            & RVroot(inc%almso,inc%almso,inc%nrootmax,19),&
            & eigwroot(inc%almso,inc%nrootmax,19),        &
            & STAT=ierr )
    if(ierr/=0) stop 'Problem allocating LVroot etc. in find_intesection_triangles'

    allocate(kpoints_vis(3,200*nkpt), eigw_vis(200*nkpt), kpoints_int(3,20*nkpt), areas_int(20*nkpt), vis2int(200*nkpt), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating kpoints_vis etc.'

    !print header of statusbar
    if(myrank==master) write(*,'("Loop over cubes: |",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100

    printstep = nkpt/50
    if(printstep==0) printstep=1

#ifdef CPP_DEBUG
    write(filename,fmt_fn_rank_ext) filename_outinfo, myrank, ext_formatted
    open(unit=iofile, file=trim(filename), action='write', form='formatted')
#endif

#ifdef CPP_DEBUG
    if(myrank==master) call cubes2VTK('testcubes.vtp',nCub3,nmarked,imarked,bounds)
    write(iofile,'(A,4I8)') 'nCub3, nmarked=', nCub3, nmarked
    write(iofile,'(A,3ES18.9)') 'bounds(:,1)=', bounds(:,1)
    write(iofile,'(A,3ES18.9)') 'bounds(:,2)=', bounds(:,2)
    write(iofile,'(A,16I8)') 'ntot_pT', ntot_pT
    write(iofile,'(A,16I8)') 'ioff_pT', ioff_pT
    write(iofile,'(A,2I8)') 'nkpt, ioff=', nkpt, ioff
    write(iofile,'(A)') 'imarked:'
    write(iofile,'(20I8)') imarked
#endif



    !**************************************!
    !************* B E G I N **************!
    !**************************************!
    !*** (parallelized) loop over cubes ***!
    !**************************************!
    if(myrank==master) write(*,FMT=190) !beginning of statusbar
    kcounter_cub=0
    areas_int = 0d0
    do icub=1,nkpt

      !update statusbar
      if(mod(icub,printstep)==0 .and. myrank==master) write(*,FMT=200)
#ifdef CPP_TIMING
      if (icub==1) call timing_start('  FS 3D - find FS intersections - time1')
      if (icub==1) call timing_start('  FS 3D - find FS intersections - time10')
#endif

!     write(iofile,'(A,I0)') 'starting cube ', imarked(icub+ioff)

      call generate_cubevertices(nCub3,imarked(icub+ioff), bounds, kverts)
      !find middle point of the cube
      kmidcube = (kverts(:,1)+kverts(:,8))/2

#ifdef CPP_DEBUG
      write(iofile,'(2X,"kverts(:,i)= ",3ES18.9)') kverts
#endif

      !==================================!
      !=== find roots along the edges ===!
      !==================================!
      do iedge=1,19
        !pick k-points
        kends(:,1) = kverts(:,cubedges(1,iedge))
        kends(:,2) = kverts(:,cubedges(2,iedge))
        call roots_along_edge( inc, lattice, cluster, tgmatrx, nsteps, kends, niter,&
                             & roottype, rooteps, nroots(iedge), lmroot(:,iedge),   &
                             & kroot_in_cube(:,:,iedge), LVroot(:,:,:,iedge),       &
!                            & RVroot(:,:,:,iedge), eigwroot(:,:,iedge), 500+iedge, &
                             & RVroot(:,:,:,iedge), eigwroot(:,:,iedge), -1,        &
                             & eigwends(:,:,iedge)                                  )
      end do!iedge

#ifdef CPP_TIMING
      if (icub==1) call timing_stop('  FS 3D - find FS intersections - time1')
      if (icub==10) call timing_stop('  FS 3D - find FS intersections - time10')
#endif

#ifdef CPP_DEBUG
      write(iofile,'(2X,"#roots found on the edges:")')
      do iedge=1,19
        write(iofile,'(4X,"iedge= ",I0,", # = ",I0)') iedge, nroots(iedge)
      end do
#endif

      !======================================================!
      !=== determine which kroots belong to the same band ===!
      !======================================================!
      lmid=0
      ! the array lmid(iroot,iedge) will contain a unique id to which band a k-point belongs.
      curid=0
      do iedge1=1,19
        do iroot1=1,nroots(iedge1)
          !pick a first root
          lm1 = lmroot(iroot1,iedge1)
          if(lmid(iroot1,iedge1)>0) cycle !this root has already an ID

          curid = curid+1
          lmid(iroot1,iedge1)=curid
#ifdef CPP_DEBUG
          write(iofile,'(4X,"Set (iedge,iroot)=( ",I0,",",I0,") to curid= ",I0)') iedge1, iroot1, curid
#endif


!         write(*,'("Next cur-id= ",I0)'), curid
          !compare to all other roots on different edges which are not yet treated
          loop3: do iedge2=iedge1+1,19
            loop4: do iroot2=1,nroots(iedge2)
              if(lmid(iroot2,iedge2)>0) cycle
              lm2 = lmroot(iroot2,iedge2)
              kends(:,1) = kroot_in_cube(:,iroot1,iedge1)
              kends(:,2) = kroot_in_cube(:,iroot2,iedge2)
              call compare_two_eigv_in_substeps( inc, lattice, cluster, tgmatrx, nsteps,&
                                               & kends, lm1, LVroot(:,:,iroot1,iedge1), &
                                               & RVroot(:,lm2,iroot2,iedge2),           &
                                               & eigwroot(lm2,iroot2,iedge2), match     )

#ifdef CPP_DEBUG
              write(iofile,'(6X,"Compare (iedge,iroot)=( ",I0,",",I0,") with (iedge,iroot)=( ",I0,",",I0,"), matching=",L1)') iedge1, iroot1, iedge2, iroot2, match
#endif

              if(match) then
                lmid(iroot2,iedge2)=curid
                cycle loop3
              end if!match
            end do loop4!iroot2
          end do loop3!iedge2

        end do!iroot1
      end do!iedge1
      nbands = curid

#ifdef CPP_DEBUG
      write(iofile,'(2X,A)') 'LM-check'
      do iedge=1,19
          write(iofile,'(4X,"iedge= ",I0," - ",20I8)') iedge, lmid(1:nroots(iedge),iedge)
      end do!iedge

      write(iofile,'(2X,A,I0)') 'nbands= ', nbands
#endif

!     !delete multiple bands on the same edge
!     do lm0=1,nbands
!       do iedge1=1,19
!         itmp=0
!         do iroot1=1,nroots(iedge1)
!           if(lmid(iroot1,iedge1)==lm0) itmp=itmp+1
!           if(lmid(iroot1,iedge1)==lm0 .and. itmp>1) lmid(iroot1,iedge1)=0
!         end do
!       end do
!     end do!lm0

!     write(*,*) 'LM-check 2'
!     do iedge=1,19
!         write(*,'(4X,"iedge= ",I0," - ",20I8)') iedge, lmid(1:nroots(iedge),iedge)
!     end do!iedge


      !==========================================!
      !=== scan the tetrahedra for each band  ===!
      !===  and find the triangles and areas  ===!
      !==========================================!
#ifdef CPP_DEBUG
          write(iofile,'(2X,A)') 'Scan all bands ands tetraheda'
#endif
      !loop over all bands
      do lm0=1,nbands
        kcounter_lm=0

        !scan all tetrahedra
        do itet=1,6
#ifdef CPP_DEBUG
          write(iofile,'(4X,2(A,I0))') 'lm0= ', lm0, ', itet= ', itet
#endif

          nfound_tet = 0

          !find all roots belonging to this this tetrahedron and this band
          edgeloop: do iedge1=1,6
            itetedge = tetedges(iedge1,itet)
!           write(*,'(20(A,I0))') 'itetedge= ', itetedge, ', #roots=', nroots(itetedge)
            do iroot1=1,nroots(itetedge)
              if(lmid(iroot1,itetedge)==lm0)then
                nfound_tet = nfound_tet+1
                if(nfound_tet>4)then
#ifdef CPP_DEBUG
                  write(iofile,'(I0,A)') imarked(icub+ioff), ' =imarked: more than 4 intersections of band with tetrahedron --> skip'
#endif
                  iproblems = iproblems+1
                  icubeproblems(1,iproblems) = imarked(icub+ioff)
                  icubeproblems(2,iproblems) = 1
                  nfound_tet=0
                  exit edgeloop
                end if
                lookup_tet(:,nfound_tet) = (/ iroot1, itetedge /)
              end if
            end do!iroot1
          end do edgeloop!iedge1

#ifdef CPP_DEBUG
          write(iofile,'(6X,A,I0)') 'nfound_tet= ', nfound_tet
#endif

          select case ( nfound_tet )
            case( 0 ); cycle
            case( 3 )

              !store the triangle k-points
              ktriangle(:,1) = kroot_in_cube(:,lookup_tet(1,1),lookup_tet(2,1))
              ktriangle(:,2) = kroot_in_cube(:,lookup_tet(1,2),lookup_tet(2,2))
              ktriangle(:,3) = kroot_in_cube(:,lookup_tet(1,3),lookup_tet(2,3))
              eigw_triangle(1) = eigwroot(lmroot(lookup_tet(1,1),lookup_tet(2,1)),lookup_tet(1,1),lookup_tet(2,1))
              eigw_triangle(2) = eigwroot(lmroot(lookup_tet(1,2),lookup_tet(2,2)),lookup_tet(1,2),lookup_tet(2,2))
              eigw_triangle(3) = eigwroot(lmroot(lookup_tet(1,3),lookup_tet(2,3)),lookup_tet(1,3),lookup_tet(2,3))
              !split triangle if BZ-face is crossed
              call split_triangle(ktriangle, eigw_triangle, kroot_lm_tet(:,:), eigw_lm_tet(:), kcounter_lm, nfaces, nvec, dscal)
!cccccccc     kroot_lm_tet(:,kcounter_lm+1:kcounter_lm+3) = ktriangle
!cccccccc     kcounter_lm = kcounter_lm+3

            case( 4 )

              !pick the eigenvalues at the four vertices of the tetrahedron
              icount = 0
              iverts_picked = 0
              eigw_picked = (0d0, 0d0)
              do iedge=1,4
                itetroot = lookup_tet(1,iedge)
                itetedge = lookup_tet(2,iedge)
                do ii=1,2
                  ivc = cubedges(ii,itetedge)
                  if(any(iverts_picked==ivc)) cycle
                  icount = icount+1
                  iverts_picked(icount) = ivc
                  eigw_picked(icount) = eigwends(ii,itetroot,itetedge)
                end do!ii
              end do
              if(icount>4) stop 'icount>4 sould not happen'

#ifdef CPP_DEBUG
              write(iofile,'(6X,"iverts_picked= ",4I8)') iverts_picked
              write(iofile,'(8X,"eigws_picked= ",2ES18.9)') eigw_picked

              write(iofile,'(6X,A)') 'Test-eigws:'
              do iedge=1,4
                itetroot = lookup_tet(1,iedge)
                itetedge = lookup_tet(2,iedge)
                do ii=1,2
                  write(iofile,'(8X,A,I0,A,I0,A,2ES18.9)') 'itetedge=',itetedge,', icorner=', cubedges(ii,itetedge),&
                                                         & ', eigw= ', eigwends(ii,itetroot,itetedge)
                end do!ii
              end do!iedge
#endif

              select case (roottype)
                case(ROOT_REAL); dtmp_ri = dble(eigw_picked)
                case(ROOT_IMAG); dtmp_ri = aimag(eigw_picked)
                case default; stop 'option not known: only real/imag'
              end select!roottype

              !sort the values
              call bubblesort(4, dtmp_ri, sorted_tmp)
              do ii=1,4
                sorted(ii) = iverts_picked(sorted_tmp(ii))
              end do
              i4edge(1) = get_edgeindex(sorted(1),sorted(3))
              i4edge(2) = get_edgeindex(sorted(1),sorted(4))
              i4edge(3) = get_edgeindex(sorted(2),sorted(3))
              i4edge(4) = get_edgeindex(sorted(2),sorted(4))

              !reorder the indices
              icount=0
              i4tet =0
              do ii=1,4
                do iedge=1,4
                  if(i4edge(ii)==lookup_tet(2,iedge))then
                    icount = icount+1
                    i4tet(icount) = iedge
                  end if
                end do!iedge
              end do!ii

              if(any(i4tet==0) .or. any(i4tet>4))then
#ifdef CPP_DEBUG
                write(iofile,'(I0,A)') imarked(icub+ioff), ' =imarked: it does not hold (1 <= i4tet(ii) <= 4) --> skip'
#endif
                iproblems = iproblems+1
                icubeproblems(1,iproblems) = imarked(icub+ioff)
                icubeproblems(2,iproblems) = 2
                cycle
              end if

              !store the first triangle k-points
              ktriangle(:,1) = kroot_in_cube(:,lookup_tet(1,i4tet(1)),lookup_tet(2,i4tet(1)))
              ktriangle(:,2) = kroot_in_cube(:,lookup_tet(1,i4tet(2)),lookup_tet(2,i4tet(2)))
              ktriangle(:,3) = kroot_in_cube(:,lookup_tet(1,i4tet(3)),lookup_tet(2,i4tet(3)))
              eigw_triangle(1) = eigwroot(lmroot(lookup_tet(1,1),lookup_tet(2,1)),lookup_tet(1,1),lookup_tet(2,1))
              eigw_triangle(2) = eigwroot(lmroot(lookup_tet(1,2),lookup_tet(2,2)),lookup_tet(1,2),lookup_tet(2,2))
              eigw_triangle(3) = eigwroot(lmroot(lookup_tet(1,3),lookup_tet(2,3)),lookup_tet(1,3),lookup_tet(2,3))
              !split triangle if BZ-face is crossed
              call split_triangle(ktriangle, eigw_triangle, kroot_lm_tet(:,:), eigw_lm_tet(:), kcounter_lm, nfaces, nvec, dscal)
!cccccccc     kroot_lm_tet(:,kcounter_lm+1:kcounter_lm+3) = ktriangle
!cccccccc     kcounter_lm = kcounter_lm+3

              !store the second triangle k-points
              ktriangle(:,1) = kroot_in_cube(:,lookup_tet(1,i4tet(2)),lookup_tet(2,i4tet(2)))
              ktriangle(:,2) = kroot_in_cube(:,lookup_tet(1,i4tet(3)),lookup_tet(2,i4tet(3)))
              ktriangle(:,3) = kroot_in_cube(:,lookup_tet(1,i4tet(4)),lookup_tet(2,i4tet(4)))
              eigw_triangle(1) = eigwroot(lmroot(lookup_tet(1,2),lookup_tet(2,2)),lookup_tet(1,2),lookup_tet(2,2))
              eigw_triangle(2) = eigwroot(lmroot(lookup_tet(1,3),lookup_tet(2,3)),lookup_tet(1,3),lookup_tet(2,3))
              eigw_triangle(3) = eigwroot(lmroot(lookup_tet(1,4),lookup_tet(2,4)),lookup_tet(1,4),lookup_tet(2,4))
              !split triangle if BZ-face is crossed
              call split_triangle(ktriangle, eigw_triangle, kroot_lm_tet(:,:), eigw_lm_tet(:), kcounter_lm, nfaces, nvec, dscal)
!cccccccc     kroot_lm_tet(:,kcounter_lm+1:kcounter_lm+3) = ktriangle
!cccccccc     kcounter_lm = kcounter_lm+3

            case default
#ifdef CPP_DEBUG
              write(iofile,'(I0,A)') imarked(icub+ioff), ' =imarked: neither 3 nor 4 intesections found for --> skip'
#endif
              iproblems = iproblems+1
              icubeproblems(1,iproblems) = imarked(icub+ioff)
              icubeproblems(2,iproblems) = 3
          end select

        end do!itet

        !if cube is empty, proceed to next cube
        if(kcounter_lm==0) cycle


        !==============================!
        !===  Find the k-point and  ===!
        !=== weight for integration ===!
        !==============================!

        !update the counter for integration k-points
        kcounter_cub = kcounter_cub+1

        !find kpoint-index representing this cube
        do ii=1,kcounter_lm
          dist(ii) = sum((kroot_lm_tet(:,ii) - kmidcube)**2)
        end do
        call findminindex(kcounter_lm, dist(1:kcounter_lm), isave)

        !save the k-point
        kpoints_int(:,kcounter_cub) = kroot_lm_tet(:,isave)

        !find the weight of the representing cube (=area of all triangles)
        weight_lm = 0d0
        do ii=1,kcounter_lm/3
          weight_lm = weight_lm + area_triangle(kroot_lm_tet(:,ii*3-2:ii*3))
        end do!ii
        areas_int(kcounter_cub) = weight_lm

        !if weight for this triangle is negligible, forget this k-point
        if(abs(weight_lm)<1d-16) kcounter_cub = kcounter_cub - 1



        !=================================!
        !===  Store the visualization  ===!
        !=== information in temp array ===!
        !=================================!

        !add the visualization-triangles to the temporary storage array
        kpoints_vis(:,npoints_vis+1:npoints_vis+kcounter_lm) = kroot_lm_tet(:,1:kcounter_lm)
        eigw_vis(npoints_vis+1:npoints_vis+kcounter_lm) = eigw_lm_tet(1:kcounter_lm)
        vis2int(npoints_vis+1:npoints_vis+kcounter_lm) = kcounter_cub
        npoints_vis = npoints_vis + kcounter_lm

      end do!lm0

    end do!icub

    if(myrank==master) write(*,*) ''
    !**************************************!
    !*** (parallelized) loop over cubes ***!
    !**************************************!
    !*************** E N D ****************!
    !**************************************!




    !*************************************!
    !************* B E G I N *************!
    !*************************************!
    !*** collect the problematic cubes ***!
    !*************************************!
#ifdef CPP_MPI
    !collect the numbers of problematic arrays
    call MPI_Allgather( iproblems,      1, MPI_INTEGER, &
                      & iproblems_tmpa, 1, MPI_INTEGER, &
                      & MPI_COMM_WORLD, ierr              )
    if(ierr/=MPI_SUCCESS) stop 'Problem in allgather npoints_vis_glob'
    iproblems_tot = sum(iproblems_tmpa)
    allocate(icubeproblems_all(2,iproblems_tot), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating icubeproblems_all'

    !calculate offsets
    ioffs_tmpa=0
    do irank=1,nranks-1
      ioffs_tmpa(irank) = ioffs_tmpa(irank-1)+iproblems_tmpa(irank-1)
    end do!irank
    ioffs_save = ioffs_tmpa

    !collect the problematic cubes
    isendarr = 2*iproblems_tmpa
    ioffsarr = 2*ioffs_tmpa
    call MPI_Allgatherv( icubeproblems, isendarr(myrank), MPI_INTEGER, &
                       & icubeproblems_all, isendarr, ioffsarr,        &
                       & MPI_INTEGER, MPI_COMM_WORLD, ierr               )
    if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv icubeproblems'
#else
    iproblems_tot = iproblems
    allocate(icubeproblems_all(2,iproblems_tot), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating icubeproblems_all'
    icubeproblems_all = icubeproblems(:,1:iproblems_tot)
#endif

#ifndef CPP_DEBUG
    if(myrank==master)then
      write(filename,fmt_fn_rank_ext) filename_outinfo, myrank, ext_formatted
      open(unit=13512, file=trim(filename), action='write', form='formatted')
      do ii=1,iproblems_tot
        select case (icubeproblems_all(2,ii))
          case(1); write(13512,'(I0,A)') icubeproblems_all(1,ii), ' =imarked: more than 4 intersections of band with tetrahedron --> skip'
          case(2); write(13512,'(I0,A)') icubeproblems_all(1,ii), ' =imarked: it does not hold (1 <= i4tet(ii) <= 4) --> skip'
          case(3); write(13512,'(I0,A)') icubeproblems_all(1,ii), ' =imarked: neither 3 nor 4 intesections found for --> skip'
          case default; write(13512,'(I0,A,I0,A)') icubeproblems_all(1,ii), ' =imarked: identifier ', icubeproblems_all(2,ii), ' not known'
        end select
      end do!ii
      close(13512)
    end if!myrank==master
#endif
    !*************************************!
    !*** collect the problematic cubes ***!
    !*************************************!
    !************** E N D ****************!
    !*************************************!





    !************************************!
    !************ B E G I N *************!
    !************************************!
    !*** output of visualization data ***!
    !************************************!

#ifdef CPP_MPI
    !collect the number of points
    call MPI_Allgather( npoints_vis,      1, MPI_INTEGER, &
                      & npoints_vis_tmpa, 1, MPI_INTEGER, &
                      & MPI_COMM_WORLD, ierr              )
    if(ierr/=MPI_SUCCESS) stop 'Problem in allgather npoints_vis_glob'
    npoints_vis_tot = sum(npoints_vis_tmpa)

    ioffs_tmpa=0
    do irank=1,nranks-1
      ioffs_tmpa(irank) = ioffs_tmpa(irank-1)+npoints_vis_tmpa(irank-1)
    end do!irank
    ioffs_save = ioffs_tmpa

    !allocate result arrays
    allocate(kpoints_vis_all(3,npoints_vis_tot), eigw_vis_all(npoints_vis_tot), vis2int_all(npoints_vis_tot), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating kpoints_vis_all etc.'

    !collect the k-points
    isendarr = 3*npoints_vis_tmpa
    ioffsarr = 3*ioffs_tmpa
    call MPI_Allgatherv( kpoints_vis, isendarr(myrank), MPI_DOUBLE_PRECISION,&
                       & kpoints_vis_all, isendarr, ioffsarr,                &
                       & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr          )
    if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv kpoints_vis_glob'

    !collect the eigenvalues
    isendarr = npoints_vis_tmpa
    ioffsarr = ioffs_tmpa
    call MPI_Allgatherv( eigw_vis, isendarr(myrank), MPI_DOUBLE_COMPLEX,&
                       & eigw_vis_all, isendarr, ioffsarr,              &
                       & MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr       )
    if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv eigv_vis_glob'

    !collect the translation array
    isendarr = npoints_vis_tmpa
    ioffsarr = ioffs_tmpa
    call MPI_Allgatherv( vis2int, isendarr(myrank), MPI_INTEGER, &
                       & vis2int_all, isendarr, ioffsarr,        &
                       & MPI_INTEGER, MPI_COMM_WORLD, ierr       )
    if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv vis2int_glob'
    !adapt the offsets for the translation array --> done below (!!!!!)
#else
    npoints_vis_tot = npoints_vis
    allocate(kpoints_vis_all(3,npoints_vis_tot), eigw_vis_all(npoints_vis_tot), vis2int_all(npoints_vis_tot), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating kpoints_vis_all etc.'
    kpoints_vis_all = kpoints_vis(:,1:npoints_vis_tot)
    eigw_vis_all = eigw_vis(1:npoints_vis_tot)
    vis2int_all = vis2int(1:npoints_vis_tot)
#endif

    deallocate( kpoints_vis, eigw_vis, STAT=ierr )
    if(ierr/=0) stop 'Problem deallocating kpoints_vis etc.'


    if(myrank==master)then

      allocate(scalardata(npoints_vis_tot,3))
      scalardata(:,1) = abs(eigw_vis_all)
      scalardata(:,2) = dble(eigw_vis_all)
      scalardata(:,3) = aimag(eigw_vis_all)
      scalarstring(1) = 'eigw_abs'
      scalarstring(2) = 'eigw_real'
      scalarstring(3) = 'eigw_imag'

      write(filename,fmt_fn_ext) filename_vtktest, ext_vtkxml
      itmparr(1) = 1
!     write(*,*) 'prepare fstest with nsym=1.. done'
      call write_pointdata_rot( trim(filename),npoints_vis_tot,kpoints_vis_all, &
                              & 3,scalardata,scalarstring,   &
                              & 0,vectordata,vectorstring,   &
                              & 1,symmetries%rotmat, itmparr )

!     write(*,*) 'write fstest with nsym=1.. done'

      write(filename,fmt_fn_sub_ext) filename_vtktest, filemode_rot, ext_vtkxml
      call write_pointdata_rot( trim(filename),npoints_vis_tot,kpoints_vis_all, &
                              & 3,scalardata,scalarstring,&
                              & 0,vectordata,vectorstring,&
                              & symmetries%nsym_used,symmetries%rotmat, &
                              & symmetries%isym_used)
!     write(*,*) 'write fstest with nsym=48.. done'

    end if!myrank==master
    !************************************!
    !*** output of visualization data ***!
    !************************************!
    !************** E N D ***************!
    !************************************!



    !*************************************!
    !************ B E G I N **************!
    !*************************************!
    !*** gathering of integration data ***!
    !*************************************!
#ifdef CPP_MPI
    !collect the number of points
    call MPI_Allgather( kcounter_cub,     1, MPI_INTEGER, &
                      & npoints_int_tmpa, 1, MPI_INTEGER, &
                      & MPI_COMM_WORLD, ierr              )
    if(ierr/=MPI_SUCCESS) stop 'Problem in allgather npoints_int_tmpa'
    npoints_int_tot = sum(npoints_int_tmpa)

    ioffs_tmpa=0
    do irank=1,nranks-1
      ioffs_tmpa(irank) = ioffs_tmpa(irank-1)+npoints_int_tmpa(irank-1)
    end do!irank

    !allocate result arrays
    allocate(kpoints_int_all(3,npoints_int_tot), areas_int_all(npoints_int_tot), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating kpoints_int_all etc.'

    !collect the k-points
    isendarr = 3*npoints_int_tmpa
    ioffsarr = 3*ioffs_tmpa
    call MPI_Allgatherv( kpoints_int, isendarr(myrank), MPI_DOUBLE_PRECISION,&
                       & kpoints_int_all, isendarr, ioffsarr,                &
                       & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr          )
    if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv kpoints_int_all'

    !collect the areas
    isendarr = npoints_int_tmpa
    ioffsarr = ioffs_tmpa
    call MPI_Allgatherv( areas_int, isendarr(myrank), MPI_DOUBLE_PRECISION,&
                       & areas_int_all, isendarr, ioffsarr,                &
                       & MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr        )
    if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv kpoints_int_all'

    !adapt the offsets for the translation array --> from above (!!!!!)
    do irank=1,nranks-1
      lb = ioffs_save(irank)+1
      ub = ioffs_save(irank)+npoints_vis_tmpa(irank)
      vis2int_all(lb:ub) = vis2int_all(lb:ub) + ioffs_tmpa(irank)
    end do!irank
#else
    npoints_int_tot = kcounter_cub
    allocate(kpoints_int_all(3,npoints_int_tot), areas_int_all(npoints_int_tot), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating kpoints_int_all etc.'
    kpoints_int_all = kpoints_int(:,1:npoints_int_tot)
    areas_int_all = areas_int(1:npoints_int_tot)
#endif


!   write(iofile,'(2(A,ES18.9))') 'Test: sum of areas on all tasks=', sum(areas_int_all)
!   write(iofile,'(2(A,ES18.9))') 'Test: sum of areas on this task=', sum(areas_int)

    close(iofile)

190 FORMAT('                 |'$)
200 FORMAT('|'$)

  end subroutine find_intesection_triangles





  subroutine cut_and_update_cubes(ncut, nCub3, nmarked, imarked)
    use mod_fermisurf_basic, only: unroll_ixyz
    implicit none

    integer,              intent(in)    :: ncut
    integer,              intent(inout) :: nCub3(3)
    integer,              intent(inout) :: nmarked
    integer, allocatable, intent(inout) :: imarked(:)

    integer :: nCub3in(3), nmarked_tmp, imarked_tmp(nmarked)
    integer :: ii3(3), ix, iy, iz, icub, icount, ierr

    !save input to temporary arrays
    nCub3in     = nCub3
    nmarked_tmp = nmarked
    imarked_tmp = imarked
    deallocate(imarked)

    !determine size of output-array
    nmarked = nmarked_tmp * ncut**3
    allocate( imarked(nmarked), STAT=ierr )
    if(ierr/=0) stop 'Problem allocating imarked in "cut_and_update_cubes"'
    nCub3 = nCub3in*ncut

    icount = 0
    do icub=1,nmarked_tmp
      call unroll_ixyz(imarked_tmp(icub), nCub3in, ii3)
      do iz=0,ncut-1
       do iy=0,ncut-1
        do ix=0,ncut-1
          icount = icount+1
          imarked(icount) = (ii3(1)*ncut+ix+1) + (ii3(2)*ncut+iy)*nCub3(1) + (ii3(3)*ncut+iz)*product(nCub3(1:2))
        end do!ix
       end do!iy
      end do!iz
    end do!icub

  end subroutine cut_and_update_cubes





  subroutine mark_cubes_FScrossold(inc, lattice, cluster, tgmatrx, nCub3, nsteps, bounds, roottype, rooteps, nmarked, imarked)

    use type_inc,       only: inc_type
    use type_data,      only: lattice_type, cluster_type, tgmatrx_type
    use mod_parutils,   only: distribute_linear_on_tasks
    use mod_mympi,      only: myrank, nranks, master
    use mod_fermisurf_basic, only: connect_eigw_in_substeps, find_roots_any_eigw, generate_cubevertices
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    type(tgmatrx_type), intent(in) :: tgmatrx
    integer,            intent(in) :: nCub3(3), nsteps
    double precision,   intent(in) :: bounds(3,2)
    integer,            intent(in) :: roottype
    double precision,   intent(in) :: rooteps
    integer,              intent(inout) :: nmarked
    integer, allocatable, intent(inout) :: imarked(:)

    integer :: ncubmarked_tmp, ioff, nkpt, ntot_pT(0:nranks-1), ioff_pT(0:nranks-1)
    integer :: ncubmarked_pT(0:nranks-1), iioff(0:nranks-1)
    double precision :: kverts(3,8)
    integer, allocatable :: icubmarked_tmp(:)

    integer          :: connection(inc%almso,nsteps+1)
    double precision :: ksub(3,nsteps+1), kends(3,2)
    double complex   :: LVeig(inc%almso,inc%almso,nsteps+1),&
                      & RVeig(inc%almso,inc%almso,nsteps+1),&
                      & eigw(inc%almso,nsteps+1)

    integer :: iedge, irank, icub, ierr
    logical :: edgeroot(19)

    !Parallelize
    call distribute_linear_on_tasks(nranks,myrank,master, nmarked, ntot_pT, ioff_pT, .false.)
    nkpt = ntot_pT(myrank)
    ioff = ioff_pT(myrank)

    !Create temp arrays
    allocate(icubmarked_tmp(nkpt), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating icubmarked_tmp in mark_cubes_FScross'
    ncubmarked_tmp = 0
    icubmarked_tmp = -1

    !Test whether cubes cross FS
    do icub=1,nkpt
      if(mod(icub,10)==0 .and. myrank==master) write(*,"(2X,F8.3,A)") dble(icub)/nkpt*100, ' percent done'
      call generate_cubevertices(nCub3,imarked(icub+ioff), bounds, kverts)
      do iedge=1,4
        kends(:,1) = kverts(:,tetdiags(1,iedge))
        kends(:,2) = kverts(:,tetdiags(2,iedge))
        call connect_eigw_in_substeps( inc, lattice, cluster, tgmatrx, nsteps, kends, &
                                     & connection, ksub, eigw, LVeig, RVeig           )
        call find_roots_any_eigw( inc%almso, nsteps, connection, eigw, &
                                & roottype, rooteps, edgeroot(iedge)   )
      end do!iedge
      if(any(edgeroot))then
        ncubmarked_tmp = ncubmarked_tmp + 1
        icubmarked_tmp(ncubmarked_tmp) = imarked(icub+ioff)
      end if!l_cut
    end do!icub

    deallocate(imarked)

#ifdef CPP_MPI
    !Combine results
    call MPI_Allgather(ncubmarked_tmp, 1, MPI_INTEGER, ncubmarked_pT, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
    if(ierr/=MPI_SUCCESS) stop 'Problem in Allgather(ncubmarked_tmp) in mark_cubes_FScross'

    nmarked = sum(ncubmarked_pT)
    allocate(imarked(nmarked), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating imarked in mark_cubes_FScross'

    iioff = 0
    do irank=1,nranks-1
      iioff(irank) = sum(ncubmarked_pT(0:irank-1))
    end do

    call MPI_Allgatherv( icubmarked_tmp, ncubmarked_tmp, MPI_INTEGER, &
                       & imarked, ncubmarked_pT, iioff, MPI_INTEGER,  &
                       & MPI_COMM_WORLD, ierr )
    if(ierr/=MPI_SUCCESS) stop 'Problem in Allgatherv(icubmarked_tmp) in mark_cubes_FScross'

#else
    nmarked=ncubmarked_tmp
    allocate(imarked(nmarked), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating imarked in mark_cubes_FScross'
    imarked = icubmarked_tmp(1:ncubmarked_tmp)
#endif

  end subroutine mark_cubes_FScrossold





  subroutine mark_cubes_in_IBZ(nCub3, nfaces, nvec, dscal, bounds, nmarked, imarked)
    use mod_parutils,   only: distribute_linear_on_tasks
    use mod_mympi,      only: myrank, nranks, master
    use mod_symmetries, only: points_in_wedge
    use mod_fermisurf_basic, only: generate_cubevertices
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    integer,              intent(in)    :: nCub3(3), nfaces
    double precision,     intent(in)    :: nvec(3,nfaces), dscal(nfaces), bounds(3,2)
    integer,              intent(inout) :: nmarked
    integer, allocatable, intent(inout) :: imarked(:)

    integer :: ncubmarked_tmp, ioff, nkpt, ntot_pT(0:nranks-1), ioff_pT(0:nranks-1)
    integer :: ncubmarked_pT(0:nranks-1), iioff(0:nranks-1)
    double precision :: kverts(3,8)
    integer, allocatable :: icubmarked_tmp(:)

    integer :: irank, icub, ierr
    logical :: cubeinwedge


    !Parallelize
    call distribute_linear_on_tasks(nranks,myrank,master, nmarked, ntot_pT, ioff_pT, .false.)
    nkpt = ntot_pT(myrank)
    ioff = ioff_pT(myrank)

    !Create temp arrays
    allocate(icubmarked_tmp(nkpt), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating icubmarked_tmp in mark_cubes_in_IBZ'
    ncubmarked_tmp = 0
    icubmarked_tmp = -1

    !Test whether cubes lie in wedge
    do icub=1,nkpt
      call generate_cubevertices(nCub3,imarked(icub+ioff),bounds,kverts)
      cubeinwedge = points_in_wedge(nfaces, nvec, dscal, 8, kverts, 'any')
      if(cubeinwedge)then
        ncubmarked_tmp = ncubmarked_tmp + 1
        icubmarked_tmp(ncubmarked_tmp) = imarked(icub+ioff)
      end if!cubeinwedge
    end do!icub

    deallocate(imarked)

#ifdef CPP_MPI
    !Combine results
    call MPI_Allgather(ncubmarked_tmp, 1, MPI_INTEGER, ncubmarked_pT, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
    if(ierr/=MPI_SUCCESS) stop 'Problem in Allgather(ncubmarked_tmp) in mark_cubes_in_IBZ'

    nmarked = sum(ncubmarked_pT)
    allocate(imarked(nmarked), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating imarked in mark_cubes_in_IBZ'

    iioff = 0
    do irank=1,nranks-1
      iioff(irank) = sum(ncubmarked_pT(0:irank-1))
    end do


    call MPI_Allgatherv( icubmarked_tmp, ncubmarked_tmp, MPI_INTEGER, &
                       & imarked, ncubmarked_pT, iioff, MPI_INTEGER,  &
                       & MPI_COMM_WORLD, ierr )
    if(ierr/=MPI_SUCCESS) stop 'Problem in Allgatherv(icubmarked_tmp) in mark_cubes_in_IBZ'

#else
    nmarked=ncubmarked_tmp
    allocate(imarked(nmarked), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating imarked in mark_cubes_in_IBZ'
    imarked = icubmarked_tmp(1:ncubmarked_tmp)
#endif

  end subroutine mark_cubes_in_IBZ








  subroutine cubes2VTK(filename,nCub3,nmarked,imarked,bounds)
    use mod_fermisurf_basic, only: generate_cubevertices
    use mod_vtkxml
    implicit none

    character(len=*), intent(in) :: filename
    integer,          intent(in) :: nCub3(3), nmarked, imarked(nmarked)
    double precision, intent(in) :: bounds(3,2)

    integer          :: faces(4,6), offsets(6)
    double precision :: kverts(3,8)
    integer :: icub, ivert, iface
    character(len=80) :: fmtstr

    open(unit=123894, file=trim(filename), action='write', form='formatted')

    !write header
    write(123894, FMT=vtkfmt_ivtkfile)
    write(123894, FMT=vtkfmt_ipolydata)
    write(123894, FMT=vtkfmt_ipiece) 8*nmarked, 6*nmarked

    !write points
    write(123894, FMT=vtkfmt_ipoints)
    write(123894, FMT=vtkfmt_idata_points)
    write(fmtstr,'("(",I0,"X,3ES14.5)")') vtkfmxXdata
    do icub=1,nmarked
      call generate_cubevertices(nCub3,imarked(icub),bounds,kverts)
      do ivert=1,8
        write(123894, FMT=fmtstr) kverts(:,ivert)
      end do!ivert
    end do!icub
    write(123894, FMT=vtkfmt_fdata)
    write(123894, FMT=vtkfmt_fpoints)

    faces(:,1) = (/ 0,1,3,2 /)
    faces(:,2) = (/ 0,1,5,4 /)
    faces(:,3) = (/ 0,2,6,4 /)
    faces(:,4) = (/ 1,3,7,5 /)
    faces(:,5) = (/ 2,3,7,6 /)
    faces(:,6) = (/ 4,5,7,6 /)

    offsets = (/ 4, 8, 12, 16, 20, 24 /)

    !write polys
    write(123894, FMT=vtkfmt_ipolys)
    write(123894, FMT=vtkfmt_idata_connectivity)
    do icub=1,nmarked
      do iface=1,6
        write(123894, FMT='(10X,4I8)' ) faces(:,iface)+(icub-1)*8
      end do!iface
    end do!icub
    write(123894, FMT=vtkfmt_fdata)
    write(123894, FMT=vtkfmt_idata_offsets)
    write(fmtstr, '("(",I0,"X,6(I0,X))")') vtkfmxXdata
    do icub=1,nmarked
      write(123894, FMT=fmtstr) offsets + (icub-1)*24
    end do!icub
    write(123894, FMT=vtkfmt_fdata)
    write(123894, FMT=vtkfmt_fpolys)

    !write tail
    write(123894, FMT=vtkfmt_fpiece)
    write(123894, FMT=vtkfmt_fpolydata)
    write(123894, FMT=vtkfmt_fvtkfile)

    close(123894)

  end subroutine cubes2VTK





  subroutine init_cube2tetralines()

    ! A cube is cut into 6 tetrahedra.
    !  - itetcorner(j,i)   is the corner-index of the cube of the jth corner of the ith tetrahedron.
    !  - lines(1,l) and lines(2,l) contains the cube-index of the
    !     start- and end-points of the l-th edge of the tetrahedron in this triangle.
    !  - itetlines(l,i) is the lth edge of the ith tetrahedron

    integer :: itet

    tetdiags = reshape( (/ 1, 8, 2, 7, 3, 6, 4, 5 /), (/ 2,4 /))

!   do ii=1,4
!     write(*,*) 't', tetdiags(:,ii)
!   end do!ii

    !define the endpoints of the edges of the tetraeda. cubedge(1,j) = k means:
    !the startpoint of the j-th edge is the k-th corner-point of the cube 

    cubedges(:,1) = (/ 1, 2 /) !\
    cubedges(:,2) = (/ 1, 3 /) ! \ = the 4 cube-edges in the base
    cubedges(:,3) = (/ 2, 4 /) ! /
    cubedges(:,4) = (/ 3, 4 /) !/

    cubedges(:,5) = (/ 1, 5 /) !\
    cubedges(:,6) = (/ 2, 6 /) ! \ = the 4 cube-edges from the base to the top
    cubedges(:,7) = (/ 3, 7 /) ! /
    cubedges(:,8) = (/ 4, 8 /) !/

    cubedges(:,9)  = (/ 5, 6 /) !\
    cubedges(:,10) = (/ 5, 7 /) ! \ = the 4 cube-edges in the top
    cubedges(:,11) = (/ 6, 8 /) ! /
    cubedges(:,12) = (/ 7, 8 /) !/

    cubedges(:,13) = (/ 1, 4 /) !\
    cubedges(:,14) = (/ 1, 6 /) ! \
    cubedges(:,15) = (/ 1, 7 /) !  \ = the diagonals along the 6 faces
    cubedges(:,16) = (/ 2, 8 /) !  /
    cubedges(:,17) = (/ 3, 8 /) ! /
    cubedges(:,18) = (/ 5, 8 /) !/

    cubedges(:,19) = (/ 1, 8 /) ! = diagonal across the cube

    tetcorners(:,1) = (/ 1, 2, 4, 8 /)
    tetcorners(:,2) = (/ 1, 3, 4, 8 /)
    tetcorners(:,3) = (/ 1, 2, 6, 8 /)
    tetcorners(:,4) = (/ 1, 5, 6, 8 /)
    tetcorners(:,5) = (/ 1, 3, 7, 8 /)
    tetcorners(:,6) = (/ 1, 5, 7, 8 /)

    do itet=1,6
      tetedges(:,itet) = get_tedges(tetcorners(:,itet))
    end do

  end subroutine init_cube2tetralines





  function get_tedges(tcorners) result(tedges)

    integer, intent(in) :: tcorners(4)
    integer :: tedges(6)

    integer :: ii, jj, icount, iedge, pair(2)

    !loop over pairs
    icount=0
    iloop: do ii=1,3
      jloop: do jj=ii+1,4
        icount = icount+1

        !order the pair
        if(tcorners(ii)<tcorners(jj))then
          pair = (/ tcorners(ii), tcorners(jj) /)
        else
          pair = (/ tcorners(jj), tcorners(ii) /)
        end if

        do iedge=1,19
          if(all( pair == cubedges(:,iedge) ))then
            tedges(icount) = iedge
            cycle jloop
          end if
        end do
        stop 'edge not found'
      end do jloop
    end do iloop

  end function get_tedges





  integer function get_edgeindex(icorner1, icorner2)

    integer, intent(in) :: icorner1, icorner2 
    integer :: ivc(2), iout, iedge

    if(icorner1<icorner2)then
      ivc(1)=icorner1
      ivc(2)=icorner2
    else
      ivc(1)=icorner2
      ivc(2)=icorner1
    end if

    iout=0

    do iedge=1,19
      if(all(ivc==cubedges(:,iedge)))then
        iout=iedge
        exit
      end if
    end do

#ifdef DEBUG
      if(iout==0) stop 'get_edgeindex: no matching case found!'
#endif

    get_edgeindex = iout
    return
  end function get_edgeindex





  double precision function area_triangle(kpoints)

    use mod_mathtools, only: crossprod

    double precision, intent(in) :: kpoints(3,3)
    double precision :: k21(3), k31(3), kcross(3)

    k21 = kpoints(:,3) - kpoints(:,1)
    k31 = kpoints(:,3) - kpoints(:,2)
    call crossprod(k21, k31, kcross)

    area_triangle = 0.5d0*sqrt(sum(kcross**2))

  end function area_triangle




  subroutine read_cubesrefine(nCub3, nmarked, imarked)

    use mod_mympi,     only: myrank, nranks, master
    use mod_ioformat,  only: filename_cubesrefine, ext_formatted, fmt_fn_ext
    use mod_iohelp,    only: get_nLines
    use mod_mathtools, only: bubblesort_int
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    integer, intent(out) :: nCub3(3), nmarked
    integer, allocatable, intent(out) :: imarked(:)
    integer :: ierr, ii, ipointer, iline, istart, istop, nlines, ncubsplit, nCub3in(3), icubtmp, nrefine, nmarkedtmp
    integer, allocatable :: imarkedtmp(:), srtidx(:), imarked2(:), imarked3(:)
    character(len=256)  :: filename
    integer, parameter  :: iounit=15369

    if(myrank==master)then
      write(filename,fmt_fn_ext) filename_cubesrefine, ext_formatted
      open(unit=iounit, file=filename, form='formatted', action='read')

      ncubsplit = 0

      nlines = get_nLines( iounit )

      !check the input data and get total number of split cubes
      rewind(iounit)
      do iline = 1,nlines
        read(iounit,*)  nCub3in, nrefine, icubtmp
!       write(*,'(5I8)') nCub3in, nrefine, icubtmp
        if(iline==1) nCub3 = nCub3in*nrefine
        if(any(nCub3in*nrefine/=nCub3)) stop 'nCub3 not consistent in refine' 
        ncubsplit = ncubsplit + nrefine**3
      end do!iline

      allocate(imarked2(ncubsplit), imarked3(ncubsplit), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating imarked in readin on master'

      rewind(iounit)
      istart=1
      istop =0
      do iline = 1,nlines
        allocate(imarkedtmp(1))
        nmarkedtmp = 1
        read(iounit,'(5I8)')  nCub3in, nrefine, imarkedtmp
        call cut_and_update_cubes(nrefine, nCub3in, nmarkedtmp, imarkedtmp)
        if(nmarkedtmp/=nrefine**3)then
          write(*,'(A,I0,A,I0)') 'Something went wrong. nmarkedtmp= ', nmarkedtmp, ', nrefine=', nrefine
          stop 'Error in read_cubesrefine'
        end if

        istop =istart-1+nrefine**3
        imarked2(istart:istop) = imarkedtmp

        deallocate(imarkedtmp)
        istart=istop+1

      end do!iline
      nmarkedtmp = istart-1
      close(iounit)

      if(nmarkedtmp/=ncubsplit) then
        write(*,'(A,I0,A,I0)') 'Something went wrong. nmarkedtmp= ', nmarkedtmp, ', ncubsplit=', ncubsplit
        stop 'Error in read_cubesrefine'
      end if

      !sort the index and filter out double indices
      allocate(srtidx(ncubsplit), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating srtidx in readin on master'
      call bubblesort_int(nmarkedtmp, imarked2, srtidx)
      imarked3 = -1
      imarked3(1) = imarked2(srtidx(1))
      ipointer = 1
      do iline=2,nmarkedtmp
        if(imarked2(srtidx(iline))>imarked3(ipointer))then
          imarked3(ipointer+1) = imarked2(srtidx(iline))
          ipointer = ipointer+1
        end if
      end do

      nmarked = ipointer
      allocate(imarked(nmarked), STAT=ierr)
      imarked = imarked3(1:nmarked)

!      write(666,'(I8)') nmarked, imarked

    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(nCub3, 3, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in Bcast(nCub3)'
    call MPI_Bcast(nmarked, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in Bcast(nmarked)'
    if(myrank/=master)then
      allocate(imarked(nmarked), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating imarked in readin on slaves'
    end if!myrank/=master
    call MPI_Bcast(imarked, nmarked, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in Bcast(imarked)'
#endif

!   write(*,'("myrank=",I0," imarked=",10I8)') myrank, imarked
!   write(2000+myrank,'(I8)')   nmarked
!   write(2000+myrank,'(10I8)') imarked
  end subroutine read_cubesrefine





  subroutine split_triangle(ktriangle, eigw_triangle, kstore, eigw_store, kcounter, nfaces, nvec, dscal)

    use mod_symmetries, only: singlepoint_in_wedge
    implicit none
    integer,            intent(in)    :: nfaces
    double precision,   intent(in)    :: ktriangle(3,3), nvec(3,nfaces), dscal(nfaces)
    double complex,     intent(in)    :: eigw_triangle(3)
    integer,            intent(inout) :: kcounter
    double precision,   intent(inout) :: kstore(3,nkpmax)
    double complex,     intent(inout) :: eigw_store(nkpmax)

    integer :: ntri, ii, iface, itri, cutcase, its(0:2), ikp(3), nLines, nLinesINP
    double precision :: dtmp, ktmp1(3,nkpmax), ktmp2(3,nkpmax), denom, rscal, ki0new(3,2)
    double complex :: etmp1(nkpmax), etmp2(nkpmax), ei0new(2)
    logical :: linIBZ(3), lplane(3)

    double precision, parameter :: eps=1d-10

    !test whether the points lie inside the IBZ
    do ii=1,3
      linIBZ(ii) = singlepoint_in_wedge(nfaces, nvec, dscal, ktriangle(:,ii))
    end do!ii

    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    if(all(linIBZ))then ! DISTINGUISH THE DIFFERENT CASES !
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


      !the triangle lies completely inside the IBZ --> just copy
      do ii=1,3
        kstore(:,ii+kcounter) = ktriangle(:,ii)
        eigw_store(ii+kcounter) = eigw_triangle(ii)
      end do
      kcounter = kcounter+3


    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    elseif(any(linIBZ))then ! DISTINGUISG THE DIFFERENT CASES !
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


      !the triangle lies partly in the IBZ --> cut it

      nLinesINP=3
      ktmp1(:,1:3) = ktriangle
      etmp1(1:3) = eigw_triangle

      !===== BEGIN =====!
      ! loop over faces !
      do iface=1,nfaces

        nLines=0!init counter for all triangles found
        ntri=nLinesINP/3

        !(updated) loop over triangles
        do itri=1,ntri

          ikp(1) = 3*(itri-1)+1
          ikp(2) = 3*(itri-1)+2
          ikp(3) = 3*itri

          !check if points are inside or outside the plane
          lplane(:) = .true.
          do ii=1,3
            dtmp = sum(nvec(:,iface)*ktmp1(:,ikp(ii)))-dscal(iface)+eps
            if(dtmp<0d0)then
              lplane(ii) = .false.
            end if!eps
          end do!ii

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
          if(all(lplane))then ! distinguish the different cases --- part 2 !
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

            !the whole triangles lies on the corect side of the plane
            ktmp2(:,nLines+1:nLines+3) = ktmp1(:,ikp)
            etmp2(nLines+1:nLines+3) = etmp1(ikp)
            nLines = nLines+3

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
          elseif(any(lplane))then ! distinguish the different cases --- part 2 !
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

            !the triangles crosses the plane --> cut it

            !=== BEGIN ===!
            !get the single point and the cutcase
            its(:)  = -1
            cutcase = 0

            if(count(lplane)==1)then

              cutcase=1
              do ii=1,3
                if(lplane(ii)) its(0)=ii
              end do!ii

            elseif(count(lplane)==2)then

              cutcase=2
              do ii=1,3
                if(.not.lplane(ii)) its(0)=ii
              end do!ii

            else
              stop 'this case should not happen'
            end if!count
            !get the single point and the cutcase
            !===  END  ===!

            !get the indices of the other points
            its(1)=mod(its(0),  3)+1 
            its(2)=mod(its(0)+1,3)+1

            !=== BEGIN ===!
            ! cut the triangle according to the case
            do ii=1,2
              denom = sum(ktmp1(:,ikp(its(ii)))*nvec(:,iface)) - sum(ktmp1(:,ikp(its(0)))*nvec(:,iface))
              rscal = ( dscal(iface) - sum(ktmp1(:,ikp(its(0)))*nvec(:,iface)) )/denom
!             if(abs(rscal-1d0)<eps)then
!               ki0new(:,ii) = ktmp1(:,ikp(its(ii)))
!             else!abs(rscal)
                ki0new(:,ii) = ktmp1(:,ikp(its(0))) + rscal*( ktmp1(:,ikp(its(ii))) - ktmp1(:,ikp(its(0))) )
                ei0new(ii)   = etmp1(ikp(its(0)))   + rscal*( etmp1(ikp(its(ii))) -   etmp1(ikp(its(0))) )
!             end if!abs(rscal)
            end do!ii

            if(cutcase==1)then
              ktmp2(:,nLines+1) = ktmp1(:,ikp(its(0)))
              ktmp2(:,nLines+2) = ki0new(:,1)
              ktmp2(:,nLines+3) = ki0new(:,2)
              etmp2(nLines+1) = etmp1(ikp(its(0)))
              etmp2(nLines+2) = ei0new(1)
              etmp2(nLines+3) = ei0new(2)
              nLines = nLines+3
            elseif(cutcase==2)then
              ktmp2(:,nLines+1) = ktmp1(:,ikp(its(1)))
              ktmp2(:,nLines+2) = ktmp1(:,ikp(its(2)))
              ktmp2(:,nLines+3) = ki0new(:,1)
              ktmp2(:,nLines+4) = ktmp1(:,ikp(its(2)))
              ktmp2(:,nLines+5) = ki0new(:,2)
              ktmp2(:,nLines+6) = ki0new(:,1)
              etmp2(nLines+1) = etmp1(ikp(its(1)))
              etmp2(nLines+2) = etmp1(ikp(its(2)))
              etmp2(nLines+3) = ei0new(1)
              etmp2(nLines+4) = etmp1(ikp(its(2)))
              etmp2(nLines+5) = ei0new(2)
              etmp2(nLines+6) = ei0new(1)
              nLines = nLines+6
            else!cutcase
              stop 'this case may not happen = 3'
            end if!cutcase
            ! cut the triangle according to the case
            !=== END ===!


          !++++++++++++++++++++++++++++++++++++++++++++++++!
          else! distinguish the different cases --- part 2 !
          !++++++++++++++++++++++++++++++++++++++++++++++++!

            !the triangle lies on the wrong side of the plane --> throw away

          !++++++++++++++++++++++++++++++++++++++++++++++++++!
          end if! distinguish the different cases --- part 2 !
          !++++++++++++++++++++++++++++++++++++++++++++++++++!

          ktmp1(:,1:nLines) = ktmp2(:,1:nLines)
          etmp1(1:nLines)   = etmp2(1:nLines)
          nLinesINP = nLines

        end do!itri

      end do!iface
      ! loop over faces !
      !=====  END  =====!

      !copy here to big array
      do ii=1,nLinesINP
        kstore(:,ii+kcounter) = ktmp1(:,ii)
        eigw_store(ii+kcounter) = etmp1(ii)
      end do
      kcounter = kcounter+nLinesINP


    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    else ! DISTINGUISG THE DIFFERENT CASES !
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !the triangle lies completely outside of the IBZ --> throw it away

    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    end if!DISTINGUISG THE DIFFERENT CASES !
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  end subroutine split_triangle








end module mod_fermisurf_3D
