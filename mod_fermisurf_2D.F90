!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module mod_fermisurf_2D

  implicit none

  private
  public :: find_fermisurface_2D

  integer, parameter :: nkpmax = 512, nedges=5 ! shall be >> 32 to contain enough triangles per cube and band

  integer :: tridiags(2,2) = -1, squedges(2,5)=-1, tricorners(3,2)=-1, triedges(3,2)=-1

contains

  subroutine find_fermisurface_2D( inc, lattice, cluster, tgmatrx, symmetries, nCub3, nFSiter, nROOTiter, nstepsconnect, &
                                 & nCut_iter, roottype, rooteps, lrefine, nrefinenew, nkpts_int, kpoints_int, areas_int  )

    use type_inc,       only: inc_type
    use type_data,      only: lattice_type, cluster_type, tgmatrx_type
    use mod_symmetries, only: symmetries_type, get_2DIBZwedge_lines
    use mod_mympi,      only: myrank, nranks, master
    use mod_ioformat,   only: filename_cubesinfo, ext_formatted, fmt_fn_ext
    use mod_iohelp,     only: file_present
    use mod_fermisurf_basic, only: get_cubesinfo_filename, read_cubesfile, save_cubesfile, mark_cubes_FScross, mark_cubes_FScross_memopt, find_kpoints_irredset, save_kpointsfile_vis
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

    integer               :: nBZlines
    double precision      :: bounds(3,2)
    double precision, allocatable :: nvec(:,:), dscal(:)

    integer :: nmarked, ntotal
    integer, allocatable :: imarked(:)

    integer :: nkpts_vis, nkpts_irr
    integer,          allocatable :: kpt2irr(:), irr2kpt(:), vis2int(:)
    double precision, allocatable :: kpoints_vis(:,:)

    integer :: isqu, iter, iter2, ierr, ii, nCub3test(3), iterstart
    double precision :: dtmp
    character(len=256) :: filename, filetest
    logical :: l_cubesfile, l_exist

    call get_2DIBZwedge_lines(lattice%recbv, symmetries%nsym_used, symmetries%rotmat, symmetries%isym_used, nBZlines, nvec, dscal, bounds)

    if(nCub3(3)/=1) then
      if(myrank==master) write(*,*) 'WARNING!!!!!!!!!!  nCub(3) =/=1 in 2D mode. set nCub(3) = 1'
      nCub3(3) = 1
    end if



    !initialize square indices arrays
    call init_square2triangles()



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
        call read_squaresrefine(nCub3test, nmarked, imarked)
        if(any(nCub3/=nCub3test)) stop 'nCub3 not consistent between cubesinfo and cubesrefine'
        if(myrank==master)then
          write(filename,'("squares_refine_step=",I0,"_",A".txt")') 0, 'read'
          call squares2TXT(filename, nCub3, nmarked, imarked, bounds)
        end if!myrank==master

        call cut_and_update_squares(nrefinenew, nCub3, nmarked, imarked)
        if(myrank==master)then
          write(filename,'("squares_refine_step=",I0,"_",A".txt")') 1, 'cut'
          call squares2TXT(filename, nCub3, nmarked, imarked, bounds)
        end if!myrank==master

        call mark_squares_in_IBZ(nCub3, nBZlines, nvec, dscal, bounds, nmarked, imarked)
        if(myrank==master)then
          write(filename,'("squares_refine_step=",I0,"_",A".txt")') 2, 'ibz'
          call squares2TXT(filename, nCub3, nmarked, imarked, bounds)
        end if!myrank==master
      end if!lrefine==1

    else!cubesfile_present

      !search for existing cubesfiles that match the grid of an iteration within the iterative refinements
      iterstart=1
      !try to take the densest file, therefore count from top to bottom:
      iter_loop: do iter=nFSiter,1,-1

        !find the density of the cubes of iteration #iter
        nCub3test = nCub3
        do iter2=2,iter
          nCub3test(1:2) = nCub3test(1:2)*nCut_iter(iter2-1)
        end do!iter2

        !search for the existing cubesfile
        filetest = get_cubesinfo_filename(nCub3test,.true.)
        l_exist = file_present(trim(filetest))

        if(l_exist)then
          iterstart=iter+1
          call read_cubesfile(nCub3, nmarked, imarked, nCub3test)
          call cut_and_update_squares(nCut_iter(iter), nCub3, nmarked, imarked)
          ntotal = nmarked
          if(myrank==master) then
            write(filename,'("squares_iter=",I0,"_step=",I0,"_",(A),".txt")') iter, 3, 'cut'
            call squares2TXT(filename, nCub3, nmarked, imarked, bounds)
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
        do isqu=1,nmarked
          imarked(isqu) = isqu
        end do!icube
      end if!iterstart==1


      !======
      != scan the cubes for intersections
      !======
      do iter=iterstart,nFSiter

        !=== mark the cubes that lie (at least with one corner) within the first BZ ===!
        call mark_squares_in_IBZ(nCub3, nBZlines, nvec, dscal, bounds, nmarked, imarked)
        if(myrank==master) then
          write(filename,'("squares_iter=",I0,"_step=",I0,"_",(A),".txt")') iter, 1, 'IBZ'
          call squares2TXT(filename, nCub3, nmarked, imarked, bounds)
          write(*,'("***** Iteration ",I0," *****")') iter
          write(*,'(2X,2(A,I0),(A,F5.1,A))') 'Cubes found within the IBZ: ', nmarked, ' of ', ntotal, ' (= ', real(nmarked*100)/ntotal,' %)'
        end if!myrank==master

        !========= mark the cubes that cross the Fermi surface =========!
        !=== (only searching across the four diagonals of the cubes) ===!
        ntotal = nmarked
        if(inc%memopt==.true.) then
          call mark_cubes_FScross_memopt( inc, lattice, cluster, tgmatrx, nCub3,          &
                               & nstepsconnect(iter), 4, 2, tridiags, bounds,    &
                               & roottype(iter), rooteps(iter), nmarked, imarked )
        else
          call mark_cubes_FScross( inc, lattice, cluster, tgmatrx, nCub3,          &
                               & nstepsconnect(iter), 4, 2, tridiags, bounds,    &
                               & roottype(iter), rooteps(iter), nmarked, imarked )
        end if!inc%memopt==.true.
        if(myrank==master) then
          write(filename,'("squares_iter=",I0,"_step=",I0,"_",(A),".txt")') iter, 2, 'mark'
          call squares2TXT(filename, nCub3, nmarked, imarked, bounds)
          write(*,'(2X,2(A,I0),(A,F5.1,A))') 'Cubes found intersecting with FS: ', nmarked, ' of ', ntotal, ' (= ', real(nmarked*100)/ntotal,' %)'
          call save_cubesfile(nCub3, nmarked, imarked, lintermediate=.true.)
        end if!myrank==master

        !=== cut the remaining cubes into smaller pieces and update the indices ===!
!       write(1000+myrank,*) iter, nCut_iter(iter)
        call cut_and_update_squares(nCut_iter(iter), nCub3, nmarked, imarked)
        ntotal = nmarked
        if(myrank==master) then
          write(filename,'("squares_iter=",I0,"_step=",I0,"_",(A),".txt")') iter, 3, 'cut'
          call squares2TXT(filename, nCub3, nmarked, imarked, bounds)
        end if!myrank==master

      end do!iter

      !=== mark the cubes that lie (maybe only partly) within the first BZ ===!
      iter = nFSiter+1
      call mark_squares_in_IBZ(nCub3, nBZlines, nvec, dscal, bounds, nmarked, imarked)
      if(myrank==master) then
        write(filename,'("squares_iter=",I0,"_step=",I0,"_",(A),".txt")') iter, 1, 'IBZ'
        call squares2TXT(filename, nCub3, nmarked, imarked, bounds)
        call save_cubesfile(nCub3, nmarked, imarked)
      end if!myrank==master

    end if!cubesfile_present
    !=====================================!
    !=== generate or read a cubesfile ====!
    !=====================================!


    iter = nFSiter+1

    if(inc%memopt==.true.)then
      call find_intesection_lines_memopt( inc, lattice, cluster, tgmatrx, symmetries,  &
                                        & nCub3, bounds, nmarked, imarked,             &
                                        & nstepsconnect(iter), nROOTiter,              &
                                        & roottype(iter), rooteps(iter), nBZlines,     &
                                        & nvec, dscal, nkpts_vis, nkpts_int,           &
                                        & kpoints_vis, kpoints_int, areas_int, vis2int )
    else
      call find_intesection_lines( inc, lattice, cluster, tgmatrx, symmetries,  &
                                 & nCub3, bounds, nmarked, imarked,             &
                                 & nstepsconnect(iter), nROOTiter,              &
                                 & roottype(iter), rooteps(iter), nBZlines,     &
                                 & nvec, dscal, nkpts_vis, nkpts_int,           &
                                 & kpoints_vis, kpoints_int, areas_int, vis2int )
    end if!

    write(1368,*) nkpts_vis
    call find_kpoints_irredset( bounds, nkpts_vis, kpoints_vis, nkpts_irr, kpt2irr, irr2kpt)

    !save the visualization k-points to a file
    if(myrank==master) call save_kpointsfile_vis(nkpts_vis, nkpts_irr, kpoints_vis, symmetries%nsym_used, symmetries%isym_used, kpt2irr, irr2kpt, vis2int)


  end subroutine find_fermisurface_2D





!*******************************************************!
!***************** BEGIN BIG SUBROUTINE ****************!
!*******************************************************!
  subroutine find_intesection_lines( inc, lattice, cluster, tgmatrx, symmetries, &
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
!   use mod_vtkxml,     only: write_pointdata_rot
    use mod_ioformat,   only: fmt_fn_ext, fmt_fn_sub_ext, filemode_rot, filename_vtktest, fmt_fn_rank_ext, filename_outinfo, ext_formatted
    use mod_fermisurf_basic, only: roots_along_edge, compare_two_eigv_in_substeps, ROOT_IMAG, ROOT_REAL, generate_squarevertices
#ifdef CPP_MPI
    use mpi
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

    integer :: curid, nbands, nroots(nedges), lmid(inc%nrootmax,nedges), lmroot(inc%nrootmax,nedges)

    double precision :: kverts(3,4),&
                      & kends(3,2),&
                      & kmidsquare(3),&
                      & kroot_in_cube(3,inc%nrootmax,nedges),&
                      & kline(3,3),&
                      & kroot_lm_tri(3,nkpmax)   ! nkpmax >> (3 intersections per triangle) * (2 triangles per tetrahedron) * (6 tetrahedra) x
                                                 !           (multiplication-factor for splitting the triangle if BZ-face is crossed)
    double complex   :: eigw_lm_tri(nkpmax),&
                      & eigw_line(3)

    double complex, allocatable :: LVroot(:,:,:,:), &
                                 & RVroot(:,:,:,:), &
                                 & eigwroot(:,:,:)

    double precision :: ksub(3,nsteps+1)
    double complex   :: eigwends(2,inc%nrootmax,nedges)

    integer :: iproblems, isquareproblems(2,10*nmarked/nranks), iproblems_tmpa(0:nranks-1), iproblems_tot
    integer, allocatable :: isquareproblems_all(:,:)

    integer :: ioff, nkpt, ntot_pT(0:nranks-1), ioff_pT(0:nranks-1), lb, ub
    integer :: ierr, ii, ivc, icub, itmp, itri, iedge, itriedge, itetroot,&
             & icount, iedge1, iedge2, iroot1, iroot2, lm0, lm1, lm2,     &
             & nfound_band, nfound_tet, lookup_tri(2,2), iverts_picked(4),&
             & sorted(4), sorted_tmp(4), i4edge(4), i4tet(4), isave,      &
             & kcounter_lm, kcounter_cub, printstep, itmparr(1)
    double precision :: dtmp_ri(4), dist(nkpmax), weight_lm, kdiff(3), dtmp
    double complex   :: eigw_picked(4)
    logical :: match
    character(len=256) :: filename
    character(len=1024) :: errormessage
    integer, parameter :: iofile=69823

!    write(*,*) 'in find_intesection_lines ok'

    !Parallelize
    call distribute_linear_on_tasks(nranks,myrank,master,nmarked,ntot_pT,ioff_pT,.true.)
    nkpt = ntot_pT(myrank)
    ioff = ioff_pT(myrank)

    !initialize
    isquareproblems = 0
    iproblems = 0

    !allocate arrays
    allocate( LVroot(inc%almso,inc%almso,inc%nrootmax,nedges),&
            & RVroot(inc%almso,inc%almso,inc%nrootmax,nedges),&
            & eigwroot(inc%almso,inc%nrootmax,nedges),        &
            & STAT=ierr, ERRMSG=errormessage )
!   dtmp = dble(inc%almso)**2*inc%nrootmax*nedges*2 + dble(inc%almso)*inc%nrootmax*nedges
!   if(myrank==master) write(*,*) 'Trying to allocate ', dtmp*16d0/1024**3, ' GB'
!   write(13200+myrank,*) 'Trying to allocate ', dtmp*16d0/1024**3, ' GB'
!   write(13200+myrank,*) inc%almso, inc%nrootmax, nedges, dtmp
    if(ierr/=0)then
       write(*,*) 'Problem allocating LVroot etc. in find_intesection_lines'
       dtmp = dble(inc%almso)**2*inc%nrootmax*nedges*2 + dble(inc%almso)*inc%nrootmax*nedges
       write(*,*) 'Error Status=', ierr
       write(*,*) 'Error message=', trim(errormessage)
       write(*,*) 'Trying to allocate ', dtmp*16d0/1024**3, ' GB'
       stop 'Problem allocating LVroot etc. in find_intesection_lines'
    end if

    npoints_vis = 0
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

#ifdef CPP_DEBUG
      write(iofile,'(A,I8)') 'starting cube ', imarked(icub+ioff)
#endif

      call generate_squarevertices(nCub3,imarked(icub+ioff), bounds, kverts)
      !find middle point of the cube
      kmidsquare = (kverts(:,1)+kverts(:,4))/2

#ifdef CPP_DEBUG
      write(iofile,'(2X,"kverts(:,i)= ",3ES18.9)') kverts
#endif


      !==================================!
      !=== find roots along the edges ===!
      !==================================!
      do iedge=1,nedges
        !pick k-points
        kends(:,1) = kverts(:,squedges(1,iedge))
        kends(:,2) = kverts(:,squedges(2,iedge))
        call roots_along_edge( inc, lattice, cluster, tgmatrx, nsteps, kends, niter,&
                             & roottype, rooteps, nroots(iedge), lmroot(:,iedge),   &
                             & kroot_in_cube(:,:,iedge), LVroot(:,:,:,iedge),       &
                             & RVroot(:,:,:,iedge), eigwroot(:,:,iedge), -1,        &
                             & eigwends(:,:,iedge)                                  )
      end do!iedge

#ifdef CPP_DEBUG
      write(iofile,'(2X,"#roots found on the edges:")')
      do iedge=1,nedges
        write(iofile,'(4X,"iedge= ",I0,", # = ",I0)') iedge, nroots(iedge)
      end do
#endif


      !======================================================!
      !=== determine which kroots belong to the same band ===!
      !======================================================!
      lmid=0
      ! the array lmid(iroot,iedge) will contain a unique id to which band a k-point belongs.
      curid=0
      do iedge1=1,nedges
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
          loop3: do iedge2=iedge1+1,nedges
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
      do iedge=1,nedges
          write(iofile,'(4X,"iedge= ",I0," - ",20I8)') iedge, lmid(1:nroots(iedge),iedge)
      end do!iedge

      write(iofile,'(2X,A,I0)') 'nbands= ', nbands
#endif



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
        do itri=1,2
#ifdef CPP_DEBUG
          write(iofile,'(4X,2(A,I0))') 'lm0= ', lm0, ', itri= ', itri
#endif

          nfound_tet = 0

          !find all roots belonging to this this tetrahedron and this band
          edgeloop: do iedge1=1,3
            itriedge = triedges(iedge1,itri)
!           write(*,'(20(A,I0))') 'itetedge= ', itetedge, ', #roots=', nroots(itetedge)
            do iroot1=1,nroots(itriedge)
              if(lmid(iroot1,itriedge)==lm0)then
                nfound_tet = nfound_tet+1
                if(nfound_tet>2)then
#ifdef CPP_DEBUG
                  write(iofile,'(I0,A)') imarked(icub+ioff), ' =imarked: more than 2 intersections of band with triangle --> skip'
#endif
                  iproblems = iproblems+1
                  isquareproblems(1,iproblems) = imarked(icub+ioff)
                  isquareproblems(2,iproblems) = 1
                  nfound_tet=0
                  exit edgeloop
                end if
                lookup_tri(:,nfound_tet) = (/ iroot1, itriedge /)
              end if
            end do!iroot1
          end do edgeloop!iedge1

#ifdef CPP_DEBUG
          write(iofile,'(6X,A,I0)') 'nfound_tet= ', nfound_tet
#endif


          select case ( nfound_tet )
            case( 0 ); cycle
            case( 2 )
              !store the triangle k-points
              kline(:,1) = kroot_in_cube(:,lookup_tri(1,1),lookup_tri(2,1))
              kline(:,2) = kroot_in_cube(:,lookup_tri(1,2),lookup_tri(2,2))
              eigw_line(1) = eigwroot(lmroot(lookup_tri(1,1),lookup_tri(2,1)),lookup_tri(1,1),lookup_tri(2,1))
              eigw_line(2) = eigwroot(lmroot(lookup_tri(1,2),lookup_tri(2,2)),lookup_tri(1,2),lookup_tri(2,2))

#ifdef CPP_DEBUG
              write(iofile,'(8X,A,3ES25.16," | ",2ES25.16)') 'k1=', kline(:,1), eigw_line(1)
              write(iofile,'(8X,A,3ES25.16," | ",2ES25.16)') 'k2=', kline(:,2), eigw_line(2)
#endif

              !truncate line if BZ-face is crossed
              call split_line(kline, eigw_line, kroot_lm_tri(:,:), eigw_lm_tri(:), kcounter_lm, nfaces, nvec, dscal)

            case default
#ifdef CPP_DEBUG
              write(iofile,'(I0,A)') imarked(icub+ioff), ' =imarked: neither 0 nor 2 intesections found for --> skip'
#endif
              iproblems = iproblems+1
              isquareproblems(1,iproblems) = imarked(icub+ioff)
              isquareproblems(2,iproblems) = 2
          end select



        end do!itri

#ifdef CPP_DEBUG
        write(iofile,'(6X,A,I0)') 'cut lines; kcounter_lm = ', kcounter_lm
        do ii=1,kcounter_lm
          write(iofile,'(8X,A,3ES25.16," | ",2ES25.16)') 'k=', kroot_lm_tri(:,ii), eigw_lm_tri(ii)
        end do!ii 
#endif

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
          dist(ii) = sum((kroot_lm_tri(:,ii) - kmidsquare)**2)
        end do
        call findminindex(kcounter_lm, dist(1:kcounter_lm), isave)

        !save the k-point
        !kpoints_int(:,kcounter_cub) = kroot_lm_tri(:,isave) 
        ! shift the kpoint toward the center of the cube by a tiny change
        ! to be able to identify the cube later on
        kpoints_int(:,kcounter_cub) = kroot_lm_tri(:,isave) +(kmidsquare-kroot_lm_tri(:,isave))*1d-6

        !find the weight of the representing cube (=area of all triangles)
        weight_lm = 0d0
        do ii=1,kcounter_lm/2
          kdiff = kroot_lm_tri(:,ii*2) - kroot_lm_tri(:,ii*2-1)
          dtmp = sqrt(sum(kdiff**2))
          weight_lm = weight_lm + dtmp
        end do!ii
        areas_int(kcounter_cub) = weight_lm

        !if weight for this triangle is negligible, forget this k-point
        if(abs(weight_lm)<1d-16) kcounter_cub = kcounter_cub - 1


        !=================================!
        !===  Store the visualization  ===!
        !=== information in temp array ===!
        !=================================!

        !add the visualization-triangles to the temporary storage array
        kpoints_vis(:,npoints_vis+1:npoints_vis+kcounter_lm) = kroot_lm_tri(:,1:kcounter_lm)
        eigw_vis(npoints_vis+1:npoints_vis+kcounter_lm) = eigw_lm_tri(1:kcounter_lm)
        vis2int(npoints_vis+1:npoints_vis+kcounter_lm) = kcounter_cub
        npoints_vis = npoints_vis + kcounter_lm


      end do!lm0

    end do!icub

    if(myrank==master) write(*,*)
    !**************************************!
    !*** (parallelized) loop over cubes ***!
    !**************************************!
    !*************** E N D ****************!
    !**************************************!




    !***************************************!
    !************* B E G I N ***************!
    !***************************************!
    !*** collect the problematic squares ***!
    !***************************************!
#ifdef CPP_MPI
    !collect the numbers of problematic arrays
    call MPI_Allgather( iproblems,      1, MPI_INTEGER, &
                      & iproblems_tmpa, 1, MPI_INTEGER, &
                      & MPI_COMM_WORLD, ierr              )
    if(ierr/=MPI_SUCCESS) stop 'Problem in allgather npoints_vis_glob'
    iproblems_tot = sum(iproblems_tmpa)
    allocate(isquareproblems_all(2,iproblems_tot), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating isquareproblems_all'

    !calculate offsets
    ioffs_tmpa=0
    do irank=1,nranks-1
      ioffs_tmpa(irank) = ioffs_tmpa(irank-1)+iproblems_tmpa(irank-1)
    end do!irank
    ioffs_save = ioffs_tmpa

    !collect the problematic squares
    isendarr = 2*iproblems_tmpa
    ioffsarr = 2*ioffs_tmpa
    call MPI_Allgatherv( isquareproblems, isendarr(myrank), MPI_INTEGER, &
                       & isquareproblems_all, isendarr, ioffsarr,        &
                       & MPI_INTEGER, MPI_COMM_WORLD, ierr               )
    if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv isquareproblems'
#else
    iproblems_tot = iproblems
    allocate(isquareproblems_all(2,iproblems_tot), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating isquareproblems_all'
    isquareproblems_all = isquareproblems(:,1:iproblems_tot)
#endif

#ifndef CPP_DEBUG
    if(myrank==master)then
      write(filename,fmt_fn_rank_ext) filename_outinfo, myrank, ext_formatted
      open(unit=13512, file=trim(filename), action='write', form='formatted')
      do ii=1,iproblems_tot
        select case (isquareproblems_all(2,ii))
          case(1); write(13512,'(I0,A)') isquareproblems_all(1,ii), ' =imarked: more than 2 intersections of band with triangle --> skip'
          case(2); write(13512,'(I0,A)') isquareproblems_all(1,ii), ' =imarked: neither 0 nor 2 intesections found for --> skip'
          case default; write(13512,'(I0,A,I0,A)') isquareproblems_all(1,ii), ' =imarked: identifier ', isquareproblems_all(2,ii), ' not known'
        end select
      end do!ii
      close(13512)
    end if!myrank==master
#endif
    !***************************************!
    !*** collect the problematic squares ***!
    !***************************************!
    !*************** E N D *****************!
    !***************************************!


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

      write(filename,fmt_fn_ext) 'fstest_eigw0', ext_formatted
      open(unit=3654,file=filename,form='formatted',action='write')
      write(3654,'(2ES25.16)') eigw_vis_all
      close(3654)


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



#ifdef CPP_DEBUG
    close(iofile)
#endif

190 FORMAT('                 |'$)
200 FORMAT('|'$)

  end subroutine find_intesection_lines

!*******************************************************!
!*****************   END BIG SUBROUTINE ****************!
!*******************************************************!


  subroutine find_intesection_lines_memopt( inc, lattice, cluster, tgmatrx, symmetries, &
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
!   use mod_vtkxml,     only: write_pointdata_rot
    use mod_ioformat,   only: fmt_fn_ext, fmt_fn_sub_ext, filemode_rot, filename_vtktest, fmt_fn_rank_ext, filename_outinfo, ext_formatted
    use mod_iohelp,     only: file_present2
    use mod_fermisurf_basic, only: roots_along_edge, compare_two_eigv_in_substeps, ROOT_IMAG, ROOT_REAL, generate_squarevertices, &
                                 & roots_along_edge_memopt, compare_two_eigv_in_substeps_memopt
#ifdef CPP_MPI
    use mpi
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

    integer :: curid, nbands, nroots(nedges), lmid(inc%nrootmax,nedges), lmroot(inc%nrootmax,nedges)

    double precision :: kverts(3,4),&
                      & kends(3,2),&
                      & kmidsquare(3),&
                      & kroot_in_cube(3,inc%nrootmax,nedges),&
                      & kline(3,3),&
                      & kroot_lm_tri(3,nkpmax)   ! nkpmax >> (3 intersections per triangle) * (2 triangles per tetrahedron) * (6 tetrahedra) x
                                                 !           (multiplication-factor for splitting the triangle if BZ-face is crossed)
    double complex   :: eigw_lm_tri(nkpmax),&
                      & eigw_line(3)

    double complex, allocatable :: LVroot(:,:,:,:), &
                                 & RVroot(:,:,:,:), &
                                 & eigwroot(:,:,:)

    double precision :: ksub(3,nsteps+1)
    double complex   :: eigwends(2,inc%nrootmax,nedges)

    integer :: iproblems, isquareproblems(2,10*nmarked/nranks), iproblems_tmpa(0:nranks-1), iproblems_tot
    integer, allocatable :: isquareproblems_all(:,:)

    integer :: ioff, nkpt, ntot_pT(0:nranks-1), ioff_pT(0:nranks-1), lb, ub, i
    integer :: ierr, ii, ivc, icub, itmp, itri, iedge, itriedge, itetroot,&
             & icount, iedge1, iedge2, iroot1, iroot2, lm0, lm1, lm2,     &
             & nfound_band, nfound_tet, lookup_tri(2,2), iverts_picked(4),&
             & sorted(4), sorted_tmp(4), i4edge(4), i4tet(4), isave,      &
             & kcounter_lm, kcounter_cub, printstep, itmparr(1)
    integer :: kcounter_cub_file, kcounter_lm_file
    double precision :: dtmp_ri(4), dist(nkpmax), weight_lm, kdiff(3), dtmp
    double complex   :: eigw_picked(4)
    logical :: match
    character(len=256) :: filename
    character(len=256) :: filename_cube
    character(len=1024) :: errormessage
    integer, parameter :: iofile=69823
    integer, parameter :: iofile_cube=69824

!    write(*,*) 'in find_intesection_lines ok'

    !Parallelize
    call distribute_linear_on_tasks(nranks,myrank,master,nmarked,ntot_pT,ioff_pT,.true.)
    nkpt = ntot_pT(myrank)
    ioff = ioff_pT(myrank)

    !initialize
    isquareproblems = 0
    iproblems = 0

    !allocate arrays
    allocate( LVroot(inc%almso,inc%neig,inc%nrootmax,nedges),&
            & RVroot(inc%almso,inc%neig,inc%nrootmax,nedges),&
            & eigwroot(inc%neig,inc%nrootmax,nedges),        &
            & STAT=ierr, ERRMSG=errormessage )
!   dtmp = dble(inc%almso)**2*inc%nrootmax*nedges*2 + dble(inc%almso)*inc%nrootmax*nedges
!   if(myrank==master) write(*,*) 'Trying to allocate ', dtmp*16d0/1024**3, ' GB'
!   write(13200+myrank,*) 'Trying to allocate ', dtmp*16d0/1024**3, ' GB'
!   write(13200+myrank,*) inc%almso, inc%nrootmax, nedges, dtmp
    if(ierr/=0)then
       write(*,*) 'Problem allocating LVroot etc. in find_intesection_lines'
       dtmp = dble(inc%almso)*inc%neig*inc%nrootmax*nedges*2 + inc%neig*inc%nrootmax*nedges
       write(*,*) 'Error Status=', ierr
       write(*,*) 'Error message=', trim(errormessage)
       write(*,*) 'Trying to allocate ', dtmp*16d0/1024**3, ' GB'
       stop 'Problem allocating LVroot etc. in find_intesection_lines'
    end if

    npoints_vis = 0
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


    !**************************************!
    !************* B E G I N **************!
    !**************************************!
    !*** (parallelized) loop over cubes ***!
    !**************************************!
    if(myrank==master) write(*,FMT=190) !beginning of statusbar
    kcounter_cub=0
    areas_int = 0d0
    do icub=1,nkpt
      kcounter_cub_file=0
      kcounter_lm_file=0
      !update statusbar
      if(mod(icub,printstep)==0 .and. myrank==master) write(*,FMT=200)

#ifdef CPP_DEBUG
      write(iofile,'(A,I8)') 'starting cube ', imarked(icub+ioff)
#endif
      write(filename_cube,'(1I6,A)') imarked(icub+ioff), ".cube"
      write(*,*) "looking for ", filename_cube
      ! compute a cube only if no file was found from previous calculation
      write(*,*) "file_present2(filename) : ",file_present2(trim(filename_cube))
      if(.not.file_present2(trim(filename_cube))) then

       call generate_squarevertices(nCub3,imarked(icub+ioff), bounds, kverts)
       !find middle point of the cube
       kmidsquare = (kverts(:,1)+kverts(:,4))/2

#ifdef CPP_DEBUG
       write(iofile,'(2X,"kverts(:,i)= ",3ES18.9)') kverts
#endif


       !==================================!
       !=== find roots along the edges ===!
       !==================================!
       do iedge=1,nedges
         !pick k-points
         kends(:,1) = kverts(:,squedges(1,iedge))
         kends(:,2) = kverts(:,squedges(2,iedge))
         call roots_along_edge_memopt( inc, lattice, cluster, tgmatrx, nsteps, kends, niter,&
                                     & roottype, rooteps, nroots(iedge), lmroot(:,iedge),   &
                                     & kroot_in_cube(:,:,iedge), LVroot(:,:,:,iedge),       &
                                     & RVroot(:,:,:,iedge), eigwroot(:,:,iedge),  1,        &
                                     & eigwends(:,:,iedge)                                  )
       end do!iedge

#ifdef CPP_DEBUG
       write(iofile,'(2X,"#roots found on the edges:")')
       do iedge=1,nedges
         write(iofile,'(4X,"iedge= ",I0,", # = ",I0)') iedge, nroots(iedge)
       end do
#endif


       !======================================================!
       !=== determine which kroots belong to the same band ===!
       !======================================================!
       lmid=0
       ! the array lmid(iroot,iedge) will contain a unique id to which band a k-point belongs.
       curid=0
       do iedge1=1,nedges
         do iroot1=1,nroots(iedge1)
           !pick a first root
           lm1 = lmroot(iroot1,iedge1)
           if(lmid(iroot1,iedge1)>0) cycle !this root has already an ID

           curid = curid+1
           lmid(iroot1,iedge1)=curid
#ifdef CPP_DEBUG
           write(iofile,'(4X,"Set (iedge,iroot)=( ",I0,",",I0,") to curid= ",I0)') iedge1, iroot1, curid
#endif


!          write(*,'("Next cur-id= ",I0)'), curid
           !compare to all other roots on different edges which are not yet treated
           loop3: do iedge2=iedge1+1,nedges
             loop4: do iroot2=1,nroots(iedge2)
               if(lmid(iroot2,iedge2)>0) cycle
               lm2 = lmroot(iroot2,iedge2)
               kends(:,1) = kroot_in_cube(:,iroot1,iedge1)
               kends(:,2) = kroot_in_cube(:,iroot2,iedge2)
               call compare_two_eigv_in_substeps_memopt( inc, lattice, cluster, tgmatrx, nsteps,&
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
       do iedge=1,nedges
           write(iofile,'(4X,"iedge= ",I0," - ",20I8)') iedge, lmid(1:nroots(iedge),iedge)
       end do!iedge

       write(iofile,'(2X,A,I0)') 'nbands= ', nbands
#endif



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
         do itri=1,2
#ifdef CPP_DEBUG
           write(iofile,'(4X,2(A,I0))') 'lm0= ', lm0, ', itri= ', itri
#endif

           nfound_tet = 0

           !find all roots belonging to this this tetrahedron and this band
           edgeloop: do iedge1=1,3
             itriedge = triedges(iedge1,itri)
!            write(*,'(20(A,I0))') 'itetedge= ', itetedge, ', #roots=', nroots(itetedge)
             do iroot1=1,nroots(itriedge)
               if(lmid(iroot1,itriedge)==lm0)then
                 nfound_tet = nfound_tet+1
                 if(nfound_tet>2)then
#ifdef CPP_DEBUG
                   write(iofile,'(I0,A)') imarked(icub+ioff), ' =imarked: more than 2 intersections of band with triangle --> skip'
#endif
                   iproblems = iproblems+1
                   isquareproblems(1,iproblems) = imarked(icub+ioff)
                   isquareproblems(2,iproblems) = 1
                   nfound_tet=0
                   exit edgeloop
                 end if
                 lookup_tri(:,nfound_tet) = (/ iroot1, itriedge /)
               end if
             end do!iroot1
           end do edgeloop!iedge1

#ifdef CPP_DEBUG
           write(iofile,'(6X,A,I0)') 'nfound_tet= ', nfound_tet
#endif


           select case ( nfound_tet )
             case( 0 ); cycle
             case( 2 )
               !store the triangle k-points
               kline(:,1) = kroot_in_cube(:,lookup_tri(1,1),lookup_tri(2,1))
               kline(:,2) = kroot_in_cube(:,lookup_tri(1,2),lookup_tri(2,2))
               eigw_line(1) = eigwroot(lmroot(lookup_tri(1,1),lookup_tri(2,1)),lookup_tri(1,1),lookup_tri(2,1))
               eigw_line(2) = eigwroot(lmroot(lookup_tri(1,2),lookup_tri(2,2)),lookup_tri(1,2),lookup_tri(2,2))

#ifdef CPP_DEBUG
               write(iofile,'(8X,A,3ES25.16," | ",2ES25.16)') 'k1=', kline(:,1), eigw_line(1)
               write(iofile,'(8X,A,3ES25.16," | ",2ES25.16)') 'k2=', kline(:,2), eigw_line(2)
#endif

               !truncate line if BZ-face is crossed
               call split_line(kline, eigw_line, kroot_lm_tri(:,:), eigw_lm_tri(:), kcounter_lm, nfaces, nvec, dscal)

             case default
#ifdef CPP_DEBUG
               write(iofile,'(I0,A)') imarked(icub+ioff), ' =imarked: neither 0 nor 2 intesections found for --> skip'
#endif
               iproblems = iproblems+1
               isquareproblems(1,iproblems) = imarked(icub+ioff)
               isquareproblems(2,iproblems) = 2
           end select



         end do!itri

#ifdef CPP_DEBUG
         write(iofile,'(6X,A,I0)') 'cut lines; kcounter_lm = ', kcounter_lm
         do ii=1,kcounter_lm
           write(iofile,'(8X,A,3ES25.16," | ",2ES25.16)') 'k=', kroot_lm_tri(:,ii), eigw_lm_tri(ii)
         end do!ii 
#endif

         !if cube is empty, proceed to next cube
         if(kcounter_lm==0) cycle


         !==============================!
         !===  Find the k-point and  ===!
         !=== weight for integration ===!
         !==============================!

         !update the counter for integration k-points
         kcounter_cub = kcounter_cub+1
         kcounter_cub_file = kcounter_cub_file+1

         !find kpoint-index representing this cube
         do ii=1,kcounter_lm
           dist(ii) = sum((kroot_lm_tri(:,ii) - kmidsquare)**2)
         end do
         call findminindex(kcounter_lm, dist(1:kcounter_lm), isave)

         !save the k-point
         kpoints_int(:,kcounter_cub) = kroot_lm_tri(:,isave)

         !find the weight of the representing cube (=area of all triangles)
         weight_lm = 0d0
         do ii=1,kcounter_lm/2
           kdiff = kroot_lm_tri(:,ii*2) - kroot_lm_tri(:,ii*2-1)
           dtmp = sqrt(sum(kdiff**2))
           weight_lm = weight_lm + dtmp
         end do!ii
         areas_int(kcounter_cub) = weight_lm

         !if weight for this triangle is negligible, forget this k-point
         if(abs(weight_lm)<1d-16) then
           kcounter_cub = kcounter_cub - 1
           kcounter_cub_file = kcounter_cub_file - 1
         endif!abs(weight_lm)<1d-16


        !=================================!
        !===  Store the visualization  ===!
        !=== information in temp array ===!
        !=================================!

         !add the visualization-triangles to the temporary storage array
         kpoints_vis(:,npoints_vis+1:npoints_vis+kcounter_lm) = kroot_lm_tri(:,1:kcounter_lm)
         eigw_vis(npoints_vis+1:npoints_vis+kcounter_lm) = eigw_lm_tri(1:kcounter_lm)
         vis2int(npoints_vis+1:npoints_vis+kcounter_lm) = kcounter_cub
         npoints_vis = npoints_vis + kcounter_lm
         kcounter_lm_file = kcounter_lm_file + kcounter_lm

       end do!lm0

       ! save results for this cube in a separate cube file
       write(*,*) "writing in ", trim(filename_cube)
       open(unit=iofile_cube, file=trim(filename_cube), action='write', form='formatted')

       write(iofile_cube,'(1I2)') kcounter_cub_file
       if(kcounter_cub_file>0) then
         do i=1,kcounter_cub_file
           write(iofile_cube,'(3ES25.16)') kpoints_int(:,kcounter_cub)
           write(iofile_cube,'(1ES25.16)') areas_int(kcounter_cub)
         end do!i=1,kcounter_cub_file
       end if!kcounter_cub_file>0
       write(iofile_cube,'(1I2)') kcounter_lm_file
       if(kcounter_lm_file>0) then
         write(iofile_cube,'(3ES25.16)') kpoints_vis(:,npoints_vis+1-kcounter_lm_file:npoints_vis)
         write(iofile_cube,'(1I8)') vis2int(npoints_vis+1-kcounter_lm_file:npoints_vis)
       end if!kcounter_cub_file>0
       close(unit=iofile_cube)              
     else
       ! if file exist for this cube, read the content
       write(*,*) "reading ", trim(filename_cube)
       open(unit=iofile_cube, file=trim(filename_cube), action='read', form='formatted')

       read(iofile_cube,'(1I2)') kcounter_cub_file
       if(kcounter_cub_file>0) then
         do i=1,kcounter_cub_file
           kcounter_cub = kcounter_cub +1
           read(iofile_cube,'(3ES25.16)') kpoints_int(:,kcounter_cub)
           read(iofile_cube,'(1ES25.16)') areas_int(kcounter_cub)
         end do!i=1,kcounter_cub_file
       end if!kcounter_cub_file>0
       read(iofile_cube,'(1I2)') kcounter_lm_file
       if(kcounter_lm_file>0) then
         read(iofile_cube,'(3ES25.16)') kpoints_vis(:,npoints_vis+1:npoints_vis+kcounter_lm_file)
         read(iofile_cube,'(1I8)') vis2int(npoints_vis+1:npoints_vis+kcounter_lm_file)
         npoints_vis = npoints_vis + kcounter_lm_file
       end if!kcounter_cub_file>0
       close(unit=iofile_cube)       
     end if!file_present(filename_cube)
    end do!icub

    if(myrank==master) write(*,*)
    !**************************************!
    !*** (parallelized) loop over cubes ***!
    !**************************************!
    !*************** E N D ****************!
    !**************************************!




    !***************************************!
    !************* B E G I N ***************!
    !***************************************!
    !*** collect the problematic squares ***!
    !***************************************!
#ifdef CPP_MPI
    !collect the numbers of problematic arrays
    call MPI_Allgather( iproblems,      1, MPI_INTEGER, &
                      & iproblems_tmpa, 1, MPI_INTEGER, &
                      & MPI_COMM_WORLD, ierr              )
    if(ierr/=MPI_SUCCESS) stop 'Problem in allgather npoints_vis_glob'
    iproblems_tot = sum(iproblems_tmpa)
    allocate(isquareproblems_all(2,iproblems_tot), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating isquareproblems_all'

    !calculate offsets
    ioffs_tmpa=0
    do irank=1,nranks-1
      ioffs_tmpa(irank) = ioffs_tmpa(irank-1)+iproblems_tmpa(irank-1)
    end do!irank
    ioffs_save = ioffs_tmpa

    !collect the problematic squares
    isendarr = 2*iproblems_tmpa
    ioffsarr = 2*ioffs_tmpa
    call MPI_Allgatherv( isquareproblems, isendarr(myrank), MPI_INTEGER, &
                       & isquareproblems_all, isendarr, ioffsarr,        &
                       & MPI_INTEGER, MPI_COMM_WORLD, ierr               )
    if(ierr/=MPI_SUCCESS) stop 'Problem in allgatherv isquareproblems'
#else
    iproblems_tot = iproblems
    allocate(isquareproblems_all(2,iproblems_tot), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating isquareproblems_all'
    isquareproblems_all = isquareproblems(:,1:iproblems_tot)
#endif

#ifndef CPP_DEBUG
    if(myrank==master)then
      write(filename,fmt_fn_rank_ext) filename_outinfo, myrank, ext_formatted
      open(unit=13512, file=trim(filename), action='write', form='formatted')
      do ii=1,iproblems_tot
        select case (isquareproblems_all(2,ii))
          case(1); write(13512,'(I0,A)') isquareproblems_all(1,ii), ' =imarked: more than 2 intersections of band with triangle --> skip'
          case(2); write(13512,'(I0,A)') isquareproblems_all(1,ii), ' =imarked: neither 0 nor 2 intesections found for --> skip'
          case default; write(13512,'(I0,A,I0,A)') isquareproblems_all(1,ii), ' =imarked: identifier ', isquareproblems_all(2,ii), ' not known'
        end select
      end do!ii
      close(13512)
    end if!myrank==master
#endif
    !***************************************!
    !*** collect the problematic squares ***!
    !***************************************!
    !*************** E N D *****************!
    !***************************************!


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

      write(filename,fmt_fn_ext) 'fstest_eigw0', ext_formatted
      open(unit=3654,file=filename,form='formatted',action='write')
      write(3654,'(2ES25.16)') eigw_vis_all
      close(3654)


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



#ifdef CPP_DEBUG
    close(iofile)
#endif

190 FORMAT('                 |'$)
200 FORMAT('|'$)

  end subroutine find_intesection_lines_memopt





  subroutine split_line(kline, eigw_line, kstore, eigw_store, kcounter, nfaces, nvec, dscal)

    use mod_symmetries, only: singlepoint_in_wedge
    implicit none
    integer,            intent(in)    :: nfaces
    double precision,   intent(in)    :: kline(3,2), nvec(3,nfaces), dscal(nfaces)
    double complex,     intent(in)    :: eigw_line(2)
    integer,            intent(inout) :: kcounter
    double precision,   intent(inout) :: kstore(3,nkpmax)
    double complex,     intent(inout) :: eigw_store(nkpmax)

    integer :: ii, iface
!    integer :: ntri, ii, iface, itri, cutcase, its(0:2), ikp(3), nLines, nLinesINP
!    double precision :: dtmp, ktmp1(3,nkpmax), ktmp2(3,nkpmax), denom, rscal, ki0new(3,2)
    double complex :: etmp(2), enew, ediff
    double precision :: ktmp(3,2), kdiff(3), knew(3), rscal, dtmp
    logical :: linIBZ(3), lplane(2)

    double precision, parameter :: eps=1d-10

    !test whether the points lie inside the IBZ
    do ii=1,2
      linIBZ(ii) = singlepoint_in_wedge(nfaces, nvec, dscal, kline(:,ii))
    end do!ii


    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    if(all(linIBZ))then ! DISTINGUISH THE DIFFERENT CASES !
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      !the line lies completely inside the IBZ --> just copy
      do ii=1,2
        kstore(:,ii+kcounter) = kline(:,ii)
        eigw_store(ii+kcounter) = eigw_line(ii)
      end do
      kcounter = kcounter+2


    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    elseif(any(linIBZ))then ! DISTINGUISG THE DIFFERENT CASES !
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      !the line lies partly in the IBZ --> cut it

      ktmp = kline
      etmp = eigw_line

      !===== BEGIN =====!
      ! loop over faces !
      do iface=1,nfaces

        !check if points are inside or outside the plane
        lplane(:) = .true.
        do ii=1,2
          dtmp = sum(nvec(:,iface)*ktmp(:,ii))-dscal(iface)+eps
          if(dtmp<0d0)then
            lplane(ii) = .false.
          end if!eps
        end do!ii


        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        if(all(lplane))then ! distinguish the different cases --- part 2      !
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
            !-- do nothing
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        elseif(any(lplane))then ! distinguish the different cases --- part 2 !
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

          !find intersection point of line with BZ face
          kdiff = ktmp(:,2)-ktmp(:,1)
          rscal = (dscal(iface) - sum(ktmp(:,1)*nvec(:,iface)))/sum(kdiff*nvec(:,iface))
          knew  = ktmp(:,1) + rscal*kdiff

          !interpolate eigenvalue
          ediff = etmp(2)-etmp(1)
          enew  = etmp(1) + rscal*ediff

          if(lplane(1))then
            ktmp(:,2) = knew
            etmp(2) = enew
          elseif(lplane(2))then
            ktmp(:,1) = knew
            etmp(1) = enew
          else
            stop 'this case should not happen in "split_line": case 1 '
          end if

        !+++++++++++++++++++++++++++++++++++++++++++++++++!
        else ! distinguish the different cases --- part 2 !
        !+++++++++++++++++++++++++++++++++++++++++++++++++!

          stop 'this case should not happen in "split_line": case 2'

        !++++++++++++++++++++++++++++++++++++++++++++++++++!
        end if! distinguish the different cases --- part 2 !
        !++++++++++++++++++++++++++++++++++++++++++++++++++!


      end do!iface
      ! loop over faces !
      !=====  END  =====!

      !copy here to big array
      do ii=1,2
        kstore(:,ii+kcounter) = ktmp(:,ii)
        eigw_store(ii+kcounter) = etmp(ii)
      end do
      kcounter = kcounter+2

    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    else ! DISTINGUISG THE DIFFERENT CASES !
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !the line lies completely outside of the IBZ --> throw it away

    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    end if!DISTINGUISG THE DIFFERENT CASES !
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  end subroutine split_line


  subroutine mark_squares_in_IBZ(nCub3, nfaces, nvec, dscal, bounds, nmarked, imarked)
    use mod_parutils,   only: distribute_linear_on_tasks
    use mod_mympi,      only: myrank, nranks, master
    use mod_symmetries, only: points_in_wedge
    use mod_fermisurf_basic, only: generate_squarevertices
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    integer,              intent(in)    :: nCub3(3), nfaces
    double precision,     intent(in)    :: nvec(3,nfaces), dscal(nfaces), bounds(3,2)
    integer,              intent(inout) :: nmarked
    integer, allocatable, intent(inout) :: imarked(:)

    integer :: nsqumarked_tmp, ioff, nkpt, ntot_pT(0:nranks-1), ioff_pT(0:nranks-1)
    integer :: nsqumarked_pT(0:nranks-1), iioff(0:nranks-1)
    double precision :: kverts(3,4)
    integer, allocatable :: isqumarked_tmp(:)

    integer :: irank, isqu, ierr
    logical :: squareinwedge


    !Parallelize
    call distribute_linear_on_tasks(nranks, myrank, master, nmarked, ntot_pT, ioff_pT, .false.)
    nkpt = ntot_pT(myrank)
    ioff = ioff_pT(myrank)

    !Create temp arrays
    allocate(isqumarked_tmp(nkpt), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating isqumarked_tmp in mark_squares_in_IBZ'
    nsqumarked_tmp = 0
    isqumarked_tmp = -1

    !Test whether squares lie in wedge
    do isqu=1,nkpt
      call generate_squarevertices(nCub3,imarked(isqu+ioff), bounds, kverts)
      squareinwedge = points_in_wedge(nfaces, nvec, dscal, 4, kverts, 'any')
      if(squareinwedge)then
        nsqumarked_tmp = nsqumarked_tmp + 1
        isqumarked_tmp(nsqumarked_tmp) = imarked(isqu+ioff)
      end if!squareinwedge
    end do!isqu

    deallocate(imarked)

#ifdef CPP_MPI
    !Combine results
    call MPI_Allgather(nsqumarked_tmp, 1, MPI_INTEGER, nsqumarked_pT, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
    if(ierr/=MPI_SUCCESS) stop 'Problem in Allgather(nsqumarked_tmp) in mark_squares_in_IBZ'

    nmarked = sum(nsqumarked_pT)
    allocate(imarked(nmarked), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating imarked in mark_squares_in_IBZ'

    iioff = 0
    do irank=1,nranks-1
      iioff(irank) = sum(nsqumarked_pT(0:irank-1))
    end do


    call MPI_Allgatherv( isqumarked_tmp, nsqumarked_tmp, MPI_INTEGER, &
                       & imarked, nsqumarked_pT, iioff, MPI_INTEGER,  &
                       & MPI_COMM_WORLD, ierr )
    if(ierr/=MPI_SUCCESS) stop 'Problem in Allgatherv(isqumarked_tmp) in mark_squares_in_IBZ'

#else
    nmarked=nsqumarked_tmp
    allocate(imarked(nmarked), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating imarked in mark_squares_in_IBZ'
    imarked = isqumarked_tmp(1:nsqumarked_tmp)
#endif

  end subroutine mark_squares_in_IBZ




  subroutine cut_and_update_squares(ncut, nCub3, nmarked, imarked)
    use mod_fermisurf_basic, only: unroll_ixyz
    use mod_mympi, only: myrank, master
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
    nmarked = nmarked_tmp * ncut**2
    allocate( imarked(nmarked), STAT=ierr )
    if(ierr/=0) stop 'Problem allocating imarked in "cut_and_update_squares"'
    nCub3 = nCub3in*ncut
    nCub3(3) = 1

    icount = 0
    do icub=1,nmarked_tmp
      call unroll_ixyz(imarked_tmp(icub), nCub3in, ii3)
!     do iz=0,ncut-1
      do iy=0,ncut-1
        do ix=0,ncut-1
          icount = icount+1
          imarked(icount) = (ii3(1)*ncut+ix+1) + (ii3(2)*ncut+iy)*nCub3(1)
        end do!ix
      end do!iy
!     end do!iz
    end do!icub

    if(icount/=nmarked .and. myrank==master) write(*,*) 'WARNING!! icount =/= nmarked in cut_and_update_squares'

  end subroutine cut_and_update_squares



  subroutine read_squaresrefine(nCub3, nmarked, imarked)

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
        ncubsplit = ncubsplit + nrefine**2
      end do!iline
      nCub3(3) = 1!set 3rd component to 1 for 2D-mode

      allocate(imarked2(ncubsplit), imarked3(ncubsplit), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating imarked in readin on master'

      rewind(iounit)
      istart=1
      istop =0
      do iline = 1,nlines
        allocate(imarkedtmp(1))
        nmarkedtmp = 1
        read(iounit,'(5I8)')  nCub3in, nrefine, imarkedtmp
        call cut_and_update_squares(nrefine, nCub3in, nmarkedtmp, imarkedtmp)
        if(nmarkedtmp/=nrefine**3)then
          write(*,'(A,I0,A,I0)') 'Something went wrong. nmarkedtmp= ', nmarkedtmp, ', nrefine=', nrefine
          stop 'Error in read_squaresrefine'
        end if

        istop =istart-1+nrefine**2
        imarked2(istart:istop) = imarkedtmp

        deallocate(imarkedtmp)
        istart=istop+1

      end do!iline
      nmarkedtmp = istart-1
      close(iounit)

      if(nmarkedtmp/=ncubsplit) then
        write(*,'(A,I0,A,I0)') 'Something went wrong. nmarkedtmp= ', nmarkedtmp, ', ncubsplit=', ncubsplit
        stop 'Error in read_squaresrefine'
      end if

!      write(111,'(A,I0)') 'ncubsplit=', ncubsplit
!      write(111,"(10I8)") imarked2(1:ncubsplit)

      !sort the index and filter out double indices
      allocate(srtidx(ncubsplit), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating srtidx in readin on master'
      call bubblesort_int(nmarkedtmp, imarked2, srtidx)
!      write(111,*) 'sorted:'
!      write(111,'(I8)') imarked2(srtidx)

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
!      write(111,'(A,I0)') 'nmarked=', nmarked
!      write(111,"(10I8)") imarked
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
  end subroutine read_squaresrefine




  subroutine squares2TXT(filename,nCub3,nmarked,imarked,bounds)
    use mod_fermisurf_basic, only: generate_squarevertices
    implicit none

    character(len=*), intent(in) :: filename
    integer,          intent(in) :: nCub3(3), nmarked, imarked(nmarked)
    double precision, intent(in) :: bounds(3,2)

    double precision :: kverts(3,4)
    integer :: isqu, ivert
    character(len=80) :: fmtstr

    open(unit=123895, file=trim(filename), action='write', form='formatted')

    !write header
    write(123895,'(A)') '# four lines define one square'

    !write points
    do isqu=1,nmarked
      call generate_squarevertices(nCub3,imarked(isqu),bounds,kverts)
      do ivert=1,4
        write(123895,'(3ES25.16)') kverts(:,ivert)
      end do!ivert
    end do!icub

    close(123895)

  end subroutine squares2TXT






  subroutine init_square2triangles()

    !   The triangles, lines and corners are labeled as follows
    !     ( numbers denote corners, [x] denote lines and (x) denote triangles):
    !
    !    3____[3]____4
    !    |          /|
    !    | (2)   /   |
    ! [4]|    /      |[2]
    !    | /   (1)   |
    !    |___________|
    !    1    [1]    2

    ! A square is cut into 2 triangles

    tridiags = reshape( (/ 1, 4, 2, 3 /), (/ 2, 2 /))

    !define the endpoints of the edges of the tetraeda. cubedge(1,j) = k means:
    !the startpoint of the j-th edge is the k-th corner-point of the cube 

    squedges(:,1) = (/ 1, 2 /) !
    squedges(:,2) = (/ 2, 4 /) !
    squedges(:,3) = (/ 3, 4 /) !
    squedges(:,4) = (/ 1, 3 /) !
    squedges(:,5) = (/ 1, 4 /) !


    tricorners(:,1) = (/ 1, 2, 4 /)
    tricorners(:,2) = (/ 1, 3, 4 /)

    triedges(:,1) = (/ 1, 2, 5 /)
    triedges(:,2) = (/ 3, 4, 5 /)

  end subroutine init_square2triangles



! 
end module mod_fermisurf_2D
