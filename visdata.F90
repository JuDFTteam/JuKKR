program visdata

  use type_inc,       only: inc_type
  use type_data,      only: lattice_type, cluster_type, tgmatrx_type
  use mod_mympi,      only: mympi_init, myrank, nranks, master
  use mod_symmetries, only: symmetries_type, set_symmetries, rotate_kpoints, expand_areas, expand_visarrays, expand_spinvalues, unfold_visarrays
  use mod_ioformat,   only: filemode_int, filemode_vis, fmt_fn_sub_ext, ext_formatted, filename_spin, filename_fvel
  use mod_iohelp,     only: getBZvolume
  use mod_read,       only: read_inc, read_TBkkrdata, read_kpointsfile_vis, read_kpointsfile_int, read_weights, read_fermivelocity, read_spinvalue
! use mod_calconfs,   only: calculate_spinmixing_int, calculate_spinmixing_vis
  use mod_vtkxml,     only: write_pointdata_rot

#ifdef CPP_MPI
    use mpi
#endif

    implicit none

    logical            :: lspin, lfvel, lcond, lcond2, lpairs
    character(len=256) :: filemode, filename

    type(inc_type)      :: inc
    type(lattice_type)  :: lattice
    type(cluster_type)  :: cluster
    type(tgmatrx_type)  :: tgmatrx

    !symmetry arrays
    integer :: nsym, nsym2
    integer, allocatable :: isym(:)
    type(symmetries_type) :: symmetries

    !local k-point arrays
    integer :: nkpts, nkpts_all, nkpts_int, nkpts_int_all, nkpts_vis_all, nkpts_inp_all
    integer, allocatable :: kpt2irr(:), irr2kpt(:), vis2int(:), ipairs(:)
    double precision :: tau_avg
    double precision, allocatable :: kpoints(:,:), areas(:), weights(:), fermivel(:,:), fermivel_2(:,:), meanfreepath(:,:,:), meanfreepath_sum(:), tauk(:), tauk2(:), kpoints_2(:,:), kpoints_3(:,:)

    !spinvalue locals
    integer              :: nsqa, ndegen1
    integer, allocatable :: ispincomb(:)
    double precision, allocatable :: nvect(:,:), spinval(:,:,:), spinval1(:,:,:), spinmix(:)

    !temp k-point arrays
    integer :: nkpts1, nkpts2, nkpts_all1, nkpts_all2
    integer,          allocatable :: kpt2irr1(:), irr2kpt1(:), kpt2irr2(:), irr2kpt2(:), vis2int2(:)
    double precision, allocatable :: areas1(:), weights1(:), kpoints1(:,:), kpoints2(:,:), fermivel1(:,:)

    double precision :: pi, BZVol, dos, dtmp
    integer          :: iset_select, ierr, isqa, ikp, isy, lb, ub, ispin

    integer                       :: nscalar, nvector, iscalar, ivector
    double precision, allocatable :: scalardata(:,:),  &
                                   & vectordata(:,:,:),&
                                   & dtmparr(:),       &
                                   & dtmparr2(:),      &
                                   & fvelabs(:)
    character(len=256), allocatable :: scalarstring(:), vectorstring(:)
    character(len=256) :: dummyline

    !init
#ifdef CPP_MPI
    call MPI_Init ( ierr )
#endif
    call mympi_init()

    !Read in TBKKR-data
    call read_inc(inc)
    call read_TBkkrdata(inc, lattice, cluster, tgmatrx)


    call set_symmetries(inc, lattice, symmetries)
    BZVol = getBZvolume(lattice%recbv)

!   write(*,*) "Which set shall I visualize?"
!   write(*,*) " 1=integration set"
!   write(*,*) " 2=visualization set"
!   read(*,*)  iset_select
    iset_select=1

    if(iset_select==1) filemode = filemode_int
    if(iset_select==2) filemode = filemode_vis



    nvector=0
    nscalar=0


    !=============================!
    != Read in the k-point files =!
    if(iset_select==1)then
      call read_kpointsfile_vis(nkpts_all2, nkpts2, kpoints2, nsym, isym, kpt2irr2, irr2kpt2,vis2int=vis2int2)
      call unfold_visarrays(nkpts_all2, nkpts2, kpt2irr2, kpoints2, kpt2irr1, irr2kpt1, kpoints1)
      nkpts1    =nkpts_all2
      nkpts_all1=nkpts_all2
      deallocate(kpt2irr2, irr2kpt2, kpoints2)

      call rotate_kpoints(symmetries%rotmat, nkpts1, kpoints1, nsym, isym, nkpts, kpoints)
      call expand_visarrays(nsym, nkpts_all1, nkpts1, kpt2irr1, irr2kpt1, kpt2irr, irr2kpt)
      nkpts_vis_all = nsym*nkpts_all1
      deallocate(kpt2irr1, irr2kpt1, kpoints1, isym)

      call read_kpointsfile_int(nkpts_int, kpoints1, areas1, nsym2, isym)
      if(nsym2/=nsym) stop 'nsym2/=nsym'
      if(vis2int2(nkpts_all1)/=nkpts_int) stop 'last element from vis2int/=nkpts_int'
      nkpts_int_all = nsym*nkpts_int

      call rotate_kpoints(symmetries%rotmat, nkpts_int, kpoints1, nsym, isym, nkpts2, kpoints_2)
      if(nkpts2/=nkpts_int_all) stop 'inconsistency in number of k-points 2'

      deallocate(isym, areas1, kpoints1)

      allocate(vis2int(nkpts_all1*nsym), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating vis2int'
      do isy=1,nsym
        lb=(isy-1)*nkpts_all1+1
        ub=isy*nkpts_all1
        vis2int(lb:ub) = vis2int2 + (isy-1)*nkpts_int
      end do!isym

      nkpts_inp_all = nkpts_int_all

    else

      call read_kpointsfile_vis(nkpts_all1, nkpts1, kpoints1, nsym, isym, kpt2irr1, irr2kpt1)
      call rotate_kpoints(symmetries%rotmat, nkpts1, kpoints1, nsym, isym, nkpts, kpoints)
      call expand_visarrays(nsym, nkpts_all1, nkpts1, kpt2irr1, irr2kpt1, kpt2irr, irr2kpt)
      nkpts_vis_all = nsym*nkpts_all1
      deallocate(kpt2irr1, irr2kpt1, kpoints1, isym)

      nkpts_inp_all = nsym*nkpts1

    end if!iset_select==1


    !=========================!
    != Read in the spinvalue =!
    write(filename,fmt_fn_sub_ext) filename_spin, trim(filemode), ext_formatted
    inquire(file=filename, exist=lspin)
    if(lspin) then
        call read_spinvalue(trim(filemode), nkpts1, nsqa, ndegen1, ispincomb, nvect, spinval1, nsym, isym)
        if(nkpts1*nsym/=nkpts_inp_all) stop 'inconsistency in number of k-points'
        if(ndegen1/=inc%ndegen) stop 'inconsistency in ndegen'
        call expand_spinvalues(nsym,ndegen1,nsqa,nkpts1,spinval1,spinval)
        deallocate(spinval1, isym)
        nscalar = nscalar + nsqa
    end if!lspin


    !==============================!
    != Read in the fermi velocity =!
    write(filename,fmt_fn_sub_ext) filename_fvel, trim(filemode), ext_formatted
    inquire(file=filename, exist=lfvel)
    if(lfvel) then
        call read_fermivelocity(trim(filemode), nkpts1, fermivel1, nsym, isym)
        call rotate_kpoints(symmetries%rotmat, nkpts1, fermivel1, nsym, isym, nkpts2, fermivel)
        if(nkpts2/=nkpts_inp_all) stop 'inconsistency in number of k-points'
        deallocate(fermivel1, isym)
        nvector = nvector+1
    end if!lfvel



    !================================!
    != Read in the cond_output file =!
    inquire(file='output_lifetime', exist=lcond2)
    lcond = lcond2.and.(iset_select==1)
    if(lcond) then
      call read_condfile(inc%ndegen,nkpts_inp_all, kpoints_3, meanfreepath, fermivel_2, tauk, tauk2, tau_avg)
      nscalar = nscalar+1
      nvector = nvector+0

      allocate(ipairs(nkpts_inp_all), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating ipairs'

      lpairs=.true.
      !calculate pairs or read pairs from file
!     inquire(file='kpoint_pairs',exist=lpairs)
!     if(lpairs)then
!       write(*,*) 'read pairs from file'
!       open(95,file='kpoint_pairs',form='formatted',action='read')
!       read(95,*) dummyline
!       read(95,'(I8)') ipairs
!       close(95)
!     else!lpairs
!       write(*,*) 'start to find pairs'
!       call find_kpoint_pairs(nkpts_inp_all, kpoints_3, ipairs)
!       write(*,*) 'pairs found'
!       open(95,file='kpoint_pairs',form='formatted',action='write')
!       write(95,'(A,I0)') 'nkpts_inp_all = ', nkpts_inp_all
!       write(95,'(I8)') ipairs
!       close(95)
!     end if!lpairs


    end if!lcond



    allocate( scalarstring(nscalar), vectorstring(nvector), STAT=ierr )
    if(ierr/=0) stop 'Problem allocating scalarstring etc.'

    !reset symmetries
    allocate(isym(1))
    nsym = 1
    isym = (/ 1 /)

    if(nvector>0) then
      allocate(vectordata(3,nkpts,nvector), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating vectordata'
    end if

    if(nscalar>0) then
      allocate(scalardata(nkpts_vis_all,nscalar), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating scalardata'
    end if

    iscalar=0
    ivector=0

    if(lfvel) then
      ivector = ivector+1
      if(iset_select==1) then
        do ikp=1,nkpts_vis_all
          vectordata(:,ikp,ivector) = fermivel(:,vis2int(ikp))
        end do!ikp
      else!iset_select==1
        vectordata(:,:,ivector) = fermivel(:,:)
      end if!iset_select==1
      vectorstring(ivector) = 'fvelocity'
    end if!lfvel

    if(lspin .and. inc%ndegen==2) then
      do isqa=1,nsqa
        iscalar = iscalar+1
        if(iset_select==1) then
          do ikp=1,nkpts_vis_all
            scalardata(ikp,iscalar) = (1d0-spinval(1,isqa,vis2int(ikp)))/2
          end do!ikp
        else!iset_select==1
          scalardata(:,iscalar) = (1d0-spinval(1,isqa,:))/2
        end if!iset_select==1
        write(scalarstring(iscalar),'(A,I0)') 'Eyaf_',isqa
      end do!isqa
    elseif(lspin .and. inc%ndegen==1) then
      do isqa=1,nsqa
        iscalar = iscalar+1
        if(iset_select==1) then
          do ikp=1,nkpts_vis_all
            scalardata(ikp,iscalar) = spinval(1,isqa,vis2int(ikp))/2
          end do!ikp
        else!iset_select==1
          scalardata(:,iscalar) = spinval(1,isqa,:)/2
        end if!iset_select==1
        write(scalarstring(iscalar),'(A,I0)') 'Spin_',isqa
      end do!isqa
    end if!lspin


    if(lcond)then
      allocate(dtmparr(nkpts_inp_all), dtmparr2(nkpts_inp_all), fvelabs(nkpts_inp_all), meanfreepath_sum(nkpts_inp_all), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating dtmparr'

      do ikp=1,nkpts_inp_all
        fvelabs(ikp)= sqrt(sum(fermivel_2(:,ikp)**2))
        meanfreepath_sum(ikp) = sum(meanfreepath(1,:,ikp))
      end do!ikp

!     iscalar = iscalar+1
!     dtmparr = fermivel_2(1,:)**2*tau_avg/fvelabs
!     write(scalarstring(iscalar),'(A)') 'RTA'
!     do ikp=1,nkpts_vis_all
!       scalardata(ikp,iscalar) = dtmparr(vis2int(ikp))
!     end do!ikp

!     iscalar = iscalar+1
!     dtmparr = fermivel_2(1,:)**2*tauk/fvelabs
!     write(scalarstring(iscalar),'(A)') 'BE_noscin'
!     do ikp=1,nkpts_vis_all
!       scalardata(ikp,iscalar) = dtmparr(vis2int(ikp))
!     end do!ikp

!     iscalar = iscalar+1
!     dtmparr = fermivel_2(1,:)*meanfreepath(1,:)/fvelabs
!     write(scalarstring(iscalar),'(A)') 's_xx'
!     do ikp=1,nkpts_vis_all
!       scalardata(ikp,iscalar) = dtmparr(vis2int(ikp))
!     end do!ikp

!     iscalar = iscalar+1
!     dtmparr = fermivel_2(2,:)*meanfreepath(1,:)/fvelabs
!     write(scalarstring(iscalar),'(A)') 's_xy'
!     do ikp=1,nkpts_vis_all
!       scalardata(ikp,iscalar) = dtmparr(vis2int(ikp))
!     end do!ikp

!     iscalar  = iscalar+1
      dtmparr  = 0d0
      dtmparr2 = 0d0
      dtmparr  = fermivel_2(2,:)*meanfreepath_sum(:)/fvelabs
      dtmparr2 = fermivel_2(1,:)*meanfreepath_sum(:)/fvelabs
!     write(scalarstring(iscalar)  ,'(A)') 's_xy'
!     write(scalarstring(iscalar+1),'(A)') 'diff_s_xy'
!     write(scalarstring(iscalar+2),'(A)') 'log(diff_s_xy)'
      write(scalarstring(iscalar+1),'(A)') 'tauk'
!     write(scalarstring(iscalar+2),'(A)') 'tauk2'
!     write(scalarstring(iscalar+3),'(A)') 's_xx'
      do ikp=1,nkpts_vis_all
!       dtmp =  dtmparr(vis2int(ikp))
!       scalardata(ikp,iscalar) = dtmp

!       dtmp = (dtmparr(vis2int(ikp))+dtmparr(ipairs(vis2int(ikp))))/2
!       scalardata(ikp,iscalar+1) = dtmp
!       scalardata(ikp,iscalar+2) = sign(log10(abs(dtmp)+1),dtmp)
        scalardata(ikp,iscalar+1) = tauk(vis2int(ikp))
!       scalardata(ikp,iscalar+2) = tauk2(vis2int(ikp))
!       scalardata(ikp,iscalar+3) = dtmparr2(vis2int(ikp))

      end do!ikp

!     ivector  = ivector+1
!     vectorstring(ivector) = 'meanfreepath'
!     do ikp=1,nkpts_vis_all
!       vectordata(:,ikp,ivector) = meanfreepath(:,1,vis2int(ikp))
!     end do!ikp

!     iscalar = iscalar+1
!     dtmparr = tauk
!     write(scalarstring(iscalar),'(A)') 'tau_k'
!     do ikp=1,nkpts_vis_all
!       scalardata(ikp,iscalar) = dtmparr(vis2int(ikp))
!     end do!ikp

!     iscalar = iscalar+1
!     dtmparr = tauk2
!     write(scalarstring(iscalar),'(A)') 'tau_k2'
!     do ikp=1,nkpts_vis_all
!       scalardata(ikp,iscalar) = dtmparr(vis2int(ikp))
!     end do!ikp

!     iscalar = iscalar+1
!     dtmparr = tauk-tauk2
!     write(scalarstring(iscalar),'(A)') 'tau_k-tau_k2'
!     do ikp=1,nkpts_vis_all
!       scalardata(ikp,iscalar) = dtmparr(vis2int(ikp))
!     end do!ikp

    end if!lcond



    write(filename,'(A)') 'testfile.vtp'
    call write_pointdata_rot( trim(filename),nkpts,kpoints,   &
                            & nscalar,scalardata,scalarstring,&
                            & nvector,vectordata,vectorstring,&
                            & nsym,symmetries%rotmat,isym     )


contains


  subroutine read_condfile(ndegen,nkptin, kpoints, meanfreepath, fermivel, tauk, tauk2, tau_avg)
    implicit none

    integer, intent(in) :: ndegen, nkptin
    double precision, intent(out) :: tau_avg
    double precision, allocatable, intent(out) :: kpoints(:,:), meanfreepath(:,:,:), fermivel(:,:),tauk(:), tauk2(:)

    integer :: ierr, ikp, nkpts
    double precision :: dtmp
    character(len=80) :: dummyline
    integer, parameter :: iou=56153

    allocate( kpoints(3,nkptin), meanfreepath(3,ndegen,nkptin),fermivel(3,nkptin),tauk(nkptin), tauk2(nkptin), STAT=ierr )
    if(ierr/=0) stop 'error allocating arrays in read_condfile'
    
    open(unit=iou,file='output_lifetime',form='formatted',action='read')

    read(iou,*) dummyline, nkpts
    write(*,*) 'check:', nkpts
    if(nkpts/=nkptin) stop 'error in nkpts in file output_lifetime'
    read(iou,*) dummyline
    read(iou,'(A15,ES25.16)') dummyline, tau_avg
    write(*,*) 'check:', tau_avg
    read(iou,*) dummyline

    do ikp=1,nkpts
      read(iou,*) kpoints(:,ikp), fermivel(:,ikp), tauk(ikp), tauk2(ikp), meanfreepath(:,:,ikp)
!     read(iou,*) kpoints(:,ikp), fermivel(:,ikp), dtmp, tauk(ikp)
    end do

    close(iou)
  end subroutine read_condfile


  subroutine find_kpoint_pairs(nkpts,kpoints,ipairs)
    !find a pair of two k-points that are (nearly) related by a mirror symmetry
    implicit none

    integer, intent(in) :: nkpts
    double precision, intent(in) :: kpoints(3,nkpts)
    integer, intent(out) :: ipairs(nkpts)

    integer :: ikp1, ikp2
    double precision :: kpick(3), diff, diffmin

    !init
    ipairs=0

    do ikp1=1,nkpts
      if(mod(ikp1,1000)==0) write(*,'(A,I8,A,I8)') ' .. doing ', ikp1, ' of ', nkpts
!     if(ipairs(ikp1)>0) cycle

      !pick a k-point
      kpick=kpoints(:,ikp1)

      !apply mirror symmetry
      kpick(2) = -kpick(2)

      !find the closest k-point
      diffmin=1e12
      do ikp2=1,nkpts

        diff=sum((kpoints(:,ikp2)-kpick)**2)
        if(diff<diffmin)then
          diffmin=diff
          ipairs(ikp1)=ikp2
        end if!diff<diffmin

      end do!ikp2
      write(93,*) sqrt(diffmin)

    end do!ikp1

    if(any(ipairs==0)) stop 'something went wrong in find_kpoint_pairs'

  end subroutine find_kpoint_pairs



end program visdata
