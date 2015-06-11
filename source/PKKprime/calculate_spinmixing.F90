program calcspinmix

  use type_inc,       only: inc_type
  use type_data,      only: lattice_type, cluster_type, tgmatrx_type
  use mod_mympi,      only: mympi_init, myrank, nranks, master
  use mod_symmetries, only: symmetries_type, set_symmetries, rotate_kpoints, expand_areas, expand_visarrays, expand_spinvalues
  use mod_ioformat,   only: filemode_int, filemode_vis, fmt_fn_sub_ext, ext_formatted, filename_spin, filename_fvel, filename_torq
  use mod_iohelp,     only: getBZvolume
  use mod_read,       only: read_inc, read_TBkkrdata, read_kpointsfile_vis, read_kpointsfile_int, read_weights, read_fermivelocity, read_spinvalue, read_torqvalue
  use mod_calconfs,   only: calculate_spinmixing_int, calculate_spinmixing_vis, calculate_torkance_CRTA_int, calculate_torkance_CRTA_vis, calculate_dos_int

#ifdef CPP_MPI
    use mpi
#endif

    implicit none

    integer            :: iter
    logical            :: l_spinfile, l_fvelfile, l_torqfile
    character(len=256) :: filemode, filename

    type(inc_type)      :: inc
    type(lattice_type)  :: lattice
    type(cluster_type)  :: cluster
    type(tgmatrx_type)  :: tgmatrx

    !symmetry arrays
    integer :: nsym
    integer, allocatable :: isym(:)
    type(symmetries_type) :: symmetries

    !local k-point arrays
    integer :: nkpts, nkpts_all
    integer, allocatable :: kpt2irr(:), irr2kpt(:)
    double precision, allocatable :: kpoints(:,:), areas(:), weights(:), fermivel(:,:)

    !spinvalue locals
    integer              :: nsqa, ndegen1
    integer, allocatable :: ispincomb(:)
    double precision, allocatable :: nvect(:,:), spinval(:,:,:), spinval1(:,:,:), spinmix(:)

    !torqvalue locals
    integer                       :: itmp
    double precision, allocatable :: torqval(:,:,:), arrtmp1(:,:), arrtmp2(:,:)

    !temp k-point arrays
    integer :: nkpts1, nkpts2, nkpts_all1
    integer,          allocatable :: kpt2irr1(:), irr2kpt1(:)
    double precision, allocatable :: areas1(:), weights1(:), kpoints1(:,:), fermivel1(:,:)

    double precision :: pi, BZVol, dos
    integer             :: ierr

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

!    write(*,*) 'before loop'
!    write(*,*) 'nBZdim=', inc%nBZdim
!    write(*,*) 'BZvol=', BZVol
!    write(*,*) 'recbv(:,3)=', lattice%recbv(:,3)

    do iter=1,2

      if(iter==1) filemode = filemode_int
      if(iter==2) filemode = filemode_vis
      write(*,'(A,I0,A,A)') 'loop iter= ', iter, ' , filemode= ', filemode

      write(filename,fmt_fn_sub_ext) filename_spin, trim(filemode), ext_formatted
      inquire(file=filename, exist=l_spinfile)

      write(filename,fmt_fn_sub_ext) filename_fvel, trim(filemode), ext_formatted
      inquire(file=filename, exist=l_fvelfile)

      write(filename,fmt_fn_sub_ext) filename_torq, trim(filemode), ext_formatted
      inquire(file=filename, exist=l_torqfile)




      !=============================!
      != Read in the k-point files =!
      if(iter==1)then
        call read_kpointsfile_int(nkpts1, kpoints1, areas1, nsym, isym)
        call rotate_kpoints(symmetries%rotmat, nkpts1, kpoints1, nsym, isym, nkpts, kpoints)
        call expand_areas(nsym,nkpts1,areas1,areas)
        deallocate(isym, areas1, kpoints1)
      else
        call read_kpointsfile_vis(nkpts_all1, nkpts1, kpoints1, nsym, isym, kpt2irr1, irr2kpt1)
        call rotate_kpoints(symmetries%rotmat, nkpts1, kpoints1, nsym, isym, nkpts, kpoints)
        call expand_visarrays(nsym, nkpts_all1, nkpts1, kpt2irr1, irr2kpt1, kpt2irr, irr2kpt)
        nkpts_all = nsym*nkpts_all1
        deallocate(kpt2irr1, irr2kpt1, kpoints1, isym)
      end if!iter==1

      !==============================!
      != Read in the fermi velocity =!
      if(l_fvelfile) then
        call read_fermivelocity(trim(filemode), nkpts1, fermivel1, nsym, isym)
        call rotate_kpoints(symmetries%rotmat, nkpts1, fermivel1, nsym, isym, nkpts2, fermivel)
        if(nkpts2/=nkpts) stop 'inconsistency in number of k-points'
        deallocate(fermivel1, isym)
      end if!l_fvelfile



      if(l_spinfile) then
        !=========================!
        != Read in the spinvalue =!
        call read_spinvalue(trim(filemode), nkpts1, nsqa, ndegen1, ispincomb, nvect, spinval1, nsym, isym)
        if(nkpts1*nsym/=nkpts) stop 'inconsistency in number of k-points'
        if(ndegen1/=inc%ndegen) stop 'inconsistency in ndegen'
        call expand_spinvalues(nsym,ndegen1,nsqa,nkpts1,spinval1,spinval)
        deallocate(spinval1, isym)
      end if!l_spinfile


      !=============================!
      != Read in the torque values =!
      if(l_torqfile) then
        call read_torqvalue(trim(filemode), nkpts1, ndegen1, torqval, nsym, isym)

        !next apply symmetries to torqval; as a first step, flatten the array
        itmp = size(torqval)/3
        allocate(arrtmp1(3,itmp), STAT=ierr)
        if(ierr/=0) stop 'problem alloaction arrtmp1'
        arrtmp1 = reshape(torqval, (/3, itmp/))

        deallocate(torqval)

        !now apply the symmetries
        call rotate_kpoints(symmetries%rotmat, nkpts1*ndegen1, arrtmp1, nsym, isym, nkpts2, arrtmp2)
        if(nkpts1*ndegen1*nsym /= nkpts2) stop 'nkpts1*ndegen1*nsym /= nkpts2'

        !transform back to old shape
        allocate(torqval(3,ndegen1,nkpts1*nsym), STAT=ierr)
        if(ierr/=0) stop 'problem allocating torqval'
        torqval = reshape(arrtmp2,(/3,ndegen1,nkpts1*nsym/))

        deallocate(arrtmp1, arrtmp2, isym)
        if(nkpts1*nsym/=nkpts) stop 'inconsistency in number of k-points'

      end if!l_torqfile



      !++++++++++++++++++++++++++++++++++++++++++++
      !++++ Perform Fermi surface integrations ++++

      ! Spin-mixing parameter
      if(l_spinfile .and. l_fvelfile) then
        nsym=1
        allocate(spinmix(nsqa), STAT=ierr)
        if(ierr/=0) stop 'Problem allocating spinmix'

        if(iter==1) call calculate_spinmixing_int(nsym,inc%ndegen,nsqa,nkpts,areas,fermivel,spinval,spinmix,dos,.true.,BZVol)
        if(iter==2) call calculate_spinmixing_vis(nsym,inc%ndegen,nsqa,nkpts,nkpts_all,kpt2irr,kpoints,fermivel,spinval,spinmix,dos,.true.,BZVol,inc%nBZdim)
        deallocate(spinmix)

      end if!l_spinfile .and. l_fvelfile


      allocate(isym(1))
      nsym    = 1
      isym(1) = 1

      if(l_torqfile .and. l_fvelfile) then
        if(iter==1)  call calculate_torkance_CRTA_int(nsym,isym,symmetries%rotmat,lattice%alat,BZVol,inc%ndegen,nkpts,areas,fermivel,torqval)
        if(iter==2)  call calculate_torkance_CRTA_vis(inc%nBZdim,nsym,isym,symmetries%rotmat,lattice%alat,BZVol,inc%ndegen,nkpts,nkpts_all,kpt2irr,irr2kpt,kpoints,fermivel,torqval)
      end if!l_torqfile.and.l_fvelfile

      if(l_fvelfile) then
        if(iter==1)  call calculate_dos_int(nsym,isym,symmetries%rotmat,lattice%alat,BZVol,nkpts,areas,fermivel,inc%nBZdim)
      end if!l_fvelfile


    end do!iter

end program calcspinmix
