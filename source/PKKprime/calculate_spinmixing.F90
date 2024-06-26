!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


program calcspinmix

  use type_inc,       only: inc_type
  use type_data,      only: lattice_type, cluster_type, tgmatrx_type
  use mod_mympi,      only: mympi_init, myrank, nranks, master
  use mod_symmetries, only: symmetries_type, set_symmetries, rotate_kpoints, expand_areas, expand_visarrays, expand_spinvalues
  use mod_ioformat,   only: filemode_int, filemode_vis, fmt_fn_sub_ext, fmt_fn_atom_sub_ext, ext_formatted, filename_spin, filename_fvel, filename_torq, filename_spinflux, filename_spinvec
  use mod_iohelp,     only: getBZvolume
  use mod_read,       only: read_inc, read_TBkkrdata, read_kpointsfile_vis, read_kpointsfile_int, read_weights, read_fermivelocity, read_spinvalue, read_spinvec_atom, read_torqvalue, read_torqvalue_atom, read_spinflux_atom
  use mod_calconfs,   only: calculate_spinmixing_int, calculate_spinmixing_vis, calculate_response_functions_CRTA_int, calculate_response_functions_CRTA_vis, calculate_dos_int

#ifdef CPP_MPI
    use mpi
#endif

    implicit none

    integer            :: iter
    logical            :: l_spinfile, l_fvelfile, l_torqfile, l_torqfile_atom, l_spinvecfile_atom, l_spinfluxfile_atom
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
    integer                       :: itmp, iat
    double precision, allocatable :: torqval(:,:,:), torqval_atom(:,:,:,:), spinvec_atom(:,:,:,:), spinflux_atom(:,:,:,:)
    double precision, allocatable :: arrtmp1(:,:), arrtmp2(:,:)

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

      ! l_torqfile_atom is true only if files are found for all atoms
      l_torqfile_atom = .true.
      do iat=1,inc%natypd
        write(filename,fmt_fn_atom_sub_ext) filename_torq, iat, trim(filemode), ext_formatted
        inquire(file=filename, exist=l_torqfile_atom)
        if(l_torqfile_atom==.false.) exit
      end do!iat=1,inc%natypd

      ! l_spinvecfile_atom is true only if files are found for all atoms
      l_spinvecfile_atom = .true.
      do iat=1,inc%natypd
        write(filename,fmt_fn_atom_sub_ext) filename_spinvec, iat, trim(filemode), ext_formatted
        inquire(file=filename, exist=l_spinvecfile_atom)
        if(l_spinvecfile_atom==.false.) exit
      end do!iat=1,inc%natypd

      ! l_spinfluxfile_atom is true only if files are found for all atoms
      l_spinfluxfile_atom = .true.
      do iat=1,inc%natypd
        write(filename,fmt_fn_atom_sub_ext) filename_spinflux, iat, trim(filemode), ext_formatted
        inquire(file=filename, exist=l_spinfluxfile_atom)
        if(l_spinfluxfile_atom==.false.) exit
      end do!iat=1,inc%natypd


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

      if(l_torqfile) then
        !=============================!
        != Read in the torque values =!
        call read_torqvalue(trim(filemode), nkpts1, ndegen1, torqval, nsym, isym)

        !next apply symmetries to torqval; as a first step, flatten the array
        itmp = size(torqval)/3
        allocate(arrtmp1(3,itmp), STAT=ierr)
        if(ierr/=0) stop 'problem alloaction arrtmp1'
        arrtmp1 = reshape(torqval, (/3, itmp/))

        deallocate(torqval)

        !now apply the symmetries
        call rotate_kpoints(symmetries%rotmat, nkpts1*ndegen1, arrtmp1, nsym, isym, nkpts2, arrtmp2)
        if(nkpts2 /= nkpts) stop 'nkpts2 /= nkpts'

        !transform back to old shape
        allocate(torqval(3,ndegen1,nkpts1*nsym), STAT=ierr)
        if(ierr/=0) stop 'problem allocating torqval'
        torqval = reshape(arrtmp2,(/3,ndegen1,nkpts1*nsym/))

        deallocate(arrtmp1, arrtmp2, isym)

      end if!l_torqfile

      if(l_torqfile_atom) then
        !===========================================!
        != Read in the torque values for each atom =!
        call read_torqvalue_atom(trim(filemode), inc%natypd, nkpts1, ndegen1, torqval_atom, nsym, isym)

        !next apply symmetries to torqval_atom; as a first step, flatten the array
        itmp = size(torqval_atom)/3
        allocate(arrtmp1(3,itmp), STAT=ierr)
        if(ierr/=0) stop 'problem alloaction arrtmp1'
        arrtmp1 = reshape(torqval_atom, (/3, itmp/))

        deallocate(torqval_atom)

        !now apply the symmetries
        call rotate_kpoints(symmetries%rotmat, inc%natypd*nkpts1*ndegen1, arrtmp1, nsym, isym, nkpts2, arrtmp2)
        if(nkpts2 /= nkpts*inc%natypd) stop 'nkpts2 /= nkpts*natpyd'

        !transform back to old shape
        allocate(torqval_atom(3,inc%natypd,ndegen1,nkpts1*nsym), STAT=ierr)
        if(ierr/=0) stop 'problem allocating torqval_atom'
        torqval_atom = reshape(arrtmp2,(/3,inc%natypd,ndegen1,nkpts1*nsym/))

        deallocate(arrtmp1, arrtmp2, isym)

      end if!l_torqval_atom

      if(l_spinvecfile_atom) then
        !===========================================!
        != Read in the torque values for each atom =!
        call read_spinvec_atom(trim(filemode), inc%natypd, nkpts1, ndegen1, spinvec_atom, nsym, isym)

        !next apply symmetries to spinvec_atom; as a first step, flatten the array
        itmp = size(spinvec_atom)/3
        allocate(arrtmp1(3,itmp), STAT=ierr)
        if(ierr/=0) stop 'problem alloaction arrtmp1'
        arrtmp1 = reshape(spinvec_atom, (/3, itmp/))

        deallocate(spinvec_atom)

        !now apply the symmetries
        call rotate_kpoints(symmetries%rotmat, inc%natypd*nkpts1*ndegen1, arrtmp1, nsym, isym, nkpts2, arrtmp2)
        if(nkpts2 /= nkpts*inc%natypd) stop 'nkpts2 /= nkpts*natpyd'

        !transform back to old shape
        allocate(spinvec_atom(3,inc%natypd,ndegen1,nkpts1*nsym), STAT=ierr)
        if(ierr/=0) stop 'problem allocating spinvec_atom'
        spinvec_atom = reshape(arrtmp2,(/3,inc%natypd,ndegen1,nkpts1*nsym/))

        deallocate(arrtmp1, arrtmp2, isym)

      end if!l_spinvec_atom

      if(l_spinfluxfile_atom) then
        !=============================================!
        != Read in the spinflux values for each atom =!
        call read_spinflux_atom(trim(filemode), inc%natypd, nkpts1, ndegen1, spinflux_atom, nsym, isym)

        !next apply symmetries to spinflux_atom; as a first step, flatten the array
        itmp = size(spinflux_atom)/3
        allocate(arrtmp1(3,itmp), STAT=ierr)
        if(ierr/=0) stop 'problem alloaction arrtmp1'
        arrtmp1 = reshape(spinflux_atom, (/3, itmp/))

        deallocate(spinflux_atom)

        !now apply the symmetries
        call rotate_kpoints(symmetries%rotmat, inc%natypd*nkpts1*ndegen1, arrtmp1, nsym, isym, nkpts2, arrtmp2)
        if(nkpts2 /= nkpts*inc%natypd) stop 'nkpts2 /= nkpts*natpyd'

        !transform back to old shape
        allocate(spinflux_atom(3,inc%natypd,ndegen1,nkpts1*nsym), STAT=ierr)
        if(ierr/=0) stop 'problem allocating spinflux_atom'
        spinflux_atom = reshape(arrtmp2,(/3,inc%natypd,ndegen1,nkpts1*nsym/))

        deallocate(arrtmp1, arrtmp2, isym)

      end if!l_spinfluxfile_atom

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
        if(iter==1)  call calculate_response_functions_CRTA_int(nsym,isym,symmetries%rotmat,lattice%alat,BZVol,inc%ndegen,inc%natypd,nkpts,areas,fermivel,torqval,torqval_atom,spinvec_atom,spinflux_atom)

        if(iter==2)  call calculate_response_functions_CRTA_vis(inc%nBZdim,nsym,isym,symmetries%rotmat,lattice%alat,BZVol,inc%ndegen,inc%natypd,nkpts,nkpts_all,kpt2irr,irr2kpt,kpoints,fermivel,torqval,torqval_atom,spinvec_atom,spinflux_atom)
      end if!l_torqfile.and.l_fvelfile

      if(l_fvelfile) then
        if(iter==1)  call calculate_dos_int(nsym,isym,symmetries%rotmat,lattice%alat,BZVol,nkpts,areas,fermivel,inc%nBZdim)
      end if!l_fvelfile


    end do!iter

end program calcspinmix
