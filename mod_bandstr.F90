!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module mod_bandstr

  implicit none

  private
  public :: calc_bandstr

contains

  subroutine calc_bandstr

    use type_inc,  only: inc_type
    use type_data, only: lattice_type, cluster_type!, tgmatrx_type
    use mod_mympi, only: mympi_init, myrank, nranks, master
    use mod_read,  only: read_inc
    use mod_parutils, only: distribute_linear_on_tasks
#ifdef CPP_MPI
    use mpi
#endif

    implicit none

    type(inc_type)      :: inc
    type(lattice_type)  :: lattice
    type(cluster_type)  :: cluster
!   type(tgmatrx_type)  :: tgmatrx

    integer :: ietot, nkpts, nept
    double precision, allocatable :: kpoints(:,:)
    double complex,   allocatable :: mineigw(:,:), mineigw_all(:,:)

    integer :: irank, ie, ik, ierr
    character(len=80) :: fmtstr, filename
    integer, allocatable :: nepts(:), ioff(:), recvcnt(:)
    integer, parameter :: ifile=135


    !initialize MPI
#ifdef CPP_MPI
    call MPI_Init ( ierr )
#endif
    call mympi_init()

    !Read in TBKKR-header
    call read_inc(inc)

    !parallelize
    ietot = get_nFiles()
    allocate(nepts(0:nranks-1), recvcnt(0:nranks-1), ioff(0:nranks-1), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating nepts'

    call distribute_linear_on_tasks(nranks, myrank, master, ietot, nepts, ioff, .true.)

    !read TBkkr-data
    call read_TBkkr_data(inc, lattice, cluster)

    call read_kpath(nkpts, kpoints)

    allocate( mineigw(nkpts,nepts(myrank)), STAT=ierr )
    if(ierr/=0) stop 'Problem allocating mineigw'

    call get_mineigw(inc, nepts(myrank), ioff(myrank), lattice, cluster, nkpts, kpoints, mineigw)

    recvcnt = nkpts*nepts
    ioff = 0
    do irank=1,nranks-1
      ioff(irank) = sum(recvcnt(0:irank-1))
    end do

    allocate( mineigw_all(nkpts,ietot), STAT=ierr )
    if(ierr/=0) stop 'Problem allocating mineigw_all'

    call MPI_Gatherv( mineigw, nkpts*nepts(myrank), MPI_DOUBLE_COMPLEX,&
                    & mineigw_all, recvcnt, ioff, MPI_DOUBLE_COMPLEX,&
                    & master, MPI_COMM_WORLD, ierr                  )
    if(ierr/=MPI_SUCCESS) stop 'Problem in gather mineigw'

    if(myrank==master)then
      open(unit=ifile,file='mineigw.txt',form='formatted',action='write')
      write(fmtstr,'("(",I0,"ES25.16)")') nkpts
        write(ifile,fmtstr) abs(mineigw_all)
      close(ifile)
    end if!myrank==master

    call MPI_Finalize(ierr)

  end subroutine calc_bandstr





  subroutine get_mineigw(inc, nept, ieoff, lattice, cluster, nkpt, kpoints, mineigw)

    use type_inc,  only: inc_type
    use type_data, only: lattice_type, cluster_type!, tgmatrx_type
    use mod_kkrmat, only: compute_kkr_eigenvectors
    use mod_mathtools, only: bubblesort
    use mod_mympi, only: myrank
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    integer,            intent(in) :: nept, ieoff, nkpt
    double precision,   intent(in) :: kpoints(3,nkpt)
    double complex,     intent(out) :: mineigw(nkpt,nept)

    integer :: ie, ik, indsrt(inc%almso)
    double precision :: kpoint(3), dtmp(inc%almso)
    double complex :: LVeig(inc%almso,inc%almso),&
                    & RVeig(inc%almso,inc%almso),&
                    & eigw(inc%almso)

    double complex :: energy,&
                    & tinvll(inc%lmmaxso,inc%lmmaxso,inc%naezd),&
                    & ginp(inc%naclsd*inc%lmmax,inc%lmmax,inc%nclsd)
    do ie=1,nept
     write(*,*) '  do ie=',ie
     call read_TBkkr_tgmatrix(ie+ieoff, inc, lattice, energy, tinvll, ginp)

     do ik=1,nkpt
      kpoint = kpoints(:,ik)
      call compute_kkr_eigenvectors(inc, lattice, cluster, ginp, tinvll, kpoint, eigw, LVeig, RVeig)
      dtmp=abs(eigw)
      call bubblesort(inc%almso, dtmp, indsrt)
      mineigw(ik,ie) = eigw(indsrt(1))
     end do!ik

    end do!ie

  end subroutine get_mineigw





  subroutine read_kpath(nkpt_tot, kpoints)

    use mod_mympi, only: myrank, nranks, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    integer, intent(out) :: nkpt_tot
    double precision, allocatable, intent(out) :: kpoints(:,:)

    integer :: npaths, nkpt_left
    integer,          allocatable :: nkptpath(:)
    double precision, allocatable :: kpoints_main(:,:)


    integer :: ipath, ii, ikp, ierr
    double precision :: deltaK(3)
    

    if(myrank==master) then

      write(*,*) 'Begin to calculate bands.'
      open(123,file='kpts_in',form='formatted',status='old', iostat=ierr)
      if(ierr/=0) then
        open(124,file='kpts_in',form='formatted',status='new', iostat=ierr)
        write(124,"((I8),(A))") 2, ' <- n = number of paths, followed by (n+1) lines.'
        write(124,"((3ES18.9),(I8))") (/ 0d0, 0d0, 0d0 /), 100
        write(124,"((3ES18.9),(I8))") (/ 0d0, 0d0, 0d0 /), 100
        write(124,"((3ES18.9),(A8))") (/ 1d0, 1d0, 1d0 /), '     END'
        stop 'Error opening "kpts_in" file. Dummy file created.'
      end if

      ! read in of k-points, which shall be interpolated between
      ! array 'kpoints_main' is filled with the main k-points, the last k-point being in entry (npaths+1),
      ! meaning for the path 'i':
      !                            starting k-point: kpoints_main(i)
      !                              ending k-point: kpoints_main(i+1)
      !                            number intervals: nkptpath(i)
      read(123,"((I8))") npaths
      write(*,"((A),(I8),(A))") 'Read in ', npaths, ' paths.'

    end if!master


    call MPI_Bcast(npaths, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) STOP 'band: nkpt'

    allocate( kpoints_main(3,npaths+1), nkptpath(npaths), STAT=ierr )
    if(ierr/=0) stop 'Problem allocating kpoints_main'
    kpoints_main = 0d0
    nkptpath = 0

    if(myrank==master) then
      do ikp=1,npaths
        read (123,"((3ES18.9),(I8))") kpoints_main(1:3,ikp), nkptpath(ikp)
      end do!ikp
      read (123,"((3ES18.9))") kpoints_main(1:3,npaths+1)
    end if!myrank==master


    call MPI_Bcast(nkptpath, npaths, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) STOP 'band: nkptsub'
    call MPI_Bcast(kpoints_main, 3*(npaths+1), MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) STOP 'band: kpoints_main'

    ! create k-path
    nkpt_tot = sum(nkptpath)+1
    allocate(kpoints(3,nkpt_tot), STAT=ierr)
    if(ierr/=0) stop 'Problem alloc. kpoints'

    ii=1
    do ipath=1,npaths
      deltaK = (kpoints_main(:,ipath+1)-kpoints_main(:,ipath))/nkptpath(ipath)
      do ikp = 1,nkptpath(ipath)
        kpoints(:,ii) = kpoints_main(:,ipath) + deltaK*(ikp-1)
        ii=ii+1
      end do!ikp
    end do!ipath
    kpoints(:,ii) = kpoints_main(:,npaths+1)

    !output the k-points
    if(myrank==master) then
      open(329,file='kpts_out',form='formatted',iostat=ierr)
      write(329,"((A))") "# ikp, k_x, k_y, k_z"
      do ii=1,nkpt_tot
        write(329,"((I8),(3ES18.9))") ii, kpoints(:,ii)
      end do!ii=1,nkpt_tot
      close(329)
    end if!myrank==master



  end subroutine read_kpath





  subroutine read_TBkkr_data(inc, lattice, cluster)
  ! This subroutine reads the TBKKR-datafile,
  ! containing the main variables from the TB-KKR calculation.

    use type_inc,     only: inc_type
    use type_data,    only: lattice_type, lattice_init,&
                          & cluster_type, cluster_init
    use mod_mympi,    only: myrank, nranks, master
    use mod_kkrmat,   only: InvertT
    use mod_ioformat, only: fmt_fn_ext, ext_formatted, filename_tbkkrcontainer
#ifdef CPP_MPI
    use type_data,    only: create_mpimask_lattice, create_mpimask_cluster
    use mpi
#endif
    implicit none

      type(inc_type),     intent(in)  :: inc
      type(lattice_type), intent(out) :: lattice
      type(cluster_type), intent(out) :: cluster

      integer            :: i1, i2, ie, lm1, lm2, isigma, ierr, idummy, iversion
      integer            :: mask_lattice, mask_cluster
      character(len=1)   :: strdum
      character(len=80)  :: filename, strline
      integer, parameter :: ifile= 12, ifilegreen=62
      double complex     :: zdummy


      call lattice_init(inc,lattice)
      call cluster_init(inc,cluster)


      if(myrank==master) then

        write(filename,fmt_fn_ext) filename_tbkkrcontainer, ext_formatted
        open(unit=ifile, file=trim(filename), form='formatted', action='read')

        !=== get file version ===!
        iversion = 1
        read(ifile,'(A)') strline
        if(strline(1:13) == '#FILEVERSION=')then
          read(strline(14:80),*) iversion
        else
          rewind(ifile)
        end if!

        !=== lattice information ===!
        read(ifile,'(A)') strdum
        read(ifile,'(ES25.16)') lattice%alat

        read(ifile,'(A)') strdum
        read(ifile,'(3ES25.16)') ((lattice%bravais(i1,i2),i1=1,3),i2=1,3)
        if( all(abs(lattice%bravais(:,3))<1d-6) .and. all(abs(lattice%bravais(3,:))<1d-6) ) lattice%bravais(3,3) = 1d16

        read(ifile,'(A)') strdum
        read(ifile,'(3ES25.16)') ((lattice%recbv(i1,i2),i1=1,3),i2=1,3)
        if( all(abs(lattice%recbv(:,3))<1d-6) .and. all(abs(lattice%recbv(3,:))<1d-6) ) lattice%recbv(3,3) = 1d0

        read(ifile,'(A)') strdum
        read(ifile,'(3ES25.16)') ((lattice%rbasis(i2,i1), i2=1,3), i1=1,inc%naezd+inc%nembd)


        !=== cluster information ===!
        read(ifile,'(A)') strdum
        read(ifile,'(1I8)') (cluster%cls(i1),i1=1,inc%natypd)

        read(ifile,'(A)') strdum
        read(ifile,'(1I8)') (cluster%nacls(i1), i1=1,inc%nclsd)

        read(ifile,'(A)') strdum
        do i1=1,inc%nclsd
         do i2=1,inc%naclsd
           read(ifile,'(3ES25.16)') cluster%rcls(:,i2,i1)
         end do
        end do

        if(iversion==1)then!backwards-compatibility to old file version of <2014
          read(ifile,'(A)') strdum
!         read(ifile,'(1I8)') (cluster%eqinv(i1),i1=1,inc%naezd)
          read(ifile,'(1I8)') (idummy, i1=1,inc%naezd)
        end if!iversion==1

        read(ifile,'(A)') strdum
        read(ifile,'(1I8)') ((cluster%ezoa(i1,i2),i1=1,inc%naclsd),i2=1,inc%naezd)

        read(ifile,'(A)') strdum
        read(ifile,'(1I8)') ((cluster%atom(i1,i2),i1=1,inc%naclsd),i2=1,inc%naezd)

        if(iversion==1)then!backwards-compatibility to old file version of <2014
          read(ifile,'(A)') strdum
!         read(ifile,'(1I8)') (cluster%kaoez(i1),   i1=1,inc%naezd+inc%nembd)
          read(ifile,'(1I8)') (idummy, i1=1,inc%naezd+inc%nembd)
        end if!iversion==1

        read(ifile,'(A)') strdum
        do i1=0,inc%nrd
          read(ifile,'(3ES25.16)') cluster%rr(:,i1)
        end do

        close(ifile)

      end if!myrank==master


#ifdef CPP_MPI
      ! Brodcast lattice information
      call create_mpimask_lattice(lattice,mask_lattice)
      call MPI_Type_commit(mask_lattice, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error commiting mask for lattice'
      call MPI_Bcast(lattice%N, 1, mask_lattice, master, MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error brodcasting inc'
      call MPI_Type_free(mask_lattice, ierr)

      ! Brodcast cluster information
      call create_mpimask_cluster(cluster,mask_cluster)
      call MPI_Type_commit(mask_cluster, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error commiting mask for cluster'
      call MPI_Bcast(cluster%N, 1, mask_cluster, master, MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error brodcasting inc'
      call MPI_Type_free(mask_cluster, ierr)
#endif

  end subroutine read_TBkkr_data



  subroutine read_TBkkr_tgmatrix(ie, inc, lattice, energy, tinvll, ginp)

    use type_inc, only: inc_type
    use type_data, only: lattice_type
    use mod_kkrmat, only: InvertT
    implicit none
    type(inc_type),     intent(in)  :: inc
    type(lattice_type), intent(in)  :: lattice
    integer,            intent(in)  :: ie
    double complex,     intent(out) :: energy,&
                                     & tinvll(inc%lmmaxso,inc%lmmaxso,inc%naezd),&
                                     & ginp(inc%naclsd*inc%lmmax,inc%lmmax,inc%nclsd)

    double complex :: tmatll(inc%lmmax*inc%nspd,inc%lmmax*inc%nspd,inc%naezd)

    character(len=1)   :: strdum
    character(len=256) :: filename
    integer :: i1, lm2, lm1
    integer, parameter :: ifilegreen=1354

    write(filename,'(A,I6.6,A)') 'TBkkr_green_band_',ie,'.txt'
    open(unit=ifilegreen,file=filename,form='formatted', action='read')

    read(ifilegreen,'(A)') strdum
    read(ifilegreen,'(2ES25.16)') energy

    ! Read in (energydependent) T-Matrix
    read(ifilegreen,'(A)') strdum
    do i1=1,inc%naezd
      do lm2=1,inc%lmmax*inc%nspd
        do lm1=1,inc%lmmax*inc%nspd
          read(ifilegreen,'(2ES25.16)') tmatll(lm1,lm2,i1)
        end do!lm1
      end do!lm2
    end do!i1

    ! Read in (energydependent) Greens function
    read(ifilegreen,'(A)') strdum
    do i1=1,inc%nclsd
      do lm2=1,inc%lmmax
        do lm1=1,inc%lmmax*inc%naclsd
          read(ifilegreen,'(2ES25.16)') ginp(lm1,lm2,i1)
        end do!lm1
      end do!lm2
    end do!i1

    close(ifilegreen)

    !perform handy modifications on t-matrix-arrays for later use in the code
    call InvertT(inc, lattice%alat, tmatll, tinvll)

 end subroutine read_TBkkr_tgmatrix



  integer function get_nFiles()

    use mod_ioformat, only: filename_outinfo, ext_formatted

    implicit none

    integer :: nFiles
    logical :: fileexists
    character(len=80) :: filename
    nFiles=0
    fileexists = .TRUE.
    do while (fileexists)
      write(filename,'(A,I6.6,A)') 'TBkkr_green_band_',nFiles+1,'.txt'
      inquire(file=filename, exist=fileexists)
      nFiles = nFiles+1
    end do
    get_nFiles = nFiles-1
    write(*,'(I0,A)') get_nFiles, ' energyfiles found'
  end function get_nFiles


end module mod_bandstr
