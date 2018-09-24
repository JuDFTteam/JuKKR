!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module mod_read

  implicit none

  private
  public :: read_inc, read_TBkkrdata, read_kpointsfile_vis, read_kpointsfile_int, read_weights, read_fermivelocity, read_spinvalue, read_torqvalue, read_torqvalue_atom, read_spinvec_atom, read_spinflux_atom


contains


  subroutine read_inc(inc)
  ! This subroutine reads the inc-file,
  ! containing the main parameters from the TB-KKR calculation.

    use type_inc,     only: inc_type
    use mod_mympi,    only: myrank, nranks, master
    use mod_ioformat, only: fmt_fn_ext, ext_formatted, filename_tbkkrparams,&
                          & filename_tbkkrrhod, filename_tbkkrtorq, filename_tbkkrspinflux,&
                          & filename_tbkkralpha
    use mod_ioinput,  only: IoInput
#ifdef CPP_MPI
    use type_inc,     only: create_mpimask_inc
    use mpi
#endif

    implicit none

      type(inc_type), intent(inout) :: inc

      integer            :: mask_inc, ierr, iversion, korbit, nspind
      character(len=80)  :: filename, uio, strline
      integer, parameter :: ifile= 11


      if(myrank==master) then
        write(filename,fmt_fn_ext) filename_tbkkrparams, ext_formatted
        open(unit=ifile, file=trim(filename), form='formatted', action='read')

        iversion = 1
        read(ifile,'(A)') strline
        if(strline(1:13) == '#FILEVERSION=')then
          read(strline(14:80),*) iversion
        else
          rewind(ifile)
        end if!

        select case (iversion)
        case(1)
          read(ifile,'(2I8)') inc%lmaxd, inc%lmax
          read(ifile,'(3I8)') inc%nspd, inc%nspo, inc%nspoh
          read(ifile,'(1I8)') inc%nrd
          read(ifile,'(2I8)') inc%nembd, inc%nemb
          read(ifile,'(2I8)') inc%nclsd, inc%ncls
          read(ifile,'(2I8)') inc%natypd, inc%natyp
          read(ifile,'(2I8)') inc%naezd, inc%naez
          read(ifile,'(1I8)') inc%naclsd
          read(ifile,'(1I8)') inc%ielast
          read(ifile,'(1I8)') inc%ins
        case(2)
          read(ifile,'(1I8)') inc%lmaxd
          read(ifile,'(1I8)') inc%lmax
          read(ifile,'(1I8)') korbit
          read(ifile,'(1I8)') nspind
          read(ifile,'(1I8)') inc%nrd
          read(ifile,'(1I8)') inc%nembd
          read(ifile,'(1I8)') inc%nemb
          read(ifile,'(1I8)') inc%nclsd
          read(ifile,'(1I8)') inc%ncls
          read(ifile,'(1I8)') inc%natypd
          read(ifile,'(1I8)') inc%natyp
          read(ifile,'(1I8)') inc%naezd
          read(ifile,'(1I8)') inc%naez
          read(ifile,'(1I8)') inc%naclsd
          read(ifile,'(1I8)') inc%ielast
          read(ifile,'(1I8)') inc%ins
          inc%nspd = max(nspind,korbit+1)
          inc%nspo  = korbit+1
          inc%nspoh = korbit+1
        case default; stop 'unknown fileversion for inc-file'
        end select

        close(ifile)

        inc%lmmax   = (inc%lmaxd+1)**2
        inc%lmmaxso = inc%lmmax*inc%nspd
        inc%alm     = inc%naezd*inc%lmmax
        inc%almso   = inc%alm*inc%nspd

        call IoInput('NDEGEN    ',uio,1,7,ierr)
        read(unit=uio,fmt=*) inc%ndegen

        call IoInput('NBZDIM    ',uio,1,7,ierr)
        read(unit=uio,fmt=*) inc%nBZdim

        call IoInput('NROOTMAX  ',uio,1,7,ierr)
        read(unit=uio,fmt=*) inc%nrootmax
        if(inc%nrootmax==0) inc%nrootmax = inc%almso/inc%ndegen

        call IoInput('MEMOPT    ',uio,1,7,ierr)

        if(ierr==0) then
          read(unit=uio,fmt=*) inc%memopt
        else ! ensures compatibility of old format input files
          write(*,*) "Warning : MEMOPT set to .false. by default !"
          inc%memopt = .false.
        end if!ierr/=0

        if(inc%memopt==.true.) then
          call IoInput('NEIG      ',uio,1,7,ierr)
          read(unit=uio,fmt=*) inc%neig
          if(inc%neig==0) then
             inc%neig = 32 ! empirical
             write(*,*)"neig automatically set to ", inc%neig
          endif!inc%neig==0

          call IoInput('REIG      ',uio,1,7,ierr)
          read(unit=uio,fmt=*) inc%reig
          if(inc%reig==0.0d0) then
             inc%reig = 8.0/inc%natyp ! empirical
             write(*,*)"reig automatically set to ", inc%reig
          endif!inc%reig==0

          call IoInput('FEAST     ',uio,1,7,ierr)
          read(unit=uio,fmt=*) inc%feast

!          inc%ielast=1
!          write(*,*)"ielast automatically set to 1 to save memory (not suitable for Fermi velocity!)"
        end if!inc%memopt==1

        write(filename,fmt_fn_ext) filename_tbkkrrhod, ext_formatted
        inquire(file=trim(filename), exist=inc%lrhod)

        write(filename,fmt_fn_ext) filename_tbkkrtorq, ext_formatted
        inquire(file=trim(filename), exist=inc%ltorq)

        write(filename,fmt_fn_ext) filename_tbkkrspinflux, ext_formatted
        inquire(file=trim(filename), exist=inc%lspinflux)

        write(filename,fmt_fn_ext) filename_tbkkralpha, ext_formatted
        inquire(file=trim(filename), exist=inc%lalpha)
      end if!myrank==master

#ifdef CPP_MPI
      ! Brodcast input parameter
      call create_mpimask_inc(inc,mask_inc,ierr)
      if(ierr/=MPI_SUCCESS) stop 'error in constructing mask for inc'
      call MPI_Type_commit(mask_inc, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error commiting mask for inc'
      call MPI_Bcast(inc%N, 1, mask_inc, master, MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error brodcasting inc'
      call MPI_Type_free(mask_inc, ierr)
#endif


  end subroutine read_inc



  subroutine read_TBkkrdata(inc, lattice, cluster, tgmatrx)
  ! This subroutine reads the TBKKR-datafile,
  ! containing the main variables from the TB-KKR calculation.

    use type_inc,     only: inc_type
    use type_data,    only: lattice_type, lattice_init,&
                          & cluster_type, cluster_init,&
                          & tgmatrx_type, tgmatrx_init
    use mod_mympi,    only: myrank, nranks, master
    use mod_kkrmat,   only: InvertT, RearrangeTMAT
    use mod_ioformat, only: fmt_fn_ext, ext_formatted, filename_tbkkrcontainer,&
                          & filename_tbkkrrhod,filename_tbkkrtorq, filename_tbkkrspinflux,&
                          & filename_tbkkralpha
#ifdef CPP_MPI
    use type_data,    only: create_mpimask_lattice, create_mpimask_cluster, create_mpimask_tgmatrx
    use mpi
#endif
    implicit none

      type(inc_type),     intent(in)  :: inc
      type(lattice_type), intent(out) :: lattice
      type(cluster_type), intent(out) :: cluster
      type(tgmatrx_type), intent(out) :: tgmatrx

      integer            :: i1, i2, ie, lm1, lm2, isigma, ierr, idummy, iversion
      integer            :: mask_lattice, mask_cluster, mask_tgmatrx
      character(len=1)   :: strdum
      character(len=80)  :: filename, strline
      integer, parameter :: ifile= 12
      double complex     :: zdummy
!     double complex     :: ginp2(inc%naclsd*inc%lmmax,inc%lmmax,inc%nclsd,inc%ielast)
!     ginp2=(0d0,0d0)

      call lattice_init(inc,lattice)
      call cluster_init(inc,cluster)
      call tgmatrx_init(inc,tgmatrx)


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



        !====================================!
        !===          S T A R T           ===!
        !====================================!
        !=== energy dependent information ===!
        !====================================!

        !**********************
        !*** spin-channel 1 ***
        !**********************
        do ie=1,inc%ielast

          read(ifile,'(A)') strdum
          read(ifile,'(2ES25.16)') tgmatrx%energies(ie)

          ! Read in (energydependent) T-Matrix
          read(ifile,'(A)') strdum
          do i1=1,inc%naezd
            do lm2=1,inc%lmmax*inc%nspoh
              do lm1=1,inc%lmmax*inc%nspoh
                read(ifile,'(2ES25.16)') tgmatrx%tmatll(lm1,lm2,i1,ie)
              end do!lm1
            end do!lm2
          end do!i1

          ! Read in (energydependent) Greens function
          read(ifile,'(A)') strdum
          do i1=1,inc%nclsd
            do lm2=1,inc%lmmax
              do lm1=1,inc%lmmax*inc%naclsd
                read(ifile,'(2ES25.16)') tgmatrx%ginp(lm1,lm2,i1,ie)
              end do!lm1
            end do!lm2
          end do!i1

        end do!ie

        !**********************
        !*** spin-channel 2 ***
        !**********************
        if(inc%nspoh==1 .and. inc%nspd==2)then
          do ie=1,inc%ielast

!           write(*,*) 'spin channel 2: ie=', ie, ' of ', inc%ielast
            read(ifile,'(A)') strdum
            read(ifile,'(2ES25.16)') tgmatrx%energies(ie)

            ! Read in (energydependent) T-Matrix
            read(ifile,'(A)') strdum
            do i1=1,inc%naezd
              do lm2=1,inc%lmmax
                do lm1=1,inc%lmmax
                  read(ifile,'(2ES25.16)') tgmatrx%tmatll(inc%lmmax+lm1,inc%lmmax+lm2,i1,ie) !read second spin block
                  tgmatrx%tmatll(lm1,inc%lmmax+lm2,i1,ie) = (0d0,0d0)                        !set spin-off-diagonal parts to zero
                  tgmatrx%tmatll(inc%lmmax+lm1,lm2,i1,ie) = (0d0,0d0)                        !set spin-off-diagonal parts to zero
                end do!lm1
              end do!lm2
            end do!i1

            ! Read in (energydependent) Greens function
            read(ifile,'(A)') strdum
            do i1=1,inc%nclsd
              do lm2=1,inc%lmmax
                do lm1=1,inc%lmmax*inc%naclsd
                  read(ifile,'(2ES25.16)') zdummy
!                  read(ifile,'(2ES25.16)') ginp2(lm1,lm2,i1,ie)
                end do!lm1
              end do!lm2
            end do!i1

          end do!ie

        end if!inc%nspoh==1 .and. inc%nspd==2
        !====================================!
        !=== energy dependent information ===!
        !====================================!
        !===            E N D             ===!
        !====================================!

        close(ifile)

        !perform handy modifications on t-matrix-arrays for later use in the code
        do ie=1,inc%ielast
          call InvertT(inc, lattice%alat, tgmatrx%tmatll(:,:,:,ie), tgmatrx%tinvll(:,:,:,ie))
          call RearrangeTMAT(inc, lattice%alat, tgmatrx%tmatll(:,:,:,ie), tgmatrx%tmat(:,:,ie))
        end do!ie

        !read in rhod
        if(inc%lrhod)then
          write(filename,fmt_fn_ext) filename_tbkkrrhod, ext_formatted
          open(unit=ifile, file=trim(filename), form='formatted', action='read')
          do isigma=1,4
            do i1=1,inc%natypd
              do lm2=1,inc%lmmaxso
                do lm1=1,inc%lmmaxso
                  read(ifile,'(2ES25.16)') tgmatrx%rhod(lm1,lm2,i1,isigma)
                end do!lm1
              end do!lm2
            end do!i1
          end do!isigma
          close(ifile)
        else
          tgmatrx%rhod(:,:,:,1) = (1d0,0d0)
        end if!inc%lrhod

        !read in torq
        if(inc%ltorq)then
          write(filename,fmt_fn_ext) filename_tbkkrtorq, ext_formatted
          open(unit=ifile, file=trim(filename), form='formatted', action='read')
          do isigma=1,3
            do i1=1,inc%natypd
              do lm2=1,inc%lmmaxso
                do lm1=1,inc%lmmaxso
                  read(ifile,'(2ES25.16)') tgmatrx%torq(lm1,lm2,i1,isigma)
                end do!lm1
              end do!lm2
            end do!i1
          end do!isigma
          close(ifile)
        end if!inc%ltorq

        !read in spinflux
        if(inc%lspinflux)then
          write(filename,fmt_fn_ext) filename_tbkkrspinflux, ext_formatted
          open(unit=ifile, file=trim(filename), form='formatted', action='read')
          do isigma=1,3
            do i1=1,inc%natypd
              do lm2=1,inc%lmmaxso
                do lm1=1,inc%lmmaxso
                  read(ifile,'(2ES25.16)') tgmatrx%spinflux(lm1,lm2,i1,isigma)
                end do!lm1
              end do!lm2
            end do!i1
          end do!isigma
          close(ifile)
        end if!inc%lspinflux
        
        !read in alpha
        if(inc%lalpha)then
          write(filename,fmt_fn_ext) filename_tbkkralpha, ext_formatted
          open(unit=ifile, file=trim(filename), form='formatted', action='read')
          do isigma=1,3
            do i1=1,inc%natypd
              do lm2=1,inc%lmmaxso
                do lm1=1,inc%lmmaxso
                  read(ifile,'(2ES25.16)') tgmatrx%alpha(lm1,lm2,i1,isigma)
                end do!lm1
              end do!lm2
            end do!i1
          end do!isigma
          close(ifile)
        end if!inc%ltorq


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

      ! Brodcast tgmatrx information
      call create_mpimask_tgmatrx(tgmatrx,mask_tgmatrx)
      call MPI_Type_commit(mask_tgmatrx, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error commiting mask for tgmatrx'
      call MPI_Bcast(tgmatrx%N, 1, mask_tgmatrx, master, MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error brodcasting inc'
      call MPI_Type_free(mask_tgmatrx, ierr)
#endif
 

  end subroutine read_TBkkrdata


  subroutine read_kpointsfile_vis(nkpts, nkpts_irr, kpoints_irr, nsym, isym, kpt2irr, irr2kpt, vis2int, filenamein)

    use mod_ioformat, only: fmt_fn_sub_ext, ext_formatted, filemode_vis, filename_fsdata
    use mod_mympi,    only: myrank, nranks, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    integer,                       intent(out) :: nkpts, nkpts_irr,  nsym
    integer,          allocatable, intent(out) :: isym(:), kpt2irr(:), irr2kpt(:)
    double precision, allocatable, intent(out) :: kpoints_irr(:,:)
    integer,          allocatable, intent(out), optional :: vis2int(:)
    character(len=*),              intent(in),  optional :: filenamein

    integer :: ii, itmp(3), ierr
    character(len=256) :: filename, dummyline

    if(myrank==master)then
      if(present(filenamein))then
        filename = filenamein
      else
        write(filename,fmt_fn_sub_ext) filename_fsdata, filemode_vis, ext_formatted
      end if
      open(unit=326520, file=trim(filename), form='formatted', action='read')
      read(326520,'(3I8)') nkpts, nkpts_irr, nsym
    end if!myrank==master

#ifdef CPP_MPI
    if(myrank==master)then
      itmp(1) = nkpts
      itmp(2) = nkpts_irr
      itmp(3) = nsym
    end if!myrank==master

    call MPI_Bcast(itmp,3,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting nkpts'

    if(myrank/=master)then
      nkpts     = itmp(1)
      nkpts_irr = itmp(2)
      nsym      = itmp(3)
    end if!myrank==master
#endif

    allocate(isym(nsym), kpoints_irr(3,nkpts_irr), kpt2irr(nkpts), irr2kpt(nkpts_irr), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating isym, kpoints_irr etc.'

    if(present(vis2int)) then
      allocate( vis2int(nkpts), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating vis2int.'
    end if!present(vis2int)

    if(myrank==master)then
      read(326520,*) dummyline
      read(326520,'(12I8)') isym
      read(326520,*) dummyline
      read(326520,'(10I8)') irr2kpt
      read(326520,*) dummyline
      read(326520,'(10I8)') kpt2irr
      read(326520,*) dummyline
      read(326520,'(3ES25.16)') kpoints_irr

      if(present(vis2int)) then
        read(326520,*) dummyline
        read(326520,'(10I8)') vis2int
      end if!present(vis2int) 

      close(326520)
    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(isym,nsym,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting isym'

    call MPI_Bcast(irr2kpt,nkpts_irr,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting irr2kpt'

    call MPI_Bcast(kpt2irr,nkpts,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting kpt2irr'

    call MPI_Bcast(kpoints_irr,3*nkpts_irr,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting kpoints_irr'
#endif


  end subroutine read_kpointsfile_vis


  subroutine read_kpointsfile_int(nkpts, kpoints, areas, nsym, isym, filenamein)

    use mod_ioformat,   only: fmt_fn_sub_ext, ext_formatted, filemode_int, filename_fsdata
    use mod_mympi,    only: myrank, nranks, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    integer, intent(out) :: nkpts, nsym
    integer,          allocatable, intent(out) :: isym(:)
    double precision, allocatable, intent(out) :: kpoints(:,:), areas(:)
    character(len=*), intent(in), optional :: filenamein

    double precision, allocatable :: dtmpinp(:,:)
    integer :: ii, itmp(2), ierr
    character(len=256) :: filename

    if(myrank==master)then
      if(present(filenamein))then
        filename = filenamein
      else
        write(filename,fmt_fn_sub_ext) filename_fsdata, filemode_int, ext_formatted
      end if
      open(unit=326523, file=trim(filename), form='formatted', action='read')
      read(326523,'(2I8)') nkpts, nsym
    end if!myrank==master

#ifdef CPP_MPI
    if(myrank==master)then
      itmp(1) = nkpts
      itmp(2) = nsym
    end if!myrank==master

    call MPI_Bcast(itmp,2,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting nkpts'

    if(myrank/=master)then
      nkpts = itmp(1)
      nsym  = itmp(2)
    end if!myrank==master
#endif

    allocate(isym(nsym), kpoints(3,nkpts), areas(nkpts), dtmpinp(4,nkpts), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating isym, kpoints'

    if(myrank==master)then
      read(326523,'(12I8)') isym
      read(326523,'(4ES25.16)') dtmpinp
      close(326523)
    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(isym,nsym,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting isym'

    call MPI_Bcast(dtmpinp,4*nkpts,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting kpoints'
#endif

    kpoints = dtmpinp(1:3,:)
    areas   = dtmpinp(4,:)

    deallocate(dtmpinp)

  end subroutine read_kpointsfile_int


  

  subroutine read_weights(nkpts, weights, nsym, isym)

    use mod_ioformat, only: fmt_fn_ext, ext_formatted, filename_weights
    use mod_mympi,    only: myrank, nranks, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    integer, intent(out) :: nkpts, nsym
    integer,          allocatable, intent(out) :: isym(:)
    double precision, allocatable, intent(out) :: weights(:)

    integer :: ii, itmp(2), ierr
    character(len=256) :: filename

    if(myrank==master)then
      write(filename,fmt_fn_ext) filename_weights, ext_formatted
      open(unit=326522, file=trim(filename), form='formatted', action='read')
      read(326522,'(2I8)') nkpts, nsym
    end if!myrank==master

#ifdef CPP_MPI
    if(myrank==master)then
      itmp(1) = nkpts
      itmp(2) = nsym
    end if!myrank==master

    call MPI_Bcast(itmp,2,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting nkpts'

    if(myrank/=master)then
      nkpts = itmp(1)
      nsym  = itmp(2)
    end if!myrank==master
#endif

    allocate(isym(nsym), weights(nkpts), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating isym, kpoints'

    if(myrank==master)then
      read(326522,'(12I8)') isym
      read(326522,'(10ES25.16)') weights
      close(326522)
    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(isym,nsym,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting isym'

    call MPI_Bcast(weights,nkpts,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting kpoints'
#endif

  end subroutine read_weights




  subroutine read_fermivelocity(filemode, nkpts, fermivel, nsym, isym)

    use mod_ioformat, only: fmt_fn_sub_ext, ext_formatted, filename_fvel
    use mod_mympi,    only: myrank, nranks, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    character(len=*), intent(in)  :: filemode
    integer, intent(out) :: nkpts, nsym
    integer,          allocatable, intent(out) :: isym(:)
    double precision, allocatable, intent(out) :: fermivel(:,:)

    integer :: ii, itmp(2), ierr
    character(len=256) :: filename

    if(myrank==master)then
      write(filename,fmt_fn_sub_ext) filename_fvel, filemode, ext_formatted
      open(unit=326524, file=trim(filename), form='formatted', action='read')
      read(326524,'(2I8)') nkpts, nsym
    end if!myrank==master

#ifdef CPP_MPI
    if(myrank==master)then
      itmp(1) = nkpts
      itmp(2) = nsym
    end if!myrank==master

    call MPI_Bcast(itmp,2,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting itmp_2'

    if(myrank/=master)then
      nkpts = itmp(1)
      nsym  = itmp(2)
    end if!myrank==master
#endif

    allocate(isym(nsym), fermivel(3,nkpts), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating isym, fermivel'

    if(myrank==master)then
      read(326524,'(12I8)') isym
      read(326524,'(3ES25.16)') fermivel
      close(326524)
    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(isym,nsym,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting isym'

    call MPI_Bcast(fermivel,3*nkpts,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting fermivel'
#endif

  end subroutine read_fermivelocity





  subroutine read_spinvalue(filemode, nkpts, nsqa, ndegen, ispincomb, nvect, spinval, nsym, isym)

    use mod_ioformat, only: fmt_fn_sub_ext, ext_formatted, filename_spin
    use mod_mympi,    only: myrank, nranks, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    character(len=*), intent(in)  :: filemode
    integer, intent(out) :: nkpts, nsym, nsqa, ndegen
    integer,          allocatable, intent(out) :: isym(:), ispincomb(:)
    double precision, allocatable, intent(out) :: nvect(:,:), spinval(:,:,:)

    integer :: ii, itmp(4), ierr
    character(len=256) :: filename

    if(myrank==master)then
      write(filename,fmt_fn_sub_ext) filename_spin, filemode, ext_formatted
      open(unit=326526, file=trim(filename), form='formatted', action='read')
      read(326526,'(2I8)') nkpts, nsym, nsqa, ndegen
    end if!myrank==master

#ifdef CPP_MPI
    if(myrank==master)then
      itmp(1) = nkpts
      itmp(2) = nsym
      itmp(3) = nsqa
      itmp(4) = ndegen
    end if!myrank==master

    call MPI_Bcast(itmp,4,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting itmp_4'

    if(myrank/=master)then
      nkpts = itmp(1)
      nsym  = itmp(2)
      nsqa  = itmp(3)
      ndegen= itmp(4)
    end if!myrank==master
#endif

    allocate( isym(nsym), ispincomb(nsqa), nvect(3,nsqa), &
            & spinval(ndegen,nsqa,nkpts), STAT=ierr       )
    if(ierr/=0) stop 'Problem allocating isym, spinval'

    if(myrank==master)then
      read(326526,'(12I8)') isym
      read(326526,'(10I8)') ispincomb
      read(326526,'(3ES25.16)') nvect
      read(326526,'(ES25.16)') spinval
      close(326526)
    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(isym,nsym,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting isym'

    call MPI_Bcast(ispincomb,nsqa,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting ispincomb'

    call MPI_Bcast(nvect,3*nsqa,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting nvect'

    call MPI_Bcast(spinval,ndegen*nsqa*nkpts,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting spinval'
#endif

  end subroutine read_spinvalue


  subroutine read_torqvalue(filemode, nkpts, ndegen, torqval, nsym, isym)

    use mod_ioformat, only: fmt_fn_sub_ext, ext_formatted, filename_torq
    use mod_mympi,    only: myrank, nranks, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    character(len=*), intent(in)  :: filemode
    integer, intent(out) :: nkpts, nsym, ndegen
    integer,          allocatable, intent(out) :: isym(:)
    double precision, allocatable, intent(out) :: torqval(:,:,:)

    integer :: ii, itmp(3), ierr
    character(len=256) :: filename

    if(myrank==master)then
      write(filename,fmt_fn_sub_ext) filename_torq, filemode, ext_formatted
      open(unit=326526, file=trim(filename), form='formatted', action='read')
      read(326526,'(3I8)') nkpts, nsym, ndegen
    end if!myrank==master

#ifdef CPP_MPI
    if(myrank==master)then
      itmp(1) = nkpts
      itmp(2) = nsym
      itmp(3) = ndegen
    end if!myrank==master

    call MPI_Bcast(itmp,3,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting itmp_2'

    if(myrank/=master)then
      nkpts = itmp(1)
      nsym  = itmp(2)
      ndegen= itmp(3)
    end if!myrank==master
#endif

    allocate(isym(nsym), torqval(3,ndegen,nkpts), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating isym, torqval'

    if(myrank==master)then
      read(326526,'(12I8)') isym
      read(326526,'(3ES25.16)') torqval
      close(326526)
    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(isym,nsym,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting isym'

    call MPI_Bcast(torqval,3*ndegen*nkpts,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting torqval'
#endif

  end subroutine read_torqvalue

  subroutine read_torqvalue_atom(filemode, natyp, nkpts, ndegen, torqval_atom, nsym, isym)

    use mod_ioformat, only: fmt_fn_atom_sub_ext, ext_formatted, filename_torq
    use mod_mympi,    only: myrank, nranks, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    character(len=*), intent(in)  :: filemode
    integer, intent(in)  :: natyp
    integer, intent(out) :: nkpts, nsym, ndegen
    integer,          allocatable, intent(out) :: isym(:)
    double precision, allocatable, intent(out) :: torqval_atom(:,:,:,:) !torqval_atom(3,natyp,ndegen,nkpts)

    integer :: ii, itmp(3), ierr
    character(len=256) :: filename

    if(myrank==master)then
        write(filename,fmt_fn_atom_sub_ext) filename_torq, 1, trim(filemode), ext_formatted
        open(unit=326526, file=trim(filename), form='formatted', action='read')
        read(326526,'(3I8)') nkpts, nsym, ndegen
        close(326526)
    end if!myrank==master

#ifdef CPP_MPI
    if(myrank==master)then
      itmp(1) = nkpts
      itmp(2) = nsym
      itmp(3) = ndegen
    end if!myrank==master

    call MPI_Bcast(itmp,3,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting itmp_2'

    if(myrank/=master)then
      nkpts = itmp(1)
      nsym  = itmp(2)
      ndegen= itmp(3)
    end if!myrank==master
#endif

    allocate(isym(nsym),torqval_atom(3,natyp,ndegen,nkpts), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating isym, torqval_atom'

    if(myrank==master)then
      do ii=1,natyp
        write(filename,fmt_fn_atom_sub_ext) filename_torq, ii, trim(filemode), ext_formatted
        open(unit=326526, file=trim(filename), form='formatted', action='read')
        read(326526,'(3I8)') nkpts, nsym, ndegen ! not needed, but used to skip lines
        read(326526,'(12I8)') isym
        read(326526,'(3ES25.16)') torqval_atom(:,ii,:,:)
        close(326526)
      end do!ii
    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(isym,nsym,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting isym'

    call MPI_Bcast(torqval_atom,3*natyp*ndegen*nkpts,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting torqval_atom'
#endif

  end subroutine read_torqvalue_atom

  subroutine read_spinvec_atom(filemode, natyp, nkpts, ndegen, spinvec_atom, nsym, isym)

    use mod_ioformat, only: fmt_fn_atom_sub_ext, ext_formatted, filename_spinvec
    use mod_mympi,    only: myrank, nranks, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    character(len=*), intent(in)  :: filemode
    integer, intent(in)  :: natyp
    integer, intent(out) :: nkpts, nsym, ndegen
    integer,          allocatable, intent(out) :: isym(:)
    double precision, allocatable, intent(out) :: spinvec_atom(:,:,:,:) !spinvec_atom(3,natyp,ndegen,nkpts)

    integer :: ii, itmp(3), ierr
    character(len=256) :: filename

    if(myrank==master)then
        write(filename,fmt_fn_atom_sub_ext) filename_spinvec, 1, trim(filemode), ext_formatted
        open(unit=326526, file=trim(filename), form='formatted', action='read')
        read(326526,'(3I8)') nkpts, nsym, ndegen
        close(326526)
    end if!myrank==master

#ifdef CPP_MPI
    if(myrank==master)then
      itmp(1) = nkpts
      itmp(2) = nsym
      itmp(3) = ndegen
    end if!myrank==master

    call MPI_Bcast(itmp,3,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting itmp_2'

    if(myrank/=master)then
      nkpts = itmp(1)
      nsym  = itmp(2)
      ndegen= itmp(3)
    end if!myrank==master
#endif

    allocate(isym(nsym),spinvec_atom(3,natyp,ndegen,nkpts), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating isym, spinvec_atom'

    if(myrank==master)then
      do ii=1,natyp
        write(filename,fmt_fn_atom_sub_ext) filename_spinvec, ii, trim(filemode), ext_formatted
        open(unit=326526, file=trim(filename), form='formatted', action='read')
        read(326526,'(3I8)') nkpts, nsym, ndegen ! not needed, but used to skip lines
        read(326526,'(12I8)') isym
        read(326526,'(3ES25.16)') spinvec_atom(:,ii,:,:)
        close(326526)
      end do!ii
    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(isym,nsym,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting isym'

    call MPI_Bcast(spinvec_atom,3*natyp*ndegen*nkpts,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting spinvec_atom'
#endif

  end subroutine read_spinvec_atom

  subroutine read_spinflux_atom(filemode, natyp, nkpts, ndegen, spinflux, nsym, isym)

    use mod_ioformat, only: fmt_fn_atom_sub_ext, ext_formatted, filename_spinflux
    use mod_mympi,    only: myrank, nranks, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    character(len=*), intent(in)  :: filemode
    integer, intent(in)  :: natyp
    integer, intent(out) :: nkpts, nsym, ndegen
    integer,          allocatable, intent(out) :: isym(:)
    double precision, allocatable, intent(out) :: spinflux(:,:,:,:) !spinflux(3,natyp,ndegen,nkpts)

    integer :: ii, itmp(3), ierr
    character(len=256) :: filename

    if(myrank==master)then
        write(filename,fmt_fn_atom_sub_ext) filename_spinflux, 1, trim(filemode), ext_formatted
        open(unit=326526, file=trim(filename), form='formatted', action='read')
        read(326526,'(3I8)') nkpts, nsym, ndegen
        close(326526)
    end if!myrank==master

#ifdef CPP_MPI
    if(myrank==master)then
      itmp(1) = nkpts
      itmp(2) = nsym
      itmp(3) = ndegen
    end if!myrank==master

    call MPI_Bcast(itmp,3,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting itmp_2'

    if(myrank/=master)then
      nkpts = itmp(1)
      nsym  = itmp(2)
      ndegen= itmp(3)
    end if!myrank==master
#endif

    allocate(isym(nsym), spinflux(3,natyp,ndegen,nkpts), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating isym, spinflux'

    if(myrank==master)then
      do ii=1,natyp
        write(filename,fmt_fn_atom_sub_ext) filename_spinflux, ii, trim(filemode), ext_formatted
        open(unit=326526, file=trim(filename), form='formatted', action='read')
        read(326526,'(3I8)') nkpts, nsym, ndegen ! not needed, but used to skip lines
        read(326526,'(12I8)') isym
        read(326526,'(3ES25.16)') spinflux(:,ii,:,:)
        close(326526)
      end do!ii
    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(isym,nsym,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting isym'

    call MPI_Bcast(spinflux,3*natyp*ndegen*nkpts,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting spinflux'
#endif

  end subroutine read_spinflux_atom

  subroutine read_alphavalue(filemode, nkpts, ndegen, alphaval, nsym, isym)

    use mod_ioformat, only: fmt_fn_sub_ext, ext_formatted, filename_alpha
    use mod_mympi,    only: myrank, nranks, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    character(len=*), intent(in)  :: filemode
    integer, intent(out) :: nkpts, nsym, ndegen
    integer,          allocatable, intent(out) :: isym(:)
    double precision, allocatable, intent(out) :: alphaval(:,:,:)

    integer :: ii, itmp(3), ierr
    character(len=256) :: filename

    if(myrank==master)then
      write(filename,fmt_fn_sub_ext) filename_alpha, filemode, ext_formatted
      open(unit=326526, file=trim(filename), form='formatted', action='read')
      read(326526,'(3I8)') nkpts, nsym, ndegen
    end if!myrank==master

#ifdef CPP_MPI
    if(myrank==master)then
      itmp(1) = nkpts
      itmp(2) = nsym
      itmp(3) = ndegen
    end if!myrank==master

    call MPI_Bcast(itmp,3,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting itmp_2'

    if(myrank/=master)then
      nkpts = itmp(1)
      nsym  = itmp(2)
      ndegen= itmp(3)
    end if!myrank==master
#endif

    allocate(isym(nsym), alphaval(3,ndegen,nkpts), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating isym, alphaval'

    if(myrank==master)then
      read(326526,'(12I8)') isym
      read(326526,'(3ES25.16)') alphaval
      close(326526)
    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(isym,nsym,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting isym'

    call MPI_Bcast(alphaval,3*ndegen*nkpts,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting alphaval'
#endif

  end subroutine read_alphavalue


end module mod_read
