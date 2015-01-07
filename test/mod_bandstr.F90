  subroutine read_TBkkr_data(inc, lattice, cluster)
  ! This subroutine reads the TBKKR-datafile,
  ! containing the main variables from the TB-KKR calculation.

    use type_inc,     only: inc_type
    use type_data,    only: lattice_type, lattice_init,&
                          & cluster_type, cluster_init
    use mod_mympi,    only: myrank, nranks, master
    use mod_kkrmat,   only: InvertT, CalcTmat
    use mod_ioformat, only: fmt_fn_ext, ext_formatted, filename_tbkkrcontainer, filename_tbkkrrhod
#ifdef CPP_MPI
    use type_data,    only: create_mpimask_lattice, create_mpimask_cluster
    use mpi
#endif
    implicit none

      type(inc_type),     intent(in)  :: inc
      type(lattice_type), intent(out) :: lattice
      type(cluster_type), intent(out) :: cluster

      integer            :: i1, i2, ie, lm1, lm2, isigma, ierr, idum1, idum2
      integer            :: mask_lattice, mask_cluster
      character(len=1)   :: strdum
      character(len=80)  :: filename
      integer, parameter :: ifile= 12, ifilegreen=62


      call lattice_init(inc,lattice)
      call cluster_init(inc,cluster)


      if(myrank==master) then

        write(filename,fmt_fn_ext) filename_tbkkrcontainer, ext_formatted
        open(unit=ifile, file=trim(filename), form='formatted', action='read')

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

        read(ifile,'(A)') strdum
        read(ifile,'(1I8)') (cluster%eqinv(i1),i1=1,inc%naezd)

        read(ifile,'(A)') strdum
        read(ifile,'(1I8)') ((cluster%ezoa(i1,i2),i1=1,inc%naclsd),i2=1,inc%naezd)

        read(ifile,'(A)') strdum
        read(ifile,'(1I8)') ((cluster%atom(i1,i2),i1=1,inc%naclsd),i2=1,inc%naezd)

        read(ifile,'(A)') strdum
        read(ifile,'(1I8)') (cluster%kaoez(i1),   i1=1,inc%naezd+inc%nembd)

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
