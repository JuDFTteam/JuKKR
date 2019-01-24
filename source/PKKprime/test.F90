!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


program TEST

  use type_inc,       only: inc_type
  use type_data,      only: lattice_type, cluster_type, tgmatrx_type
  use mod_read,       only: read_inc, read_TBkkrdata
  use mod_pointgrp,   only: pointgrp
  use mod_symmetries, only: rotate_kpoints
  use mod_read,       only: read_kpointsfile_int
  use mod_symmetries, only: symmetries_type, set_symmetries, rotate_kpoints
  use mod_mympi,      only: mympi_init, myrank, nranks, master
#ifdef CPP_MPI
  use mpi
#endif

  implicit none

    type(inc_type)     :: inc
    type(lattice_type) :: lattice
    type(cluster_type) :: cluster
    type(tgmatrx_type) :: tgmatrx

    !symmetry arrays
    integer :: nsym
    integer, allocatable :: isym(:)
    type(symmetries_type) :: symmetries

    !local k-point arrays
    integer :: nkpts
    double precision, allocatable :: kpoints(:,:), areas(:), weights(:), fermivel(:,:)

    !temp k-point arrays
    integer :: nkpts1, nkpts2
    double precision, allocatable :: areas1(:), weights1(:), kpoints1(:,:), fermivel1(:,:), temparr(:,:)

    integer :: ierr, ii, itest
    integer, parameter :: nFSiter=3

    !initialize MPI
#ifdef CPP_MPI
    call MPI_Init ( ierr )
#endif
    call mympi_init()

    call read_inc(inc)
    call read_TBkkrdata(inc, lattice, cluster, tgmatrx)

    call set_symmetries(inc, lattice, symmetries)

      do itest=nFSiter,1,-1
        write(*,*) 'itest=', itest
      end do


    !=============================!
    != Read in the k-point files =!
    call read_kpointsfile_int(nkpts1, kpoints1, areas1, nsym, isym)
    call rotate_kpoints(symmetries%rotmat, nkpts1, kpoints1, nsym, isym, nkpts, kpoints)
    deallocate(isym, areas1, kpoints1)



    if(myrank==master) write(*,'(A,3ES25.16)') 'kpoints-sum in old cart. system:  ', sum(kpoints,  dim=2)
    call project_fermivel_newaxis(nkpts,kpoints)

    if(myrank==master) write(*,'(A,3ES25.16)') 'kpoints-sum in new cart. system:  ', sum(kpoints,  dim=2)

#ifdef CPP_MPI
    call MPI_Finalize( ierr )
#endif

contains

  subroutine project_fermivel_newaxis(nkpts,fermivel)
    implicit none
    integer,          intent(in)    :: nkpts
    double precision, intent(inout) :: fermivel(3,nkpts)

    double precision :: fermivel_tmp(3,nkpts), n1(3), n2(3), n3(3)
    integer :: ikp


    n1=(/ 1d0, -1d0,  0d0 /)/sqrt(2d0)
    n2=(/ 1d0,  1d0, -2d0 /)/sqrt(6d0)
    n3=(/ 1d0,  1d0,  1d0 /)/sqrt(3d0)

    fermivel_tmp = fermivel

    do ikp=1,nkpts
      fermivel(1,ikp) = sum(fermivel_tmp(:,ikp)*n1)
      fermivel(2,ikp) = sum(fermivel_tmp(:,ikp)*n2)
      fermivel(3,ikp) = sum(fermivel_tmp(:,ikp)*n3)
    end do!ikp

  end subroutine project_fermivel_newaxis

end program TEST
