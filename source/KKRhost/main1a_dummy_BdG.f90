!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------
!> Summary: Dummy version of `main1a` which is used to test the `tmat_newsolver`
!> Author: Philipp Ruessmann
!> Category: unit-test, KKRhost 
!> Deprecated: False 
!> Dummy version of `main1a` which is used to test the `tmat_newsolver`
!-------------------------------------------------------------------------------
subroutine main1a_dummy
  use :: mod_profiling, only: memocc
  use :: mod_constants, only: 
  use :: global_variables, only: nrmaxd, lmpotd, nspotd, iemxd, ncleb, ntotd, mmaxd, nspind, natypd
  use :: mod_datatypes, only: dp

  use :: mod_types, only: t_tgmat, t_inc, t_lloyd, t_dtmatjij, init_t_dtmatjij, init_t_dtmatjij_at
  use :: mod_mympi, only: nranks, master, myrank
  use :: mod_timing, only: 
  use :: mod_wunfiles, only: 
  use :: mod_jijhelp, only: set_jijcalc_flags

  use :: mod_tmatnewsolver, only: tmat_newsolver
  use :: mod_main0, only: ielast, natyp, nspin, lmax, nsra, iend, lly, deltae, idoldau, ncheb, cleb, icleb, ez, npan_tot, &
    lopt, ipan_intervall, zat, socscale, rnew, rpan_intervall, wldau
  use :: mod_ioinput, only: ioinput
  use :: mod_runoptions, only:write_BdG_tests


  implicit none

  ! .. Local variables
  integer :: i1
  integer :: ipot
  integer :: iltmp
  integer :: ispin
  integer :: itmpdir
  character (len=80) :: tmpdir
  logical :: lrefsys
  integer :: lrectmt
  integer :: lrectra
  ! .. Local arrays
  real (kind=dp), dimension (natypd) :: phi
  real (kind=dp), dimension (natypd) :: theta
  real (kind=dp), dimension (:, :, :), allocatable :: vinsnew

  integer :: i1_run, i1_start, i1_end, ierr, i_stat, i_all

  ! for data import:
  character (len=25) :: dummy
  integer :: ier
  character (len=:), allocatable :: uio           ! NCOLIO=256


  ! BdG specific:
  logical :: use_bdg
  complex (kind=dp) :: delta_bdg


  ! ================
  ! start read-in

  write (*, *) 'start reading BdG-inputs from inputcard ...'

  call ioinput('use_BdG         ', uio, 1, 7, ier)
  if (ier==0) then
    read (unit=uio, fmt=*) use_bdg
    write (*, *) 'use_BdG= ', use_bdg
  else
    stop '[main1a_dummy] error "use_BdG" not found in inputcard'
  end if

  call ioinput('delta_BdG       ', uio, 1, 7, ier)
  if (ier==0) then
    read (unit=uio, fmt=*) delta_bdg
    write (*, *) 'delta_BdG= ', delta_bdg
  else
    stop '[main1a_dummy] error "delta_BdG" not found in inputcard'
  end if

  ! end read-in
  ! ================


  write (*, *) 'now read atom-specific input for tmat_newsolver'

  allocate (vinsnew(nrmaxd,lmpotd,nspotd), stat=i_stat)
  call memocc(i_stat, product(shape(vinsnew))*kind(vinsnew), 'VINSNEW', 'main1a')
  vinsnew = 0.0d0

  i1_start = 1
  i1_end = natyp

  do i1_run = i1_start, i1_end

    if (write_BdG_tests) then
      ! read out inputs for tmat_newsolver to extract first BdG
      if (nranks>1) stop 'test option BdG_dev can only be used in serial!'
      if (i1_run==1) open (887766, file='BdG_tmat_inputs.txt', form='formatted')
      if (i1_run==1) then
        read (887766, *) dummy
        read (887766, *)
        read (887766, *) dummy, ielast
        read (887766, *) dummy, nspin
        read (887766, *) dummy, lmax
        read (887766, *) dummy, nsra
        read (887766, *) dummy, iend
        read (887766, *) dummy, lmpotd
        read (887766, *) dummy, lly
        read (887766, '(A,2ES21.9)') dummy, deltae
        read (887766, *) dummy, idoldau
        read (887766, *) dummy, ncleb
        read (887766, *) dummy, ncheb
        read (887766, *) dummy, ntotd
        read (887766, *) dummy, mmaxd
        read (887766, *) dummy, nspind
        read (887766, *) dummy, iemxd
        read (887766, *) dummy, nrmaxd
        read (887766, *) dummy, nspotd
        read (887766, *) dummy
        read (887766, *) cleb(:, 1)
        read (887766, *) dummy
        read (887766, *) icleb(:, :)
        read (887766, *) dummy
        read (887766, '(2ES21.9)') ez
        read (887766, *)
        read (887766, *) dummy
        read (887766, *)
      end if
      read (887766, *) dummy, i1
      read (887766, *) dummy, ipot
      read (887766, *) dummy, npan_tot(i1)
      read (887766, *) dummy, lopt(i1)
      read (887766, *) dummy, ipan_intervall(:, i1)
      read (887766, *) dummy, zat(i1)
      read (887766, *) dummy, phi(i1)
      read (887766, *) dummy, theta(i1)
      read (887766, *) dummy, socscale(i1)
      read (887766, *) dummy, rnew(:, i1)
      read (887766, *) dummy, rpan_intervall(:, i1)
      read (887766, *) dummy, wldau(:, :, :, i1)
      read (887766, *) dummy, vinsnew
      read (887766, *)
      if (i1==i1_end) then
        close (887766)
        write (*, *) 'done reading tmat_newsolver input of test option BdG_dev'
      end if
    end if

    write (*, *) 'start tmat_newsolver ...'

    call init_t_dtmatjij(t_inc, t_dtmatjij)

    call tmat_newsolver(ielast, nspin, lmax, zat(i1), socscale(i1), ez, nsra, cleb(:,1), icleb, iend, ncheb, npan_tot(i1), rpan_intervall(:,i1), ipan_intervall(:,i1), rnew(:,i1), &
      vinsnew, theta(i1), phi(i1), i1, ipot, lmpotd, lly, deltae, idoldau, lopt(i1), wldau(:,:,:,i1), t_dtmatjij(i1), 1)

  end do                           ! I1, atom loop

end subroutine main1a_dummy
