!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_gll13

  private
  public :: gll13

contains

  !-------------------------------------------------------------------------------
  !> Summary: Solution of the DYSON equation for a cluster
  !> Author: P. Zahn, Ph. Mavropoulos
  !> Date: Oct. 2013
  !> Category: KKRhost, structural-greensfunction, reference-system
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> solution of the DYSON equation for a cluster of potentials
  !> (TMATLL) centered at positions RATOM in free space,
  !>
  !> (modified version of GLL91 by P. Zahn, Sept. 95)
  !> (modified by Phivos Mavropoulos to apply Lloyds formula
  !> ported from KKRnano, Oct. 2013)
  !-------------------------------------------------------------------------------
  subroutine gll13(ez, cleb, icleb, loflm, iend, tmatll, dtmatll, atom, refpot, ratom, natom, tolrdif, alat, out_wr, gref0, dgdeout, naclsmax, lly_g0tr, lly)


#ifdef CPP_HYBRID
    use :: omp_lib
#endif
    use :: mod_types, only: t_inc
    use :: mod_runoptions, only: print_program_flow, print_refpot
    use :: global_variables, only: lmgf0d, nrefd, ncleb, naclsd, naezd, nembd, lm2d
    use :: mod_datatypes, only: dp
    use :: mod_gfree13, only: gfree13
    use :: mod_grefsy13, only: grefsy13
    use :: mod_constants, only: cone, czero
    implicit none
    ! ..
    ! .. Scalar Arguments ..
    complex (kind=dp) :: ez
    real (kind=dp) :: alat, tolrdif ! Set free GF to zero if R<TOLRDIF in case of virtual atoms
    integer :: iend, natom, out_wr, naclsmax
    ! ..
    ! .. Array Arguments ..
    complex (kind=dp) :: gref0(naclsmax*lmgf0d, lmgf0d), tmatll(lmgf0d, lmgf0d, nrefd)
    real (kind=dp) :: cleb(ncleb), ratom(3, naclsd)
    integer :: atom(naclsd), icleb(ncleb, 4), loflm(lm2d), refpot(naezd+nembd)
    ! ..
    ! .. Local Scalars ..
    integer :: i, lm1, lm2, m, n, n1, n2, ndim, nlm1, nlm2, info, ngd1
    ! ..
    ! .. Local Arrays ..
    integer, allocatable :: ipvt(:)
    real (kind=dp) :: rdiff(3), absrdiff
    complex (kind=dp) :: dtmatll(lmgf0d, lmgf0d, nrefd) ! Derivative of ref.-sys t-matrix
    complex (kind=dp) :: dgdeout(naclsmax*lmgf0d, lmgf0d)
    complex (kind=dp), allocatable :: gref(:, :), gll(:, :), gtref(:, :)
    complex (kind=dp), allocatable :: dgtde(:, :), dgtde0(:, :) ! LLY (1-gt)^-1 * d(1-gt)/dE (after grefsy13)
    complex (kind=dp), allocatable :: dgllde(:, :), dgde(:, :)
    complex (kind=dp) :: lly_g0tr  ! LLY Trace of  DTGLL for Lloyds formula
    integer :: lly                 ! LLY =0 : no Lloyd's formula; <>0: use Lloyd's formula
    ! ..
    ! .. External Functions ..
    logical, external :: test
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic :: abs
    ! ..
    ! .. OpenMP stuff ..
#ifdef CPP_HYBRID
    integer :: thread_id
#endif


    ngd1 = naclsmax*lmgf0d
    ndim = lmgf0d*natom

    allocate (gtref(ngd1,lmgf0d), dgtde(ngd1,lmgf0d), stat=lm1)
    if (lm1/=0) stop 'Error allocating gtref etc. <GLL95>'
    gtref(:, :) = czero
    dgtde(:, :) = czero

    allocate (gref(ngd1,ngd1), ipvt(ngd1), stat=lm1)
    if (lm1/=0) stop 'Error allocating gref etc. <GLL95>'
    gref(:, :) = czero
    ipvt(:) = 0
    if (lly/=0) then
      allocate (dgtde0(ngd1,ngd1), dgde(ngd1,ngd1), stat=lm1)
      if (lm1/=0) stop 'Error allocating dgtde0 etc. <GLL95>'
      dgtde0(:, :) = czero
      dgtde(:, :) = czero
    end if
100 format (6x, 'ERROR: failed to allocate array(s) :', a, /)
    if (print_program_flow .and. (t_inc%i_write>0)) write (1337, fmt=*) '>>> GLL95'


#ifdef CPP_HYBRID
    ! $omp parallel default(shared) &
    ! $omp& private(n1, n2, rdiff, absrdiff, lm2, lm1, nlm2, nlm1, GLL) &
    ! $omp& private(thread_id, gtref, dgtde, DGLLDE)
    thread_id = omp_get_thread_num()
#endif
    ! allocate here, inside omp parallel region
    allocate (gll(lmgf0d,lmgf0d), dgllde(lmgf0d,lmgf0d), stat=lm1)
    if (lm1/=0) stop 'Error allocating gll etc. <GLL95>'
    gll(:, :) = czero
    dgllde(:, :) = czero


    ! ---> construct free Green's function

#ifdef CPP_HYBRID
    ! $omp do
#endif
    do n1 = 1, natom
      do n2 = 1, natom
        rdiff(1:3) = -(ratom(1:3,n1)-ratom(1:3,n2))*alat
        absrdiff = sqrt(rdiff(1)**2+rdiff(2)**2+rdiff(3)**2)

        if (n1/=n2 .and. (absrdiff>tolrdif)) then
          call gfree13(rdiff, ez, gll, dgllde, cleb, icleb, loflm, iend)
          do lm2 = 1, lmgf0d
            nlm2 = (n2-1)*lmgf0d + lm2
            do lm1 = 1, lmgf0d
              nlm1 = (n1-1)*lmgf0d + lm1
              gref(nlm1, nlm2) = gll(lm1, lm2)
              if (lly/=0) dgde(nlm1, nlm2) = dgllde(lm1, lm2)
            end do
          end do
        else
          do lm2 = 1, lmgf0d
            nlm2 = (n2-1)*lmgf0d + lm2
            do lm1 = 1, lmgf0d
              nlm1 = (n1-1)*lmgf0d + lm1
              gref(nlm1, nlm2) = czero
              if (lly/=0) dgde(nlm1, nlm2) = czero
            end do
          end do
        end if

      end do
    end do
#ifdef CPP_HYBRID
    ! $omp end do
    ! deallocate in omp parallel region
#endif
    deallocate (gll, dgllde, stat=lm1)
    if (lm1/=0) stop ' [gll13] dealloc'
#ifdef CPP_HYBRID
    ! $omp end parallel
#endif
    if (print_program_flow .and. (t_inc%i_write>0)) write (1337, fmt=*) 'GFREE o.k.'
    ! ----------------------------------------------------------------------

    ! GREF0 = g:= gfree
    call zcopy(ngd1*lmgf0d, gref, 1, gref0, 1)

    ! ----------------------------------------------------------------------
    ! LLY Lloyd
    ! Prepare source term -dg/dE * t - g * dt/dE
    if (lly/=0) then

      do n2 = 1, natom
        nlm2 = (n2-1)*lmgf0d + 1
        ! GTREF = -DGDE*t = -dg/dE * t
        call zgemm('N', 'N', ndim, lmgf0d, lmgf0d, -cone, dgde(1,nlm2), ngd1, tmatll(1,1,refpot(abs(atom(n2)))), lmgf0d, czero, gtref, ngd1)
        ! GTREF = GTREF - GREF*DTMATLL = -dg/dE * t - g * dt/dE  (here, GREF=g:=gfree)
        call zgemm('N', 'N', ndim, lmgf0d, lmgf0d, -cone, gref(1,nlm2), ngd1, dtmatll(1,1,refpot(abs(atom(n2)))), lmgf0d, cone, gtref, ngd1)
        call zcopy(ngd1*lmgf0d, gtref, 1, dgtde0(1,nlm2), 1)
      end do
      do n2 = 1, lmgf0d
        do n1 = 1, ngd1
          dgtde(n1, n2) = dgtde0(n1, n2)
        end do
      end do
      ! Now DGTDE = DGTDE0 = -dg/dE * t - g * dt/dE
      ! (DGTDE is reduced matrix; DGTDE0 is full matrix)

    end if                         ! (LLY.NE.0)
    ! LLY Lloyd
    ! ----------------------------------------------------------------------
    do n2 = 1, natom
      nlm2 = (n2-1)*lmgf0d + 1
      ! GTREF =  -g*t
      call zgemm('N', 'N', ndim, lmgf0d, lmgf0d, -cone, gref(1,nlm2), ngd1, tmatll(1,1,refpot(abs(atom(n2)))), lmgf0d, czero, gtref, ngd1)
      call zcopy(ngd1*lmgf0d, gtref, 1, gref(1,nlm2), 1)
      ! Now GREF =  -g*t
      if (print_refpot .and. (t_inc%i_write>0)) write (1337, fmt=*) n2, refpot(abs(atom(n2))), atom(n2)
    end do

    call grefsy13(gref, gref0, dgtde, lly_g0tr, ipvt, ndim, lly, lmgf0d, ngd1)
    ! Now GREF contains LU(1-gt) (full matrix NGD1xNGD1)
    ! DGTDE contains (1-gt)^-1 * d(1-gt)/dE (Thiess PhD Eq.5.28)
    ! between atoms 0 and n (n running through the cluster)
    ! and LLY_G0TR contains -Trace[ (1-gt)^-1 * d(1-gt)/dE ] (Thiess PhD Eq.5.38)


    ! ----------------------------------------------------------------------
    ! LLY Lloyd
    dgdeout(:, :) = czero
    if (lly/=0) then

      ! Prepare dg/de + (dg/dE * t + g * dt/dE)*Gref (Thiess PhD Eq.5.42)

      ! DGDE = DGDE - DGTDE0*GREF0 = dg/de + (dg/dE * t + g * dt/dE)*Gref
      call zgemm('N', 'N', ndim, lmgf0d, ndim, -cone, dgtde0, ngd1, gref0, ngd1, cone, dgde, ngd1)

      ! Solve linear system: (remember GREF contains LU(1-gt))
      ! (1-gt)*DGDE = dg/de + (dg/dE * t + g * dt/dE)*Gref
      call zgetrs('N', ndim, lmgf0d, gref, ngd1, ipvt, dgde, ngd1, info)
      ! Result is DGDE = dGref/dE

      do n2 = 1, lmgf0d
        do n1 = 1, ngd1
          dgdeout(n1, n2) = dgde(n1, n2)
        end do
      end do

    end if
    ! LLY Lloyd
    ! ----------------------------------------------------------------------


    if (print_program_flow .and. (t_inc%i_write>0)) write (1337, fmt=*) 'GREFSY o.k.'

    if (out_wr>0) write (out_wr)((gref0(n,m),m=1,lmgf0d), n=1, ngd1)

    ! deallocate arrays
    deallocate (gtref, dgtde, gref, ipvt, stat=lm1)
    if (lm1/=0) stop ' [gll13] dealloc'
    if (lly/=0) then
      deallocate (dgtde0, dgde, stat=lm1)
      if (lm1/=0) stop ' [gll13] dealloc'
    end if
  end subroutine gll13

end module mod_gll13
