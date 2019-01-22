!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Wrapper for the calculation of the irregular solutions
!> Author: 
!> Wrapper for the calculation of the irregular solutions
!------------------------------------------------------------------------------------
module mod_sll_global_solutions

contains

  !-------------------------------------------------------------------------------
  !> Summary: Wrapper for the calculation of the irregular solutions
  !> Author: 
  !> Category: single-site, KKRhost
  !> Deprecated: False 
  !> Wrapper for the calculation of the irregular solutions for the impurity code `KKRimp`
  !-------------------------------------------------------------------------------
  subroutine sll_global_solutions(rpanbound,rmesh,vll,sll,ncheb,npan,lmsize,lmsize2,&
    lbessel,nrmax, nvec,jlk_index,hlk,jlk,hlk2,jlk2,gmatprefactor,cmodesll,         &
    use_sratrick1)               ! LLY
    ! ************************************************************************
    ! for description see rllsll routine
    ! ************************************************************************
    use :: mod_timing, only: timing_start, timing_stop  ! timing routine
#ifdef CPP_HYBRID
    use :: omp_lib               ! omp functions
#endif

    use :: mod_constants, only: czero, cone
    use :: mod_datatypes, only: dp
    use :: mod_chebint, only: chebint
    use :: mod_sll_local_solutions, only: sll_local_solutions

    implicit none
    integer :: ncheb             ! number of chebyshev nodes
    integer :: npan              ! number of panels
    integer :: lmsize            ! lm-components * nspin
    integer :: lmsize2           ! lmsize * nvec
    integer :: nvec              ! spinor integer
    ! nvec=1 non-rel, nvec=2 for sra and dirac
    integer :: nrmax             ! total number of rad. mesh points
    integer :: lbessel, use_sratrick1 ! dimensions etc., needed only for host code interface

    ! running indices
    integer :: lm1, lm2
    integer :: icheb, ipan, mn

    ! source terms
    complex (kind=dp) :: gmatprefactor ! prefactor of green function
    ! non-rel: = kappa = sqrt e
    complex (kind=dp) :: hlk(lbessel, nrmax), jlk(lbessel, nrmax), hlk2(lbessel, nrmax), jlk2(lbessel, nrmax)

    integer :: jlk_index(2*lmsize)

    character (len=1) :: cmodesll                           ! These define the op(V(r))

    ! cmodesll ="1" : op( )=identity
    ! cmodesll ="T" : op( )=transpose in L

    complex (kind=dp) :: sll(lmsize2, lmsize, nrmax), & ! irr. volterra sol.
      vll(lmsize*nvec, lmsize*nvec, nrmax) ! potential term in 5.7
    ! bauer, phd

    complex (kind=dp), allocatable :: work(:, :), cllp(:, :, :), dllp(:, :, :), mihvy(:, :, :), mihvz(:, :, :), &
      mijvy(:, :, :), mijvz(:, :, :) ! ***************
    complex (kind=dp), allocatable :: yif(:, :, :, :), & ! source terms (different array
      zif(:, :, :, :)
    ! chebyshev arrays
    real (kind=dp) :: c1(0:ncheb, 0:ncheb), rpanbound(0:npan), drpan2
    real (kind=dp) :: cslc1(0:ncheb, 0:ncheb), & ! Integration matrix from left ( C*S_L*C^-1 in eq. 5.53)
      csrc1(0:ncheb, 0:ncheb), & ! Same from right ( C*S_R*C^-1 in eq. 5.54)
      tau(0:ncheb, 0:npan), &    ! Radial mesh point
      slc1sum(0:ncheb), rmesh(nrmax)

    integer :: ierror, use_sratrick
    integer :: idotime
    integer, parameter :: directsolv = 1

#ifdef CPP_HYBRID
    ! openMP variable --sacin 23/04/2015
    integer :: thread_id
#endif

    ! ***********************************************************************
    ! SRA trick
    ! ***********************************************************************
    ! on page 68 of Bauer, PhD, a method is described how to speed up the
    ! calculations in case of the SRA. A similar approach is implemented
    ! here by using Eq. 4.132 and substituting DV from 4.133, and discretising
    ! the radial mesh of the Lippmann-Schwinger eq. according to 5.68.
    ! The Lippmann-Schwinger equation leads to a matrix inversion
    ! problem. The matrix M which needs to be inverted has a special form
    ! if the SRA approximation is used:

    ! matrix A ( C 1)     (same as in eq. 5.68)
    ! ( B 0)
    ! (C, B are matricies here)

    ! inverse of A is   (C^-1    0 )
    ! (-B C^-1 1 )
    ! Thus, it is sufficient to only inverse the matrix C which saves computational
    ! time. This is refered to as the SRA trick.
    ! ***********************************************************************
    ! in future implementation equation 4.134 is supposed to be
    ! implemented which should lead to an additional speed-up.
    ! ***********************************************************************

    if (lmsize==1) then
      use_sratrick = 0
    else
      use_sratrick = use_sratrick1
    end if

    ! turn timing output off if in the host code
    idotime = 0
#ifdef test_run
    idotime = 1
#endif
    if (idotime==1) call timing_start('sll-glob')


    do ipan = 1, npan
      do icheb = 0, ncheb
        mn = ipan*ncheb + ipan - icheb
        tau(icheb, ipan) = rmesh(mn)
      end do
    end do

    call chebint(cslc1, csrc1, slc1sum, c1, ncheb)

    allocate (mihvy(lmsize,lmsize,npan), mihvz(lmsize,lmsize,npan))
    allocate (mijvy(lmsize,lmsize,npan), mijvz(lmsize,lmsize,npan))
    allocate (yif(lmsize2,lmsize,0:ncheb,npan))
    allocate (zif(lmsize2,lmsize,0:ncheb,npan))

    allocate (work(lmsize,lmsize))
    allocate (cllp(lmsize,lmsize,0:npan), dllp(lmsize,lmsize,0:npan))

    if (idotime==1) call timing_start('local')

#ifdef CPP_HYBRID
    ! $OMP PARALLEL DEFAULT (PRIVATE) &
    ! $OMP& SHARED(tau,npan,drpan2,rpanbound,mihvy,mihvz,mijvy,mijvz,yif) &
    ! $OMP& SHARED(zif,nvec,lmsize,lmsize2,ncheb,jlk,jlk2,jlk_index) &
    ! $OMP& SHARED(vll,gmatprefactor,hlk,hlk2,cslc1,csrc1,slc1sum) &
    ! $OMP& SHARED(cmodesll,use_sratrick, rmesh)

    thread_id = omp_get_thread_num()
#endif

    ! loop over subintervals
#ifdef CPP_HYBRID
    ! openMP pragmas added sachin, parallel region starts earlier to get allocations of arrays right

    ! $OMP DO
#endif
    do ipan = 1, npan

      drpan2 = (rpanbound(ipan)-rpanbound(ipan-1))/2.d0 ! *(b-a)/2 in eq. 5.53, 5.54
      call sll_local_solutions(vll,tau(0,ipan),drpan2,csrc1, slc1sum,               &
        mihvy(1,1,ipan),mihvz(1,1,ipan),mijvy(1,1,ipan),mijvz(1,1,ipan),            &
        yif(1,1,0,ipan),zif(1,1,0,ipan),ncheb,ipan,lmsize,lmsize2,nrmax,nvec,       &
        jlk_index,hlk,jlk,hlk2,jlk2,gmatprefactor,cmodesll,lbessel,use_sratrick1)

    end do                       ! ipan
#ifdef CPP_HYBRID
    ! $OMP END DO
    ! $OMP END PARALLEL
#endif
    ! end the big loop over the subintervals

    if (idotime==1) call timing_stop('local')
    if (idotime==1) call timing_start('afterlocal')

    ! ***********************************************************************
    ! calculate A(M), B(M), C(M), D(M)
    ! according to 5.17-5.18 (regular solution) of Bauer PhD
    ! C,D are calculated accordingly for the irregular solution
    ! (starting condition: A(0) = 1, B(0) = 0, C(MMAX) = 0 and D(MMAX) = 1)
    ! ***********************************************************************

    ! irregular
    do lm2 = 1, lmsize
      do lm1 = 1, lmsize
        dllp(lm1, lm2, npan) = czero
        cllp(lm1, lm2, npan) = czero
      end do
    end do

    do lm1 = 1, lmsize
      dllp(lm1, lm1, npan) = cone
    end do

    do ipan = npan, 1, -1
      call zcopy(lmsize*lmsize, cllp(1,1,ipan), 1, cllp(1,1,ipan-1), 1)
      call zcopy(lmsize*lmsize, dllp(1,1,ipan), 1, dllp(1,1,ipan-1), 1)
      call zgemm('n', 'n', lmsize, lmsize, lmsize, cone, mihvz(1,1,ipan), lmsize, cllp(1,1,ipan), lmsize, cone, cllp(1,1,ipan-1), lmsize)
      call zgemm('n', 'n', lmsize, lmsize, lmsize, cone, mihvy(1,1,ipan), lmsize, dllp(1,1,ipan), lmsize, cone, cllp(1,1,ipan-1), lmsize)
      call zgemm('n', 'n', lmsize, lmsize, lmsize, -cone, mijvz(1,1,ipan), lmsize, cllp(1,1,ipan), lmsize, cone, dllp(1,1,ipan-1), lmsize)
      call zgemm('n', 'n', lmsize, lmsize, lmsize, -cone, mijvy(1,1,ipan), lmsize, dllp(1,1,ipan), lmsize, cone, dllp(1,1,ipan-1), lmsize)
    end do

    ! ***********************************************************************
    ! determine the irregular solution sll by using 5.14
    ! ***********************************************************************
    do ipan = 1, npan
      do icheb = 0, ncheb
        mn = ipan*ncheb + ipan - icheb
        call zgemm('n', 'n', lmsize2, lmsize, lmsize, cone, zif(1,1,icheb,ipan), lmsize2, cllp(1,1,ipan), lmsize, czero, sll(1,1,mn), lmsize2)
        call zgemm('n', 'n', lmsize2, lmsize, lmsize, cone, yif(1,1,icheb,ipan), lmsize2, dllp(1,1,ipan), lmsize, cone, sll(1,1,mn), lmsize2)
      end do
    end do

    if (idotime==1) call timing_stop('afterlocal')
    if (idotime==1) call timing_start('endstuff')

    if (idotime==1) call timing_stop('endstuff')
    if (idotime==1) call timing_start('checknan')
    if (idotime==1) call timing_stop('checknan')
    if (idotime==1) call timing_stop('sll-glob')

    deallocate (work, cllp, dllp, mihvy, mihvz, mijvy, mijvz, yif, zif, stat=ierror)
    if (ierror/=0) stop '[sll-glob] ERROR in deallocating arrays'
  end subroutine sll_global_solutions

end module mod_sll_global_solutions
