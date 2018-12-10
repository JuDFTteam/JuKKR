!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Wrapper for the calculation of the regular solutions
!> Author: 
!> Wrapper for the calculation of the regular solutions
!------------------------------------------------------------------------------------
!> @note preprocessor options: change this definition if used in host/impurity code
!> this is commented out, since then the logical hostcode is not defined
!> and thus `#indef hostcode` returns true and `#ifdef hostcode` false 
!>
!> * **Host code**: leave the line `#define hostcode` **uncommented**.
!> * **Impurity code**: **comment** the line`#define hostcode`.
!>
!> This allows one to choose between interface for impurity and host code (different calling lists)
!> @endnote
!------------------------------------------------------------------------------------
module mod_rll_global_solutions

contains

#define hostcode 

#ifndef hostcode
  !-------------------------------------------------------------------------------
  !> Summary: Wrapper for the calculation of the regular solutions
  !> Author: 
  !> Category: single-site, KKRhost
  !> Deprecated: False 
  !> Wrapper for the calculation of the regular solutions for the host code `KKRhost`
  !-------------------------------------------------------------------------------
  !> @note One can comment out the line `#define hostcode` to instead choose the
  !> calculation of the regular solutions for the impuirty code.
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine rll_global_solutions(rpanbound,rmesh,vll,rll,tllp,ncheb,npan,lmsize,   &
    lmsize2,nrmax,nvec,jlk_index,hlk,jlk,hlk2,jlk2,gmatprefactor,cmoderll,idotime)
#else
  !-------------------------------------------------------------------------------
  !> Summary: Wrapper for the calculation of the regular solutions
  !> Author: 
  !> Category: single-site, KKRhost
  !> Deprecated: False 
  !> Wrapper for the calculation of the regular solutions for the impurity code `KKRimp`
  !-------------------------------------------------------------------------------
  !> @note One can uncomment out the line `#define hostcode` to instead choose the
  !> calculation of the regular solutions for the host code.
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine rll_global_solutions(rpanbound,rmesh,vll,rll,tllp,ncheb,npan,lmsize,   &
    lmsize2,lbessel,nrmax,nvec,jlk_index,hlk,jlk,hlk2,jlk2,gmatprefactor,cmoderll,  &
    use_sratrick1,alphaget)     ! LLY
#endif
    ! ************************************************************************
    ! for description see rllsll routine
    ! ************************************************************************
#ifndef hostcode
    use :: mod_beshank           ! calculates bessel and hankel func.
    use :: mod_chebint           ! chebyshev integration routines
    use :: mod_config, only: config_testflag ! reads if testflags are present
    use :: mod_rllslltools, only: inverse ! inverse matrix routine
    use :: mod_physic_params, only: cvlight ! speed of light
    use :: sourceterms
    use :: mod_chebyshev
#endif
    use :: mod_timing            ! timing routine
#ifdef CPP_HYBRID
    use :: omp_lib               ! omp functions
#endif

    use :: mod_constants
    use :: mod_datatypes, only: dp
    use :: mod_chebint
    use :: mod_rll_local_solutions

    implicit none
    integer :: ncheb             ! number of chebyshev nodes
    integer :: npan              ! number of panels
    integer :: lmsize            ! lm-components * nspin
    integer :: lmsize2           ! lmsize * nvec
    integer :: nvec              ! spinor integer
    ! nvec=1 non-rel, nvec=2 for sra and dirac
    integer :: nrmax             ! total number of rad. mesh points
#ifdef hostcode
    integer :: lbessel, use_sratrick1 ! dimensions etc., needed only for host code interface
#endif

    ! running indices
    integer :: lm1, lm2
    integer :: info, icheb, ipan, mn, nm

    ! source terms
    complex (kind=dp) :: gmatprefactor ! prefactor of green function
    ! non-rel: = kappa = sqrt e
#ifndef hostcode
    complex (kind=dp) :: hlk(:, :), jlk(:, :), & ! right sol. source terms
      hlk2(:, :), jlk2(:, :)     ! left sol. source terms
    ! (typically bessel and hankel fn)
#else
    complex (kind=dp) :: hlk(lbessel, nrmax), jlk(lbessel, nrmax), hlk2(lbessel, nrmax), jlk2(lbessel, nrmax)
#endif

#ifndef hostcode
    integer :: jlk_index(:)      ! mapping array l = jlk_index(lm)
    ! in: lm-index
    ! corresponding l-index used hlk,..
    ! hlk(l) = jlk_index(lm)
#else
    integer :: jlk_index(2*lmsize)
#endif

    character (len=1) :: cmoderll                           ! These define the op(V(r))

    ! cmoderll ="1" : op( )=identity
    ! cmoderll ="T" : op( )=transpose in L

    complex (kind=dp) :: rll(lmsize2, lmsize, nrmax), & ! reg. fredholm sol.
      tllp(lmsize, lmsize), &    ! t-matrix
      vll(lmsize*nvec, lmsize*nvec, nrmax) ! potential term in 5.7
    ! bauer, phd
    complex (kind=dp), allocatable :: ull(:, :, :) ! reg. volterra sol.

    complex (kind=dp), allocatable :: work(:, :), allp(:, :, :), bllp(:, :, :), & ! eq. 5.9, 5.10 for reg. sol
      mrnvy(:, :, :), mrnvz(:, :, :), & ! ***************
      mrjvy(:, :, :), mrjvz(:, :, :) ! eq. 5.19-5.22
    complex (kind=dp), allocatable :: yrf(:, :, :, :), & ! source terms (different array
      zrf(:, :, :, :)            ! ordering)
    ! chebyshev arrays
    real (kind=dp) :: c1(0:ncheb, 0:ncheb), rpanbound(0:npan), drpan2
    real (kind=dp) :: cslc1(0:ncheb, 0:ncheb), & ! Integration matrix from left ( C*S_L*C^-1 in eq. 5.53)
      csrc1(0:ncheb, 0:ncheb), & ! Same from right ( C*S_R*C^-1 in eq. 5.54)
      tau(0:ncheb, 0:npan), &    ! Radial mesh point
      slc1sum(0:ncheb), rmesh(nrmax)

    integer :: ipiv(0:ncheb, lmsize2)
    integer :: ierror, use_sratrick
    integer :: idotime
    integer, parameter :: directsolv = 1
#ifdef hostcode
    complex (kind=dp) :: alphaget(lmsize, lmsize) ! LLY
#endif

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

    ! matrix A ( C 0)     (same as in eq. 5.68)
    ! ( B 1)
    ! (C, B are matricies here)

    ! inverse of A is   (C^-1    0 )
    ! (-B C^-1 1 )
    ! Thus, it is sufficient to only inverse the matrix C which saves computational
    ! time. This is refered to as the SRA trick.
    ! ***********************************************************************
    ! in future implementation equation 4.134 is supposed to be
    ! implemented which should lead to an additional speed-up.
    ! ***********************************************************************

#ifndef hostcode
    if (config_testflag('nosph') .or. lmsize==1) then
      use_sratrick = 0
    else if (.not. config_testflag('nosph')) then
      use_sratrick = 1
    else
      stop '[rll-glob] use_sratrick error'
    end if
#else
    if (lmsize==1) then
      use_sratrick = 0
    else
      use_sratrick = use_sratrick1
    end if
#endif

#ifdef hostcode
    ! turn timing output off if in the host code
    idotime = 0
#endif
#ifdef test_run
    idotime = 1
#endif
    if (idotime==1) call timing_start('rll-glob')


    do ipan = 1, npan
      do icheb = 0, ncheb
        mn = ipan*ncheb + ipan - icheb
        tau(icheb, ipan) = rmesh(mn)
      end do
    end do

    call chebint(cslc1, csrc1, slc1sum, c1, ncheb)

    allocate (mrnvy(lmsize,lmsize,npan), mrnvz(lmsize,lmsize,npan))
    allocate (mrjvy(lmsize,lmsize,npan), mrjvz(lmsize,lmsize,npan))
    allocate (yrf(lmsize2,lmsize,0:ncheb,npan))
    allocate (zrf(lmsize2,lmsize,0:ncheb,npan))

    allocate (work(lmsize,lmsize))
    allocate (allp(lmsize,lmsize,0:npan), bllp(lmsize,lmsize,0:npan))

    allocate (ull(lmsize2,lmsize,nrmax))

    if (idotime==1) call timing_start('local')

#ifdef CPP_HYBRID
    ! $OMP PARALLEL DEFAULT (PRIVATE) &
    ! $OMP& SHARED(tau,npan,drpan2,rpanbound,mrnvy,mrnvz,mrjvy,mrjvz,yrf) &
    ! $OMP& SHARED(zrf,nvec,lmsize,lmsize2,ncheb,jlk,jlk2,jlk_index,vll) &
    ! $OMP& SHARED(gmatprefactor,hlk,hlk2,cslc1,csrc1,slc1sum) &
    ! $OMP& SHARED(cmoderll,use_sratrick, rmesh)

    thread_id = omp_get_thread_num()
#endif

    ! loop over subintervals
#ifdef CPP_HYBRID
    ! openMP pragmas added sachin, parallel region starts earlier to get allocations of arrays right

    ! $OMP DO
#endif
    do ipan = 1, npan

      drpan2 = (rpanbound(ipan)-rpanbound(ipan-1))/2.d0 ! *(b-a)/2 in eq. 5.53, 5.54
      call rll_local_solutions(vll,tau(0,ipan),drpan2,cslc1,slc1sum,mrnvy(1,1,ipan),& 
        mrnvz(1,1,ipan),mrjvy(1,1,ipan),mrjvz(1,1,ipan),yrf(1,1,0,ipan),            &
        zrf(1,1,0,ipan),ncheb,ipan,lmsize,lmsize2,nrmax,nvec,jlk_index,hlk,jlk,hlk2,& 
        jlk2,gmatprefactor,cmoderll,lbessel,use_sratrick1)

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

    ! regular
    do lm2 = 1, lmsize
      do lm1 = 1, lmsize
        bllp(lm1, lm2, 0) = czero
        allp(lm1, lm2, 0) = czero
      end do
    end do

    do lm1 = 1, lmsize
      allp(lm1, lm1, 0) = cone
    end do

    do ipan = 1, npan
      call zcopy(lmsize*lmsize, allp(1,1,ipan-1), 1, allp(1,1,ipan), 1)
      call zcopy(lmsize*lmsize, bllp(1,1,ipan-1), 1, bllp(1,1,ipan), 1)
      call zgemm('n','n',lmsize,lmsize,lmsize,-cone,mrnvy(1,1,ipan),lmsize,         &
        allp(1,1,ipan-1),lmsize,cone,allp(1,1,ipan),lmsize)
      call zgemm('n','n',lmsize,lmsize,lmsize,-cone,mrnvz(1,1,ipan),lmsize,         &
        bllp(1,1,ipan-1),lmsize,cone,allp(1,1,ipan),lmsize)
      call zgemm('n','n',lmsize,lmsize,lmsize,cone,mrjvy(1,1,ipan),lmsize,          &
        allp(1,1,ipan-1),lmsize,cone,bllp(1,1,ipan),lmsize)
      call zgemm('n','n',lmsize,lmsize,lmsize,cone,mrjvz(1,1,ipan),lmsize,          &
        bllp(1,1,ipan-1),lmsize,cone,bllp(1,1,ipan),lmsize)
    end do

    ! ***********************************************************************
    ! determine the regular solution ull by using 5.14
    ! ***********************************************************************
    do ipan = 1, npan
      do icheb = 0, ncheb
        mn = ipan*ncheb + ipan - icheb
        call zgemm('n','n',lmsize2,lmsize,lmsize,cone,yrf(1,1,icheb,ipan),lmsize2,  &
          allp(1,1,ipan-1),lmsize,czero,ull(1,1,mn),lmsize2)
        call zgemm('n','n',lmsize2,lmsize,lmsize,cone,zrf(1,1,icheb,ipan),lmsize2,  &
          bllp(1,1,ipan-1),lmsize,cone,ull(1,1,mn),lmsize2)
      end do
    end do

    if (idotime==1) call timing_stop('afterlocal')
    if (idotime==1) call timing_start('endstuff')

    ! ***********************************************************************
    ! next part converts the volterra solution u of equation (5.7) to
    ! the fredholm solution r by employing eq. 4.122 and 4.120 of bauer, phd
    ! and the t-matrix is calculated
    ! ***********************************************************************

    call zgetrf(lmsize, lmsize, allp(1,1,npan), lmsize, ipiv, info) ! invert alpha
    call zgetri(lmsize, allp(1,1,npan), lmsize, ipiv, work, lmsize*lmsize, info) ! invert alpha -> transformation matrix rll=alpha^-1*rll
#ifdef hostcode
    ! get alpha matrix
    do lm1 = 1, lmsize           ! LLY
      do lm2 = 1, lmsize         ! LLY
        alphaget(lm1, lm2) = allp(lm1, lm2, npan) ! LLY
      end do                     ! LLY
    end do                       ! LLY
#endif
    ! calculation of the t-matrix ! calc t-matrix tll = bll*alpha^-1
    call zgemm('n','n',lmsize,lmsize,lmsize,cone/gmatprefactor,bllp(1,1,npan),      & 
      lmsize,allp(1,1,npan),lmsize,czero,tllp,lmsize)

    do nm = 1, nrmax
      call zgemm('n','n',lmsize2,lmsize,lmsize,cone,ull(1,1,nm),lmsize2,            &
        allp(1,1,npan),lmsize,czero,rll(1,1,nm),lmsize2)
    end do

    if (idotime==1) call timing_stop('endstuff')
    if (idotime==1) call timing_start('checknan')
    if (idotime==1) call timing_stop('checknan')
    if (idotime==1) call timing_stop('rll-glob')

    deallocate (work, allp, bllp, mrnvy, mrnvz, mrjvy, mrjvz, yrf, zrf, stat=ierror)
    if (ierror/=0) stop '[rll-glob] ERROR in deallocating arrays'

  end subroutine rll_global_solutions

end module mod_rll_global_solutions
