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
    use :: mod_runoptions, only: use_cheby_quadprec
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
    integer :: icheb, ipan, mn, nm
    integer :: info, ipiv(lmsize)
    integer :: iter_beta, niter_beta

    ! source terms
    complex (kind=dp) :: gmatprefactor ! prefactor of green function
    ! non-rel: = kappa = sqrt e
    complex (kind=dp) :: hlk(lbessel, nrmax), jlk(lbessel, nrmax), hlk2(lbessel, nrmax), jlk2(lbessel, nrmax)

    integer :: jlk_index(2*lmsize)

    character (len=1) :: cmodesll                           ! These define the op(V(r))

    ! cmodesll ="1" : op( )=identity
    ! cmodesll ="T" : op( )=transpose in L
 
    complex*32, allocatable :: qcllp(:, :, :), qdllp(:, :, :)
    complex*32, allocatable :: qmihvy(:, :), qmihvz(:, :), qmijvy(:, :), qmijvz(:, :)
    complex*32, allocatable :: qyif(:, :, :)
    complex*32, allocatable :: qcllptemp(:, :), qdllptemp(:, :)
    complex*32, allocatable :: qsll(:, :)
    complex*32, allocatable :: qcone, qczero
    complex*32, allocatable :: qbetainv(:, :), qbetainv_save(:, :)

    complex (kind=dp) :: sll(lmsize2, lmsize, nrmax), & ! irr. volterra sol.
      vll(lmsize*nvec, lmsize*nvec, nrmax) ! potential term in 5.7
    ! bauer, phd

    real (kind=dp) :: dllpmax,dllpval
    complex (kind=dp), allocatable :: cllptemp(:, :), dllptemp(:, :)
    complex (kind=dp), allocatable :: work(:, :), cllp(:, :, :), dllp(:, :, :), mihvy(:, :, :), mihvz(:, :, :), &
      mijvy(:, :, :), mijvz(:, :, :) ! ***************
    complex (kind=dp), allocatable :: yif(:, :, :, :), zif(:, :, :, :)
    complex (kind=dp), allocatable :: betainv(:, :), betainv_save(:, :) 

    ! chebyshev arrays
    real (kind=dp) :: c1(0:ncheb, 0:ncheb), rpanbound(0:npan), drpan2
    real (kind=dp) :: cslc1(0:ncheb, 0:ncheb), & ! Integration matrix from left ( C*S_L*C^-1 in eq. 5.53)
      csrc1(0:ncheb, 0:ncheb), & ! Same from right ( C*S_R*C^-1 in eq. 5.54)
      tau(0:ncheb, 0:npan), &    ! Radial mesh point
      slc1sum(0:ncheb), rmesh(nrmax)

    integer :: ierror
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

    ! turn timing output off if in the host code
    idotime = 0
#ifdef test_run
    idotime = 1
#endif
    if (idotime==1) call timing_start('sll-glob')

    allocate (mihvy(lmsize,lmsize,npan), mihvz(lmsize,lmsize,npan))
    allocate (mijvy(lmsize,lmsize,npan), mijvz(lmsize,lmsize,npan))
    allocate (yif(lmsize2,lmsize,0:ncheb,npan))
    allocate (zif(lmsize2,lmsize,0:ncheb,npan))
    allocate (betainv(lmsize,lmsize), betainv_save(lmsize,lmsize))
    allocate (cllp(lmsize,lmsize,0:npan), dllp(lmsize,lmsize,0:npan))
    allocate (cllptemp(lmsize,lmsize), dllptemp(lmsize,lmsize))
    allocate (work(lmsize,lmsize))

    do ipan = 1, npan
      do icheb = 0, ncheb
        mn = ipan*ncheb + ipan - icheb
        tau(icheb, ipan) = rmesh(mn)
      end do
    end do

    call chebint(cslc1, csrc1, slc1sum, c1, ncheb)

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
    ! calculate C(M), D(M)
    ! (starting condition: C(MMAX) = 0 and D(MMAX) = 1)
    ! ***********************************************************************

    niter_beta = 3

    cllp(:, :, npan) = czero
    dllp(:, :, npan) = czero
    do lm1 = 1, lmsize
      dllp(lm1, lm1, npan) = cone
    end do

    do ipan = npan, 1, -1

      cllp(:, :, ipan-1) = cllp(:, :, ipan)
      dllp(:, :, ipan-1) = dllp(:, :, ipan)

      call zgemm('n', 'n', lmsize, lmsize, lmsize, cone, mihvz(1,1,ipan), lmsize, cllp(1,1,ipan), lmsize, cone, cllp(1,1,ipan-1), lmsize)
      call zgemm('n', 'n', lmsize, lmsize, lmsize, cone, mihvy(1,1,ipan), lmsize, dllp(1,1,ipan), lmsize, cone, cllp(1,1,ipan-1), lmsize)
      call zgemm('n', 'n', lmsize, lmsize, lmsize, -cone, mijvz(1,1,ipan), lmsize, cllp(1,1,ipan), lmsize, cone, dllp(1,1,ipan-1), lmsize)
      call zgemm('n', 'n', lmsize, lmsize, lmsize, -cone, mijvy(1,1,ipan), lmsize, dllp(1,1,ipan), lmsize, cone, dllp(1,1,ipan-1), lmsize)

    end do

    ! ***********************************************************************
    ! invert beta = dllp(:, :, 0)
    ! ***********************************************************************

    betainv = dllp(:, :, 0)
    call zgetrf(lmsize, lmsize, betainv, lmsize, ipiv, info)
    call zgetri(lmsize, betainv, lmsize, ipiv, work, lmsize*lmsize, info)
 
    if(use_cheby_quadprec) then
      allocate (qcone, qczero)
      qcone = (1.q0,0.0q0)
      qczero = (0.q0,0.0q0)
      allocate (qmihvy(lmsize,lmsize), qmihvz(lmsize,lmsize))
      allocate (qmijvy(lmsize,lmsize), qmijvz(lmsize,lmsize))
      allocate (qyif(lmsize2,lmsize,0:ncheb))
      allocate (qbetainv(lmsize,lmsize), qbetainv_save(lmsize,lmsize))
      allocate (qsll(lmsize2,lmsize))
      allocate (qcllp(lmsize,lmsize,0:npan), qdllp(lmsize,lmsize,0:npan))
      allocate (qcllptemp(lmsize,lmsize), qdllptemp(lmsize,lmsize))
      qbetainv = betainv
    end if

    do iter_beta = 1, niter_beta

    if(.not.use_cheby_quadprec) then

    dllp(:, :, npan) = betainv
    cllp(:, :, npan) = czero

    do lm2 = 1, lmsize
        dllp(lm2, lm2, npan) = betainv(lm2,lm2) - cone
    end do

    do ipan = npan, 1, -1

      cllp(:, :, ipan-1) = cllp(:, :, ipan) + mihvy(:, :, ipan)
      dllp(:, :, ipan-1) = dllp(:, :, ipan) - mijvy(:, :, ipan)

      call zgemm('n', 'n', lmsize, lmsize, lmsize, cone, mihvz(1,1,ipan), lmsize, cllp(1,1,ipan), lmsize, cone, cllp(1,1,ipan-1), lmsize)
      call zgemm('n', 'n', lmsize, lmsize, lmsize, cone, mihvy(1,1,ipan), lmsize, dllp(1,1,ipan), lmsize, cone, cllp(1,1,ipan-1), lmsize)
      call zgemm('n', 'n', lmsize, lmsize, lmsize, -cone, mijvz(1,1,ipan), lmsize, cllp(1,1,ipan), lmsize, cone, dllp(1,1,ipan-1), lmsize)
      call zgemm('n', 'n', lmsize, lmsize, lmsize, -cone, mijvy(1,1,ipan), lmsize, dllp(1,1,ipan), lmsize, cone, dllp(1,1,ipan-1), lmsize)

    end do

    betainv_save = betainv

    call zgemm('n', 'n', lmsize, lmsize, lmsize, -cone, betainv_save, lmsize, dllp(1,1,0), lmsize, cone, betainv, lmsize)

    dllpmax = 0.d0
    do lm1 = 1,lmsize
      do lm2 = 1,lmsize
      dllpval = dllp(lm1,lm2,0)
      if(lm1.ne.lm2.and.abs(dllpval).gt.dllpmax) dllpmax = abs(dllpval)
      if(lm1.eq.lm2.and.abs(dllpval-cone).gt.dllpmax) dllpmax = abs(dllpval-cone)
      end do
    end do

    else

    qdllp(:, :, npan) = qbetainv
    qcllp(:, :, npan) = qczero

    do lm2 = 1, lmsize
        qdllp(lm2, lm2, npan) = qbetainv(lm2,lm2) - qcone
    end do

    do ipan = npan, 1, -1

          qmihvz(:, :) = mihvz(:, :, ipan) 
          qmihvy(:, :) = mihvy(:, :, ipan) 
          qmijvz(:, :) = mijvz(:, :, ipan) 
          qmijvy(:, :) = mijvy(:, :, ipan) 

      qcllp(:, :, ipan-1) = qcllp(:, :, ipan) + qmihvy(:, :)
      qdllp(:, :, ipan-1) = qdllp(:, :, ipan) - qmijvy(:, :)

      call cqgemm(lmsize,lmsize,lmsize,qcone,qmihvz,lmsize,qcllp(1,1,ipan),lmsize,qcone,qcllp(1,1,ipan-1),lmsize)
      call cqgemm(lmsize,lmsize,lmsize,qcone,qmihvy,lmsize,qdllp(1,1,ipan),lmsize,qcone,qcllp(1,1,ipan-1),lmsize)
      call cqgemm(lmsize,lmsize,lmsize,-qcone,qmijvz,lmsize,qcllp(1,1,ipan),lmsize,qcone,qdllp(1,1,ipan-1),lmsize)
      call cqgemm(lmsize,lmsize,lmsize,-qcone,qmijvy,lmsize,qdllp(1,1,ipan),lmsize,qcone,qdllp(1,1,ipan-1),lmsize)

    end do
    
    qbetainv_save = qbetainv

    call cqgemm(lmsize, lmsize, lmsize, -qcone, qbetainv_save, lmsize, qdllp(1,1,0), lmsize, qcone, qbetainv, lmsize)
    dllpmax = 0.0d0
    do lm1 = 1,lmsize
      do lm2 = 1,lmsize
      dllpval = qdllp(lm1,lm2,0)
      if(abs(dllpval).gt.dllpmax) dllpmax = abs(dllpval)
      end do
    end do
    end if
 
    ! test writeout
    !write(6,*) 'dllpmax',dllpmax,'iter_beta',iter_beta

    end do ! niter_beta

    ! ***********************************************************************
    ! determine the irregular solution sll by using 5.14
    ! ***********************************************************************

    if(.not.use_cheby_quadprec) then

    do ipan = 0, npan
       do lm1 = 1, lmsize
        dllp(lm1,lm1,ipan) = dllp(lm1,lm1,ipan) + cone
       end do
    end do

      do ipan = 1, npan
      cllptemp(:, :) = cllp(:, :, ipan)
      dllptemp(:, :) = dllp(:, :, ipan)
      cllp(:, :,ipan) = cllptemp(:, :)*(cone+cone)
      dllp(:, :,ipan) = dllptemp(:, :)*(cone+cone)
      call zgemm('n', 'n', lmsize, lmsize, lmsize, -cone, cllptemp, lmsize, dllp(1,1,0), lmsize, cone, cllp(1,1,ipan), lmsize)
      call zgemm('n', 'n', lmsize, lmsize, lmsize, -cone, dllptemp, lmsize, dllp(1,1,0), lmsize, cone, dllp(1,1,ipan), lmsize)
    end do

    do ipan = 1, npan
      do icheb = 0, ncheb
        mn = ipan*ncheb + ipan - icheb
        call zgemm('n', 'n', lmsize2, lmsize, lmsize, cone, zif(1,1,icheb,ipan), lmsize2, cllp(1,1,ipan), lmsize, czero, sll(1,1,mn), lmsize2)
        call zgemm('n', 'n', lmsize2, lmsize, lmsize, cone, yif(1,1,icheb,ipan), lmsize2, dllp(1,1,ipan), lmsize, cone, sll(1,1,mn), lmsize2)
      end do
    end do
   
    else

    do ipan = 0, npan
       do lm1 = 1, lmsize
        qdllp(lm1,lm1,ipan) = qdllp(lm1,lm1,ipan) + qcone
       end do
    end do

    do ipan = 1, npan
      qcllptemp(:, :) = qcllp(:, :, ipan)
      qdllptemp(:, :) = qdllp(:, :, ipan)
      qcllp(:, :,ipan) = qcllptemp(:, :)*(qcone + qcone)
      qdllp(:, :,ipan) = qdllptemp(:, :)*(qcone + qcone)
      call cqgemm(lmsize, lmsize, lmsize, -qcone, qcllptemp, lmsize, qdllp(1,1,0), lmsize, qcone, qcllp(1,1,ipan), lmsize)
      call cqgemm(lmsize, lmsize, lmsize, -qcone, qdllptemp, lmsize, qdllp(1,1,0), lmsize, qcone, qdllp(1,1,ipan), lmsize)
    end do

      cllp = qcllp
    do ipan = 1, npan
      qyif(:,:,:) = yif(:,:,:,ipan)
      do icheb = 0, ncheb
        mn = ipan*ncheb + ipan - icheb
        call zgemm('n', 'n', lmsize2, lmsize, lmsize, cone, zif(1,1,icheb,ipan), lmsize2, cllp(1,1,ipan), lmsize, czero, sll(1,1,mn), lmsize2)
        call cqgemm(lmsize2, lmsize, lmsize, qcone, qyif(1,1,icheb), lmsize2, qdllp(1,1,ipan), lmsize, qczero, qsll, lmsize2)

      sll(:,:,mn) = sll(:,:,mn) + qsll(:,:) 

      end do
    end do
    end if
 
    if (idotime==1) call timing_stop('afterlocal')
    if (idotime==1) call timing_start('endstuff')

    if (idotime==1) call timing_stop('endstuff')
    if (idotime==1) call timing_start('checknan')
    if (idotime==1) call timing_stop('checknan')
    if (idotime==1) call timing_stop('sll-glob')

    deallocate (work, betainv, betainv_save, cllp, dllp, cllptemp, dllptemp, mihvy, mihvz, mijvy, mijvz, yif, zif, stat=ierror)
    if (ierror/=0) stop '[sll-glob] ERROR in deallocating arrays'
    if(use_cheby_quadprec) deallocate (qbetainv, qbetainv_save, qcllp, qdllp, qcllptemp, qdllptemp, qmihvy, qmihvz, qmijvy, qmijvz, qyif, qsll, stat=ierror)
    if (ierror/=0) stop '[sll-glob] ERROR in deallocating arrays'
  end subroutine sll_global_solutions

end module mod_sll_global_solutions

      SUBROUTINE CQGEMM (M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
      IMPLICIT NONE
!     .. Scalar Arguments ..
      INTEGER            M, N, K, LDA, LDB, LDC
      COMPLEX*32         ALPHA, BETA
!     .. Array Arguments ..
      COMPLEX*32         A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!     .. Local Scalars ..
      COMPLEX*32         TEMP
      INTEGER            I, J, L
!     .. Parameters ..
      COMPLEX*32         ONE
      PARAMETER        ( ONE  = ( 1.0Q+0, 0.0Q+0 ) )
      COMPLEX*32         ZERO
      PARAMETER        ( ZERO = ( 0.0Q+0, 0.0Q+0 ) )
!     ..
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
      RETURN
      END
