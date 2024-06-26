!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_inversion

contains

  !-------------------------------------------------------------------------------
  !> Summary: Driver of matrix inversion used in Dyson equation 
  !> Author: 
  !> Category: KKRhost, structural-greensfunction
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This subroutine calculates the inversion of a matrix
  !> in 4 different ways depending on the form of the matrix
  !>
  !> INVMOD = 0  ----> total inversion scheme
  !> INVMOD = 1  ----> band matrix inversion scheme
  !> INVMOD = 2  ----> corner band matrix inversion scheme
  !> INVMOD = 3  ----> godfrin module
  !-------------------------------------------------------------------------------
  subroutine inversion(gllke, invmod, icheck)

    use :: mod_datatypes, only: dp
    use :: global_variables, only: alm, ndim_slabinv, nlayerd, naezd, nprincd, lmmaxd
    use :: godfrin, only: sparse_inverse, t_godfrin
    use :: mod_invslab, only: invslab
    use :: mod_invsupercell, only: invsupercell
    use :: mod_constants, only: czero, cone
    implicit none

    complex (kind=dp) :: gllke(alm, alm), gdi(ndim_slabinv, ndim_slabinv, nlayerd), gup(ndim_slabinv, ndim_slabinv, nlayerd), gdow(ndim_slabinv, ndim_slabinv, nlayerd)
    complex (kind=dp), allocatable :: gtemp(:, :)
    integer :: i, i1, ip1, ii1, il1, ldi1, ip2, ii2, il2, ldi2, j, invmod
    integer :: lm1, lm2, info, ipvt(alm)
    integer :: icheck(naezd/nprincd, naezd/nprincd)
    ! total matrix inversion

    allocate (gtemp(alm,alm))

    ! ---------------------------------------------------------------------
    ! full matrix inversion
    ! ---------------------------------------------------------------------
    if (invmod==0) then

      ! initialize unit matrix
      do i = 1, alm
        do j = 1, alm
          gtemp(i, j) = czero
          if (i==j) then
            gtemp(i, j) = cone
          end if
        end do
      end do

      call zgetrf(alm, alm, gllke, alm, ipvt, info)
      call zgetrs('N', alm, alm, gllke, alm, ipvt, gtemp, alm, info)
      call zcopy(alm*alm, gtemp, 1, gllke, 1)

      ! ---------------------------------------------------------------------
      ! block tridiagonal inversion (slab or supercell)
      ! ---------------------------------------------------------------------
    else if ((invmod==1) .or. (invmod==2)) then
      ! copy blocks of assumed tridiagonal matrix
      do i1 = 1, nlayerd
        do ip1 = 1, nprincd
          do ip2 = 1, nprincd
            ii1 = (i1-1)*nprincd + ip1
            ii2 = (i1-1)*nprincd + ip2
            do lm1 = 1, lmmaxd
              do lm2 = 1, lmmaxd
                ldi1 = lmmaxd*(ip1-1) + lm1
                il1 = lmmaxd*(ii1-1) + lm1
                ldi2 = lmmaxd*(ip2-1) + lm2
                il2 = lmmaxd*(ii2-1) + lm2
                gdi(ldi1, ldi2, i1) = gllke(il1, il2)
              end do
            end do
          end do
        end do
      end do
      ! this part now is correct also for    ! changes 20/10/99
      ! supercell geometry : 20/10/99
      ! --> upper linear part
      do i1 = 1, nlayerd
        do ip1 = 1, nprincd
          do ip2 = 1, nprincd
            do lm1 = 1, lmmaxd
              do lm2 = 1, lmmaxd
                ldi1 = lmmaxd*(ip1-1) + lm1
                ldi2 = lmmaxd*(ip2-1) + lm2
                if (i1<=(nlayerd-1)) then
                  ii1 = (i1-1)*nprincd + ip1
                  ii2 = i1*nprincd + ip2
                  il1 = lmmaxd*(ii1-1) + lm1
                  il2 = lmmaxd*(ii2-1) + lm2
                  gup(ldi1, ldi2, i1) = gllke(il1, il2)
                else
                  ii1 = ip1
                  ii2 = (nlayerd-1)*nprincd + ip2
                  il1 = lmmaxd*(ii1-1) + lm1
                  il2 = lmmaxd*(ii2-1) + lm2
                  gdow(ldi1, ldi2, i1) = gllke(il1, il2)
                end if
              end do
            end do
          end do
        end do
      end do
      ! --> lower linear part
      do i1 = 1, nlayerd
        do ip1 = 1, nprincd
          do ip2 = 1, nprincd
            do lm1 = 1, lmmaxd
              do lm2 = 1, lmmaxd
                ldi1 = lmmaxd*(ip1-1) + lm1
                ldi2 = lmmaxd*(ip2-1) + lm2
                if (i1<=(nlayerd-1)) then
                  ii1 = i1*nprincd + ip1
                  ii2 = (i1-1)*nprincd + ip2
                  il1 = lmmaxd*(ii1-1) + lm1
                  il2 = lmmaxd*(ii2-1) + lm2
                  gdow(ldi1, ldi2, i1) = gllke(il1, il2)
                else
                  ii1 = (nlayerd-1)*nprincd + ip1
                  ii2 = ip2
                  il1 = lmmaxd*(ii1-1) + lm1
                  il2 = lmmaxd*(ii2-1) + lm2
                  gup(ldi1, ldi2, i1) = gllke(il1, il2)
                end if
              end do
            end do
          end do
        end do
      end do
      ! end of the corrected part  20/10/99

      ! slab: matrix is block tridiagonal
      if (invmod==1) then
        ! write (6,*) '-------slab calculation--------'
        call invslab(gdi, gup, gdow, gllke, icheck)
      else if (invmod==2) then
        ! write (6,*) '-------supercell calculation--------'
        call invsupercell(gdi, gup, gdow, gllke, icheck)
      end if

      ! -----------------------------------------------------------------
      ! godfrin module
      ! -----------------------------------------------------------------
    else if (invmod==3) then
      ! write (*, '("na=",i8,"  nb=",i8,"  ldiag=",l2,"  lper=",l2,"  lpardiso=",l2)') t_godfrin%na, t_godfrin%nb, t_godfrin%ldiag, t_godfrin%lper, t_godfrin%lpardiso
      ! write (*, '("bdims(1:nb)=",100i8)') t_godfrin%bdims(:)
      call sparse_inverse(gllke, t_godfrin%na, t_godfrin%nb, t_godfrin%bdims, t_godfrin%ldiag, t_godfrin%lper, t_godfrin%lpardiso) ! GODFRIN Flaviano
      ! ------------------------------------------------------------------
    else
      ! if it gets here, did you have a coffee before running the code?
      stop 'UNKNOWN INVERSION MODE !!!'
    end if
    ! ------------------------------------------------------------------

    deallocate (gtemp)

    return

  end subroutine inversion

end module mod_inversion
