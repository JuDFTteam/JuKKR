!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_dlke1

contains

  !-------------------------------------------------------------------------------
  !> Summary: Fourier transformation of the cluster Greens function GINP
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Fourier transformation of the cluster Greens function GINP
  !-------------------------------------------------------------------------------
  subroutine dlke1(gllke, alat, nacls, naclsmax, rr, ezoa, atom, bzkp, ic, ginp, rcls)

    use :: mod_datatypes, only: dp
    use :: global_variables, only: almgf0, nrd, lmgf0d, naezd 
    use :: mod_constants, only: ci
    use :: mod_types, only: t_inc
    use :: mod_cinit, only: cinit
    implicit none

    ! ..
    ! .. Local Scalars ..
    real (kind=dp) :: alat
    integer :: ic, naclsmax
    ! ..
    ! .. Local Arrays ..
    integer :: atom(*), ezoa(*), nacls(*)
    complex (kind=dp) :: gllke(almgf0, *), ginp(lmgf0d*naclsmax, *)
    real (kind=dp) :: bzkp(*), rr(3, 0:nrd), rcls(3, *)
    ! ..
    real (kind=dp) :: convpu, tpi
    integer :: am, i, ii, im, lm2, m
    complex (kind=dp) :: eikr, tt
    ! ..
    complex (kind=dp) :: arg(3)
    ! ..
    logical :: opt, test
    external :: test, opt

    ii = 3
    if (opt('COMPLEX ')) ii = 6
    if (test('BZKP    ') .and. (t_inc%i_write>0)) write (1337, fmt='(6f12.6)')(bzkp(i), i=1, ii)

    tpi = 8.0d0*atan(1.0d0)
    convpu = alat/tpi

    call cinit(lmgf0d*naezd*lmgf0d, gllke)

    do m = 1, nacls(ic)

      ! --->  for option 'WIRE': avoid artifical couplings in the structural
      ! Greens Function in in-plane-direction (perp. to c-axis)

      
      if (atom(m)<0) cycle

      if (opt('ONEBULK ')) then    ! added 1.02.2000, corrected on 25.02.2000

        !     if the phase factor exp(ik(r-r')) is included      ~
        !     in the G...so if we resolve the dyson part for the G
        !     and not for the G (see Peter Lang Ph.D thesis)
        !     
        !     Here we do   --                           nn'
        !                  \                            ii'          ii'
        !                  /  exp( ik(X  -X + R  -R  ))G   (E)  =   G   (k,E)
        !                  --          n'  n   i'  i    LL'          LL'
        !                  n'
        !                   
        ! In this case rcls is always (by constraction symmetric around each
        ! atom this means that a minus sign will not affect the result of the
        ! summation

        arg(1) = -ci*tpi*rcls(1, m)
        arg(2) = -ci*tpi*rcls(2, m)
        arg(3) = -ci*tpi*rcls(3, m)

      else ! opt('ONEBULK ')
        
        !     Here we do   --                  nn'
        !                  \                   ii'          ii'
        !                  /  exp(+ik(x  -x ))G   (E)  =   G   (k,E)
        !                  --          n'  n   LL'          LL'
        !                  n'
        !  Be carefull a minus sign must be included here. RR is not
        !  symmetric around each atom. The minus comes from the fact that
        !  the repulsive potential GF is calculated for 0n and not n0!                   
        !  and that is why we nead a minus sign extra!

        arg(1) = -ci*tpi*rr(1, ezoa(m))
        arg(2) = -ci*tpi*rr(2, ezoa(m))
        arg(3) = -ci*tpi*rr(3, ezoa(m))

      end if ! opt('ONEBULK ')

      tt = bzkp(1)*arg(1) + bzkp(2)*arg(2) + bzkp(3)*arg(3)

      if (opt('COMPLEX ')) then
        tt = tt + ci*(bzkp(4)*arg(1)+bzkp(5)*arg(2)+bzkp(6)*arg(3))
      end if

      eikr = exp(tt)*convpu  ! convert to p.u.

      im = 1 + (m-1)*lmgf0d
      am = 1 + (atom(m)-1)*lmgf0d
      do lm2 = 1, lmgf0d
        call zaxpy(lmgf0d, eikr, ginp(im,lm2), 1, gllke(am,lm2), 1)
      end do

    end do  ! M = 1,NACLS(IC)

    return
100 format (3f12.4)
110 format (2f18.10)
  end subroutine dlke1

end module mod_dlke1
