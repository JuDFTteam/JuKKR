! 01.06.99 *************************************************************
subroutine dlke1(gllke, alat, nacls, naclsmax, rr, ezoa, atom, bzkp, ic, ginp, &
  rcls)
  ! **********************************************************************

  ! Fourier transformation of the cluster Greens function GINP

  ! ----------------------------------------------------------------------
  use :: mod_types, only: t_inc
  use :: mod_datatypes, only: dp
  use global_variables
  implicit none

  complex (kind=dp) :: ci
  parameter (ci=(0.0d0,1.0d0))
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
  ! .. External Subroutines ..
  real (kind=dp) :: convpu, tpi
  integer :: am, i, ii, im, lm2, m
  complex (kind=dp) :: eikr, tt
  logical :: opt, test
  ! ..
  ! .. Intrinsic Functions ..
  complex (kind=dp) :: arg(3)
  ! ..
  ! .. Save statement ..
  external :: cinit, test, opt, zaxpy
  ! ..

  intrinsic :: atan, exp
  ! = 2*PI

  save

  ii = 3
  if (opt('COMPLEX ')) ii = 6
  if (test('BZKP    ') .and. (t_inc%i_write>0)) write (1337, fmt='(6f12.6)') &
    (bzkp(i), i=1, ii)

  tpi = 8.0d0*atan(1.0d0)
  convpu = alat/tpi
  ! --->  for option 'WIRE': avoid artifical couplings in the structural
  call cinit(lmgf0d*naezd*lmgf0d, gllke)
  ! Greens Function in in-plane-direction (perp. to c-axis)
  do m = 1, nacls(ic)


    ! added 1.02.2000
    ! corrected on 25.02.2000
    ! if the phase factor exp(ik(r-r')) is included      ~
    if (atom(m)<0) cycle
    ! in the G...so if we resolve the dyson part for the G
    if (opt('ONEBULK ')) then      ! and not for the G (see Peter Lang Ph.D
                                   ! thesis)

      ! Here we do   --                           nn'
      ! \                            ii'          ii'
      ! /  exp( ik(X  -X + R  -R  ))G   (E)  =   G   (k,E)
      ! --          n'  n   i'  i    LL'          LL'
      ! n'

      ! In this case rcls is always (by constraction symmetric around each
      ! atom this means that a minus sign will not affect the result of the
      ! summation


      ! Here we do   --                  nn'
      ! \                   ii'          ii'
      ! /  exp(+ik(x  -x ))G   (E)  =   G   (k,E)
      arg(1) = -ci*tpi*rcls(1, m)
      arg(2) = -ci*tpi*rcls(2, m)
      arg(3) = -ci*tpi*rcls(3, m)
    else
      ! --          n'  n   LL'          LL'
      ! n'
      ! Be carefull a minus sign must be included here. RR is not
      ! symmetric around each atom. The minus comes from the fact that
      ! the repulsive potential GF is calculated for 0n and not n0!
      ! and that is why we nead a minus sign extra!





      ! write(6,*) BZKP(1),ARG(1),BZKP(2),ARG(2),BZKP(3),ARG(3)
      arg(1) = -ci*tpi*rr(1, ezoa(m))
      arg(2) = -ci*tpi*rr(2, ezoa(m))
      arg(3) = -ci*tpi*rr(3, ezoa(m))
    end if
    ! write(6,*) BZKP(4),BZKP(5),BZKP(6)
    tt = bzkp(1)*arg(1) + bzkp(2)*arg(2) + bzkp(3)*arg(3)
    ! write(6,*) 'm,atom(m),tt',m,atom(m),tt
    if (opt('COMPLEX ')) then
      tt = tt + ci*(bzkp(4)*arg(1)+bzkp(5)*arg(2)+bzkp(6)*arg(3))
    end if

    ! convert to p.u.

    ! write(6,*) 'eikr',eikr

    eikr = exp(tt)*convpu

    ! 01.06.99 *************************************************************
    ! **********************************************************************
    im = 1 + (m-1)*lmgf0d
    am = 1 + (atom(m)-1)*lmgf0d
    do lm2 = 1, lmgf0d
      call zaxpy(lmgf0d, eikr, ginp(im,lm2), 1, gllke(am,lm2), 1)
    end do

  end do
  ! Fourier transformation of the cluster Greens function GINP
  return
100 format (3f12.4)
110 format (2f18.10)
end subroutine dlke1
