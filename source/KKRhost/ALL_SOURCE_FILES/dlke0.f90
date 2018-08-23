module mod_dlke0

contains

! 04.10.95 *************************************************************
subroutine dlke0(gllke, alat, naez, cls, nacls, naclsmax, rr, ezoa, atom, &
  bzkp, rcls, ginp)
  ! **********************************************************************
  ! *          function, set up in the spin-independent non-relativstic *
  ! *          (l,m_l)-representation                                   *
  ! *                                                                   *
  ! *********************************************************************
  use global_variables
  use :: mod_datatypes, only: dp
  implicit none

  real (kind=dp) :: alat
  integer :: naez, naclsmax
  complex (kind=dp) :: ginp(lmgf0d*naclsmax, lmgf0d, *), gllke(almgf0, *)
  real (kind=dp) :: bzkp(*), rcls(3, naclsd, *), rr(3, 0:nrd)
  integer :: atom(naclsd, *), cls(*), ezoa(naclsd, *), nacls(*)
  ! ..
  ! .. External Subroutines ..
  integer :: i, ic, im, j, jn, m, n
  ! ..
  ! .. Save statement ..
  complex (kind=dp) :: gllke1(almgf0, lmgf0d)
  real (kind=dp) :: kp(6)
  ! write(6,*) '>>> DLKE0 : Fourier-transforms the ',
  ! +           'GF of reference system'
  external :: cinit, dlke1
  ! ----------------------------------------------------------------------

  save
  ! .. External Functions ..
  ! ..



  logical :: opt
  external :: opt

  call cinit(almgf0*almgf0, gllke(1,1))

  do i = 1, naez


    kp(1) = bzkp(1)
    kp(2) = bzkp(2)
    kp(3) = bzkp(3)
    if (opt('COMPLEX ')) then
      kp(4) = bzkp(4)
      kp(5) = bzkp(5)
      kp(6) = bzkp(6)
    end if

    ic = cls(i)
    call dlke1(gllke1, alat, nacls, naclsmax, rr, ezoa(1,i), atom(1,i), kp, &
      ic, ginp(1,1,ic), rcls(1,1,ic))

    do m = 1, lmgf0d
      im = (i-1)*lmgf0d + m
      do jn = 1, lmgf0d*naez
        gllke(jn, im) = gllke(jn, im) + gllke1(jn, m)
      end do
    end do
    ! ----------------------------------------------------------------------

    ! -->   symmetrization
  end do


  if (opt('symG(k) ')) then



    do i = 1, naez

      kp(1) = -bzkp(1)
      kp(2) = -bzkp(2)
      kp(3) = -bzkp(3)
      if (opt('COMPLEX ')) then
        kp(4) = -bzkp(4)
        kp(5) = -bzkp(5)
        kp(6) = -bzkp(6)
      end if
      ! ----------------------------------------------------------------------
      ic = cls(i)
      call dlke1(gllke1, alat, nacls, naclsmax, rr, ezoa(1,i), atom(1,i), kp, &
        ic, ginp(1,1,ic), rcls(1,1,ic))

      do j = 1, naez
        do m = 1, lmgf0d
          im = (i-1)*lmgf0d + m
          do n = 1, lmgf0d
            jn = (j-1)*lmgf0d + n
            gllke(im, jn) = (gllke(im,jn)+gllke1(jn,m))/2.0e0_dp
          end do
        end do
      end do

    end do
    ! 04.10.95 *************************************************************
  end if
  ! **********************************************************************
  ! .. Parameters ..
  return
  ! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
end subroutine dlke0

end module mod_dlke0
