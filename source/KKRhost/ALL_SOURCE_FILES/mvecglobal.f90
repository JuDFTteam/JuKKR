!------------------------------------------------------------------------------------
!> Summary: Transforms the magnetization to the cartesian global frame of reference 
!> Author: 
!>Transforms the magnetization to the cartesian global frame of reference, first
!> by transforming from the local $$\pm z$$ coordinate system, to the local 
!> cartesian coordinate system, and then to the gloabl reference frame.
!------------------------------------------------------------------------------------
!> @note This routine has been build up from the last part of the original
!> Munich `CALCMVEC` routine
!> @endnote
!------------------------------------------------------------------------------------
module mod_mvecglobal
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Transforms the magnetization to the cartesian global frame of reference 
  !> Author: 
  !> Category: numerical-tools, physical-observables, KKRhost 
  !> Deprecated: False 
  !> Transforms the magnetization to the cartesian global frame of reference, first
  !> by transforming from the local $$\pm z$$ coordinate system, to the local 
  !> cartesian coordinate system, and then to the gloabl reference frame.
  !-------------------------------------------------------------------------------
  !> @note This routine has been build up from the last part of the original
  !> Munich `CALCMVEC` routine
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine mvecglobal(it,iq,natyp,qmphi,qmtet,mvevi,mvevil,mvevief,natypd,lmaxd,  &
    nmvecmax)
    ! ********************************************************************
    ! *                                                                  *
    ! *  this routine has been build up from the last part of the        *
    ! *  original Munich CALCMVEC routine.                               *
    ! *  on exit, MVEVI,MVEVIL,MVEVIEF are in the CARTESIAN GLOBAL       *
    ! *                                           frame of reference     *
    ! *                                                                  *
    ! ********************************************************************
    use :: mod_datatypes
    use :: mod_calcrotmat
    use :: constants, only : ci,czero,pi
    implicit none

    ! .. Parameter definitions
    integer :: lmaxdloc
    parameter (lmaxdloc=8)

    ! .. Input variables
    integer, intent(in) :: it !! Index of the current atom type
    integer, intent(in) :: iq !! Chemical type of the current atom
    integer, intent(in) :: lmaxd  !! Maximum l component in wave function expansion
    integer, intent(in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent(in) :: natypd !! Number of kinds of atoms in unit cell
    integer, intent(in) :: nmvecmax !! 4
    real (kind=dp), intent(in) :: qmphi !! $$ \phi $$ angle of the agnetization with respect to the z-axis
    real (kind=dp), intent(in) :: qmtet !! $$ \theta $$ angle of the agnetization with respect to the z-axis
    !.. In/Out variables
    complex (kind=dp), dimension(natypd, 3, nmvecmax), intent(inout) :: mvevi
    complex (kind=dp), dimension(natypd, 3, nmvecmax), intent(inout) :: mvevief
    complex (kind=dp), dimension(0:lmaxd, natypd, 3, nmvecmax), intent(inout) :: mvevil

    ! Local Scalars
    integer :: icall, i, j, k, l, imv, nmvec
    complex (kind=dp) :: cs
    complex (kind=dp) :: amin, apls
    real (kind=dp) :: mv, mvx, mvxy, mvy, mvz, wsq2

    ! Local Arrays
    complex (kind=dp), dimension(3,3) :: usc
    complex (kind=dp), dimension(3,3) :: w3x3
    complex (kind=dp), dimension(4,4) :: drot4
    complex (kind=dp), dimension(3,nmvecmax) :: mvg
    complex (kind=dp), dimension(3,nmvecmax) :: mvgef
    complex (kind=dp), dimension(0:lmaxd, 3, nmvecmax) :: mvgl
    real (kind=dp), dimension(nmvecmax) :: mvphi
    real (kind=dp), dimension(nmvecmax) :: mvtet
    real (kind=dp), dimension(0:100)    :: fact
    real (kind=dp), dimension(3,3)          :: mrot
    real (kind=dp), dimension(3, nmvecmax)  :: mvglo
    real (kind=dp), dimension(0:lmaxd, 3, nmvecmax) :: mvglol
    character (len=1), dimension(0:lmaxdloc) :: txtl

    ! Data Statements
    data icall/0/

    ! Save Statements
    save :: icall, nmvec, usc, fact, txtl

    icall = icall + 1
    ! =======================================================================
    if (icall==1) then

      if (lmaxd>lmaxdloc) then
        write (6, *)
        write (6, *) ' Please increase parameter LMAXDLOC to ', lmaxd
        write (6, *) ' in the < MVECGLOBAL > routine.'
        stop ' < TBKKR2 > '
      end if

      txtl(0) = 's'
      txtl(1) = 'p'
      txtl(2) = 'd'
      if (lmaxd>=3) then
        do l = 3, lmaxd
          txtl(l) = char(ichar('f')+l-3)
        end do
      end if

      write (1337, '(78("#"))')
      write (1337, 100)
      write (1337, '(78("#"))')
      write (1337, *)
      write (1337, 110)

      nmvec = 2

      fact(0) = 1.0e0_dp
      do i = 1, 100
        fact(i) = fact(i-1)*real(i, kind=dp)
      end do
      ! -----------------------------------------------------------------------
      ! create transformation matrix   U  cartesian/sperical ccordinates
      ! -----------------------------------------------------------------------
      ! RC,RCP  vectors in cartesian coordinates
      ! RS,RSP  vectors in spherical coordinates
      ! RS  = USC * R!..                          (4.40)
      ! RSP = MS  * RS                                 (4.37)
      ! MS(i,j) = D(j,i)                                   (4.42)
      ! D  rotation matrix for complex spherical harmonics

      ! ordering of: m=-1,0,+1 >>> row 1 and 3 interchanged compared to (4.44)

      wsq2 = 1.0e0_dp/sqrt(2.0e0_dp)

      usc(1, 1) = wsq2
      usc(1, 2) = -ci*wsq2
      usc(1, 3) = 0.0e0_dp
      usc(2, 1) = 0.0e0_dp
      usc(2, 2) = 0.0e0_dp
      usc(2, 3) = 1.0e0_dp
      usc(3, 1) = -wsq2
      usc(3, 2) = -ci*wsq2
      usc(3, 3) = 0.0e0_dp
      ! -----------------------------------------------------------------------
    end if
    ! =======================================================================

    ! -----------------------------------------------------------------------
    ! create the rotation matrices  DROT4 for complex spherical harmonics
    ! -----------------------------------------------------------------------
    call calcrotmat(2, 1, qmphi, qmtet, 0.0e0_dp, drot4, fact, 4)
    ! -----------------------------------------------------------------------
    ! create the rotation matrix  MROT for vectors in cartesian coordinates
    ! NOTE:  U^+ D^T U gives the inverse of the real matrix  M
    ! for that reason  the transposed matrix is stored as  MROT(J,I)
    ! -----------------------------------------------------------------------
    do i = 1, 3
      do j = 1, 3
        cs = 0.0e0_dp
        do k = 1, 3
          cs = cs + drot4(k+1, i+1)*usc(k, j)
        end do
        w3x3(i, j) = cs
      end do
    end do

    do i = 1, 3
      do j = 1, 3
        cs = 0.0e0_dp
        do k = 1, 3
          cs = cs + conjg(usc(k,i))*w3x3(k, j)
        end do
        if (aimag(cs)>1e-8_dp) write (*, *) ' MROT', i, j, cs, ' ???????????'
        ! see above >> MROT(I,J) = DREAL(CS)
        mrot(j, i) = real(cs)
      end do
    end do
    ! -----------------------------------------------------------------------

    ! **********************************************************************
    do imv = 1, nmvec
      ! -----------------------------------------------------------------------
      ! transform from (+,-,z) to cartesian coordinates  (x,y,z)
      ! note the convention
      ! -----------------------------------------------------------------------
      apls = mvevi(it, 1, imv)
      amin = mvevi(it, 2, imv)
      mvevi(it, 1, imv) = (amin+apls)*0.5e0_dp
      mvevi(it, 2, imv) = (amin-apls)*0.5e0_dp*ci

      apls = mvevief(it, 1, imv)
      amin = mvevief(it, 2, imv)
      mvevief(it, 1, imv) = (amin+apls)*0.5e0_dp
      mvevief(it, 2, imv) = (amin-apls)*0.5e0_dp*ci

      do l = 0, lmaxd
        apls = mvevil(l, it, 1, imv)
        amin = mvevil(l, it, 2, imv)
        mvevil(l, it, 1, imv) = (amin+apls)*0.5e0_dp
        mvevil(l, it, 2, imv) = (amin-apls)*0.5e0_dp*ci
      end do
      ! -----------------------------------------------------------------------
      ! transform from LOCAL cartesian coordinates (x,y,z)
      ! to GLOBAL cartesian coordinates
      ! -----------------------------------------------------------------------
      do i = 1, 3
        mvg(i, imv) = czero
        mvgef(i, imv) = czero
        do j = 1, 3
          mvg(i, imv) = mvg(i, imv) + mrot(i, j)*mvevi(it, j, imv)
          mvgef(i, imv) = mvgef(i, imv) + mrot(i, j)*mvevief(it, j, imv)
        end do
        mvglo(i, imv) = aimag(mvg(i,imv))

        do l = 0, lmaxd
          mvgl(l, i, imv) = czero
          do j = 1, 3
            mvgl(l, i, imv) = mvgl(l, i, imv) + mrot(i, j)*mvevil(l, it, j, imv)
          end do
          mvglol(l, i, imv) = aimag(mvgl(l,i,imv))
        end do

      end do
      ! ......................................................................
      do i = 1, 3
        mvevi(it, i, imv) = mvg(i, imv)
        mvevief(it, i, imv) = mvgef(i, imv)

        do l = 0, lmaxd
          mvevil(l, it, i, imv) = mvgl(l, i, imv)
        end do
      end do
      ! -----------------------------------------------------------------------
      ! calculate the angles
      ! -----------------------------------------------------------------------
      mvx = mvglo(1, imv)
      mvy = mvglo(2, imv)
      mvz = mvglo(3, imv)

      mv = sqrt(mvx**2+mvy**2+mvz**2)
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (mv<1e-8_dp) then
        mvphi(imv) = 0e0_dp
        mvtet(imv) = 0e0_dp
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      else
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        mvxy = sqrt(mvx**2+mvy**2)
        ! ======================================================================
        if (abs(mvxy)<1e-8_dp) then
          mvphi(imv) = 0e0_dp
          ! ======================================================================
        else
          ! ======================================================================
          if (mvy>=0e0_dp) then
            mvphi(imv) = acos(mvx/mvxy)
          else if (mvx<0e0_dp) then
            mvphi(imv) = pi + acos(-mvx/mvxy)
          else
            mvphi(imv) = 2*pi - acos(mvx/mvxy)
          end if
          mvphi(imv) = mvphi(imv)*180e0_dp/pi
          if (abs(mvphi(imv)-360.0e0_dp)<1e-8_dp) mvphi(imv) = 0e0_dp
        end if
        ! ======================================================================
        if (mvphi(imv)>=345.e0_dp) mvphi(imv) = 360.e0_dp - mvphi(imv)
        mvtet(imv) = acos(mvz/mv)*180e0_dp/pi
      end if
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! -----------------------------------------------------------------------
    end do
    ! **********************************************************************
    ! output vector components, in and out angles
    ! ----------------------------------------------------------------------
    l = 0
    write (1337, 120) it, iq, txtl(l), ((mvglol(l,i,imv),i=1,3), imv=1, 2)
    write (1337, 130)(txtl(l), ((mvglol(l,i,imv),i=1,3),imv=1,2), l=1, lmaxd)
    write (1337, 140)((mvglo(i,imv),i=1,3), imv=1, 2)
    write (1337, 150) qmphi, qmtet, (mvphi(imv), mvtet(imv), imv=1, 2)
    ! ----------------------------------------------------------------------
    if (it<natyp) then
      write (1337, '(3X,75("="))')
    else
      write (1337, *)
      write (1337, '(78("#"))')
    end if
    ! ----------------------------------------------------------------------

100 format (15x, 'vectorial magnetic properties given with respect', /, 15x, '   to the GLOBAL (crystal) frame of reference')
110 format (29x, 'm_spin', 27x, 'm_orb', /, 3x, 'ATOM/SITE     ', '    x         y         z', 8x, '    x         y         z', /, 3x, 75('='))
120 format (3x, i3, '/', i3, 2x, a1, ' =', 3f10.5, 3x, 3f10.5)
130 format (12x, a1, ' =', 3f10.5, 3x, 3f10.5)
140 format (12x, 66('-'), /, 12x, 'sum', 3f10.5, 3x, 3f10.5, /)
150 format (3x, 'angles (IN)   TET =', f9.4, ' PHI =', f9.4, /, 3x, 'angles (calc) TET =', f9.4, ' PHI =', f9.4, '   TET =', f9.4, ' PHI =', f9.4)
  end subroutine mvecglobal

end module mod_mvecglobal
