module mod_vinterface

contains

! -------------------------------------------------------------------------------
! SUBROUTINE: VINTERFACE
! > @brief This is calculating the intra-atomic contibution of the potential
! in
! >  the case of an interface taking into account the bulk potential on
! >  the two sides.

! > @details It uses the structure dependent matrices AVMAD which are
! calculated
! >  once in the subroutine MADELUNG2D() and saved in the DA-file
! >  avmad.unformatted ( May 2004)
! >
! >  For each site in a layer the summation in all other layers is split
! >  into three parts: within the slab, over the NLEFT*NLBASIS left host
! >  sites and over the NRIGHT*NRBASIS right host sites, the last two
! >  steps only in case of decimation run

! -------------------------------------------------------------------------------
! > @note
! > - Adapted for the case of more atoms on the same site, summation is
! >  done over the occupants of that site, the charge is weighted with
! >  the appropriate concentration of the occupant  V. Popescu feb. 2002
! -------------------------------------------------------------------------------
! >
! > - Impurity-program adopted feb. 2004 (according to N. Papanikalou)
! >
! > - Jonathan Chico Feb. 2018: Removed inc.p dependencies and rewrote to
! Fortran90
! -------------------------------------------------------------------------------
subroutine vinterface(cmom, cminst, lpot, nspin, nlayers, natyp, v, zat, r, &
  irws, ircut, ipan, kshape, noq, kaoez, iqat, conc, catom, icc, hostimp, &
  nlbasis, nleft, nrbasis, nright, cmomhost, chrgnt, vinters, naez, lmpot)

  use :: constants
  use :: global_variables
  use :: mod_datatypes, only: dp

  implicit none

  ! .. Input variables ..
  integer, intent (in) :: icc      ! < Enables the calculation of off-diagonal
                                   ! elements of the GF.(0=SCF/DOS; 1=cluster;
                                   ! -1=custom)
  integer, intent (in) :: lpot     ! < Maximum l component in potential
                                   ! expansion
  integer, intent (in) :: naez     ! < Number of atoms in unit cell
  integer, intent (in) :: lmpot    ! < (LPOT+1)**2
  integer, intent (in) :: nspin    ! < Counter for spin directions
  integer, intent (in) :: natyp    ! < Number of kinds of atoms in unit cell
  integer, intent (in) :: nleft    ! < Number of repeated basis for left host
                                   ! to get converged electrostatic potentials
  integer, intent (in) :: nright   ! < Number of repeated basis for right host
                                   ! to get converged electrostatic potentials
  integer, intent (in) :: kshape   ! < Exact treatment of WS cell
  integer, intent (in) :: nlayers
  integer, intent (in) :: nlbasis  ! < Number of basis layers of left host
                                   ! (repeated units)
  integer, intent (in) :: nrbasis  ! < Number of basis layers of right host
                                   ! (repeated units)
  real (kind=dp), intent (in) :: chrgnt
  integer, dimension (naez), intent (in) :: noq ! < Number of diff. atom types
                                                ! located
  integer, dimension (natyp), intent (in) :: irws ! < R point at WS radius
  integer, dimension (natyp), intent (in) :: ipan ! < Number of panels in
                                                  ! non-MT-region
  integer, dimension (natyp), intent (in) :: iqat ! < The site on which an
                                                  ! atom is located on a given
                                                  ! site
  integer, dimension (0:natyp), intent (in) :: hostimp
  integer, dimension (0:ipand, natyp), intent (in) :: ircut ! < R points of
                                                            ! panel borders
  integer, dimension (natyp, naez+nembd1-1), intent (in) :: kaoez ! < Kind of
                                                                  ! atom at
                                                                  ! site in
                                                                  ! elem. cell
  real (kind=dp), dimension (natyp), intent (in) :: zat ! < Nuclear charge
  real (kind=dp), dimension (natyp), intent (in) :: conc ! < Concentration of
                                                         ! a given atom
  real (kind=dp), dimension (natyp), intent (in) :: catom
  real (kind=dp), dimension (irmd, natyp), intent (in) :: r ! < Radial mesh (
                                                            ! in units a Bohr)
  real (kind=dp), dimension (lmpot, natyp), intent (in) :: cmom ! < LM moment
                                                                ! of total
                                                                ! charge
  real (kind=dp), dimension (lmpot, natyp), intent (in) :: cminst ! < charge
                                                                  ! moment of
                                                                  ! interstitial
  real (kind=dp), dimension (lmpot, nembd1), intent (in) :: cmomhost ! <
                                                                     ! Charge
                                                                     ! moments
                                                                     ! of each
                                                                     ! atom of
                                                                     ! the
                                                                     ! (left/right)
                                                                     ! host
  ! .. In/out variables
  real (kind=dp), dimension (irmd, lmpot, npotd), intent (inout) :: v

  ! .. Local variables
  integer :: ileft, iright
  integer :: i, iatom, ib, ih1, ilay1, ilay2, io2
  integer :: ipot, irs1, ispin, it1, it2, l, lm, lm2, m
  integer :: lrecamad, irec, nleftoff, nrightoff, nleftall, nrightall
  real (kind=dp) :: cm1, fpi
  logical :: opt, test, lread
  real (kind=dp), dimension (lmpot) :: ac
  real (kind=dp), dimension (lmpot) :: cm
  real (kind=dp), dimension (2) :: charge
  real (kind=dp), dimension (naez) :: monopol
  real (kind=dp), dimension (lmpot, lmpot) :: avmad
  real (kind=dp), dimension (lmpot, naez) :: vinters

  external :: opt, test

  if (test('flow    ')) write (1337, *) '>>>>>> Vinterface'

  inquire (file='avmad.unformatted', exist=lread) ! ewald2d

  if (lread) then
    lrecamad = wlength*2*lmpot*lmpot
    open (69, access='direct', recl=lrecamad, file='avmad.unformatted', &
      form='unformatted')
  else
    lrecamad = wlength*2*lmpot*lmpot + wlength*2*lmpot
    open (69, access='direct', recl=lrecamad, file='abvmad.unformatted', &
      form='unformatted')
  end if

  write (1337, fmt=100)
  write (1337, fmt=110)

  fpi = 4.e0_dp*pi

  if (opt('DECIMATE')) then
    ! -------------------------------------------------------------------------
    ! Setup the charges to put in the ghost layers in the case of
    ! decimation technique to achieve charge neutrality
    ! -------------------------------------------------------------------------
    charge(1) = -chrgnt/(2.e0_dp*sqrt(fpi))
    charge(2) = -chrgnt/(2.e0_dp*sqrt(fpi))

    nleftoff = nlayers*nlayers     ! record offsets
    nrightoff = nleftoff + nlayers*nleft*nlbasis ! left and right
    nleftall = nleft*nlbasis
    nrightall = nright*nrbasis
  end if
  ! ----------------------------------------------------------------------------
  ! START CALCULATION IN THE LAYERS
  ! ----------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------
  ! Loop over atoms in slab
  ! ----------------------------------------------------------------------------

  do it1 = 1, natyp
    ! -------------------------------------------------------------------------
    ! Take a site occupied by IT1
    ! -------------------------------------------------------------------------
    ilay1 = iqat(it1)

    if (kshape/=0) then
      irs1 = ircut(ipan(it1), it1)
    else
      irs1 = irws(it1)
    end if

    do lm = 1, lmpot
      ac(lm) = 0.e0_dp
    end do
    ! -------------------------------------------------------------------------
    ! 1.  Summation in all layers in the slab
    ! -------------------------------------------------------------------------
    do ilay2 = 1, nlayers
      irec = ilay2 + nlayers*(ilay1-1)
      read (69, rec=irec) avmad

      ! ----------------------------------------------------------------------
      ! Keep the monopole term -- Hoshino is doing (SMONOPOL(I) -SMONOPOL(0))
      ! ----------------------------------------------------------------------
      if (ilay1==ilay2) monopol(ilay1) = avmad(1, 1)
      ! ----------------------------------------------------------------------
      ! Loop over all occupants of site ILAY2
      ! ----------------------------------------------------------------------
      do io2 = 1, noq(ilay2)
        it2 = kaoez(io2, ilay2)

        do lm = 1, lmpot
          cm(lm) = cmom(lm, it2)
          ! ----------------------------------------------------------------
          ! Add contribution of interstial in case of shapes
          ! ----------------------------------------------------------------
          if (kshape/=0) cm(lm) = cm(lm) + cminst(lm, it2)
        end do
        cm(1) = cm(1) - zat(it2)/sqrt(fpi)

        do lm = 1, lmpot
          do lm2 = 1, lmpot
            ac(lm) = ac(lm) + avmad(lm, lm2)*cm(lm2)*conc(it2)
          end do
        end do
      end do
      ! ----------------------------------------------------------------------
      ! Loop over all occupants of site ILAY2
      ! ----------------------------------------------------------------------
    end do                         ! ILAY2 loop in all interface planes
    ! -------------------------------------------------------------------------
    do ilay2 = 1, nlayers
      ! ----------------------------------------------------------------------
      ! Loop over all occupants of site ILAY2
      ! ----------------------------------------------------------------------
      do io2 = 1, noq(ilay2)
        it2 = kaoez(io2, ilay2)

        cm1 = cmom(1, it2)
        if (kshape/=0) cm1 = cm1 + cminst(1, it2)

        cm1 = cm1 - zat(it2)/sqrt(fpi)
        ac(1) = ac(1) - monopol(ilay1)*cm1*conc(it2)

      end do
      ! ----------------------------------------------------------------------
      ! Loop over all occupants of site ILAY2
      ! ----------------------------------------------------------------------
    end do
    ! -------------------------------------------------------------------------
    ! Correction: charge neutrality is imposed (see P. Lang)
    ! -------------------------------------------------------------------------
    if (opt('DECIMATE')) then
      ! ----------------------------------------------------------------------
      ! 2.  Summation in the LEFT bulk side
      ! ----------------------------------------------------------------------
      ! ----------------------------------------------------------------------
      ! Loop over all occupants of LEFT host
      ! ----------------------------------------------------------------------
      ileft = 0
      do ih1 = 1, nleft
        do ib = 1, nlbasis
          ileft = ileft + 1
          irec = ileft + nleftall*(ilay1-1) + nleftoff
          read (69, rec=irec) avmad

          iatom = ib
          do lm = 1, lmpot
            do lm2 = 1, lmpot
              ac(lm) = ac(lm) + avmad(lm, lm2)*cmomhost(lm2, iatom)
            end do
          end do

          if ((ih1==1) .and. (ib==1)) then
            ac(1) = ac(1) + (avmad(1,1)-monopol(ilay1))*charge(1)
          end if
        end do
      end do
      ! ----------------------------------------------------------------------
      if (ileft/=nleftall) then
        write (6, *) ' < VINTERFACE > : index error ', &
          'ILEFT <> NLEFT*NLBASIS'
        stop
      end if
      ! ----------------------------------------------------------------------
      ! 3.  Summation in the RIGHT bulk side
      ! ----------------------------------------------------------------------
      ! ----------------------------------------------------------------------
      ! Loop over all occupants of RIGHT host
      ! ----------------------------------------------------------------------
      iright = 0
      do ih1 = 1, nright
        do ib = 1, nrbasis
          iright = iright + 1
          irec = iright + nrightall*(ilay1-1) + nrightoff
          read (69, rec=irec) avmad

          iatom = nlbasis + ib
          do lm = 1, lmpot
            do lm2 = 1, lmpot
              ac(lm) = ac(lm) + avmad(lm, lm2)*cmomhost(lm2, iatom)
            end do
          end do

          if ((ih1==1) .and. (ib==1)) then
            ac(1) = ac(1) + (avmad(1,1)-monopol(ilay1))*charge(2)
          end if
        end do
      end do
      ! ----------------------------------------------------------------------
      if (iright/=nrightall) then
        write (6, *) ' < VINTERFACE > : index error ', &
          'IRIGHT <> NRIGHT*NRBASIS'
        stop
      end if
    end if                         ! (OPT(DECIMATE)
    ! -------------------------------------------------------------------------
    write (1337, fmt=120) it1, (catom(it1)-zat(it1)), &
      (ac(1)/sqrt(4.e0_dp*pi)), (ac(3)/sqrt(4.e0_dp*pi))
    ! -------------------------------------------------------------------------
    ! Loop over spins of atom IT1
    ! -------------------------------------------------------------------------
    do ispin = 1, nspin
      ! ----------------------------------------------------------------------
      ! Determine the right potential number
      ! ----------------------------------------------------------------------
      ipot = nspin*(it1-1) + ispin
      ! ----------------------------------------------------------------------
      ! In the case of l=0 : r(1)**l is not defined
      ! ----------------------------------------------------------------------
      v(1, 1, ipot) = v(1, 1, ipot) + ac(1)

      do l = 0, lpot
        do m = -l, l
          lm = l*l + l + m + 1
          do i = 2, irs1
            v(i, lm, ipot) = v(i, lm, ipot) + (-r(i,it1))**l*ac(lm)
          end do
        end do
      end do
    end do
    ! -------------------------------------------------------------------------
    ! This part (ICC.GT.0) should be presumably reconsidered for impurity
    ! calculation in host-CPA case
    ! -------------------------------------------------------------------------
    if (icc>0 .or. opt('KKRFLEX ')) then
      do l = 0, lpot
        do m = -l, l
          lm = l*l + l + m + 1
          vinters(lm, ilay1) = ac(lm)
        end do
      end do
    end if
    ! -------------------------------------------------------------------------
  end do
  ! ----------------------------------------------------------------------------
  close (69)
  write (1337, '(15X,45("-"),/)')
  write (1337, '(79("="))')
  if ((icc==0) .and. (.not. opt('KKRFLEX '))) return
  ! ----------------------------------------------------------------------------
  ! Now Prepare output for Impurity calculation
  ! ----------------------------------------------------------------------------
  open (91, file='intercell_ref', status='unknown', form='formatted')
  write (1337, *)
  write (1337, *) '                     ', &
    'Writing intercell potential for impurity'
  write (1337, '(/,20X,55("-"))')
  write (1337, 130) hostimp(0), lmpot
  write (1337, '(20X,55("-"),/,35X,"  i host lm  Vint")')
  do i = 1, hostimp(0)
    write (1337, *)
    lm = 1
    write (1337, '(35X,I4,I4,I3,1X,F10.6)') i, hostimp(i), lm, &
      vinters(lm, hostimp(i))
    do lm = 2, 9
      write (1337, '(43X,I3,1X,F10.6)') lm, vinters(lm, hostimp(i))
    end do
    write (1337, '(20X,55("-"))')
  end do
  write (1337, '(79("="),/)')

  write (91, 140) hostimp(0), lmpot
  do i = 1, hostimp(0)
    write (91, 150)(vinters(lm,hostimp(i)), lm=1, lmpot)
  end do
  close (91)

  return

100 format (79('='), /, 25x, ' INTERFACE MADELUNG POTENTIALS ')
110 format (/, 15x, ' ATOM ', '  Delta_Q  ', '   MONOPOLE       DIPOLE', /, &
    15x, 45('-'))
120 format (15x, i4, 2x, f10.6, 1x, 1p, d13.6, 1x, 1p, d13.6)
130 format (22x, i4, ' host atoms, LMPOT = ', i2, ' output up to LM = 9')
140 format (3i6)
150 format (4d20.10)
end subroutine vinterface

end module mod_vinterface
