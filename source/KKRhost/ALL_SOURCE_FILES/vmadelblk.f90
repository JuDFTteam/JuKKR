!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculate the madelung potentials and add these to the potential \(V\) (in he spin-polarized case for each spin-direction this is the same)
!> Author: B. Drittler
!> It uses the structure dependent matrices `AVMAD` and `BVMAD` which
!> are calculated once in the subroutine `MADELUNG3D()` and saved in
!> the DA-file abvmad.unformatted (May 2004)
!> The charge-moments are calculated in the subroutine `vintras`,
!> therefore vintras has to be called first.
!> The madelung-potential is expanded into spherical harmonics.
!> The lm-term of the potential \(V\) of the atom \(i\) is given by
!> \begin{equation}
!> V(r,lm,i) = \sum_{i2}^{N} \sum_{l'm'} (-r)^l \left\{avmad(i,i2,lm,l'm')cmom(i2,l'm') +bvmad(i,i2,lm)z(i2)\right\}
!> \end{equation}
!> where \(N\) is the number of atoms
!------------------------------------------------------------------------------------
!> @note
!>  
!> - V. Popescu Feb. 2002: Adopted for the case of more atoms on the same site, 
!> summation is done over the occupants of that site, the charge is weighted with 
!> the appropriate concentration of the occupant
!>
!> - Impurity-program adopted feb. 2004 (according to N. Papanikalou)
!> @endnote
!------------------------------------------------------------------------------------
module mod_vmadelblk

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculate the madelung potentials and add these to the potential \(V\) (in he spin-polarized case for each spin-direction this is the same)
  !> Author: B. Drittler
  !> Category: electrostatics, potential, KKRhost
  !> Deprecated: False 
  !> It uses the structure dependent matrices `AVMAD` and `BVMAD` which
  !> are calculated once in the subroutine `MADELUNG3D()` and saved in
  !> the DA-file abvmad.unformatted (May 2004)
  !> The charge-moments are calculated in the subroutine `vintras`,
  !> therefore vintras has to be called first.
  !> The madelung-potential is expanded into spherical harmonics.
  !> The lm-term of the potential \(V\) of the atom \(i\) is given by
  !> \begin{equation}
  !> V(r,lm,i) = \sum_{i2}^{N} \sum_{l'm'} (-r)^l \left\{avmad(i,i2,lm,l'm')cmom(i2,l'm') +bvmad(i,i2,lm)z(i2)\right\}
  !> \end{equation}
  !> where \(N\) is the number of atoms
  !-------------------------------------------------------------------------------
  !> @note
  !> 
  !> - V. Popescu Feb. 2002: Adopted for the case of more atoms on the same site, 
  !> summation is done over the occupants of that site, the charge is weighted with 
  !> the appropriate concentration of the occupant
  !>
  !> - Impurity-program adopted feb. 2004 (according to N. Papanikalou)
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine vmadelblk(cmom,cminst,lmax,nspin,naez,v,zat,r,irws,ircut,ipan,kshape,  &
    noq,kaoez,conc,catom,icc,hostimp,vinters,nemb,lmpot,natyp)

    use :: mod_constants, only: pi, czero
    use :: global_variables, only: wlength, ipand, irmd, npotd
    use :: mod_main0, only: npol
    use :: mod_datatypes, only: dp
    use :: mod_types, only: t_madel

    implicit none

    ! .. Input variables
    integer, intent (in) :: icc    !! Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
    integer, intent (in) :: naez   !! Number of atoms in unit cell
    integer, intent (in) :: lmax   !! Maximum l component in wave function expansion
    integer, intent (in) :: nemb   !! Number of 'embedding' positions
    integer, intent (in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: nspin  !! Counter for spin directions
    integer, intent (in) :: lmpot  !! (LPOT+1)**2
    integer, intent (in) :: kshape !! Exact treatment of WS cell
    ! .. Array Arguments
    integer, dimension (naez), intent (in) :: noq !! Number of diff. atom types located
    integer, dimension (natyp), intent (in) :: irws !! Position of atoms in the unit cell in units of bravais vectors
    integer, dimension (natyp), intent (in) :: ipan !! Number of panels in non-MT-region
    integer, dimension (0:natyp), intent (in) :: hostimp
    integer, dimension (0:ipand, natyp), intent (in) :: ircut !! R points of panel borders
    integer, dimension (natyp, naez+nemb), intent (in) :: kaoez !! Kind of atom at site in elem. cell
    real (kind=dp), dimension (natyp), intent (in) :: zat !! Nuclear charge
    real (kind=dp), dimension (natyp), intent (in) :: conc !! Concentration of a given atom
    real (kind=dp), dimension (natyp), intent (in) :: catom
    real (kind=dp), dimension (irmd, natyp), intent (in) :: r !! Radial mesh ( in units a Bohr)
    real (kind=dp), dimension (lmpot, natyp), intent (in) :: cmom !! LM moment of total charge
    real (kind=dp), dimension (lmpot, natyp), intent (in) :: cminst !! charge moment of interstitial
    ! .. Input/Ouput variables
    real (kind=dp), dimension (irmd, lmpot, npotd), intent (inout) :: v
    ! .. Output variables
    real (kind=dp), dimension (lmpot, naez), intent (out) :: vinters
    ! .. Local Scalars
    integer :: lrecabmad, irec
    integer :: i, l, lm, lm2, lmmax, m, io1, io2, ipot, iq1, iq2
    integer :: irs1, ispin, it1, it2, noqval
    real (kind=dp) :: ac
    ! .. Local Arrays
    real (kind=dp), dimension (lmpot) :: bvmad !! Structure dependent matrix
    real (kind=dp), dimension (lmpot, lmpot) :: avmad !! Structure dependent matrix
    logical :: opt
    ! .. Intrinsic Functions ..
    intrinsic :: sqrt
    logical, external :: test
    ! ----------------------------------------------------------------------------
    write (1337, fmt=100)
    write (1337, fmt=110)

    lrecabmad = wlength*2*lmpot*lmpot + wlength*2*lmpot
    if (test('madelfil')) open (69, access='direct', recl=lrecabmad, file='abvmad.unformatted', form='unformatted')

    lmmax = (lmax+1)*(lmax+1)

    if (icc/=0) then
      do iq1 = 1, naez
        do lm = 1, lmpot
          vinters(lm, iq1) = 0e0_dp
        end do
      end do
    end if
    ! ----------------------------------------------------------------------------
    ! Loop over all types in unit cell
    ! ----------------------------------------------------------------------------
    do iq1 = 1, naez               ! added bauer 2/7/2012
      noqval = noq(iq1)            ! added bauer 2/7/2012
      if (noqval<1) noqval = 1     ! added bauer 2/7/2012
      do io1 = 1, noqval           ! added bauer 2/7/2012
        it1 = kaoez(io1, iq1)      ! added bauer 2/7/2012

        ! ----------------------------------------------------------------------
        ! Take a site occupied by atom IT1
        ! ----------------------------------------------------------------------
        if (it1/=-1) then          ! added bauer 2/7/2012
          if (kshape/=0) then
            irs1 = ircut(ipan(it1), it1)
          else
            irs1 = irws(it1)
          end if
        end if                     ! added bauer 2/7/2012
        ! ----------------------------------------------------------------------
        do l = 0, lmax
          ! -------------------------------------------------------------------
          do m = -l, l
            lm = l*l + l + m + 1
            ac = 0.0e0_dp
            ! ----------------------------------------------------------------
            if (naez==1) then
              irec = iq1 + naez*(iq1-1)
              if (npol==0) then
                avmad(:,:) = czero
                bvmad(:) = czero
              elseif (test('madelfil')) then
                read (69, rec=irec) avmad, bvmad
              else
                 avmad(:,:) = t_madel%avmad(irec,:,:)
                 bvmad(:) = t_madel%bvmad(irec,:)
              end if
              ! -------------------------------------------------------------
              ! Loop over all occupants of site IQ2=IQ1
              ! -------------------------------------------------------------
              do io2 = 1, noq(iq1)
                it2 = kaoez(io2, iq1)
                ! ----------------------------------------------------------
                ! lm = 1 component disappears if there is only one host atom
                ! take moments of sphere
                ! ----------------------------------------------------------
                do lm2 = 2, lmmax
                  ac = ac + avmad(lm, lm2)*cmom(lm2, it2)*conc(it2)
                end do
                ! ----------------------------------------------------------
                ! Add contribution of interstial in case of shapes
                ! ----------------------------------------------------------
                if (kshape/=0) then
                  do lm2 = 2, lmmax
                    ac = ac + avmad(lm, lm2)*cminst(lm2, it2)*conc(it2)
                  end do
                end if
              end do
              ! -------------------------------------------------------------
            else
              ! -------------------------------------------------------------
              ! Loop over all sites
              ! -------------------------------------------------------------
              do iq2 = 1, naez
                irec = iq2 + naez*(iq1-1)
                if (npol==0) then
                  avmad(:,:) = czero
                  bvmad(:) = czero
                elseif (test('madelfil')) then
                  read (69, rec=irec) avmad, bvmad
                else
                  avmad(:,:) = t_madel%avmad(irec,:,:)
                  bvmad(:) = t_madel%bvmad(irec,:)
                end if
                ! ----------------------------------------------------------
                ! Loop over all occupants of site IQ2
                ! ----------------------------------------------------------
                do io2 = 1, noq(iq2)

                  it2 = kaoez(io2, iq2)
                  ac = ac + bvmad(lm)*zat(it2)*conc(it2)
                  ! -------------------------------------------------------
                  ! Take moments of sphere
                  ! -------------------------------------------------------
                  do lm2 = 1, lmmax
                    ac = ac + avmad(lm, lm2)*cmom(lm2, it2)*conc(it2)
                  end do
                  ! -------------------------------------------------------
                  ! Add contribution of interstial in case of shapes
                  ! -------------------------------------------------------
                  if (kshape/=0) then
                    do lm2 = 1, lmmax
                      ac = ac + avmad(lm, lm2)*cminst(lm2, it2)*conc(it2)
                    end do
                  end if
                end do             ! IO2 = 1, NOQ(IQ2)
                ! ----------------------------------------------------------
              end do               ! IQ2 = 1, NAEZ
              ! -------------------------------------------------------------
            end if                 ! NAEZ.GT.1
            ! ----------------------------------------------------------------
            if (lm==1) then
              write (1337, fmt=120) it1, (catom(it1)-zat(it1)), (ac/sqrt(4.e0_dp*pi))
            end if
            ! ----------------------------------------------------------------
            ! Add to v the intercell-potential
            ! ----------------------------------------------------------------
            ! ----------------------------------------------------------------
            ! SPIN
            ! ----------------------------------------------------------------
            do ispin = 1, nspin
              ! -------------------------------------------------------------
              ! Determine the right potential number
              ! -------------------------------------------------------------
              ipot = nspin*(it1-1) + ispin
              ! -------------------------------------------------------------
              ! In the case of l=0 : r(1)**l is not defined
              ! -------------------------------------------------------------
              if (it1/=-1) then    ! added bauer 2/7/2012
                if (l==0) v(1, 1, ipot) = v(1, 1, ipot) + ac
                do i = 2, irs1
                  v(i, lm, ipot) = v(i, lm, ipot) + (-r(i,it1))**l*ac
                end do
              end if
            end do                 ! added bauer 2/7/2012
            ! ----------------------------------------------------------------
            ! SPIN
            ! ----------------------------------------------------------------
            if (icc/=0 .or. opt('KKRFLEX ')) then
              lm = l*l + l + m + 1
              write (1337, *) 'ac', iq1, lm, ac
              vinters(lm, iq1) = ac
            end if

          end do
          ! -------------------------------------------------------------------
        end do
        ! ----------------------------------------------------------------------
      end do
    end do
    ! ----------------------------------------------------------------------------
    if (test('madelfil')) close(69)
    ! ----------------------------------------------------------------------------
    write (1337, *) 'ICC in VMADELBLK', icc
    write (1337, '(25X,30("-"),/)')
    write (1337, '(79("="))')

    if ((icc==0) .and. (.not. opt('KKRFLEX '))) return
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now Prepare output for Impurity calculation
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open (91, file='intercell_ref', status='unknown', form='formatted')
    write (1337, *)
    write (1337, *) '                     ', 'Writing intercell potential for impurity'
    write (1337, '(/,20X,55("-"))')
    write (1337, 130) hostimp(0), lmmax
    write (1337, '(20X,55("-"),/,35X,"  i host lm  Vint")')
    do i = 1, hostimp(0)
      write (1337, *)
      lm = 1
      write (1337, '(35X,I4,I4,I3,1X,F10.6)') i, hostimp(i), lm, vinters(lm, hostimp(i))
      do lm = 2, 9
        write (1337, '(43X,I3,1X,F10.6)') lm, vinters(lm, hostimp(i))
      end do
      write (1337, '(20X,55("-"))')
    end do
    write (1337, '(79("="),/)')

    write (91, 140) hostimp(0), lmmax
    do i = 1, hostimp(0)
      write (91, 150)(vinters(lm,hostimp(i)), lm=1, lmmax)
    end do
    close (91)

    return

100 format (79('='), /, 18x, ' MADELUNG POTENTIALS ', '(spherically averaged) ')
110 format (/, 25x, ' ATOM ', '  Delta_Q  ', '     VMAD', /, 25x, 30('-'))
120 format (25x, i4, 2x, f10.6, 1x, f12.6)
130 format (22x, i4, ' host atoms, LMPOT = ', i2, ' output up to LM = 9')
140 format (3i6)
150 format (4d20.10)
  end subroutine vmadelblk

end module mod_vmadelblk
