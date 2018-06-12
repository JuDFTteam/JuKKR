! 13.10.95 ***************************************************************
subroutine mixstr(rmsavq, rmsavm, ins, lpot, lmpot, natref, nshell, nstart, &
  nend, conc, nspin, itc, rfpi, fpi, ipf, mixing, fcm, irc, irmin, r, drdi, &
  vons, visp, vins, vspsmo, vspsme, lsmear)
  use :: mod_datatypes, only: dp
  ! ************************************************************************
  implicit none
  ! .. Parameters ..
  include 'inc.p'
  integer :: lmpotd, irmind
  parameter (lmpotd=(lpotd+1)**2, irmind=irmd-irnsd)
  ! ..
  ! .. Local Scalars ..
  real (kind=dp) :: fcm, fpi, mixing, rfpi, rmsavm, rmsavq
  integer :: ins, ipf, itc, lmpot, lpot, natref, nend, nspin, nstart
  integer :: lsmear
  ! ..
  ! .. Intrinsic Functions ..
  real (kind=dp) :: drdi(irmd, natypd), r(irmd, natypd), &
    vins(irmind:irmd, lmpotd, *), visp(irmd, *), vons(irmd, lmpotd, *), &
    conc(natypd), vspsmo(irmd, nspotd), vspsme(irmd, nspotd)
  integer :: irc(natypd), irmin(natypd), nshell(0:nsheld)
  ! ..

  real (kind=dp) :: fac, rmserm, rmserq, vmn, vnm, vnp, voldm, voldp, vpn
  real (kind=dp) :: natom
  integer :: i, ih, ihp1, irc1, irmin1, j, lm, np
  ! ---> final construction of the potentials
  ! attention : the spherical averaged potential is the lm=1
  intrinsic :: mod, real, sqrt
  ! component of vons times sqrt(4 pi).
  rmsavq = 0.0e0_dp
  rmsavm = 0.0e0_dp

  ! first mixing scheme : straight mixing
  ! ---> determination of the root mean sqare error





  natom = 0.0e0_dp

  do np = nstart, nend

    i = np - natref
    natom = natom + real(nshell(i), kind=dp)*conc(i)
    ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
    if (nspin==2) then
      ih = 2*np - 1
      ihp1 = ih + 1
    else
      ih = np
      ihp1 = ih
    end if

    irc1 = irc(np)
    rmserq = 0.0e0_dp
    rmserm = 0.0e0_dp
    fac = 0.5e0_dp/rfpi
    ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
    do j = 1, irc1
      vnp = fac*(vons(j,1,ih)+vons(j,1,ihp1))
      vnm = fac*(vons(j,1,ih)-vons(j,1,ihp1))
      voldp = 0.5e0_dp*(visp(j,ih)+visp(j,ihp1))
      voldm = 0.5e0_dp*(visp(j,ih)-visp(j,ihp1))
      rmserq = rmserq + 2.0e0_dp*real(1+mod(j,2), kind=dp)*r(j, np)*r(j, np)* &
        drdi(j, np)*(vnp-voldp)*(vnp-voldp)
      rmserm = rmserm + 2.0e0_dp*real(1+mod(j,2), kind=dp)*r(j, np)*r(j, np)* &
        drdi(j, np)*(vnm-voldm)*(vnm-voldm)
      vpn = voldp + mixing*(vnp-voldp)
      vmn = voldm + fcm*mixing*(vnm-voldm)
      vons(j, 1, ihp1) = vpn - vmn
      vons(j, 1, ih) = vpn + vmn
    end do


    if (lsmear>=3) then
      do j = 1, irc1
        vnp = 0.5e0_dp*(vspsmo(j,ih)+vspsmo(j,ihp1))
        vnm = 0.5e0_dp*(vspsmo(j,ih)-vspsmo(j,ihp1))
        voldp = 0.5e0_dp*(vspsme(j,ih)+vspsme(j,ihp1))
        voldm = 0.5e0_dp*(vspsme(j,ih)-vspsme(j,ihp1))
        vpn = voldp + mixing*(vnp-voldp)
        vmn = voldm + fcm*mixing*(vnm-voldm)
        vspsmo(j, ihp1) = vpn - vmn
        vspsmo(j, ih) = vpn + vmn
      end do
    end if

    if ((lsmear==1) .or. (lsmear==2)) then
      do j = 1, irc1
        vspsme(j, ihp1) = vspsmo(j, ihp1)
        vspsme(j, ih) = vspsmo(j, ih)
      end do
    end if


    rmserq = rmserq/(r(irc1,np)**3)
    rmserm = rmserm/(r(irc1,np)**3)
    rmsavq = rmsavq + rmserq*nshell(i)*conc(i)
    rmsavm = rmsavm + rmserm*nshell(i)*conc(i)

    if (nspin==2) then
      write (ipf, fmt=100) i, sqrt(rmserq), sqrt(rmserm)
    else
      write (ipf, fmt=120) i, sqrt(rmserq)
    end if

    if (ins/=0 .and. lpot>0) then

      rmserq = 0.0e0_dp
      rmserm = 0.0e0_dp
      irmin1 = irmin(np)
      do lm = 2, lmpot
        do j = irmin1, irc1
          vnp = 0.5e0_dp*(vons(j,lm,ih)+vons(j,lm,ihp1))
          vnm = 0.5e0_dp*(vons(j,lm,ih)-vons(j,lm,ihp1))
          voldp = 0.5e0_dp*(vins(j,lm,ih)+vins(j,lm,ihp1))
          voldm = 0.5e0_dp*(vins(j,lm,ih)-vins(j,lm,ihp1))
          rmserq = rmserq + 2.0e0_dp*real(1+mod(j,2), kind=dp)*r(j, np)*r(j, &
            np)*drdi(j, np)*(vnp-voldp)*(vnp-voldp)
          rmserm = rmserm + 2.0e0_dp*real(1+mod(j,2), kind=dp)*r(j, np)*r(j, &
            np)*drdi(j, np)*(vnm-voldm)*(vnm-voldm)
          vpn = voldp + mixing*(vnp-voldp)
          vmn = voldm + fcm*mixing*(vnm-voldm)
          vons(j, lm, ihp1) = vpn - vmn
          vons(j, lm, ih) = vpn + vmn
        end do
      end do
      rmserq = rmserq/(r(irc1,np)**3)/fpi
      rmserm = rmserm/(r(irc1,np)**3)/fpi
      rmsavq = rmsavq + rmserq*nshell(i)*conc(i)
      rmsavm = rmsavm + rmserm*nshell(i)*conc(i)

      if (nspin==2) then
        write (ipf, fmt=110) i, sqrt(rmserq), sqrt(rmserm)
      else
        write (ipf, fmt=130) i, sqrt(rmserq)
      end if

    end if

  end do
  ! 13.10.95 ***************************************************************
  ! ************************************************************************
  rmsavq = sqrt(rmsavq/natom)
  rmsavm = sqrt(rmsavm/natom)
  ! .. Parameters ..
  write (1337, '(79("-"),/)')
  if (nspin==2) then
    write (ipf, fmt=140) itc, rmsavq, rmsavm
    write (6, fmt=140) itc, rmsavq, rmsavm
  else
    write (ipf, fmt=150) itc, rmsavq
    write (6, fmt=150) itc, rmsavq
  end if
  write (1337, '(79("-"))')
  ! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
100 format (5x, ' rms-error for atom', i3, 1x, ':', 'v+ + v- = ', 1p, d11.4, &
    2x, ',', 2x, 'v+ - v- = ', 1p, d11.4)
110 format (5x, ' rms-error non spherical contribution for atom ', i3, 1x, &
    ':', 'v+ + v- = ', 1p, d11.4, 02x, ',', 2x, 'v+ - v- = ', 1p, d11.4)
120 format (5x, ' rms-error for atom', i3, 1x, ':', 'v+ + v- = ', 1p, d11.4)
130 format (5x, ' rms-error non spherical contribution for atom ', i3, 1x, &
    ':', 'v+ + v- = ', 1p, d11.4)
140 format ('      ITERATION', i4, ' average rms-error : v+ + v- = ', 1p, &
    d11.4, /, 39x, ' v+ - v- = ', 1p, d11.4)
150 format ('      ITERATION', i4, ' average rms-error : v+ + v- = ', 1p, &
    d11.4)
end subroutine mixstr
