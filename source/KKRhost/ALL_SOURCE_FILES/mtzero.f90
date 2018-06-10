! ************************************************************************
subroutine mtzero(lmpot, natyp, conc, nspin, v, vbc, z, r, drdi, imt, ircut, &
  ipan, ntcell, lmsp, ifunm, thetas, irws, eshift, ishift, nshell, lsurf)
! ************************************************************************

!     determine muffin tin zero and shift potential to muffin tin zero

!     for spin polarized calculations muffin tin zero is related to
!         the average of the 2 spins

!                                            may,2000 (new version)

!-----------------------------------------------------------------------
  implicit none
!.. Parameters ..
  include 'inc.p'
  integer :: lmpotd
  parameter (lmpotd=(lpotd+1)**2)
!..
!.. Local Scalars ..
  double precision :: eshift, vbc(*)
  integer :: ishift, lmpot, natyp, nspin
!..
!.. Local Arrays ..
  double precision :: drdi(irmd, *), conc(natypd), r(irmd, *), &
    thetas(irid, nfund, *), v(irmd, lmpotd, *), z(*)
  integer :: ifunm(natypd, *), imt(*), ipan(*), ircut(0:ipand, *), irws(*), &
    lmsp(natypd, *), ntcell(*), nshell(0:nsheld)
  logical :: lsurf
!..
!.. External Subroutines ..
  double precision :: fpi, rfpi, vav0, vol0, zzor
  integer :: icell, ifun, ih, imt1, ipan1, ipot, ir, irc1, irh, is, lm
!..
!.. Intrinsic Functions ..
  double precision :: v1(irmd), v2(irmd), vav1(2), vol1(2)
!..

  logical :: test, opt
  external :: simp3, simpk, test


  intrinsic :: atan, sqrt

  fpi = 16.0d0*atan(1.0d0)
  rfpi = sqrt(fpi)
!---  >     muffin tin or atomic sphere calculation
  vav0 = 0.0d0
  vol0 = 0.0d0
  vav1(1) = 0.d0
  vav1(2) = 0.d0
  vol1(1) = 0.d0
  vol1(2) = 0.d0
  do ih = 1, natyp

    do ir = 1, irmd
      v1(ir) = 0.0d0
      v2(ir) = 0.0d0
    end do
    do is = 1, nspin
      ipot = nspin*(ih-1) + is
      ipan1 = ipan(ih)
      imt1 = imt(ih)

      if (ipan1==1) then

! (IPAN1.EQ.1)

        irc1 = irws(ih)
        do ir = imt1, irc1
          v2(ir) = fpi*r(ir, ih)**2
          zzor = 2.0d0*z(ih)/r(ir, ih)
          v1(ir) = (v(ir,1,ipot)/rfpi-zzor)*v2(ir)
        end do
!---  >     full potential calculation
        call simp3(v1, vav1(is), imt1, irc1, drdi(1,ih))
        call simp3(v2, vol1(is), imt1, irc1, drdi(1,ih))

      else 



        irc1 = ircut(ipan1, ih)
        icell = ntcell(ih)
        imt1 = imt(ih)
        do ir = imt1 + 1, irc1
          v2(ir) = r(ir, ih)**2*thetas(ir-imt1, 1, icell)*rfpi
          zzor = 2.0d0*z(ih)/r(ir, ih)
          v1(ir) = (v(ir,1,ipot)/rfpi-zzor)*v2(ir)
        end do
        do lm = 2, lmpot
          if (lmsp(icell,lm)>0) then
            ifun = ifunm(icell, lm)

            do ir = imt1 + 1, irc1
              irh = ir - imt1
              v1(ir) = v1(ir) + r(ir, ih)**2*v(ir, lm, ipot)*thetas(irh, ifun, &
                icell)
            end do
! (IPAN1.EQ.1)
          end if

        end do
! SPIN LOOP
        call simpk(v1, vav1(is), ipan1, ircut(0,ih), drdi(1,ih))
        call simpk(v2, vol1(is), ipan1, ircut(0,ih), drdi(1,ih))

      end if 
!     19.5.99   Nikos
    end do !     This way it is compatible with old kkr and tb-kkr
    if (nspin==1) then
      vav1(2) = vav1(1)
      vol1(2) = vol1(1)
    end if

! added 10.11.99 to fix vbc


    if (lsurf .and. (ih==1)) write (1337, *) 'Vacancies are ignored for VBC'
!---  > shift potential to muffin tin zero
    if (lsurf .and. (z(ih)<1.d0)) cycle
    vav0 = vav0 + conc(ih)*nshell(ih)*(vav1(1)+vav1(2))/2.d0
    vol0 = vol0 + conc(ih)*nshell(ih)*(vol1(1)+vol1(2))/2.d0
  end do
  if (.not. (opt('DECIMATE'))) then 
    vbc(1) = 0.0d0
    if (abs(vav0)>1d-10) vbc(1) = -vav0/vol0
    if (ishift>0) vbc(1) = vbc(1) + eshift
  end if

  write (1337, fmt=100) vol0, vav0, vbc(1)
  vbc(2) = vbc(1)


! ************************************************************************
  do is = 1, nspin
    do ih = 1, natyp
      ipot = nspin*(ih-1) + is
      do ir = 1, ircut(ipan(ih), ih)
        v(ir, 1, ipot) = v(ir, 1, ipot) + rfpi*vbc(is)
      end do
! ************************************************************************
    end do
  end do

  return
!     determine muffin tin zero and shift potential to muffin tin zero
100 format ('  VOL INT.', f16.9, '  VAV INT.', f16.9, '  VMT ZERO', f16.9)
110 format ('  ATOM ', i4, ' VMT ZERO :', f16.9)
end subroutine
