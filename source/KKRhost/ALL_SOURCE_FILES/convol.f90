! ************************************************************************
subroutine convol(imt1, irc1, icell, imaxsh, ilm_map, ifunm, lmpot, gsh, &
  thetas, thesme, z, rfpi, r, vons, vspsmo, lmsp)
! ************************************************************************
!.. Parameters ..
  include 'inc.p'
  integer :: lmpotd
  parameter (lmpotd=(lpotd+1)**2)
!..
!.. Local Scalars ..
  double precision :: rfpi, z
  integer :: icell, imaxsh, imt1, irc1, lmpot
!..
!.. Local Arrays ..
  double precision :: gsh(*), r(*), thetas(irid, nfund, *), vons(irmd, *), &
    thesme(irid, nfund, *), vspsmo(irmd)
  integer :: ifunm(natypd, *), ilm_map(ngshd, 3), lmsp(natypd, *)
!..

  double precision :: zzor
  integer :: i, ifun, ir, irh, lm, lm1, lm2, lm3


  double precision :: vstore(irid, lmpotd), vstsme(irid, lmpotd)


  do lm = 1, lmpot
    do ir = 1, irc1 - imt1
      vstore(ir, lm) = 0.0d0
      vstsme(ir, lm) = 0.0d0
    end do
  end do
!     COPY THE PART INSIDE THE MT SPHERE
  do ir = imt1 + 1, irc1
    zzor = 2.0d0*z/r(ir)*rfpi
    vons(ir, 1) = vons(ir, 1) - zzor
  end do

  do i = 1, imaxsh
    lm1 = ilm_map(i, 1)
    lm2 = ilm_map(i, 2)
    lm3 = ilm_map(i, 3)
    if (lmsp(icell,lm3)>0) then
      ifun = ifunm(icell, lm3)
      do ir = imt1 + 1, irc1
        irh = ir - imt1
        vstore(irh, lm1) = vstore(irh, lm1) + gsh(i)*vons(ir, lm2)*thetas(irh, &
          ifun, icell)
        vstsme(irh, lm1) = vstsme(irh, lm1) + gsh(i)*vons(ir, lm2)*thesme(irh, &
          ifun, icell)
      end do
    end if
  end do

  do ir = imt1 + 1, irc1
    irh = ir - imt1
    zzor = 2.0d0*z/r(ir)*rfpi
    vons(ir, 1) = vstore(irh, 1) + zzor
    vspsmo(ir) = (vstsme(irh,1)+zzor)/rfpi
  end do

! ************************************************************************
  do ir = 1, imt1
    vspsmo(ir) = vons(ir, 1)/rfpi
  end do
! ************************************************************************
  do lm = 2, lmpot
    do ir = imt1 + 1, irc1
      irh = ir - imt1
      vons(ir, lm) = vstore(irh, lm)
    end do
  end do
!.. Parameters ..
  return
! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
end subroutine
