! 13.10.95 ***************************************************************
subroutine epotinb(epotin, nspin, natyp, rho2ns, vm2z, r, drdi, ins, irmin, &
  irws, lpot, vins, ircut, ipan, z)
! ************************************************************************

!     attention : energy zero ---> electro static zero

!                 since input potential and single particle energies
!                 are using muffin tin zero as zero the energy shift
!                 is cancelled in the kinetic energy contribution !


!     calculate the energy of the input potential
!     the energy for the representive atom i is given by

!                               rws
!       epotin(i) = - sqrt(4 pi) {  dr' vm2z(r',i)*rho2ns(r',1,i)
!                                0

!     in case of non spherical input potential one has to add

!                 rirt
!            {  -  {  dr' vins(r',lm,i)rho2ns(r',lm,i)   }
!                 rmin
!                                        (summed over lm)

!     remember : the non spherical part of the input potential is
!                different from zero only between r(irmin) and r(irt)

!             (see notes by b.drittler)

!     attention: vm2z is the spherically averaged input potential ,
!                vins contains the non spherical contribution of the
!                potential and rho2ns(...,1) is the  real charge density
!                times r**2. vins and rho2ns are expanded into spherical
!                harmonics . (see deck rholm or rhons)

!     remember :  in case of shape corrections  the contribution of
!                 the nuclear potential - 2*Z/r has to be explicitly
!                 taken into account between muffin tin sphere and
!                 circum scribed sphere .
!                 only within the muffin tin sphere this term is
!                 analytically cancelled wtih the contribution of
!                 the coulomb potential - see deck ecoulom


!                 modified for non spherical potential and shape correc-
!                  tions

!                               b.drittler   oct. 1989
!-----------------------------------------------------------------------
!.. Parameters ..
  include 'inc.p'
!.. Array Arguments ..
  integer :: lmpotd
  parameter (lmpotd=(lpotd+1)**2)
  integer :: irmind
  parameter (irmind=irmd-irnsd)
!..
!.. Local Scalars ..
  integer :: ins, lpot, natyp, nspin
!..
!.. Local Arrays ..
  double precision :: drdi(irmd, *), epotin(*), r(irmd, *), &
    rho2ns(irmd, lmpotd, natypd, *), vins(irmind:irmd, lmpotd, *), &
    vm2z(irmd, *), z(*)
  integer :: ipan(*), ircut(0:ipand, *), irmin(*), irws(*)
!..
!.. External Subroutines ..
  double precision :: pi, r2rhod, r2rhou, rfpi, temp, zzor
  integer :: i, iatyp, ic, ipan1, ipotd, ipotu, irc1, irmin1, irs1, l1, lm, m1
!..
!.. Intrinsic Functions ..
  double precision :: ens(0:lpotd, natypd), er(irmd)
  integer :: ircutm(0:ipand)
!     ..

  external :: simp3, simpk


  intrinsic :: atan, sqrt

  pi = 4.0d0*atan(1.0d0)
  rfpi = sqrt(4.0d0*pi)

  do iatyp = 1, natyp

    ipan1 = ipan(iatyp)
    irc1 = ircut(ipan1, iatyp)
!---> calculate charge density times input potential
    if (ipan1>1) then
      irs1 = ircut(1, iatyp)
    else
      irs1 = irws(iatyp)
    end if

    if (nspin==1) then
      ipotu = iatyp
      ipotd = iatyp
    else
      ipotu = 2*iatyp - 1
      ipotd = 2*iatyp
    end if

    do i = 1, irs1
!--->  remember the form of vm2z between mt sphere and rirc


      r2rhou = (rho2ns(i,1,iatyp,1)-rho2ns(i,1,iatyp,nspin))/2.0d0
      r2rhod = (rho2ns(i,1,iatyp,1)+rho2ns(i,1,iatyp,nspin))/2.0d0
      er(i) = -(r2rhou*vm2z(i,ipotu)+r2rhod*vm2z(i,ipotd))*rfpi
    end do
!--->   now integrate er to get epotin


    if (ipan1>1) then
      do i = irs1 + 1, irc1
        r2rhou = (rho2ns(i,1,iatyp,1)-rho2ns(i,1,iatyp,nspin))/2.0d0
        r2rhod = (rho2ns(i,1,iatyp,1)+rho2ns(i,1,iatyp,nspin))/2.0d0
        zzor = 2.0d0*z(iatyp)/r(i, iatyp)
        er(i) = -(r2rhou*(vm2z(i,ipotu)-zzor)+r2rhod*(vm2z(i, &
          ipotd)-zzor))*rfpi
      end do
    end if

!--->   add non spher. contribution in case of non spher. input potential

    if (ipan1>1) then
      call simpk(er, temp, ipan(iatyp), ircut(0,iatyp), drdi(1,iatyp))
    else
      call simp3(er, temp, 1, irs1, drdi(1,iatyp))
    end if

    epotin(iatyp) = temp
    ens(0, iatyp) = temp



    do l1 = 1, lpot
      ens(l1, iatyp) = 0.0d0
    end do

    if (ins/=0) then

      irmin1 = irmin(iatyp)
      if (irmin1<=irs1) then

        ircutm(0) = irmin1 - 1
        do ic = 1, ipan1
          ircutm(ic) = ircut(ic, iatyp)
        end do
!---> calculate charge density times potential
        do l1 = 1, lpot

          do i = 1, irmd
            er(i) = 0.0d0
          end do

          do m1 = -l1, l1
            lm = l1*(l1+1) + m1 + 1
            do i = irmin1, irc1

! (IRMIN1.LE.IRS1)

              r2rhou = (rho2ns(i,lm,iatyp,1)-rho2ns(i,lm,iatyp,nspin))/2.0d0
              r2rhod = (rho2ns(i,lm,iatyp,1)+rho2ns(i,lm,iatyp,nspin))/2.0d0
              er(i) = er(i) - r2rhou*vins(i, lm, ipotu) - &
                r2rhod*vins(i, lm, ipotd)
            end do
          end do
          call simpk(er, temp, ipan1, ircutm, drdi(1,iatyp))
! (INS.NE.0)
          epotin(iatyp) = epotin(iatyp) + temp
          ens(l1, iatyp) = temp
        end do

      end if 
! 13.10.95 ***************************************************************
    end if ! ************************************************************************

  end do
!     attention : energy zero ---> electro static zero
  return
end subroutine
