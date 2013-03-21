!>    @param[out] CATOM  CATOM(1) charge, CATOM(2) magn. moment
!>    @param[in,out] RHO2NS is modified on output! - core charge added

subroutine rhototb_new(nspin,rho2ns,rhoc, &
drdi, &
ircut,lpot,nfu,llmsp,thetas,ipan, &
catom, &
irmd, irid, ipand, nfund)

  implicit none
  ! ************************************************************************
  !     add core and valence density expanded in spherical harmonics
  !         ( convention see subroutine rholm )
  !     in the paramagnetic case (nspin=1) the core valence charge times
  !         r**2 is add to the valence charge density times r**2
  !         then only rho2ns(irmd,lmxtsq,natypd,1) is used .
  !     in the spin-polarized case (nspin=2) the spin-splitted core
  !         charge density times r**2 is converted into core charge
  !         density times r**2 and core spin density times r**2 .
  !         then these parts are added to corresponding parts of
  !         the valence densities times r**2 , that are rho2ns(...,1)
  !         which contains the charge density  and rho2ns(...,2) which
  !         contains in that case the spin density .
  !             (see notes by b.drittler)
  !
  !     attention : the core density is spherically averaged and multi-
  !                 plied by 4 pi. therefore the core density is only
  !                 added to l=0 part .
  !
  !                               b.drittler   nov. 1989
  !
  !-----------------------------------------------------------------------
  !     .. Parameters ..

  integer irmd
  integer irid
  integer ipand
  integer nfund

  !     INTEGER LMPOTD
  !     PARAMETER (LMPOTD= (LPOTD+1)**2) ! = (2*LMAX+1)**2
  !     ..
  !     .. Scalar Arguments ..
  integer lpot,nspin

  double precision drdi(irmd)
  double precision rho2ns(irmd,(lpot+1)**2,2)
  double precision rhoc(irmd,2)
  double precision thetas(irid,nfund)
  double precision catom(nspin)

  integer ipan,ircut(0:ipand),llmsp(nfu),nfu

  !     ..
  !     .. Local Scalars ..
  double precision rfpi
  integer i,ifun,ipan1,irc1,irs1,ispin, &
  lm,lmpot
  !     ..
  !     .. Local Arrays ..
  double precision rho(irmd)
  !     ..
  !     .. External Subroutines ..
  external simpk
  !     ..
  !     .. Intrinsic Functions ..
  intrinsic atan,sqrt

  !     ..
  rfpi = sqrt(16.0d0*atan(1.0d0))
  lmpot = (lpot+1)**2

  ! ======================================================================
  !
  ipan1 = ipan
  irs1 = ircut(1)
  irc1 = ircut(ipan1)

  !-----------------------------------------------------------------------
  if(nspin.eq.2) then
    do i = 2,irs1
      rho2ns(i,1,1) = rho2ns(i,1,1) + (rhoc(i,1)+rhoc(i,2))/rfpi
      rho2ns(i,1,2) = rho2ns(i,1,2) + (rhoc(i,2)-rhoc(i,1))/rfpi
    end do
  else
    do i = 2,irs1
      rho2ns(i,1,1) = rho2ns(i,1,1) + rhoc(i,1)/rfpi
    end do
  end if
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  !
  !--->   calculate  charge and moment of the atom
  !
  do ispin = 1,nspin
    !
    !--->       convolute charge density with shape function to get the
    !           charge in the exact cell
    !
    do i = 1,irs1
      rho(i) = rho2ns(i,1,ispin)*rfpi
    end do
    
    do i = irs1 + 1,irc1
      rho(i) = 0.0d0
    end do
    
    do ifun = 1,nfu
      lm = llmsp(ifun)
      if (lm.le.lmpot) then
        do i = irs1 + 1,irc1
          rho(i) = rho(i) + rho2ns(i,lm,ispin)* &
          thetas(i-irs1,ifun)
        end do
      end if
    end do
    !
    !--->       integrate over circumscribed sphere
    !
    call simpk(rho,catom(ispin),ipan1, &
    ircut, drdi)

  end do                      ! ISPIN = 1,NSPIN
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end
