!     initialise VAV0 and VOL0 to zero before calling!!!

! ************************************************************************
subroutine mtzero_new(lmpot,nspin,vons,z,r,drdi,imt1,ircut, &
ipan1,lmsp,ifunm,thetas,irws,vav0,vol0, &
irmd, irid, nfund, ipand)

  implicit none
  ! ************************************************************************
  !
  !     determine muffin tin zero - average of potential in interstitial
  !
  !     for spin polarized calculations muffin tin zero is related to
  !         the average of the 2 spins
  !
  !                                            may,2000 (new version)
  !
  !-----------------------------------------------------------------------

  integer irmd
  integer irid
  integer nfund
  integer ipand

  !     .. Scalar Arguments ..
  integer lmpot,nspin
  !     ..
  !     .. Array Arguments ..
  double precision &
  drdi(irmd), &
  r(irmd), &
  thetas(irid,nfund), &
  vons(irmd,lmpot,2), &
  z
  integer &
  ifunm(*),ircut(0:ipand), &
  lmsp(*), irws
  !     ..
  !     .. Local Scalars ..
  double precision fpi,rfpi,vav0,vol0,zzor
  integer ifun,imt1,ipan1,ir,irc1,irh,is,lm
  !     ..
  !     .. Local Arrays ..
  double precision v1(irmd),v2(irmd),vav1(2),vol1(2)
  !     ..
  !     .. External Subroutines ..
  logical test
  external simp3,simpk,test
  !     ..
  !     .. Intrinsic Functions ..
  intrinsic atan,sqrt
  !     ..
  fpi = 16.0d0*atan(1.0d0)
  rfpi = sqrt(fpi)

  do 10 ir = 1,irmd
    v1(ir) = 0.0d0
    v2(ir) = 0.0d0
10 continue


   do is = 1,nspin
     
     if (ipan1.eq.1) then
       !
       !---  >     muffin tin or atomic sphere calculation
       !
       irc1 = irws
       do 20 ir = imt1,irc1
         v2(ir) = fpi*r(ir)**2
         zzor = 2.0d0* z / r(ir)
         v1(ir) = (vons(ir,1,is)/rfpi-zzor)*v2(ir)
20     continue
        
       call simp3(v1,vav1(is),imt1,irc1,drdi)
       call simp3(v2,vol1(is),imt1,irc1,drdi)
        
     else                 ! (IPAN1.EQ.1)
       !
       !---  >     full potential calculation
       !
       irc1 = ircut(ipan1)

       do 30 ir = imt1 + 1,irc1
         v2(ir) = r(ir)**2*thetas(ir-imt1, 1) * rfpi
         zzor = 2.0d0 * z / r(ir)
         v1(ir) = (vons(ir,1,is)/rfpi-zzor)*v2(ir)
30     continue

       do 50 lm = 2,lmpot

         if (lmsp(lm).gt.0) then
           ifun = ifunm(lm)
              
           do 40 ir = imt1 + 1,irc1
             irh = ir - imt1
             v1(ir) = v1(ir) + r(ir)**2*vons(ir,lm,is)* &
             thetas(irh,ifun)
40         continue
              
         end if
           
50     continue          ! LM = 2,LMPOT
        
       call simpk(v1,vav1(is),ipan1,ircut(0),drdi)
       call simpk(v2,vol1(is),ipan1,ircut(0),drdi)
        
     end if               ! (IPAN1.EQ.1)
   !
   !
   !
   end do                  ! SPIN LOOP
   !
   !
   !
   if (nspin.eq.1) then
     vav1(2) = vav1(1)
     vol1(2) = vol1(1)
   end if
  
   vav0 = vav0 + (vav1(1)+vav1(2))/2.d0
   vol0 = vol0 + (vol1(1)+vol1(2))/2.d0

   return

 end
