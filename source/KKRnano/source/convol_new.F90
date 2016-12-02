!>    Convolute the potential with the shapefunctions.
!>    Operates on one spin channel and one site.
!>    @param GSH shape-Gaunt-coefficients
!>    @param Z    nuclear charge
!>    @param R    radial mesh
!>    @param VONS non-spherical part of potential (only for one spin
!>                channel)
subroutine convol_new(imt1,irc1,imaxsh,ilm,ifunm,lmpot,gsh, thetas,z,r,vons,lmsp, irid, nfund, irmd, ngshd)
  implicit none
  integer, intent(in) :: irid !< number of radial mesh points for shape functions
  integer, intent(in) :: nfund !< maximal number of shape functions
  integer, intent(in) :: irmd
  integer, intent(in) :: ngshd
  integer, intent(in) :: imaxsh,imt1,irc1,lmpot
  double precision, intent(in) :: gsh(*), r(*), thetas(irid,nfund)
  double precision, intent(inout) :: vons(irmd,*)
  integer, intent(in) :: ifunm(*), ilm(ngshd,3), lmsp(*)

  double precision :: zzor
  integer :: i,ifun,ir,irh,lm,lm1,lm2,lm3
  double precision :: vstore(irid,lmpot)
  double precision :: rfpi, z

    rfpi = sqrt(16.d0 * atan(1.d0))

    do lm = 1, lmpot
      do ir = 1, irc1 - imt1
        vstore(ir,lm) = 0.d0
      enddo ! ir
    enddo ! lm
    
    do ir = imt1 + 1, irc1
      zzor = 2.d0*z/r(ir)*rfpi
      vons(ir,1) = vons(ir,1) - zzor
    enddo ! ir

    do i = 1, imaxsh
      lm1 = ilm(i,1)
      lm2 = ilm(i,2)
      lm3 = ilm(i,3)
      if (lmsp(lm3) > 0) then
        ifun = ifunm(lm3)
        do ir = imt1 + 1, irc1
          irh = ir - imt1
          vstore(irh,lm1) = vstore(irh,lm1) + gsh(i)*vons(ir,lm2)*thetas(irh,ifun)
        enddo ! ir
      endif
    enddo ! i

    do ir = imt1 + 1, irc1
      irh = ir - imt1
      zzor = 2.d0*z/r(ir)*rfpi
      vons(ir,1) = vstore(irh,1) + zzor
    enddo ! ir

    do lm = 2,lmpot
      do ir = imt1 + 1, irc1
        irh = ir - imt1
        vons(ir,lm) = vstore(irh,lm)
      enddo ! ir
    enddo ! lm

endsubroutine ! convol_new
