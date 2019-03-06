!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_calrmt
  use :: mod_datatypes, only: dp
  private :: dp

contains


  !-------------------------------------------------------------------------------
  !> Summary: Calculates the muffin-tin radius
  !> Author: 
  !> Category: KKRhost, radial-grid
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This subroutine calculates imt and rmt(cal-rmt)
  !> and prints some informations about the used meshes
  !> imt  = number of meshpoint generating a new mt-radius closer to
  !> mt-radius than every other meshpoint
  !-------------------------------------------------------------------------------  
  subroutine calrmt(ipf, ipfe, ipe, imt, z, rmt, rws, rmtnew, alat, drdi, a, b, irws, r, ifile, kshape)
    ! ***********************************************************************
    ! 
    ! ***********************************************************************
    use :: mod_rcstop, only: rcstop
    implicit none
    ! .. Scalar Arguments ..
    real (kind=dp) :: a, alat, b, rmt, rmtnew, rws, z
    integer :: ifile, imt, ipe, ipf, ipfe, irws, kshape
    ! ..
    ! .. Array Arguments ..
    real (kind=dp) :: drdi(*), r(*)
    ! ..
    ! .. Local Scalars ..
    real (kind=dp) :: drd1, drdws, rimt, rimtm1, rnuc
    integer :: ih, irwsm2
    ! ..
    if (kshape==0) then
      rimt = log(rmt/b+1.e0_dp)/a + 1.e0_dp
      imt = nint(rimt)
      irwsm2 = irws - 2
      rimtm1 = real(imt-1, kind=dp)
      rmtnew = b*exp(a*rimtm1) - b

      if (imt>irwsm2) then
        write (ipf, fmt=100)
        call rcstop('calrmt  ')

      end if

    else

      if (mod(imt,2)==0) then
        write (ipf, fmt=*) ' error stop in calrmt - imt = ', imt, ' has to be odd to get proper core charge  '
        call rcstop('29      ')

      end if

    end if

    ih = irws/2
    drd1 = drdi(1)
    drdws = drdi(irws)
    ! ----- nucleus radius rnuc in bohr's radii
    rnuc = 2.2677022e-5_dp*(2.e0_dp*z)**(1.0e0_dp/3.0e0_dp)
    ! -----
    if (ifile/=0) then
      write (ipf, fmt=110) z, a, b, rnuc, r(2), ih, r(ih), drd1, drdws
      write (ipf, fmt=120) irws, imt, rws, rmt, rmtnew, alat
      if (ipe==1) write (ipfe, fmt=110) z, a, b, rnuc, r(2), ih, r(ih), drd1, drdws
      if (ipe==1) write (ipfe, fmt=120) irws, imt, rws, rmt, rmtnew, alat
    end if



100 format (1x, 'potentials need more meshpoints', /, 50('*'))
110 format (' rmesh  z=', f5.2, '  a=', f7.4, '  b=', f9.6, '  rnuc=', f11.8, '  r(2)=', f11.8, /, ' r(', i3, ')=', f7.4, '   drdi(1)=', f11.8, '   drdi(irws)=', f9.6)
120 format (' irws=', i6, ' imt=', i6, /, ' rws=', f12.8, ' rmt=', f12.8, ' rmtnew=', f12.8, ' alat=', f12.8)
  end subroutine calrmt

end module mod_calrmt
