! ************************************************************************
subroutine calrmt(ipf, ipfe, ipe, imt, z, rmt, rws, rmtnew, alat, drdi, a, b, &
  irws, r, ifile, kshape)
!***********************************************************************
!     this subroutine calculates imt and rmt(cal-rmt)
!                     and prints some informations about the used meshes
!        imtl = maximumnumber of meshpoints generating a radius
!               less or equal than rmt
!        imt  = number of meshpoint generating a new mt-radius closer th
!               mt-radius than every ather meshpoint
!***********************************************************************
!.. Scalar Arguments ..
  double precision :: a, alat, b, rmt, rmtnew, rws, z
  integer :: ifile, imt, ipe, ipf, ipfe, irws, kshape
!..
!.. Array Arguments ..
  double precision :: drdi(*), r(*)
!..
!.. Local Scalars ..
  double precision :: drd1, drdws, rimt, rimtm1, rnuc
  integer :: idelta, ih, imtl, irwsm2
!..
!.. Intrinsic Functions ..
  intrinsic :: exp, log, mod, real
!..
!.. External Subroutines ..
  external :: rcstop
!     ..
  if (kshape==0) then
    rimt = log(rmt/b+1.d0)/a + 1.d0
    imtl = rimt
    irwsm2 = irws - 2
    idelta = (rimt-imtl)*2
    if (idelta==0) imt = imtl
    if (idelta>0) imt = imtl + 1
    rimtm1 = real(imt-1)
    rmtnew = b*exp(a*rimtm1) - b

    if (imt>irwsm2) then
      write (ipf, fmt=100)
      call rcstop('calrmt  ')

    end if

  else

    if (mod(imt,2)==0) then
      write (ipf, fmt=*) ' error stop in calrmt - imt = ', imt, &
        ' has to be odd to get proper core charge  '
      call rcstop('29      ')

    end if

  end if

  ih = irws/2
  drd1 = drdi(1)
  drdws = drdi(irws)
!----- nucleus radius rnuc in bohr's radii
  rnuc = 2.2677022d-5*(2.d0*z)**(1.0d0/3.0d0)
!-----
  if (ifile/=0) then
    write (ipf, fmt=110) z, a, b, rnuc, r(2), ih, r(ih), drd1, drdws
    write (ipf, fmt=120) irws, imt, rws, rmt, rmtnew, alat
    if (ipe==1) write (ipfe, fmt=110) z, a, b, rnuc, r(2), ih, r(ih), drd1, &
      drdws
    if (ipe==1) write (ipfe, fmt=120) irws, imt, rws, rmt, rmtnew, alat
  end if



100 format (1x, 'potentials need more meshpoints', /, 50('*'))
110 format (' rmesh  z=', f5.2, '  a=', f7.4, '  b=', f9.6, '  rnuc=', f11.8, &
    '  r(2)=', f11.8, /, ' r(', i3, ')=', f7.4, '   drdi(1)=', f11.8, &
    '   drdi(irws)=', f9.6)
120 format (' irws=', i6, ' imt=', i6, /, ' rws=', f12.8, ' rmt=', f12.8, &
    ' rmtnew=', f12.8, ' alat=', f12.8)
end subroutine
