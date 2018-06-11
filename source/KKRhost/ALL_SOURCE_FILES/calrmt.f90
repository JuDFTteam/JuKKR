! ************************************************************************
    Subroutine calrmt(ipf, ipfe, ipe, imt, z, rmt, rws, rmtnew, alat, drdi, a, &
      b, irws, r, ifile, kshape)
      Use mod_datatypes, Only: dp
!***********************************************************************
!     this subroutine calculates imt and rmt(cal-rmt)
!                     and prints some informations about the used meshes
!        imtl = maximumnumber of meshpoints generating a radius
!               less or equal than rmt
!        imt  = number of meshpoint generating a new mt-radius closer th
!               mt-radius than every ather meshpoint
!***********************************************************************
!.. Scalar Arguments ..
      Real (Kind=dp) :: a, alat, b, rmt, rmtnew, rws, z
      Integer :: ifile, imt, ipe, ipf, ipfe, irws, kshape
!..
!.. Array Arguments ..
      Real (Kind=dp) :: drdi(*), r(*)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: drd1, drdws, rimt, rimtm1, rnuc
      Integer :: idelta, ih, imtl, irwsm2
!..
!.. Intrinsic Functions ..
      Intrinsic :: exp, log, mod, real
!..
!.. External Subroutines ..
      External :: rcstop
!     ..
      If (kshape==0) Then
        rimt = log(rmt/b+1.E0_dp)/a + 1.E0_dp
        imtl = rimt
        irwsm2 = irws - 2
        idelta = (rimt-imtl)*2
        If (idelta==0) imt = imtl
        If (idelta>0) imt = imtl + 1
        rimtm1 = real(imt-1, kind=dp)
        rmtnew = b*exp(a*rimtm1) - b

        If (imt>irwsm2) Then
          Write (ipf, Fmt=100)
          Call rcstop('calrmt  ')

        End If

      Else

        If (mod(imt,2)==0) Then
          Write (ipf, Fmt=*) ' error stop in calrmt - imt = ', imt, &
            ' has to be odd to get proper core charge  '
          Call rcstop('29      ')

        End If

      End If

      ih = irws/2
      drd1 = drdi(1)
      drdws = drdi(irws)
!----- nucleus radius rnuc in bohr's radii
      rnuc = 2.2677022E-5_dp*(2.E0_dp*z)**(1.0E0_dp/3.0E0_dp)
!-----
      If (ifile/=0) Then
        Write (ipf, Fmt=110) z, a, b, rnuc, r(2), ih, r(ih), drd1, drdws
        Write (ipf, Fmt=120) irws, imt, rws, rmt, rmtnew, alat
        If (ipe==1) Write (ipfe, Fmt=110) z, a, b, rnuc, r(2), ih, r(ih), &
          drd1, drdws
        If (ipe==1) Write (ipfe, Fmt=120) irws, imt, rws, rmt, rmtnew, alat
      End If



100   Format (1X, 'potentials need more meshpoints', /, 50('*'))
110   Format (' rmesh  z=', F5.2, '  a=', F7.4, '  b=', F9.6, '  rnuc=', &
        F11.8, '  r(2)=', F11.8, /, ' r(', I3, ')=', F7.4, '   drdi(1)=', &
        F11.8, '   drdi(irws)=', F9.6)
120   Format (' irws=', I6, ' imt=', I6, /, ' rws=', F12.8, ' rmt=', F12.8, &
        ' rmtnew=', F12.8, ' alat=', F12.8)
    End Subroutine
