    Subroutine irwns(cr, dr, efac, qns, vnspll, icst, ipan, ircut, nsra, pzlm, &
      qzlm, pzekdr, qzekdr, cder, cmat, dder, dmat, irmind, irmd, irmin, &
      irmax, ipand, lmmaxd) ! Added IRMIN,IRMAX 1.7.2014
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------------
!     determines the irregular non spherical wavefunctions in the n-th.
!       born approximation ( n given by input parameter icst ) .


!     using the wave functions pz and qz ( regular and irregular
!       solution ) of the spherically averaged potential , the ir-
!       regular wavefunction qns is determined by

!           qns(ir,lm1,lm2) = cr(ir,lm1,lm2)*pz(ir,l1)

!                                   + dr(ir,lm1,lm2)*qz(ir,l1)

!      the matrices cr and dr are determined by integral equations
!        containing qns and only the non spherical contributions of
!        the potential , stored in vinspll . these integral equations
!        are solved iteratively with born approximation up to given n.

!     the original way of writing the cr and dr matrices in the equa-
!        tion above caused numerical troubles . therefore here are used
!        rescaled cr and dr matrices (compare subroutine wftsca):

!              ~
!              cr(ir,lm1,lm2) = sqrt(e)**(l1+l2)
!                             * cr(ir,lm1,lm2)/((2*l1-1)!!*(2*l2-1)!!)

!              ~
!              dr(ir,lm1,lm2) = sqrt(e)**(l2-l1)
!                             * dr(ir,lm1,lm2)*((2*l1-1)!!/(2*l2-1)!!)

!     attention :  the sign of the dr matrix is changed to reduce the
!     ===========  number of floating point operations

!     modified for the use of shape functions

!                              (see notes by b.drittler)

!                                b.drittler   mar.  1989
!-----------------------------------------------------------------------
!     modified by R. Zeller      Aug. 1994
!-----------------------------------------------------------------------
      Implicit None
!.. Parameters ..
      Complex (Kind=dp) :: cone
      Parameter (cone=(1.E0_dp,0.E0_dp))
!..
!.. Scalar Arguments ..
      Integer :: icst, ipan, ipand, irmd, irmind, lmmaxd, nsra, irmin, irmax
!..
!.. Array Arguments ..
      Complex (Kind=dp) :: cder(lmmaxd, lmmaxd, irmind:irmd), &
        cmat(lmmaxd, lmmaxd, irmind:irmd), cr(lmmaxd, lmmaxd), &
        dder(lmmaxd, lmmaxd, irmind:irmd), dmat(lmmaxd, lmmaxd, irmind:irmd), &
        dr(lmmaxd, lmmaxd), efac(lmmaxd), pzekdr(lmmaxd, irmind:irmd, 2), &
        pzlm(lmmaxd, irmind:irmd, 2), qns(lmmaxd, lmmaxd, irmind:irmd, 2), &
        qzekdr(lmmaxd, irmind:irmd, 2), qzlm(lmmaxd, irmind:irmd, 2)
      Real (Kind=dp) :: vnspll(lmmaxd, lmmaxd, irmind:irmd)
      Integer :: ircut(0:ipand)
!..
!.. Local Scalars ..
      Complex (Kind=dp) :: efac2
      Integer :: i, ir, irc1, j, lm1, lm2
!..
!.. External Subroutines ..
      External :: csinwd, wfint, wfint0
!..
      irc1 = ircut(ipan)
      Do i = 0, icst
!---> set up integrands for i-th born approximation
        If (i==0) Then
          Call wfint0(cder, dder, qzlm, qzekdr, pzekdr, vnspll, nsra, irmind, &
            irmd, lmmaxd, irmin, irmax) ! Added IRMIN,IRMAX 1.7.2014
        Else
          Call wfint(qns, cder, dder, qzekdr, pzekdr, vnspll, nsra, irmind, &
            irmd, lmmaxd, irmin, irmax) ! Added IRMIN,IRMAX 1.7.2014
        End If
!---> call integration subroutines
        Call csinwd(cder, cmat, lmmaxd**2, irmind, irmd, irmin, ipan, ircut) ! Added IRMIN 1.7.2014
        Call csinwd(dder, dmat, lmmaxd**2, irmind, irmd, irmin, ipan, ircut) ! Added IRMIN 1.7.2014
        Do ir = irmin, irc1
          Do lm2 = 1, lmmaxd
            dmat(lm2, lm2, ir) = dmat(lm2, lm2, ir) - cone
          End Do
        End Do
!---> calculate non sph. wft. in i-th born approximation
        Do j = 1, nsra
          Do ir = irmin, irc1
            Do lm1 = 1, lmmaxd
              Do lm2 = 1, lmmaxd
                qns(lm1, lm2, ir, j) = cmat(lm1, lm2, ir)*pzlm(lm1, ir, j) - &
                  dmat(lm1, lm2, ir)*qzlm(lm1, ir, j)
              End Do
            End Do
          End Do
        End Do
      End Do
      Do lm2 = 1, lmmaxd
!---> store c - and d - matrix
        Do lm1 = 1, lmmaxd
          cr(lm1, lm2) = cmat(lm1, lm2, irmin)
          dr(lm1, lm2) = -dmat(lm1, lm2, irmin)
        End Do
      End Do
!---> rescale with efac
      Do j = 1, nsra
        Do lm2 = 1, lmmaxd
          efac2 = 1.E0_dp/efac(lm2)
          Do ir = irmin, irc1
            Do lm1 = 1, lmmaxd
              qns(lm1, lm2, ir, j) = qns(lm1, lm2, ir, j)*efac2
            End Do
          End Do
        End Do
      End Do
    End Subroutine
