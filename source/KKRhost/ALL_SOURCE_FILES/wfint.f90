    Subroutine wfint(qns, cder, dder, qzekdr, pzekdr, vnspll, nsra, irmind, &
      irmd, lmmaxd, irmin, irmax) ! Added IRMIN,IRMAX 1.7.2014
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------------
!      determines the integrands CDER, DDER or ADER, BDER in the
!        integral equations for the non-spherical wavefunctions from
!        the non-spherical contributions of the potential vinsPLL.

!      R. Zeller      Aug. 1994
!-----------------------------------------------------------------------
      Implicit None
!.. Scalar Arguments ..
      Integer :: irmd, irmind, lmmaxd, nsra, irmin, irmax
!..
!.. Array Arguments ..
      Complex (Kind=dp) :: cder(lmmaxd, lmmaxd, irmind:irmd), &
        dder(lmmaxd, lmmaxd, irmind:irmd), pzekdr(lmmaxd, irmind:irmd, 2), &
        qns(lmmaxd, lmmaxd, irmind:irmd, 2), qzekdr(lmmaxd, irmind:irmd, 2)
      Real (Kind=dp) :: vnspll(lmmaxd, lmmaxd, irmind:irmd)
!..
!.. Local Scalars ..
      Integer :: ir, lm1, lm2
!..
!.. External Subroutines ..
      External :: dgemm
!..
!.. Local Arrays ..
      Real (Kind=dp) :: qnsi(lmmaxd, lmmaxd), qnsr(lmmaxd, lmmaxd), &
        vtqnsi(lmmaxd, lmmaxd), vtqnsr(lmmaxd, lmmaxd)
!..
!.. Intrinsic Functions ..
      Intrinsic :: cmplx, aimag, real
!     ..

      Do ir = irmin, irmax
        Do lm2 = 1, lmmaxd
          Do lm1 = 1, lmmaxd
            qnsr(lm1, lm2) = real(qns(lm1,lm2,ir,1))
            qnsi(lm1, lm2) = aimag(qns(lm1,lm2,ir,1))
          End Do
        End Do
        Call dgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, 1.E0_dp, vnspll(1,1,ir), &
          lmmaxd, qnsr, lmmaxd, 0.E0_dp, vtqnsr, lmmaxd)
        Call dgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, 1.E0_dp, vnspll(1,1,ir), &
          lmmaxd, qnsi, lmmaxd, 0.E0_dp, vtqnsi, lmmaxd)
        Do lm1 = 1, lmmaxd
          Do lm2 = 1, lmmaxd
            cder(lm1, lm2, ir) = qzekdr(lm1, ir, 1)* &
              cmplx(vtqnsr(lm1,lm2), vtqnsi(lm1,lm2), kind=dp)
            dder(lm1, lm2, ir) = pzekdr(lm1, ir, 1)* &
              cmplx(vtqnsr(lm1,lm2), vtqnsi(lm1,lm2), kind=dp)
          End Do
        End Do
        If (nsra==2) Then
          Do lm2 = 1, lmmaxd
            Do lm1 = 1, lmmaxd
              qnsr(lm1, lm2) = real(qns(lm1,lm2,ir,2))
              qnsi(lm1, lm2) = aimag(qns(lm1,lm2,ir,2))
            End Do
          End Do
          Call dgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, 1.E0_dp, &
            vnspll(1,1,ir), lmmaxd, qnsr, lmmaxd, 0.E0_dp, vtqnsr, lmmaxd)
          Call dgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, 1.E0_dp, &
            vnspll(1,1,ir), lmmaxd, qnsi, lmmaxd, 0.E0_dp, vtqnsi, lmmaxd)
          Do lm2 = 1, lmmaxd
            Do lm1 = 1, lmmaxd
              cder(lm1, lm2, ir) = cder(lm1, lm2, ir) + &
                qzekdr(lm1, ir, 2)*cmplx(vtqnsr(lm1,lm2), vtqnsi(lm1,lm2), &
                kind=dp)
              dder(lm1, lm2, ir) = dder(lm1, lm2, ir) + &
                pzekdr(lm1, ir, 2)*cmplx(vtqnsr(lm1,lm2), vtqnsi(lm1,lm2), &
                kind=dp)
            End Do
          End Do
        End If

      End Do
    End Subroutine
