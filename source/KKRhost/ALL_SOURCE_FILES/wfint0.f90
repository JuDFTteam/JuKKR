    Subroutine wfint0(cder, dder, qzlm, qzekdr, pzekdr, vnspll, nsra, irmind, &
      irmd, lmmaxd, irmin, irmax) ! Added IRMIN,IRMAX 1.7.2014
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------------
!      determines the integrands CDER, DDER or ADER, BDER in the
!        integral equations for the non-spherical wavefunctions from
!        the non-spherical contributions of the potential vinsPLL.
!        (This subroutine is used in zeroth order Born approximation,
!         otherwise subroutine WFINT must be used)
!      R. Zeller      Aug. 1994
!-----------------------------------------------------------------------
      Implicit None
!.. Scalar Arguments ..
      Integer :: irmd, irmind, lmmaxd, nsra, irmin, irmax
!..
!.. Array Arguments ..
      Complex (Kind=dp) :: cder(lmmaxd, lmmaxd, irmind:irmd), &
        dder(lmmaxd, lmmaxd, irmind:irmd), pzekdr(lmmaxd, irmind:irmd, 2), &
        qzekdr(lmmaxd, irmind:irmd, 2), qzlm(lmmaxd, irmind:irmd, 2)
      Real (Kind=dp) :: vnspll(lmmaxd, lmmaxd, irmind:irmd)
!..
!.. Local Scalars ..
      Complex (Kind=dp) :: v1
      Integer :: ir, lm1, lm2

      Do ir = irmin, irmax
        Do lm2 = 1, lmmaxd
          Do lm1 = 1, lmmaxd
            v1 = vnspll(lm1, lm2, ir)*qzlm(lm2, ir, 1)
            cder(lm1, lm2, ir) = qzekdr(lm1, ir, 1)*v1
            dder(lm1, lm2, ir) = pzekdr(lm1, ir, 1)*v1
          End Do
        End Do
        If (nsra==2) Then
          Do lm2 = 1, lmmaxd
            Do lm1 = 1, lmmaxd
              v1 = vnspll(lm1, lm2, ir)*qzlm(lm2, ir, 2)
              cder(lm1, lm2, ir) = cder(lm1, lm2, ir) + qzekdr(lm1, ir, 2)*v1
              dder(lm1, lm2, ir) = dder(lm1, lm2, ir) + pzekdr(lm1, ir, 2)*v1
            End Do
          End Do
        End If

      End Do
    End Subroutine
