! Added IRMIN,IRMAX 1.7.2014  &
    Subroutine wftsca(drdi, efac, pz, qz, fz, sz, nsra, pzlm, qzlm, pzekdr, &
      qzekdr, ek, loflm, irmind, irmd, irmin, irmax, lmaxd, lmmaxd)
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------------
!                 R. Zeller      Oct. 1993
!-----------------------------------------------------------------------
      Implicit None
!.. Parameters ..
      Complex (Kind=dp) :: cone
      Parameter (cone=(1.E0_dp,0.E0_dp))
!..
!.. Scalar Arguments ..
      Complex (Kind=dp) :: ek
      Integer :: irmd, irmind, lmaxd, lmmaxd, nsra, irmin, irmax
!..
!.. Array Arguments ..
      Complex (Kind=dp) :: efac(lmmaxd), fz(irmd, 0:lmaxd), pz(irmd, 0:lmaxd), &
        pzekdr(lmmaxd, irmind:irmd, 2), pzlm(lmmaxd, irmind:irmd, 2), &
        qz(irmd, 0:lmaxd), qzekdr(lmmaxd, irmind:irmd, 2), &
        qzlm(lmmaxd, irmind:irmd, 2), sz(irmd, 0:lmaxd)
      Real (Kind=dp) :: drdi(*)
      Integer :: loflm(*)
!..
!.. Local Scalars ..
      Complex (Kind=dp) :: efac1, v1
      Integer :: ir, j, l, l1, lm, lm1, m
!..
!.. Intrinsic Functions ..
      Intrinsic :: real
!..


!---> set up array efac : efac(lm) = sqrt(e)**l/(2l - 1)!!

      efac(1) = cone
      v1 = cone
      Do l = 1, lmaxd
        v1 = v1*ek/real(2*l-1, kind=dp)
        Do m = -l, l
          lm = l*(l+1) + m + 1
          efac(lm) = v1
        End Do
      End Do


!---> get wfts of same magnitude by scaling with efac

      Do lm1 = 1, lmmaxd
        l1 = loflm(lm1)
        efac1 = efac(lm1)
        Do ir = irmin, irmax
          pzlm(lm1, ir, 1) = pz(ir, l1)/efac1
          qzlm(lm1, ir, 1) = qz(ir, l1)*efac1
        End Do
        If (nsra==2) Then
          Do ir = irmin, irmax
            pzlm(lm1, ir, nsra) = fz(ir, l1)/efac1
            qzlm(lm1, ir, nsra) = sz(ir, l1)*efac1
          End Do
        End If

        Do j = 1, nsra
          Do ir = irmin, irmax
            pzekdr(lm1, ir, j) = pzlm(lm1, ir, j)*ek*drdi(ir)
            qzekdr(lm1, ir, j) = qzlm(lm1, ir, j)*ek*drdi(ir)
          End Do
        End Do
      End Do


    End Subroutine
