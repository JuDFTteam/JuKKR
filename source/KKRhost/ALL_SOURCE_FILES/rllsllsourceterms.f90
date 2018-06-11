!-------------------------------------------------------------------------------
!> @brief Calculates the source terms J,H and the left solution J2, H2 for:
!> - non-relativistic
!> - scalar-relativistic
!> - full-relativistic
!> calculations
!-------------------------------------------------------------------------------
    Subroutine rllsllsourceterms(nsra, nvec, eryd, rmesh, nrmax, nrmaxd, lmax, &
      lmsize, use_fullgmat, jlk_index, hlk, jlk, hlk2, jlk2, gmatprefactor)

      Use constants
      Use mod_datatypes, Only: dp

      Implicit None

      Integer :: nsra, lmax, nrmax, nrmaxd, nvec
      Integer :: lmsize
      Complex (Kind=dp) :: eryd
      Real (Kind=dp), Dimension (nrmaxd) :: rmesh
      Integer, Dimension (2*lmsize) :: jlk_index
      Integer :: l1, lm1, m1, ivec, ispinfullgmat, ir
      Integer :: use_fullgmat
      Complex (Kind=dp) :: ek, ek2, gmatprefactor
      Complex (Kind=dp), Dimension (1:4*(lmax+1), nrmax) :: hlk, jlk
      Complex (Kind=dp), Dimension (1:4*(lmax+1), nrmax) :: hlk2, jlk2

      If (nsra==2) Then
        nvec = 2
      Else If (nsra==1) Then
        nvec = 1
      End If

      lm1 = 1
      Do ivec = 1, nvec
        Do ispinfullgmat = 0, use_fullgmat
          Do l1 = 0, lmax
            Do m1 = -l1, l1
              jlk_index(lm1) = l1 + (ivec-1)*(lmax+1) + 1
              lm1 = lm1 + 1
            End Do
          End Do
        End Do !ispinorbit=0,use_fullgmat
      End Do !nvec

      If (nsra==1) Then
        ek = sqrt(eryd)
        ek2 = sqrt(eryd)
      Else If (nsra==2) Then
        ek = sqrt(eryd+(eryd/cvlight)**2)
        ek2 = sqrt(eryd+(eryd/cvlight)**2)*(1.0E0_dp+eryd/cvlight**2)
      End If

      Do ir = 1, nrmax

        Call beshank(hlk(:,ir), jlk(:,ir), ek*rmesh(ir), lmax)
        If (nsra==2) Then
          Call beshank_smallcomp(hlk(:,ir), jlk(:,ir), ek*rmesh(ir), &
            rmesh(ir), eryd, lmax)
        End If

        Do l1 = 1, nvec*(lmax+1)
          hlk(l1, ir) = -ci*hlk(l1, ir)
        End Do

        If (nsra==1) Then
          Do l1 = 1, nvec*(lmax+1)
            jlk2(l1, ir) = jlk(l1, ir)
            hlk2(l1, ir) = hlk(l1, ir)
          End Do
        Else If (nsra==2) Then
          Do l1 = 1, lmax + 1
            jlk2(l1, ir) = jlk(l1, ir)
            hlk2(l1, ir) = hlk(l1, ir)
          End Do
          Do l1 = lmax + 2, 2*(lmax+1)
            jlk2(l1, ir) = -jlk(l1, ir)
            hlk2(l1, ir) = -hlk(l1, ir)
          End Do
        End If

      End Do
      gmatprefactor = ek2
    End Subroutine
