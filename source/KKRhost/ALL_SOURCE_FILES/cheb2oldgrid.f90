    Subroutine cheb2oldgrid(nrmax, nrmaxnew, lmmaxpot, rmesh, ncheb, npan_tot, &
      rpan_intervall, ipan_intervall, arrayin, arrayout, irmd)
      Use mod_datatypes, Only: dp

!use mod_cheb, only: getCCmatrix, getCinvmatrix
      Implicit None

! Interpolate from Chebyshev mesh to the old radial mesh
! Programmed by David Bauer, 2011-2013
      Integer :: ncheb, npan_tot, nrmax, nrmaxnew, lmmaxpot, irmd
      Real (Kind=dp) :: rmesh(irmd)
      Real (Kind=dp) :: rpan_intervall(0:npan_tot)
      Integer :: ipan_intervall(0:npan_tot)
      Complex (Kind=dp) :: arrayin(nrmaxnew, lmmaxpot)
      Complex (Kind=dp) :: arrayout(nrmax, lmmaxpot)
!local
      Integer :: in, ir, ir2, ir3, it, ilm
      Integer :: intsub(2, nrmax)
      Real (Kind=dp) :: cinvmatrix(0:ncheb, 0:ncheb)
      Real (Kind=dp), Allocatable :: ccmatrix(:, :)
      Complex (Kind=dp) :: alphaparams(0:ncheb, lmmaxpot)
      Real (Kind=dp) :: rmeshnorm(nrmax)
      Real (Kind=dp) :: tol
      Parameter (tol=1.E-09_dp)

! divide the mesh into subintervals
      intsub(1, :) = 0
      intsub(2, :) = -1
      in = 1

      ir = 1
      it = 1
      Do While (ir<=nrmax .And. in<=npan_tot)
        If (abs(rmesh(ir)-rpan_intervall(in))<tol) Then
          intsub(it, in) = ir
          intsub(2, in) = ir
          it = 1
          ir = ir + 1
          in = in + 1
        Else If (((rmesh(ir)-rpan_intervall(in)<-tol)) .And. ((rmesh(ir)- &
            rpan_intervall(in-1))>-tol)) Then
          intsub(it, in) = ir
          intsub(2, in) = ir
          it = 2
          ir = ir + 1
        Else If ((rmesh(ir)-rpan_intervall(in-1))<-tol) Then
          intsub(1, in) = 0
          intsub(2, in) = -1
          it = 1
          ir = ir + 1
        Else
          it = 1
          in = in + 1
        End If
      End Do
      If (abs(rmesh(nrmax)-rpan_intervall(npan_tot))>tol) Then
        Write (*, *) rmesh(nrmax), rpan_intervall(npan_tot), &
          rmesh(nrmax) - rpan_intervall(npan_tot)
        Stop 'error with the new and old mesh'
      End If


      Call getcinvmatrix(ncheb, cinvmatrix)

      in = 0
      Do While (in<=npan_tot)
        in = in + 1
        If (intsub(2,in)<intsub(1,in)) Cycle
        alphaparams = 0.0E0_dp
        Do ilm = 1, lmmaxpot
          Do ir2 = 0, ncheb
            Do ir = 0, ncheb
              alphaparams(ir2, ilm) = alphaparams(ir2, ilm) + &
                cinvmatrix(ir2, ir)*arrayin(ipan_intervall(in-1)+1+ir, ilm)
            End Do
          End Do
        End Do

! Transform to normalized coordinates between [-1,1]
! Shift upper & lower end to be centered around zero
        Do ir = intsub(1, in), intsub(2, in)
          ir2 = ir + 1 - intsub(1, in)
          rmeshnorm(ir2) = (2E0_dp*rmesh(ir)-(rpan_intervall( &
            in)+rpan_intervall(in-1)))/(rpan_intervall(in)-rpan_intervall(in-1 &
            ))
          If (abs(rmeshnorm(ir2))>1.0E0_dp) Then
            If (abs(abs(rmeshnorm(ir2))-1.0E0_dp)<tol) Then
              rmeshnorm(ir2) = sign(1.0E0_dp, rmeshnorm(ir2))
            Else
              Write (*, *) 'ir, rmeshnorm(ir2)', ir, rmeshnorm(ir2)
              Write (*, *) 'rmeshnorm not in interval [-1,1] in cheb2oldgrid'
              Write (*, *) 'Check radial mesh (panels) or increase tolerance.'
              Stop 'rmeshnorm not in interval [-1,1] in cheb2oldgrid'
            End If
          End If
        End Do

        Allocate (ccmatrix(intsub(2,in)-intsub(1,in)+1,0:ncheb))

        Call getccmatrix(ncheb, rmeshnorm(1:intsub(2,in)-intsub(1, &
          in)+1), intsub(2,in)-intsub(1,in)+1, ccmatrix)
        Do ilm = 1, lmmaxpot
          Do ir2 = intsub(1, in), intsub(2, in)
            Do ir = 0, ncheb
              ir3 = ir2 - intsub(1, in) + 1
              arrayout(ir2, ilm) = arrayout(ir2, ilm) + &
                ccmatrix(ir3, ir)*alphaparams(ir, ilm)
            End Do
          End Do
        End Do !ilm=1,lmmaxpot
        Deallocate (ccmatrix) !(CCmatrix(intsub(2,in)-intsub(1,in)+1,0:ncheb))

      End Do !in

    End Subroutine
