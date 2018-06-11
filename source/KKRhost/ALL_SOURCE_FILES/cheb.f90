!module mod_cheb

!contains 

    Subroutine getcmatrix(ncheb, cmatrix)
! calculates the C matrix according to:
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
      Use mod_datatypes
      Implicit None
      Integer, Intent (In) :: ncheb
      Real (Kind=dp), Intent (Out) :: cmatrix(0:ncheb, 0:ncheb)
      Real (Kind=dp) :: pi
!local
      Integer :: icheb1, icheb2

      pi = 4E0_dp*atan(1E0_dp)
      Do icheb1 = 0, ncheb
        Do icheb2 = 0, ncheb
! maybe incorrect
          cmatrix(icheb2, icheb1) = cos(icheb1*pi*((ncheb- &
            icheb2)+0.5E0_dp)/(ncheb+1))
        End Do
      End Do
    End Subroutine


    Subroutine getcinvmatrix(ncheb, cinvmatrix)
! calculates the C**-1 matrix according to:
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
      Use mod_datatypes
      Implicit None
      Integer, Intent (In) :: ncheb
      Real (Kind=dp), Intent (Out) :: cinvmatrix(0:ncheb, 0:ncheb)
!local
      Real (Kind=dp) :: pi
      Integer :: icheb1, icheb2
      Real (Kind=dp) :: fac

      pi = 4E0_dp*atan(1E0_dp)
      fac = 1.0E0_dp/(ncheb+1)
      Do icheb1 = 0, ncheb
        Do icheb2 = 0, ncheb
          cinvmatrix(icheb1, icheb2) = fac*cos(icheb1*pi*((ncheb- &
            icheb2)+0.5E0_dp)/(ncheb+1))
        End Do
        fac = 2.0E0_dp/(ncheb+1)
      End Do
    End Subroutine


    Subroutine getccmatrix(ncheb, rmesh, nrmesh, cmatrix)
! calculates the C matrix according to:
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
      Use mod_datatypes
      Implicit None
      Integer, Intent (In) :: ncheb, nrmesh
      Real (Kind=dp), Intent (In) :: rmesh(nrmesh)
      Real (Kind=dp), Intent (Out) :: cmatrix(1:nrmesh, 0:ncheb)
      Integer :: icheb, ir

      Do ir = 1, nrmesh
        Do icheb = 0, ncheb
          cmatrix(ir, icheb) = cos(real(icheb,kind=dp)*acos(rmesh(ir)))
        End Do
      End Do
    End Subroutine


    Subroutine getlambda(ncheb, lambda)
! set up the Lambda matrix which differentiates the coefficients of an
! Chebyshev expansion 
      Use mod_datatypes
      Implicit None
      Integer, Intent (In) :: ncheb
      Real (Kind=dp), Intent (Out) :: lambda(0:ncheb, 0:ncheb)
!local
      Integer :: icheb, icheb2

      Do icheb2 = 1, ncheb, 2
        lambda(0, icheb2) = icheb2
      End Do
      Do icheb = 1, ncheb
        Do icheb2 = icheb + 1, ncheb, 2
          lambda(icheb, icheb2) = icheb2*2
        End Do
      End Do
    End Subroutine


    Subroutine getclambdacinv(ncheb, clambdacinv)
      Use mod_datatypes
      Implicit None
! set up the Lambda matrix which differentiates the coefficients of an
! Chebyshev expansion
      Integer :: ncheb
      Real (Kind=dp) :: clambdacinv(0:ncheb, 0:ncheb)
!local
      Real (Kind=dp) :: lambda(0:ncheb, 0:ncheb)
      Real (Kind=dp) :: cmatrix(0:ncheb, 0:ncheb)
      Real (Kind=dp) :: cinvmatrix(0:ncheb, 0:ncheb)
      Real (Kind=dp) :: temp1(0:ncheb, 0:ncheb)
      Integer :: n

      lambda = (0.0E0_dp, 0.0E0_dp)
      cmatrix = (0.0E0_dp, 0.0E0_dp)
      cinvmatrix = (0.0E0_dp, 0.0E0_dp)
      lambda = (0.0E0_dp, 0.0E0_dp)
      temp1 = (0.0E0_dp, 0.0E0_dp)

      Call getlambda(ncheb, lambda)
      Call getcinvmatrix(ncheb, cinvmatrix)
      Call getcmatrix(ncheb, cmatrix)
      n = ncheb + 1
      Call dgemm('N', 'N', n, n, n, 1E0_dp, lambda, n, cinvmatrix, n, 0E0_dp, &
        temp1, n)
      Call dgemm('N', 'N', n, n, n, 1E0_dp, cmatrix, n, temp1, n, 0E0_dp, &
        clambdacinv, n)
    End Subroutine


    Subroutine getclambda2cinv(ncheb, clambda2cinv)
      Use mod_datatypes
      Implicit None
! set up the Lambda matrix which differentiates the coefficients of an
! Chebyshev expansion
      Integer :: ncheb
      Real (Kind=dp) :: clambda2cinv(0:ncheb, 0:ncheb)
!local
      Real (Kind=dp) :: lambda(0:ncheb, 0:ncheb)
      Real (Kind=dp) :: cmatrix(0:ncheb, 0:ncheb)
      Real (Kind=dp) :: cinvmatrix(0:ncheb, 0:ncheb)
      Real (Kind=dp) :: temp1(0:ncheb, 0:ncheb)
      Real (Kind=dp) :: temp2(0:ncheb, 0:ncheb)

      lambda = (0.0E0_dp, 0.0E0_dp)
      cmatrix = (0.0E0_dp, 0.0E0_dp)
      cinvmatrix = (0.0E0_dp, 0.0E0_dp)
      lambda = (0.0E0_dp, 0.0E0_dp)
      temp1 = (0.0E0_dp, 0.0E0_dp)

      Call getlambda(ncheb, lambda)
      Call getcinvmatrix(ncheb, cinvmatrix)
      Call getcmatrix(ncheb, cmatrix)

      Call matmat_dmdm(lambda, lambda, ncheb, temp1)
      Call matmat_dmdm(temp1, cinvmatrix, ncheb, temp2)
      Call matmat_dmdm(cmatrix, temp2, ncheb, clambda2cinv)
    End Subroutine


    Subroutine diffcheb(fn, ncheb, dfndr)
      Use mod_datatypes
      Implicit None
      Integer :: ncheb
      Real (Kind=dp) :: fn(0:ncheb)
      Real (Kind=dp) :: dfndr(0:ncheb)
      Real (Kind=dp) :: clambdacinv(0:ncheb, 0:ncheb)

!needs to be checked!!!!!!1
      Call getclambdacinv(ncheb, clambdacinv(0:ncheb,0:ncheb))
      Call matvec_dmdm(ncheb, clambdacinv(0:ncheb,0:ncheb), fn(0:ncheb), &
        dfndr(0:ncheb))
    End Subroutine


! helper functions

    Subroutine matvec_dmdm(ncheb, mat1, vec1, outvec)
      Use mod_datatypes
      Implicit None
      Integer, Intent (In) :: ncheb
      Real (Kind=dp), Intent (In) :: mat1(0:ncheb, 0:ncheb), vec1(0:ncheb)
      Real (Kind=dp), Intent (Out) :: outvec(0:ncheb)
      Integer :: n, m

      m = size(mat1, 1)
      n = size(mat1, 2)
      If (size(vec1,1)/=n) Stop &
        'matvec_dmdm: dimensions of first input array differ.'
      Call dgemv('N', m, n, 1.0E0_dp, mat1, m, vec1, 1, 0.0E0_dp, outvec, 1)
    End Subroutine

    Subroutine matmat_dmdm(mat1, mat2, ncheb, outmat)
      Use mod_datatypes
      Implicit None
      Integer, Intent (In) :: ncheb
      Real (Kind=dp), Intent (In) :: mat1(0:ncheb, 0:ncheb), &
        mat2(0:ncheb, 0:ncheb)
      Real (Kind=dp), Intent (Out) :: outmat(0:ncheb, 0:ncheb)

      Integer :: n

      n = ncheb + 1
      Call dgemm('N', 'N', n, n, n, 1E0_dp, mat1, n, mat2, n, 0E0_dp, outmat, &
        n)
    End Subroutine
