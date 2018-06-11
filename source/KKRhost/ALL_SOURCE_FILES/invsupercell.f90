! ************************************************************************
    Subroutine invsupercell(m2, m1, m3, gin, icheck)
      Use mod_datatypes, Only: dp
! ************************************************************************
!
! ---> ALGORITM FOR SUPERCELL GEOMETRY
!
! ------------------------------------------------------------------------
!
! ---> factorization D ^-1 = (prod L) * M * (prod U)
!
!      see notes R. Zeller
!
! ------------------------------------------------------------------------

      Implicit None

!.. Parameters ..
      Include 'inc.p'
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************
!

!.. array arguments

!.. local scalars
      Integer, Parameter :: lmmaxd = (krel+korbit+1)*(lmaxd+1)**2
      Integer, Parameter :: ndim = max(nprincd*lmmaxd, 1)
      Integer, Parameter :: almd = naezd*lmmaxd
      Complex (Kind=dp), Parameter :: czero = (0.E0_dp, 0.E0_dp)
      Complex (Kind=dp), Parameter :: cone = (1.E0_dp, 0.E0_dp)

!.. local arrays
      Complex (Kind=dp) :: m1(ndim, ndim, nlayerd), m2(ndim, ndim, nlayerd), &
        m3(ndim, ndim, nlayerd), gin(almd, almd)
! --------------------------------------------------------------------
!.. External Subroutines ..
      Integer :: n, lm, info, irow, icol, nl


      Integer :: ipvt(ndim), icheck(nlayerd, nlayerd)
      Complex (Kind=dp) :: a(ndim, ndim), b(ndim, ndim, nlayerd), &
        c(ndim, ndim, nlayerd), d(ndim, ndim, nlayerd), e(ndim, ndim), &
        f(ndim, ndim), g(ndim, ndim), cunit(ndim, ndim)
! ---> START OF THE FACTORIZATION L * M * U

      External :: cinit, zcopy, zgemm, zgetrf, zgetrs, btom
      Intrinsic :: abs, aimag

! ---> N =1

!     initialize all the matricex


! ------------------------------------------------------------------------

! ---> cunit = complex unity matrix of order NDIM
      Call cinit(ndim*ndim, a)
      Call cinit(ndim*ndim*nlayerd, b)
      Call cinit(ndim*ndim*nlayerd, c)
      Call cinit(ndim*ndim*nlayerd, d)
      Call cinit(ndim*ndim, e)
      Call cinit(ndim*ndim, f)
      Call cinit(ndim*ndim, g)
      Call cinit(ndim*ndim, cunit)





      Do n = 1, ndim
        cunit(n, n) = cone
      End Do


      Call zcopy(ndim*ndim, m2(1,1,1), 1, e(1,1), 1)
      Call zcopy(ndim*ndim, cunit, 1, d(1,1,1), 1)
      Call zgetrf(ndim, ndim, e(1,1), ndim, ipvt, info)
      Call zgetrs('N', ndim, ndim, e(1,1), ndim, ipvt, d(1,1,1), ndim, info)


      nl = nlayerd

! ------------------------------------------------------------------------
      If (nl==1) Go To 120

      Call zcopy(ndim*ndim, m2(1,1,nl), 1, a(1,1), 1)
      Call zcopy(ndim*ndim, m1(1,1,nl), 1, b(1,1,1), 1)
      Call zcopy(ndim*ndim, m3(1,1,nl), 1, c(1,1,1), 1)
! ---> 2 <= N < NL-1





      If (nl==2) Go To 110
! ---> E = D(N-1) * C(N-1)
      If (nl==3) Go To 100

      Do n = 2, nl - 2



! ---> F = D(N-1) * M1(N-1)

        Call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,n-1), ndim, &
          c(1,1,n-1), ndim, czero, e(1,1), ndim)



! ---> A = A - B(N-1)*D(N-1)*C(N-1)

        Call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,n-1), ndim, &
          m1(1,1,n-1), ndim, czero, f(1,1), ndim)



! ---> B(N) = - B(N-1)*D(N-1)*M1(N-1)

        Call zgemm('N', 'N', ndim, ndim, ndim, -cone, b(1,1,n-1), ndim, &
          e(1,1), ndim, cone, a(1,1), ndim)


! ---> C(N) = - M3(N-1)*D(N-1)*C(N-1)

        Call zgemm('N', 'N', ndim, ndim, ndim, -cone, b(1,1,n-1), ndim, &
          f(1,1), ndim, czero, b(1,1,n), ndim)



! ---> D(N) = [ M2(N) - M3(N-1)*D(N-1)*M1(N-1) ]^-1

        Call zgemm('N', 'N', ndim, ndim, ndim, -cone, m3(1,1,n-1), ndim, &
          e(1,1), ndim, czero, c(1,1,n), ndim)




! ------------------------------------------------------------------------
        Call zcopy(ndim*ndim, m2(1,1,n), 1, e(1,1), 1)
        Call zgemm('N', 'N', ndim, ndim, ndim, -cone, m3(1,1,n-1), ndim, &
          f(1,1), ndim, cone, e(1,1), ndim)
        Call zcopy(ndim*ndim, cunit(1,1), 1, d(1,1,n), 1)
        Call zgetrf(ndim, ndim, e(1,1), ndim, ipvt, info)
        Call zgetrs('N', ndim, ndim, e(1,1), ndim, ipvt, d(1,1,n), ndim, info)

! ---> N = NL - 1
      End Do





! ---> E = D(N-1) * C(N-1)

100   n = nl - 1



! ---> F = D(N-1) * M1(N-1)

      Call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,n-1), ndim, &
        c(1,1,n-1), ndim, czero, e(1,1), ndim)



! ---> A = A - B(N-1)*D(N-1)*C(N-1)

      Call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,n-1), ndim, &
        m1(1,1,n-1), ndim, czero, f(1,1), ndim)



! ---> B(N) = - B(N-1)*D(N-1)*M1(N-1) + M3(N)

      Call zgemm('N', 'N', ndim, ndim, ndim, -cone, b(1,1,n-1), ndim, e(1,1), &
        ndim, cone, a(1,1), ndim)



! ---> C(N) = - M3(N-1)*D(N-1)*C(N-1) + M1(N)

      Call zcopy(ndim*ndim, m3(1,1,n), 1, b(1,1,n), 1)
      Call zgemm('N', 'N', ndim, ndim, ndim, -cone, b(1,1,n-1), ndim, f(1,1), &
        ndim, cone, b(1,1,n), ndim)



! ---> D(N) = [ M2(N) - M3(N-1)*D(N-1)*M1(N-1) ]^-1

      Call zcopy(ndim*ndim, m1(1,1,n), 1, c(1,1,n), 1)
      Call zgemm('N', 'N', ndim, ndim, ndim, -cone, m3(1,1,n-1), ndim, e(1,1), &
        ndim, cone, c(1,1,n), ndim)



! ------------------------------------------------------------------------

      Call zcopy(ndim*ndim, m2(1,1,n), 1, e(1,1), 1)
      Call zgemm('N', 'N', ndim, ndim, ndim, -cone, m3(1,1,n-1), ndim, f(1,1), &
        ndim, cone, e(1,1), ndim)
      Call zcopy(ndim*ndim, cunit(1,1), 1, d(1,1,n), 1)
      Call zgetrf(ndim, ndim, e(1,1), ndim, ipvt, info)
      Call zgetrs('N', ndim, ndim, e(1,1), ndim, ipvt, d(1,1,n), ndim, info)
! ---> N = NL




! ---> D(NL) = (A - B(NL-1)*D(NL-1)*C(NL-1))^-1

110   n = nl

! ---> E = D(NL-1) * C(NL-1)




! ---> A = A - B(NL-1) * E

      Call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,nl-1), ndim, &
        c(1,1,nl-1), ndim, czero, e(1,1), ndim)



! ---> D(NL) = (A)^-1

      Call zgemm('N', 'N', ndim, ndim, ndim, -cone, b(1,1,nl-1), ndim, e(1,1), &
        ndim, cone, a(1,1), ndim)



! jump label for NL=1

      Call zcopy(ndim*ndim, cunit(1,1), 1, d(1,1,nl), 1)
      Call zgetrf(ndim, ndim, a(1,1), ndim, ipvt, info)
      Call zgetrs('N', ndim, ndim, a(1,1), ndim, ipvt, d(1,1,nl), ndim, info)
! ------------------------------------------------------------------------

120   Continue ! --->  END OF FACTORIZATION

! ------------------------------------------------------------------------


! ---> HERE IT STARTS LOOP OVER THE DIAGONAL ELEMENT.
! ---> THE PROGRAM CHECKS ID ICHECK(N,N) = 1 AND, IF SO,
! ---> IT CALCULATES THE DIAGONAL BLOCK (N,N)
! ---> THEN IT MAKES TWO LOOPS, ONE OVER A ROW INDEX `IROW`
! ---> AND THE OTHER OVER A COLUMN INDEX `ICOL`, AND CHECKS
! ---> WHICH ARE THE ELEMENTS THAT HAS TO BE CALCULATED.
! ---> (THE ONES FOR WHICH ICHECK = 1)


! ---> IT STARTS THE LOOP OVER N


! START OF THE LOOP OVER THE DIAGONAL


!     write (6,*) 'it calculates the element ','(',nl,',',nl,')'
      Do n = nl, 1, (-1)

        If (n==nl) Then
! ---> GTOT(NL,NL) = D(NL)






          Call zcopy(ndim*ndim, d(1,1,nl), 1, e(1,1), 1)

          Call btom(nl, nl, e, ndim, gin, almd, .False.)
! ---> GTOT(NL-1,NL-1) = D(NL-1) + D(NL-1)*C(NL-1)*D(NL)*B(NL-1)*D(NL-1)

        Else If (n==nl-1) Then

          If (icheck(nl-1,nl-1)==1) Then




            Call zgemm('N', 'N', ndim, ndim, ndim, cone, b(1,1,nl-1), ndim, &
              d(1,1,nl-1), ndim, czero, e(1,1), ndim)

            Call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,nl), ndim, &
              e(1,1), ndim, czero, f(1,1), ndim)

            Call zgemm('N', 'N', ndim, ndim, ndim, cone, c(1,1,nl-1), ndim, &
              f(1,1), ndim, czero, e(1,1), ndim)
!     write (6,*) 'it calculates the element ','(',nl-1,',',nl-1,')'
            Do lm = 1, ndim
              e(lm, lm) = cone + e(lm, lm)
            End Do

            Call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,nl-1), ndim, &
              e(1,1), ndim, czero, f(1,1), ndim)

            Call btom(nl-1, nl-1, f, ndim, gin, almd, .False.)



          End If
! ---> GTOT(N,N) = D(N) + D(N)*( M(N,N+1)*GTOT(N+1,N+1) +
        Else
!                  + C(N)*Z(NL,N+1) )*M(N+1,N)*D(N) -
          If (icheck(n,n)==1) Then
!                  - Z(N,NL)*B(N)*D(N)





            Call bofm(nl, n+1, f, ndim, gin, almd)

            Call zgemm('N', 'N', ndim, ndim, ndim, cone, c(1,1,n), ndim, &
              f(1,1), ndim, czero, e(1,1), ndim)

            Call bofm(n+1, n+1, f, ndim, gin, almd)

            Call zgemm('N', 'N', ndim, ndim, ndim, cone, m1(1,1,n), ndim, &
              f(1,1), ndim, cone, e(1,1), ndim)

            Call zgemm('N', 'N', ndim, ndim, ndim, cone, m3(1,1,n), ndim, &
              d(1,1,n), ndim, czero, f(1,1), ndim)

            Call zgemm('N', 'N', ndim, ndim, ndim, cone, e(1,1), ndim, f(1,1), &
              ndim, czero, g(1,1), ndim)

            Do lm = 1, ndim
              g(lm, lm) = cone + g(lm, lm)
            End Do

            Call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,n), ndim, &
              g(1,1), ndim, czero, e(1,1), ndim)

            Call zgemm('N', 'N', ndim, ndim, ndim, cone, b(1,1,n), ndim, &
              d(1,1,n), ndim, czero, f(1,1), ndim)
!     write (6,*) 'it calculates the element ','(',n,',',n,')'
            Call bofm(n, nl, g, ndim, gin, almd)

            Call zgemm('N', 'N', ndim, ndim, ndim, -cone, g(1,1), ndim, &
              f(1,1), ndim, cone, e(1,1), ndim)

            Call btom(n, n, e, ndim, gin, almd, .False.)


! LOOP OVER THE ROW FOR THE COLUMN N
          End If

        End If


        Do irow = (n-1), 1, (-1)

          If (icheck(irow,n)==1) Then

            If ((n==nl) .And. (irow==(nl-1))) Then

              Call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,nl-1), ndim, &
                c(1,1,nl-1), ndim, czero, f(1,1), ndim)
!     M(I,I+1) * Z(I+1,J)
              Call zgemm('N', 'N', ndim, ndim, ndim, -cone, f(1,1), ndim, &
                d(1,1,nl), ndim, czero, g(1,1), ndim)

              Call btom(nl-1, nl, g, ndim, gin, almd, .False.)

            Else

!     M(I,I+1) * Z(I+1,J) + C(I) * Z(N,J)

              Call bofm(irow+1, n, e, ndim, gin, almd)

              Call zgemm('N', 'N', ndim, ndim, ndim, cone, m1(1,1,irow), ndim, &
                e(1,1), ndim, czero, g(1,1), ndim)

!     -D(I) * ( M(I,I+1)*Z(I+1,J)+C(I)*Z(N,J) )

              Call bofm(nl, n, e, ndim, gin, almd)

              Call zgemm('N', 'N', ndim, ndim, ndim, cone, c(1,1,irow), ndim, &
                e(1,1), ndim, cone, g(1,1), ndim)


!     write (6,*) 'it calculates the element ','(',irow,',',n,')'
              Call zgemm('N', 'N', ndim, ndim, ndim, -cone, d(1,1,irow), ndim, &
                g(1,1), ndim, czero, e(1,1), ndim)

              Call btom(irow, n, e, ndim, gin, almd, .False.)

            End If
! LOOP OVER THE ROW FOR THE COLUMN N


          End If
! LOOP OVER THE COLUMN FOR THE ROW N
        End Do


        Do icol = (n-1), 1, (-1)

          If (icheck(n,icol)==1) Then

            If ((n==nl) .And. (icol==(nl-1))) Then

              Call zgemm('N', 'N', ndim, ndim, ndim, cone, b(1,1,nl-1), ndim, &
                d(1,1,nl-1), ndim, czero, e(1,1), ndim)
!     Z(I,J+1) * M(J+1,J)
              Call zgemm('N', 'N', ndim, ndim, ndim, -cone, d(1,1,nl), ndim, &
                e(1,1), ndim, czero, g(1,1), ndim)

              Call btom(nl, nl-1, g, ndim, gin, almd, .False.)

            Else

!     Z(I,J+1) * M(J+1,J) + Z(I,N) * B(J)

              Call bofm(n, icol+1, e, ndim, gin, almd)

              Call zgemm('N', 'N', ndim, ndim, ndim, cone, e(1,1), ndim, &
                m3(1,1,icol), ndim, czero, g(1,1), ndim)

!     -( Z(I,J+1) * M(J+1,J)+Z(I,N) * B(J) ) * D(J)

              Call bofm(n, nl, e, ndim, gin, almd)

              Call zgemm('N', 'N', ndim, ndim, ndim, cone, e(1,1), ndim, &
                b(1,1,icol), ndim, cone, g(1,1), ndim)


!     write (6,*) 'it calculates the element ','(',n,',',icol,')'
              Call zgemm('N', 'N', ndim, ndim, ndim, -cone, g(1,1), ndim, &
                d(1,1,icol), ndim, czero, e(1,1), ndim)

              Call btom(n, icol, e, ndim, gin, almd, .False.)

            End If
! LOOP OVER THE COLUMN FOR THE ROW N


          End If

        End Do


      End Do
! ************************************************************************
! ************************************************************************
!
      Return
! ---> ALGORITM FOR SUPERCELL GEOMETRY
    End Subroutine
