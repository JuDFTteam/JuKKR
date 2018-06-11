! ************************************************************************
    Subroutine invslab(gdi, gup, gdow, gin, icheck)
      Use mod_datatypes, Only: dp
! ************************************************************************
! ************************************************************************
!
! ---> ALGORITM FOR SLAB GEOMETRY
!
! ------------------------------------------------------------------------
!
! ---> factorization D ^-1 = (prod L) * M * (prod U)
!
!      see notes R. Zeller
!
! ------------------------------------------------------------------------
      Implicit None

!     .. parameters ..
      Include 'inc.p'
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************

!.. local scalars ..
!.. data statements ..

!.. external subroutines ..
      Integer, Parameter :: lmmaxd = (krel+korbit+1)*(lmaxd+1)**2
      Integer, Parameter :: ndim = nprincd*lmmaxd
      Integer, Parameter :: almd = naezd*lmmaxd
      Complex (Kind=dp), Parameter :: ci = (0.E0_dp, 1.E0_dp)
      Complex (Kind=dp), Parameter :: czero = (0.E0_dp, 0.E0_dp)
      Complex (Kind=dp), Parameter :: cone = (1.E0_dp, 0.E0_dp)
      Integer :: ipvt(ndim)
      Complex (Kind=dp) :: gup(ndim, ndim, nlayerd), &
        gdow(ndim, ndim, nlayerd), gdi(ndim, ndim, nlayerd), &
        dmat(ndim, ndim, nlayerd), dinver(ndim, ndim, nlayerd), f(ndim, ndim), &
        e(ndim, ndim), g(ndim, ndim), cunit(ndim, ndim), gin(almd, almd), &
        gdiold(ndim, ndim, nlayerd)

      Integer :: n, lm, info, i, j, irow
!---> to be changed: set up the triangular matrix.
      Integer :: icheck(nlayerd, nlayerd)
!     first version:

      External :: cinit, zcopy, zgemm, zgetrf, zgetrs, btom





      Call cinit(ndim*ndim, e)
      Call cinit(ndim*ndim, f)
      Call cinit(ndim*ndim, g)
      Call cinit(ndim*ndim, cunit)
      Call cinit(ndim*ndim*nlayerd, dmat)
!---> calculate D_1 = (M_11)**(-1)
      Do n = 1, ndim
        cunit(n, n) = cone
      End Do

      Do n = 1, nlayerd
        Call zcopy(ndim*ndim, gdi(1,1,n), 1, gdiold(1,1,n), 1)
      End Do



!---> claculate D_N (2 <= N <= NLAYERD)

      Call zcopy(ndim*ndim, gdi(1,1,1), 1, e(1,1), 1)
      Call zcopy(ndim*ndim, cunit, 1, dmat(1,1,1), 1)

      Call zgetrf(ndim, ndim, e(1,1), ndim, ipvt, info)
      Call zgetrs('N', ndim, ndim, e(1,1), ndim, ipvt, dmat(1,1,1), ndim, &
        info)
!---> F = D(N-1) * M1(N-1)


      Do n = 2, nlayerd


!---> D(N) = [MDI(N) - MUP(N-1)*DMAT(N-1)*DOWM(N-1) ]^(-1)

        Call zgemm('N', 'N', ndim, ndim, ndim, cone, dmat(1,1,n-1), ndim, &
          gup(1,1,n-1), ndim, czero, f(1,1), ndim)




        Call zcopy(ndim*ndim, gdi(1,1,n), 1, e(1,1), 1)
        Call zcopy(ndim*ndim, cunit(1,1), 1, dmat(1,1,n), 1)

        Call zgemm('N', 'N', ndim, ndim, ndim, -cone, gdow(1,1,n-1), ndim, &
          f(1,1), ndim, cone, e(1,1), ndim)
!     At this point the matrix DMAT(ndim,ndim,nlayerd) contains the
        Call zcopy(ndim*ndim, e(1,1), 1, dinver(1,1,n), 1)
!     matrices [of dimension (ndim,ndim)]  D^n, n=1,..,nlayerd
        Call zgetrf(ndim, ndim, e(1,1), ndim, ipvt, info)
        Call zgetrs('N', ndim, ndim, e(1,1), ndim, ipvt, dmat(1,1,n), ndim, &
          info)

      End Do
!---> calculate Z_n for 1 =< n <= n-1





      Do n = nlayerd, 1, (-1)

        If (n==nlayerd) Then

          Call zcopy(ndim*ndim, dmat(1,1,nlayerd), 1, e(1,1), 1)

          Call btom(nlayerd, nlayerd, e, ndim, gin, almd, .False.)

          Call zcopy(ndim*ndim, dmat(1,1,nlayerd), 1, gdi(1,1,nlayerd), 1)
        Else

          Call zgemm('N', 'N', ndim, ndim, ndim, cone, gdow(1,1,n), ndim, &
            dmat(1,1,n), ndim, czero, f(1,1), ndim)

          Call bofm(n+1, n+1, e, ndim, gin, almd)

          Call zgemm('N', 'N', ndim, ndim, ndim, cone, e(1,1), ndim, f(1,1), &
            ndim, czero, g(1,1), ndim)
          Call zgemm('N', 'N', ndim, ndim, ndim, cone, gup(1,1,n), ndim, &
            g(1,1), ndim, czero, f(1,1), ndim)

          Do lm = 1, ndim
            f(lm, lm) = cone + f(lm, lm)
          End Do

          Call zgemm('N', 'N', ndim, ndim, ndim, cone, dmat(1,1,n), ndim, &
            f(1,1), ndim, czero, e(1,1), ndim)
!     here start the two loops on the row index,
          Call btom(n, n, e, ndim, gin, almd, .False.)
!     in order to span all the matrix
          Call zcopy(ndim*ndim, e(1,1), 1, gdi(1,1,n), 1)
!     and to calculate just the blocks that are needed
        End If
!     for the construction of the cluster of green's function



!     this is the loop for element G_ij with i<j

        If (icheck(n,n)==0) Go To 110

        If (n==1) Go To 100



        Do irow = (n-1), 1, (-1)

          If (icheck(irow,n)==1) Then

            Call bofm(irow+1, n, e, ndim, gin, almd)

            Call zgemm('N', 'N', ndim, ndim, ndim, cone, gup(1,1,irow), ndim, &
              e(1,1), ndim, czero, f(1,1), ndim)

            Call zgemm('N', 'N', ndim, ndim, ndim, -cone, dmat(1,1,irow), &
              ndim, f(1,1), ndim, czero, e(1,1), ndim)

            Call btom(irow, n, e, ndim, gin, almd, .False.)

          End If

!     this is the loop for element G_ij with i>j
        End Do

100     Continue

        If (n==nlayerd) Go To 110



        Do irow = n + 1, nlayerd, 1

          If (icheck(irow,n)==1) Then

            Call zcopy(ndim*ndim, cunit(1,1), 1, e(1,1), 1)

            Call bofm(irow, irow, f, ndim, gin, almd)

            Call zgetrf(ndim, ndim, f(1,1), ndim, ipvt, info)
            Call zgetrs('N', ndim, ndim, f(1,1), ndim, ipvt, e(1,1), ndim, &
              info)

            Do i = 1, ndim
              Do j = 1, ndim
                f(i, j) = gdiold(i, j, irow) - (dinver(i,j,irow)-e(i,j))
              End Do
            End Do

            Call zcopy(ndim*ndim, cunit(1,1), 1, e(1,1), 1)

            Call zgetrf(ndim, ndim, f(1,1), ndim, ipvt, info)
            Call zgetrs('N', ndim, ndim, f(1,1), ndim, ipvt, e(1,1), ndim, &
              info)

            Call zgemm('N', 'N', ndim, ndim, ndim, -cone, e(1,1), ndim, &
              gdow(1,1,irow-1), ndim, czero, f(1,1), ndim)
!     corrected 15.3.2000
            Call bofm(irow-1, n, e, ndim, gin, almd)


            Call zgemm('N', 'N', ndim, ndim, ndim, cone, f(1,1), ndim, e(1,1), &
              ndim, czero, g(1,1), ndim)



            Call btom(irow, n, g, ndim, gin, almd, .False.)


          End If


        End Do

110     Continue

      End Do


! ************************************************************************
! ************************************************************************
! ************************************************************************
      Return
!
    End Subroutine
