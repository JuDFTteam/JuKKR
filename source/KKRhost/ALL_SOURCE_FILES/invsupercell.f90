! ************************************************************************
subroutine invsupercell(m2, m1, m3, gin, icheck)
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

  implicit none

!.. Parameters ..
  include 'inc.p'
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************
!


!.. array arguments

!.. local scalars
  integer, parameter :: lmmaxd = (krel+korbit+1)*(lmaxd+1)**2
  integer, parameter :: ndim = max(nprincd*lmmaxd, 1)
  integer, parameter :: almd = naezd*lmmaxd
  double complex, parameter :: czero = (0.d0, 0.d0)
  double complex, parameter :: cone = (1.d0, 0.d0)

!.. local arrays
  double complex :: m1(ndim, ndim, nlayerd), m2(ndim, ndim, nlayerd), &
    m3(ndim, ndim, nlayerd), gin(almd, almd)
! --------------------------------------------------------------------
!.. External Subroutines ..
  integer :: n, lm, info, irow, icol, nl


  integer :: ipvt(ndim), icheck(nlayerd, nlayerd)
  double complex :: a(ndim, ndim), b(ndim, ndim, nlayerd), &
    c(ndim, ndim, nlayerd), d(ndim, ndim, nlayerd), e(ndim, ndim), &
    f(ndim, ndim), g(ndim, ndim), cunit(ndim, ndim)
! ---> START OF THE FACTORIZATION L * M * U

  external :: cinit, zcopy, zgemm, zgetrf, zgetrs, btom
  intrinsic :: abs, dimag

! ---> N =1

!     initialize all the matricex


! ------------------------------------------------------------------------

! ---> cunit = complex unity matrix of order NDIM
  call cinit(ndim*ndim, a)
  call cinit(ndim*ndim*nlayerd, b)
  call cinit(ndim*ndim*nlayerd, c)
  call cinit(ndim*ndim*nlayerd, d)
  call cinit(ndim*ndim, e)
  call cinit(ndim*ndim, f)
  call cinit(ndim*ndim, g)
  call cinit(ndim*ndim, cunit)





  do n = 1, ndim
    cunit(n, n) = cone
  end do


  call zcopy(ndim*ndim, m2(1,1,1), 1, e(1,1), 1)
  call zcopy(ndim*ndim, cunit, 1, d(1,1,1), 1)
  call zgetrf(ndim, ndim, e(1,1), ndim, ipvt, info)
  call zgetrs('N', ndim, ndim, e(1,1), ndim, ipvt, d(1,1,1), ndim, info)


  nl = nlayerd

! ------------------------------------------------------------------------
  if (nl==1) go to 120

  call zcopy(ndim*ndim, m2(1,1,nl), 1, a(1,1), 1)
  call zcopy(ndim*ndim, m1(1,1,nl), 1, b(1,1,1), 1)
  call zcopy(ndim*ndim, m3(1,1,nl), 1, c(1,1,1), 1)
! ---> 2 <= N < NL-1





  if (nl==2) go to 110
! ---> E = D(N-1) * C(N-1)
  if (nl==3) go to 100

  do n = 2, nl - 2



! ---> F = D(N-1) * M1(N-1)

    call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,n-1), ndim, c(1,1,n-1), &
      ndim, czero, e(1,1), ndim)



! ---> A = A - B(N-1)*D(N-1)*C(N-1)

    call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,n-1), ndim, &
      m1(1,1,n-1), ndim, czero, f(1,1), ndim)



! ---> B(N) = - B(N-1)*D(N-1)*M1(N-1)

    call zgemm('N', 'N', ndim, ndim, ndim, -cone, b(1,1,n-1), ndim, e(1,1), &
      ndim, cone, a(1,1), ndim)


! ---> C(N) = - M3(N-1)*D(N-1)*C(N-1)

    call zgemm('N', 'N', ndim, ndim, ndim, -cone, b(1,1,n-1), ndim, f(1,1), &
      ndim, czero, b(1,1,n), ndim)



! ---> D(N) = [ M2(N) - M3(N-1)*D(N-1)*M1(N-1) ]^-1

    call zgemm('N', 'N', ndim, ndim, ndim, -cone, m3(1,1,n-1), ndim, e(1,1), &
      ndim, czero, c(1,1,n), ndim)




! ------------------------------------------------------------------------
    call zcopy(ndim*ndim, m2(1,1,n), 1, e(1,1), 1)
    call zgemm('N', 'N', ndim, ndim, ndim, -cone, m3(1,1,n-1), ndim, f(1,1), &
      ndim, cone, e(1,1), ndim)
    call zcopy(ndim*ndim, cunit(1,1), 1, d(1,1,n), 1)
    call zgetrf(ndim, ndim, e(1,1), ndim, ipvt, info)
    call zgetrs('N', ndim, ndim, e(1,1), ndim, ipvt, d(1,1,n), ndim, info)

! ---> N = NL - 1
  end do





! ---> E = D(N-1) * C(N-1)

100 n = nl - 1



! ---> F = D(N-1) * M1(N-1)

  call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,n-1), ndim, c(1,1,n-1), &
    ndim, czero, e(1,1), ndim)



! ---> A = A - B(N-1)*D(N-1)*C(N-1)

  call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,n-1), ndim, m1(1,1,n-1), &
    ndim, czero, f(1,1), ndim)



! ---> B(N) = - B(N-1)*D(N-1)*M1(N-1) + M3(N)

  call zgemm('N', 'N', ndim, ndim, ndim, -cone, b(1,1,n-1), ndim, e(1,1), &
    ndim, cone, a(1,1), ndim)



! ---> C(N) = - M3(N-1)*D(N-1)*C(N-1) + M1(N)

  call zcopy(ndim*ndim, m3(1,1,n), 1, b(1,1,n), 1)
  call zgemm('N', 'N', ndim, ndim, ndim, -cone, b(1,1,n-1), ndim, f(1,1), &
    ndim, cone, b(1,1,n), ndim)



! ---> D(N) = [ M2(N) - M3(N-1)*D(N-1)*M1(N-1) ]^-1

  call zcopy(ndim*ndim, m1(1,1,n), 1, c(1,1,n), 1)
  call zgemm('N', 'N', ndim, ndim, ndim, -cone, m3(1,1,n-1), ndim, e(1,1), &
    ndim, cone, c(1,1,n), ndim)



! ------------------------------------------------------------------------

  call zcopy(ndim*ndim, m2(1,1,n), 1, e(1,1), 1)
  call zgemm('N', 'N', ndim, ndim, ndim, -cone, m3(1,1,n-1), ndim, f(1,1), &
    ndim, cone, e(1,1), ndim)
  call zcopy(ndim*ndim, cunit(1,1), 1, d(1,1,n), 1)
  call zgetrf(ndim, ndim, e(1,1), ndim, ipvt, info)
  call zgetrs('N', ndim, ndim, e(1,1), ndim, ipvt, d(1,1,n), ndim, info)
! ---> N = NL




! ---> D(NL) = (A - B(NL-1)*D(NL-1)*C(NL-1))^-1

110 n = nl

! ---> E = D(NL-1) * C(NL-1)




! ---> A = A - B(NL-1) * E

  call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,nl-1), ndim, c(1,1,nl-1), &
    ndim, czero, e(1,1), ndim)



! ---> D(NL) = (A)^-1

  call zgemm('N', 'N', ndim, ndim, ndim, -cone, b(1,1,nl-1), ndim, e(1,1), &
    ndim, cone, a(1,1), ndim)



! jump label for NL=1

  call zcopy(ndim*ndim, cunit(1,1), 1, d(1,1,nl), 1)
  call zgetrf(ndim, ndim, a(1,1), ndim, ipvt, info)
  call zgetrs('N', ndim, ndim, a(1,1), ndim, ipvt, d(1,1,nl), ndim, info)
! ------------------------------------------------------------------------

120 continue ! --->  END OF FACTORIZATION

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
  do n = nl, 1, (-1) 

    if (n==nl) then
! ---> GTOT(NL,NL) = D(NL)






      call zcopy(ndim*ndim, d(1,1,nl), 1, e(1,1), 1)

      call btom(nl, nl, e, ndim, gin, almd, .false.)
! ---> GTOT(NL-1,NL-1) = D(NL-1) + D(NL-1)*C(NL-1)*D(NL)*B(NL-1)*D(NL-1)

    else if (n==nl-1) then

      if (icheck(nl-1,nl-1)==1) then




        call zgemm('N', 'N', ndim, ndim, ndim, cone, b(1,1,nl-1), ndim, &
          d(1,1,nl-1), ndim, czero, e(1,1), ndim)

        call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,nl), ndim, e(1,1), &
          ndim, czero, f(1,1), ndim)

        call zgemm('N', 'N', ndim, ndim, ndim, cone, c(1,1,nl-1), ndim, &
          f(1,1), ndim, czero, e(1,1), ndim)
!     write (6,*) 'it calculates the element ','(',nl-1,',',nl-1,')'
        do lm = 1, ndim
          e(lm, lm) = cone + e(lm, lm)
        end do

        call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,nl-1), ndim, &
          e(1,1), ndim, czero, f(1,1), ndim)

        call btom(nl-1, nl-1, f, ndim, gin, almd, .false.)



      end if
! ---> GTOT(N,N) = D(N) + D(N)*( M(N,N+1)*GTOT(N+1,N+1) +
    else
!                  + C(N)*Z(NL,N+1) )*M(N+1,N)*D(N) -
      if (icheck(n,n)==1) then
!                  - Z(N,NL)*B(N)*D(N)





        call bofm(nl, n+1, f, ndim, gin, almd)

        call zgemm('N', 'N', ndim, ndim, ndim, cone, c(1,1,n), ndim, f(1,1), &
          ndim, czero, e(1,1), ndim)

        call bofm(n+1, n+1, f, ndim, gin, almd)

        call zgemm('N', 'N', ndim, ndim, ndim, cone, m1(1,1,n), ndim, f(1,1), &
          ndim, cone, e(1,1), ndim)

        call zgemm('N', 'N', ndim, ndim, ndim, cone, m3(1,1,n), ndim, &
          d(1,1,n), ndim, czero, f(1,1), ndim)

        call zgemm('N', 'N', ndim, ndim, ndim, cone, e(1,1), ndim, f(1,1), &
          ndim, czero, g(1,1), ndim)

        do lm = 1, ndim
          g(lm, lm) = cone + g(lm, lm)
        end do

        call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,n), ndim, g(1,1), &
          ndim, czero, e(1,1), ndim)

        call zgemm('N', 'N', ndim, ndim, ndim, cone, b(1,1,n), ndim, d(1,1,n), &
          ndim, czero, f(1,1), ndim)
!     write (6,*) 'it calculates the element ','(',n,',',n,')'
        call bofm(n, nl, g, ndim, gin, almd)

        call zgemm('N', 'N', ndim, ndim, ndim, -cone, g(1,1), ndim, f(1,1), &
          ndim, cone, e(1,1), ndim)

        call btom(n, n, e, ndim, gin, almd, .false.)


! LOOP OVER THE ROW FOR THE COLUMN N
      end if

    end if


    do irow = (n-1), 1, (-1) 

      if (icheck(irow,n)==1) then

        if ((n==nl) .and. (irow==(nl-1))) then

          call zgemm('N', 'N', ndim, ndim, ndim, cone, d(1,1,nl-1), ndim, &
            c(1,1,nl-1), ndim, czero, f(1,1), ndim)
!     M(I,I+1) * Z(I+1,J)
          call zgemm('N', 'N', ndim, ndim, ndim, -cone, f(1,1), ndim, &
            d(1,1,nl), ndim, czero, g(1,1), ndim)

          call btom(nl-1, nl, g, ndim, gin, almd, .false.)

        else

!     M(I,I+1) * Z(I+1,J) + C(I) * Z(N,J)

          call bofm(irow+1, n, e, ndim, gin, almd)

          call zgemm('N', 'N', ndim, ndim, ndim, cone, m1(1,1,irow), ndim, &
            e(1,1), ndim, czero, g(1,1), ndim)

!     -D(I) * ( M(I,I+1)*Z(I+1,J)+C(I)*Z(N,J) )

          call bofm(nl, n, e, ndim, gin, almd)

          call zgemm('N', 'N', ndim, ndim, ndim, cone, c(1,1,irow), ndim, &
            e(1,1), ndim, cone, g(1,1), ndim)


!     write (6,*) 'it calculates the element ','(',irow,',',n,')'
          call zgemm('N', 'N', ndim, ndim, ndim, -cone, d(1,1,irow), ndim, &
            g(1,1), ndim, czero, e(1,1), ndim)

          call btom(irow, n, e, ndim, gin, almd, .false.)

        end if
! LOOP OVER THE ROW FOR THE COLUMN N


      end if
! LOOP OVER THE COLUMN FOR THE ROW N
    end do 


    do icol = (n-1), 1, (-1) 

      if (icheck(n,icol)==1) then

        if ((n==nl) .and. (icol==(nl-1))) then

          call zgemm('N', 'N', ndim, ndim, ndim, cone, b(1,1,nl-1), ndim, &
            d(1,1,nl-1), ndim, czero, e(1,1), ndim)
!     Z(I,J+1) * M(J+1,J)
          call zgemm('N', 'N', ndim, ndim, ndim, -cone, d(1,1,nl), ndim, &
            e(1,1), ndim, czero, g(1,1), ndim)

          call btom(nl, nl-1, g, ndim, gin, almd, .false.)

        else

!     Z(I,J+1) * M(J+1,J) + Z(I,N) * B(J)

          call bofm(n, icol+1, e, ndim, gin, almd)

          call zgemm('N', 'N', ndim, ndim, ndim, cone, e(1,1), ndim, &
            m3(1,1,icol), ndim, czero, g(1,1), ndim)

!     -( Z(I,J+1) * M(J+1,J)+Z(I,N) * B(J) ) * D(J)

          call bofm(n, nl, e, ndim, gin, almd)

          call zgemm('N', 'N', ndim, ndim, ndim, cone, e(1,1), ndim, &
            b(1,1,icol), ndim, cone, g(1,1), ndim)


!     write (6,*) 'it calculates the element ','(',n,',',icol,')'
          call zgemm('N', 'N', ndim, ndim, ndim, -cone, g(1,1), ndim, &
            d(1,1,icol), ndim, czero, e(1,1), ndim)

          call btom(n, icol, e, ndim, gin, almd, .false.)

        end if
! LOOP OVER THE COLUMN FOR THE ROW N


      end if

    end do 


  end do
! ************************************************************************
! ************************************************************************
!
  return
! ---> ALGORITM FOR SUPERCELL GEOMETRY
end subroutine
