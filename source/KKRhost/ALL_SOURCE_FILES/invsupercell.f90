module mod_invsupercell

contains

  !-------------------------------------------------------------------------------
  !> Summary: Matrix inversion for supercell geometry
  !> Author: 
  !> Category: KKRhost, structural-greensfunction
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> ---> ALGORITM FOR SUPERCELL GEOMETRY
  !>
  !> ---> factorization D ^-1 = (prod L) * M * (prod U)
  !>
  !> see notes R. Zeller
  !-------------------------------------------------------------------------------
  subroutine invsupercell(m2, m1, m3, gin, icheck)

    use :: global_variables, only: ndim_slabinv, nlayerd, alm
    use :: mod_datatypes, only: dp
    use :: mod_bofm, only: bofm
    use :: mod_btom, only: btom
    use :: mod_cinit, only: cinit
    use :: mod_constants, only: czero, cone
    implicit none

    complex (kind=dp) :: m1(ndim_slabinv, ndim_slabinv, nlayerd), m2(ndim_slabinv, ndim_slabinv, nlayerd), m3(ndim_slabinv, ndim_slabinv, nlayerd), gin(alm, alm)
    integer :: n, lm, info, irow, icol, nl
    integer :: ipvt(ndim_slabinv), icheck(nlayerd, nlayerd)
    complex (kind=dp) :: a(ndim_slabinv, ndim_slabinv), b(ndim_slabinv, ndim_slabinv, nlayerd), c(ndim_slabinv, ndim_slabinv, nlayerd), d(ndim_slabinv, ndim_slabinv, nlayerd), &
      e(ndim_slabinv, ndim_slabinv), f(ndim_slabinv, ndim_slabinv), g(ndim_slabinv, ndim_slabinv), cunit(ndim_slabinv, ndim_slabinv)

    ! ---> START OF THE FACTORIZATION L * M * U
    !
    ! ---> N =1
    !
    !     initialize all the matrices
    call cinit(ndim_slabinv*ndim_slabinv, a)
    call cinit(ndim_slabinv*ndim_slabinv*nlayerd, b)
    call cinit(ndim_slabinv*ndim_slabinv*nlayerd, c)
    call cinit(ndim_slabinv*ndim_slabinv*nlayerd, d)
    call cinit(ndim_slabinv*ndim_slabinv, e)
    call cinit(ndim_slabinv*ndim_slabinv, f)
    call cinit(ndim_slabinv*ndim_slabinv, g)
    call cinit(ndim_slabinv*ndim_slabinv, cunit)

    ! ------------------------------------------------------------------------
    !
    ! ---> cunit = complex unity matrix of order NDIM

    do n = 1, ndim_slabinv
      cunit(n, n) = cone
    end do

    call zcopy(ndim_slabinv*ndim_slabinv, m2(1,1,1), 1, e(1,1), 1)
    call zcopy(ndim_slabinv*ndim_slabinv, cunit, 1, d(1,1,1), 1)
    call zgetrf(ndim_slabinv, ndim_slabinv, e(1,1), ndim_slabinv, ipvt, info)
    call zgetrs('N', ndim_slabinv, ndim_slabinv, e(1,1), ndim_slabinv, ipvt, d(1,1,1), ndim_slabinv, info)

    nl = nlayerd

    if (nl==1) go to 120

    call zcopy(ndim_slabinv*ndim_slabinv, m2(1,1,nl), 1, a(1,1), 1)
    call zcopy(ndim_slabinv*ndim_slabinv, m1(1,1,nl), 1, b(1,1,1), 1)
    call zcopy(ndim_slabinv*ndim_slabinv, m3(1,1,nl), 1, c(1,1,1), 1)

    ! ------------------------------------------------------------------------
    !
    ! ---> 2 <= N < NL-1

    if (nl==2) go to 110

    if (nl==3) go to 100

    do n = 2, nl - 2
      !---> E = D(N-1) * C(N-1)
      call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, d(1,1,n-1), ndim_slabinv, c(1,1,n-1), ndim_slabinv, czero, e(1,1), ndim_slabinv)
      !---> F = D(N-1) * M1(N-1)
      call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, d(1,1,n-1), ndim_slabinv, m1(1,1,n-1), ndim_slabinv, czero, f(1,1), ndim_slabinv)
      !---> A = A - B(N-1)*D(N-1)*C(N-1)
      call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, b(1,1,n-1), ndim_slabinv, e(1,1), ndim_slabinv, cone, a(1,1), ndim_slabinv)
      !---> B(N) = - B(N-1)*D(N-1)*M1(N-1)
      call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, b(1,1,n-1), ndim_slabinv, f(1,1), ndim_slabinv, czero, b(1,1,n), ndim_slabinv)
      !---> C(N) = - M3(N-1)*D(N-1)*C(N-1)
      call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, m3(1,1,n-1), ndim_slabinv, e(1,1), ndim_slabinv, czero, c(1,1,n), ndim_slabinv)
      !---> D(N) = [ M2(N) - M3(N-1)*D(N-1)*M1(N-1) ]^-1
      call zcopy(ndim_slabinv*ndim_slabinv, m2(1,1,n), 1, e(1,1), 1)
      call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, m3(1,1,n-1), ndim_slabinv, f(1,1), ndim_slabinv, cone, e(1,1), ndim_slabinv)
      call zcopy(ndim_slabinv*ndim_slabinv, cunit(1,1), 1, d(1,1,n), 1)
      call zgetrf(ndim_slabinv, ndim_slabinv, e(1,1), ndim_slabinv, ipvt, info)
      call zgetrs('N', ndim_slabinv, ndim_slabinv, e(1,1), ndim_slabinv, ipvt, d(1,1,n), ndim_slabinv, info)

    end do

    !------------------------------------------------------------------------
    !
    ! ---> N = NL - 1

100 n = nl - 1

    !---> E = D(N-1) * C(N-1)
    call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, d(1,1,n-1), ndim_slabinv, c(1,1,n-1), ndim_slabinv, czero, e(1,1), ndim_slabinv)
    !---> F = D(N-1) * M1(N-1)
    call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, d(1,1,n-1), ndim_slabinv, m1(1,1,n-1), ndim_slabinv, czero, f(1,1), ndim_slabinv)
    !---> A = A - B(N-1)*D(N-1)*C(N-1)
    call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, b(1,1,n-1), ndim_slabinv, e(1,1), ndim_slabinv, cone, a(1,1), ndim_slabinv)
    !---> B(N) = - B(N-1)*D(N-1)*M1(N-1) + M3(N)
    call zcopy(ndim_slabinv*ndim_slabinv, m3(1,1,n), 1, b(1,1,n), 1)
    call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, b(1,1,n-1), ndim_slabinv, f(1,1), ndim_slabinv, cone, b(1,1,n), ndim_slabinv)
    !---> C(N) = - M3(N-1)*D(N-1)*C(N-1) + M1(N)
    call zcopy(ndim_slabinv*ndim_slabinv, m1(1,1,n), 1, c(1,1,n), 1)
    call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, m3(1,1,n-1), ndim_slabinv, e(1,1), ndim_slabinv, cone, c(1,1,n), ndim_slabinv)
    !---> D(N) = [ M2(N) - M3(N-1)*D(N-1)*M1(N-1) ]^-1
    call zcopy(ndim_slabinv*ndim_slabinv, m2(1,1,n), 1, e(1,1), 1)
    call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, m3(1,1,n-1), ndim_slabinv, f(1,1), ndim_slabinv, cone, e(1,1), ndim_slabinv)
    call zcopy(ndim_slabinv*ndim_slabinv, cunit(1,1), 1, d(1,1,n), 1)
    call zgetrf(ndim_slabinv, ndim_slabinv, e(1,1), ndim_slabinv, ipvt, info)
    call zgetrs('N', ndim_slabinv, ndim_slabinv, e(1,1), ndim_slabinv, ipvt, d(1,1,n), ndim_slabinv, info)

    ! ------------------------------------------------------------------------
    !
    ! ---> N = NL

110 n = nl

    ! ---> D(NL) = (A - B(NL-1)*D(NL-1)*C(NL-1))^-1
    !
    !
    ! ---> E = D(NL-1) * C(NL-1)
    call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, d(1,1,nl-1), ndim_slabinv, c(1,1,nl-1), ndim_slabinv, czero, e(1,1), ndim_slabinv)
    !---> A = A - B(NL-1) * E
    call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, b(1,1,nl-1), ndim_slabinv, e(1,1), ndim_slabinv, cone, a(1,1), ndim_slabinv)
    !---> D(NL) = (A)^-1
    call zcopy(ndim_slabinv*ndim_slabinv, cunit(1,1), 1, d(1,1,nl), 1)
    call zgetrf(ndim_slabinv, ndim_slabinv, a(1,1), ndim_slabinv, ipvt, info)
    call zgetrs('N', ndim_slabinv, ndim_slabinv, a(1,1), ndim_slabinv, ipvt, d(1,1,nl), ndim_slabinv, info)


120 continue ! jump label for NL=1

    ! ------------------------------------------------------------------------
    !
    ! --->  END OF FACTORIZATION
    !
    ! ------------------------------------------------------------------------

    !
    ! ---> HERE IT STARTS LOOP OVER THE DIAGONAL ELEMENT.
    ! ---> THE PROGRAM CHECKS ID ICHECK(N,N) = 1 AND, IF SO,
    ! ---> IT CALCULATES THE DIAGONAL BLOCK (N,N)
    ! ---> THEN IT MAKES TWO LOOPS, ONE OVER A ROW INDEX `IROW`
    ! ---> AND THE OTHER OVER A COLUMN INDEX `ICOL`, AND CHECKS
    ! ---> WHICH ARE THE ELEMENTS THAT HAS TO BE CALCULATED.
    ! ---> (THE ONES FOR WHICH ICHECK = 1)

    !
    ! ---> IT STARTS THE LOOP OVER N

    do n = nl, 1, (-1) ! START OF THE LOOP OVER THE DIAGONAL

      if (n==nl) then
        !---> GTOT(NL,NL) = D(NL)
        call zcopy(ndim_slabinv*ndim_slabinv, d(1,1,nl), 1, e(1,1), 1)
        call btom(nl, nl, e, ndim_slabinv, gin, alm, .false.)

      else if (n==nl-1) then

        if (icheck(nl-1,nl-1)==1) then

          !---> GTOT(NL-1,NL-1) = D(NL-1) + D(NL-1)*C(NL-1)*D(NL)*B(NL-1)*D(NL-1)
          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, b(1,1,nl-1), ndim_slabinv, d(1,1,nl-1), ndim_slabinv, czero, e(1,1), ndim_slabinv)
          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, d(1,1,nl), ndim_slabinv, e(1,1), ndim_slabinv, czero, f(1,1), ndim_slabinv)
          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, c(1,1,nl-1), ndim_slabinv, f(1,1), ndim_slabinv, czero, e(1,1), ndim_slabinv)

          do lm = 1, ndim_slabinv
            e(lm, lm) = cone + e(lm, lm)
          end do

          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, d(1,1,nl-1), ndim_slabinv, e(1,1), ndim_slabinv, czero, f(1,1), ndim_slabinv)

          call btom(nl-1, nl-1, f, ndim_slabinv, gin, alm, .false.)

        end if

      else
        if (icheck(n,n)==1) then

          ! ---> GTOT(N,N) = D(N) + D(N)*( M(N,N+1)*GTOT(N+1,N+1) +
          !                  + C(N)*Z(NL,N+1) )*M(N+1,N)*D(N) -
          !                  - Z(N,NL)*B(N)*D(N)
          call bofm(nl, n+1, f, ndim_slabinv, gin, alm)

          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, c(1,1,n), ndim_slabinv, f(1,1), ndim_slabinv, czero, e(1,1), ndim_slabinv)

          call bofm(n+1, n+1, f, ndim_slabinv, gin, alm)

          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, m1(1,1,n), ndim_slabinv, f(1,1), ndim_slabinv, cone, e(1,1), ndim_slabinv)

          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, m3(1,1,n), ndim_slabinv, d(1,1,n), ndim_slabinv, czero, f(1,1), ndim_slabinv)

          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, e(1,1), ndim_slabinv, f(1,1), ndim_slabinv, czero, g(1,1), ndim_slabinv)

          do lm = 1, ndim_slabinv
            g(lm, lm) = cone + g(lm, lm)
          end do

          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, d(1,1,n), ndim_slabinv, g(1,1), ndim_slabinv, czero, e(1,1), ndim_slabinv)

          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, b(1,1,n), ndim_slabinv, d(1,1,n), ndim_slabinv, czero, f(1,1), ndim_slabinv)

          call bofm(n, nl, g, ndim_slabinv, gin, alm)

          call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, g(1,1), ndim_slabinv, f(1,1), ndim_slabinv, cone, e(1,1), ndim_slabinv)

          call btom(n, n, e, ndim_slabinv, gin, alm, .false.)

        end if

      end if


      do irow = (n-1), 1, (-1) ! LOOP OVER THE ROW FOR THE COLUMN N

        if (icheck(irow,n)==1) then

          if ((n==nl) .and. (irow==(nl-1))) then

            call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, d(1,1,nl-1), ndim_slabinv, c(1,1,nl-1), ndim_slabinv, czero, f(1,1), ndim_slabinv)
            ! M(I,I+1) * Z(I+1,J)
            call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, f(1,1), ndim_slabinv, d(1,1,nl), ndim_slabinv, czero, g(1,1), ndim_slabinv)

            call btom(nl-1, nl, g, ndim_slabinv, gin, alm, .false.)

          else

            ! M(I,I+1) * Z(I+1,J)
            call bofm(irow+1, n, e, ndim_slabinv, gin, alm)
            call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, m1(1,1,irow), ndim_slabinv, e(1,1), ndim_slabinv, czero, g(1,1), ndim_slabinv)

            ! M(I,I+1) * Z(I+1,J) + C(I) * Z(N,J)
            call bofm(nl, n, e, ndim_slabinv, gin, alm)
            call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, c(1,1,irow), ndim_slabinv, e(1,1), ndim_slabinv, cone, g(1,1), ndim_slabinv)

            ! -D(I) * ( M(I,I+1)*Z(I+1,J)+C(I)*Z(N,J) )
            call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, d(1,1,irow), ndim_slabinv, g(1,1), ndim_slabinv, czero, e(1,1), ndim_slabinv)
            call btom(irow, n, e, ndim_slabinv, gin, alm, .false.)

          end if

        end if
      end do ! LOOP OVER THE ROW FOR THE COLUMN N


      do icol = (n-1), 1, (-1) ! LOOP OVER THE COLUMN FOR THE ROW N

        if (icheck(n,icol)==1) then

          if ((n==nl) .and. (icol==(nl-1))) then

            call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, b(1,1,nl-1), ndim_slabinv, d(1,1,nl-1), ndim_slabinv, czero, e(1,1), ndim_slabinv)
            call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, d(1,1,nl), ndim_slabinv, e(1,1), ndim_slabinv, czero, g(1,1), ndim_slabinv)
            call btom(nl, nl-1, g, ndim_slabinv, gin, alm, .false.)

          else

            ! Z(I,J+1) * M(J+1,J)
            call bofm(n, icol+1, e, ndim_slabinv, gin, alm)
            call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, e(1,1), ndim_slabinv, m3(1,1,icol), ndim_slabinv, czero, g(1,1), ndim_slabinv)

            ! Z(I,J+1) * M(J+1,J) + Z(I,N) * B(J)
            call bofm(n, nl, e, ndim_slabinv, gin, alm)
            call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, cone, e(1,1), ndim_slabinv, b(1,1,icol), ndim_slabinv, cone, g(1,1), ndim_slabinv)

            ! -( Z(I,J+1) * M(J+1,J)+Z(I,N) * B(J) ) * D(J)
            call zgemm('N', 'N', ndim_slabinv, ndim_slabinv, ndim_slabinv, -cone, g(1,1), ndim_slabinv, d(1,1,icol), ndim_slabinv, czero, e(1,1), ndim_slabinv)
            call btom(n, icol, e, ndim_slabinv, gin, alm, .false.)

          end if


        end if

      end do ! LOOP OVER THE COLUMN FOR THE ROW N


    end do ! END OF THE LOOP OVER THE DIAGONAL

    return

  end subroutine invsupercell

end module mod_invsupercell
