subroutine projtau(icpaflag, cpachng, kmrot, wrtau, wrtaumq, ifiltau, eryd, &
  nt, nq, nkmq, msst, mssq, nlinq, iqat, conc, tauq, taut, tautlin, ikm1lin, &
  ikm2lin, drotq, ntmax, nqmax, nkmmax, linmax)
!   ********************************************************************
!   *                                                                  *
!   *   calculate the component projected TAU - matrices               *
!   *                                                                  *
!   *      TAU(IT) =  TAU(IQ) * ( 1 + (m(t)-m(c))*TAU(IQ) )**(-1)      *
!   *                                                                  *
!   *   NOTE: it is assumed that all equivalent sites  IQ  have the    *
!   *   same TAU-matrix  TAUQ(IQ). To get  TAU(IT)  the first site IQ  *
!   *   occupied by type IT is taken to be representative for          *
!   *   all other (NAT(IT)-1) sites occupied by IT                     *
!   *                                                                  *
!   *   allows an atom type IT to have different orientation of        *
!   *   its moment on different but equivalent sites  IQ               *
!   *                                                                  *
!   * 01/11/00                                                         *
!   ********************************************************************
  implicit none

! PARAMETER definitions
  double precision :: tol
  parameter (tol=1.0d-6)
  double complex :: c0, c1
  parameter (c0=(0.0d0,0.0d0), c1=(1.0d0,0.0d0))

! Dummy arguments
  double precision :: cpachng
  double complex :: eryd
  integer :: icpaflag, ifiltau, kmrot, linmax, nkmmax, nq, nqmax, nt, ntmax
  logical :: wrtau, wrtaumq
  double precision :: conc(ntmax)
  double complex :: drotq(nkmmax, nkmmax, nqmax), mssq(nkmmax, nkmmax, nqmax), &
    msst(nkmmax, nkmmax, ntmax), tauq(nkmmax, nkmmax, nqmax), &
    taut(nkmmax, nkmmax, ntmax), tautlin(linmax, ntmax)
  integer :: ikm1lin(linmax), ikm2lin(linmax), iqat(ntmax), nkmq(nqmax), &
    nlinq(nqmax)

! Local variables
  double precision :: cpac
  double complex :: dmamc(nkmmax, nkmmax), dmattg(nkmmax, nkmmax), &
    dtiltg(nkmmax, nkmmax), rmss, rtau, w1(nkmmax, nkmmax)
  integer :: i, icpaf, iq, it, j, lin, m, n

  do it = 1, nt

! ---------- pick first site IQ occupied by type IT to be representative
! ----------- all other (NAT(IT)-1) occupied sites have to be equivalent

    iq = iqat(it)
    m = nkmmax
    n = nkmq(iq)

    if (conc(it)<0.995) then

! ------------------------- rotate the single site m-matrix if necessary
      if (kmrot/=0) then

        call rotate(msst(1,1,it), 'L->G', w1, n, drotq(1,1,iq), m)

        call getdmat(tauq(1,1,iq), dmattg, dtiltg, dmamc, n, mssq(1,1,iq), w1, &
          m)

      else

        call getdmat(tauq(1,1,iq), dmattg, dtiltg, dmamc, n, mssq(1,1,iq), &
          msst(1,1,it), m)

      end if

!     -------------------------------------------
!              TAU(t) = TAU * D~(t)
!     ----------------------------------------
      call zgemm('N', 'N', n, n, n, c1, tauq(1,1,iq), m, dtiltg, m, c0, &
        taut(1,1,it), m)

      icpaf = icpaflag
      cpac = cpachng

    else

!     CONC > 0.995:  COPY TAU TO TAUTLIN

      do j = 1, n
        call zcopy(n, tauq(1,j,iq), 1, taut(1,j,it), 1)
      end do

      icpaf = 0
      cpac = 0d0

    end if

!     -------------------------------------------
!            rotate  TAU(t)  if required
!     -------------------------------------------

    if (kmrot/=0) then

      do j = 1, n
        call zcopy(n, taut(1,j,it), 1, w1(1,j), 1)
      end do

      call rotate(w1, 'G->L', taut(1,1,it), n, drotq(1,1,iq), m)

    end if

!     -------------------------------------------
!        STORE TAU(t) IN LINEAR ARRAY TAUTLIN
!     -------------------------------------------

    do lin = 1, nlinq(iq)
      tautlin(lin, it) = taut(ikm1lin(lin), ikm2lin(lin), it)
    end do

    if (wrtau) then
      write (ifiltau, 100) eryd, it, iq, icpaf, cpac
      do i = 1, n
        do j = 1, n
          if (i==j) then
            write (ifiltau, 120) i, j, taut(i, j, it)
          else
            if (cdabs(taut(i,j,it)/taut(i,i,it))>tol) write (ifiltau, 120) i, &
              j, taut(i, j, it)
          end if
        end do
      end do
    end if

  end do
!================================================================= IT ==

  if (wrtaumq) then
    do iq = 1, nq
      write (ifiltau, 110) eryd, iq, icpaflag, cpachng
      do i = 1, n
        do j = 1, n
          if (i==j) then
            write (ifiltau, 120) i, j, tauq(i, j, iq), mssq(i, j, iq)
          else
            rtau = tauq(i, j, iq)/tauq(i, i, iq)
            rmss = mssq(i, j, iq)/mssq(i, i, iq)
            if ((cdabs(rtau)>tol) .or. (cdabs(rmss)>tol)) write (ifiltau, 120) &
              i, j, tauq(i, j, iq), mssq(i, j, iq)
          end if

        end do
      end do

    end do
  end if
!--------------------------------------------------------- FORMAT IFMT=2
100 format (/, 80('*'), /, 2f21.15, ' RYD   TAU FOR IT=', i2, '  IQ=', i2, :, &
    '  CPA:', i2, f15.6)
110 format (/, 80('*'), /, 2f21.15, ' RYD   TAU-C M-C  FOR IQ=', i2, :, &
    '  CPA:', i2, f15.6)
120 format (2i5, 1p, 4e22.14)

end subroutine
