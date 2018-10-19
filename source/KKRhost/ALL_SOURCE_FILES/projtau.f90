!------------------------------------------------------------------------------------
!> Summary: Calculate the component projected TAU - matrices
!> Author: 
!> Calculate the component projected TAU - matrices
!> \begin{equation}
!> TAU(IT) =  TAU(IQ) \left( 1 + (m(t)-m(c)) TAU(IQ) \right)^{-1}
!> \end{equation}
!------------------------------------------------------------------------------------
!> @note It is assumed that all equivalent sites IQ have the same TAU-matrix 
!> TAUQ(IQ). To get TAU(IT) the first site IQ occupied by type IT is taken to be 
!> representative for all other (NAT(IT)-1) sites occupied by IT.
!>
!> Allows an atom type IT to have different orientation of its moment on different 
!> but equivalent sites IQ
!> @endnote
!------------------------------------------------------------------------------------
module mod_projtau
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculate the component projected TAU - matrices
  !> Author: 
  !> Category: single-site, coherent-potential-approx, KKRhost 
  !> Deprecated: False 
  !> Calculate the component projected TAU - matrices
  !> \begin{equation}
  !> TAU(IT) = TAU(IQ) \left( 1 + (m(t)-m(c)) TAU(IQ) \right)^{-1}
  !> \end{equation}
  !-------------------------------------------------------------------------------
  !> @note It is assumed that all equivalent sites IQ have the same TAU-matrix 
  !> TAUQ(IQ). To get TAU(IT) the first site IQ occupied by type IT is taken to be 
  !> representative for all other (NAT(IT)-1) sites occupied by IT.
  !>
  !> Allows an atom type IT to have different orientation of its moment on different 
  !> but equivalent sites IQ
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine projtau(icpaflag,cpachng,kmrot,wrtau,wrtaumq,ifiltau,eryd,nt,nq,nkmq,  &
    msst,mssq,nlinq,iqat,conc,tauq,taut,tautlin,ikm1lin,ikm2lin,drotq,ntmax,nqmax,  &
    nkmmax,linmax)

    use :: mod_getdmat
    use :: mod_rotate
    use :: constants, only: cone,czero
    implicit none

    ! PARAMETER definitions
    real (kind=dp) :: tol
    parameter (tol=1.0e-6_dp)

    ! Dummy arguments
    real (kind=dp) :: cpachng
    complex (kind=dp) :: eryd
    integer :: icpaflag, ifiltau, kmrot, linmax, nkmmax, nq, nqmax, nt, ntmax
    logical :: wrtau, wrtaumq
    real (kind=dp) :: conc(ntmax)
    complex (kind=dp) :: drotq(nkmmax, nkmmax, nqmax), mssq(nkmmax, nkmmax, nqmax), msst(nkmmax, nkmmax, ntmax), tauq(nkmmax, nkmmax, nqmax), taut(nkmmax, nkmmax, ntmax), &
      tautlin(linmax, ntmax)
    integer :: ikm1lin(linmax), ikm2lin(linmax), iqat(ntmax), nkmq(nqmax), nlinq(nqmax)

    ! Local variables
    real (kind=dp) :: cpac
    complex (kind=dp) :: dmamc(nkmmax, nkmmax), dmattg(nkmmax, nkmmax), dtiltg(nkmmax, nkmmax), rmss, rtau, w1(nkmmax, nkmmax)
    integer :: i, icpaf, iq, it, j, lin, m, n

    do it = 1, nt

      ! ---------- pick first site IQ occupied by type IT to be representative
      ! ----------- all other (NAT(IT)-1) occupied sites have to be equivalent

      iq = iqat(it)
      m = nkmmax
      n = nkmq(iq)

      if (conc(it)<0.995_dp) then
        ! ------------------------- rotate the single site m-matrix if necessary
        if (kmrot/=0) then
          call rotate(msst(1,1,it), 'L->G', w1, n, drotq(1,1,iq), m)
          call getdmat(tauq(1,1,iq), dmattg, dtiltg, dmamc, n, mssq(1,1,iq), w1, m)
        else
          call getdmat(tauq(1,1,iq), dmattg, dtiltg, dmamc, n, mssq(1,1,iq), msst(1,1,it), m)
        end if
        ! -------------------------------------------
        ! TAU(t) = TAU * D~(t)
        ! ----------------------------------------
        call zgemm('N', 'N', n, n, n, cone, tauq(1,1,iq), m, dtiltg, m, czero, taut(1,1,it), m)
        icpaf = icpaflag
        cpac = cpachng
      else
        ! CONC > 0.995:  COPY TAU TO TAUTLIN
        do j = 1, n
          call zcopy(n, tauq(1,j,iq), 1, taut(1,j,it), 1)
        end do
        icpaf = 0
        cpac = 0e0_dp
      end if

      ! -------------------------------------------
      ! rotate  TAU(t)  if required
      ! -------------------------------------------
      if (kmrot/=0) then
        do j = 1, n
          call zcopy(n, taut(1,j,it), 1, w1(1,j), 1)
        end do
        call rotate(w1, 'G->L', taut(1,1,it), n, drotq(1,1,iq), m)
      end if
      ! -------------------------------------------
      ! STORE TAU(t) IN LINEAR ARRAY TAUTLIN
      ! -------------------------------------------

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
              if (abs(taut(i,j,it)/taut(i,i,it))>tol) write (ifiltau, 120) i, j, taut(i, j, it)
            end if
          end do
        end do
      end if

    end do
    ! ================================================================= IT ==

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
              if ((abs(rtau)>tol) .or. (abs(rmss)>tol)) write (ifiltau, 120) i, j, tauq(i, j, iq), mssq(i, j, iq)
            end if

          end do
        end do

      end do
    end if
    ! --------------------------------------------------------- FORMAT IFMT=2
100 format (/, 80('*'), /, 2f21.15, ' RYD   TAU FOR IT=', i2, '  IQ=', i2, :, '  CPA:', i2, f15.6)
110 format (/, 80('*'), /, 2f21.15, ' RYD   TAU-C M-C  FOR IQ=', i2, :, '  CPA:', i2, f15.6)
120 format (2i5, 1p, 4e22.14)

  end subroutine projtau

end module mod_projtau
