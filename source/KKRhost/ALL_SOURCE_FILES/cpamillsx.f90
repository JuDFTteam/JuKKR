subroutine cpamillsx(itcpa, cpaerr, cpacorr, cpachng, iprint, icpa, nq, nkmq, &
  noq, itoq, conc, mssq, msst, tauq, dmssq, kmrot, drotq, ntmax, nqmax, &
  nkmmax)
!   ********************************************************************
!   *                                                                  *
!   *   perform  CPA-iteration according    MILLS's  algorithm         *
!   *                                                                  *
!   *   the CPA - iteration step for site IQ is omitted if             *
!   *   ICPA(IQ) = 0 ( set in <INITALL> )                              *
!   *                                                                  *
!   *   only the projection matrix  DTILT(G) is needed                 *
!   *   this is set up with respect to the global frame                *
!   *   for this reason MSST has to be rotated prior calling <GETDMAT> *
!   *                                                                  *
!   *   allows an atom type IT to have different orientation of        *
!   *   its moment on different but equivalent sites  IQ               *
!   *                                                                  *
!   * 15/12/03                                                         *
!   ********************************************************************
  use :: mod_types, only: t_inc
      Use mod_datatypes, Only: dp
  implicit complex *16(a-h, o-z)

! PARAMETER definitions
  real (kind=dp) :: tol, sclstd
  parameter (tol=10d0, sclstd=1d0)
  complex (kind=dp) :: c0, c1
  parameter (c0=(0.0d0,0.0d0), c1=(1.0d0,0.0d0))

! Dummy arguments
  real (kind=dp) :: cpachng, cpacorr, cpaerr
  integer :: iprint, itcpa, kmrot, nkmmax, nq, nqmax, ntmax
  real (kind=dp) :: conc(ntmax)
  complex (kind=dp) :: drotq(nkmmax, nkmmax, nqmax), mssq(nkmmax, nkmmax, nqmax), &
    dmssq(nkmmax, nkmmax, nqmax), msst(nkmmax, nkmmax, ntmax), &
    tauq(nkmmax, nkmmax, nqmax)
  integer :: icpa(nqmax), itoq(ntmax, nqmax), nkmq(nqmax), noq(nqmax)

! Local variables
  logical :: check
  real (kind=dp) :: cpachngl, cpacorrl, scl
  complex (kind=dp) :: csum
  complex (kind=dp) :: dmamc(nkmmax, nkmmax), dmattg(nkmmax, nkmmax), &
    dq(nkmmax, nqmax), dtiltg(nkmmax, nkmmax), err(nkmmax, nkmmax), &
    w1(nkmmax, nkmmax), w2(nkmmax, nkmmax)
  integer :: i, icparun, info, io
  integer :: ipiv(nkmmax), iq, it, iw, iw0, j, m, n
  real (kind=dp) :: p1, p2
  save :: cpachngl, cpacorrl, scl

  data icparun/0/

  cpaerr = 0.0d0
  cpacorr = 0.0d0
  cpachng = 0.0d0
  check = .true.
  check = .false.

  if (itcpa==1) then
    scl = sclstd
    cpachngl = 1d+20
    cpacorrl = 1d+20
    icparun = icparun + 1
  end if

  do iq = 1, nq
    if (icpa(iq)/=0) then

      m = nkmmax
      n = nkmq(iq)

      do j = 1, n
        dq(j, iq) = -c1
        call cinit(n, err(1,j))
      end do

!================================================================= IT ==
      do io = 1, noq(iq)
        it = itoq(io, iq)

! ------------------------- rotate the single site m-matrix if necessary
        if (kmrot/=0) then

          call rotate(msst(1,1,it), 'L->G', w1, n, drotq(1,1,iq), m)

          call getdmat(tauq(1,1,iq), dmattg, dtiltg, dmamc, n, mssq(1,1,iq), &
            w1, m)

        else

          call getdmat(tauq(1,1,iq), dmattg, dtiltg, dmamc, n, mssq(1,1,iq), &
            msst(1,1,it), m)

        end if

        do i = 1, n
          dq(i, iq) = dq(i, iq) + conc(it)*dtiltg(i, i)
        end do

!     -------------------------------------------
!            - E[a] = D~[a] * ( m[a] - m[c] )
!     -------------------------------------------
        call zgemm('N', 'N', n, n, n, c1, dtiltg(1,1), m, dmamc, m, c0, w1, m)

!     -------------------------------------------
!            E = SUM[a]  c[a] *  E[a]
!     -------------------------------------------
        do j = 1, n
          do i = 1, n
            err(i, j) = err(i, j) - conc(it)*w1(i, j)
          end do
        end do

      end do
!================================================================= IT ==

!     -------------------------------------------
!                   E * TAU
!     -------------------------------------------

      call zgemm('N', 'N', n, n, n, c1, err, m, tauq(1,1,iq), m, c0, w2, m)

!     -------------------------------------------
!                1 + E * TAU
!     -------------------------------------------
      do i = 1, n
        w2(i, i) = c1 + w2(i, i)
      end do

!     -------------------------------------------
!               ( 1 + E * TAU )**(-1)
!     -------------------------------------------

      call zgetrf(n, n, w2, m, ipiv, info)
      call zgetri(n, w2, m, ipiv, w1, m*m, info)

!     -------------------------------------------
!           ( 1 + E * TAU )**(-1) * E
!     -------------------------------------------

      call zgemm('N', 'N', n, n, n, c1, w2, m, err, m, c0, w1, m)

!     -------------------------------------------
!           APPLY CORRECTION  TO  MEFF
!     m{n+1} = m{n} -  ( 1 + E * TAU )**(-1) * E
!     -------------------------------------------
      do j = 1, n
        cpaerr = cpaerr + abs(dreal(dq(j,iq))) + abs(dimag(dq(j,iq)))
        cpacorr = cpacorr + abs(dreal(w1(j,j))) + abs(dimag(w1(j,j)))
        cpachng = max(cpachng, abs(w1(j,j)/mssq(j,j,iq)))
      end do
      cpachng = scl*cpachng
      cpacorr = scl*cpacorr

      if (cpachng>tol*cpachngl .or. cpacorr>cpacorrl) then
        write (*, *) '############### CPA step back'

        p1 = 0.5d0 ! P1 = 0.05D0
        p2 = 1d0 - p1
        cpachng = p1*cpachngl
        cpacorr = p1*cpacorrl
        cpachng = cpachngl
        cpacorr = cpacorrl
        scl = p1
        do j = 1, n
          do i = 1, n
            mssq(i, j, iq) = mssq(i, j, iq) + p2*dmssq(i, j, iq)
            dmssq(i, j, iq) = dmssq(i, j, iq)*p1
          end do
        end do
      else
        do j = 1, n
          do i = 1, n
            w1(i, j) = scl*w1(i, j)
            mssq(i, j, iq) = mssq(i, j, iq) - w1(i, j)
            dmssq(i, j, iq) = w1(i, j)
          end do
        end do
      end if

      cpaerr = cpachng

      if (iprint>=2 .or. check) then
        csum = c0
        do i = 1, n
          csum = csum + mssq(i, i, iq)
        end do
        if (t_inc%i_write>0) then
          write (1337, 100) iq, cpaerr, cpacorr, csum
        end if
      end if
!-----------------------------------------------------------------------
      if (check) then
        if (icparun==2) stop 'CPA-iter written to for...'
        iw0 = 100*iq
        do i = 1, n
          iw = iw0 + i
          write (iw, '(I4,4E14.4)') itcpa, mssq(i, i, iq), w1(i, i)
        end do
      end if
!-----------------------------------------------------------------------
    end if

  end do
!================================================================= IQ ==

  cpachngl = cpachng
  cpacorrl = cpacorr

100 format (' CPA:  IQ', i3, '  ERR', f12.5, '  CORR', f13.5, '  M', 18(1x,2( &
    1p,e14.6)))
end subroutine
