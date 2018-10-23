!------------------------------------------------------------------------------------
!> Summary: Calculation of lattice sums for \(l \leq 2l_{pot}\)
!> Author: B. Drittler
!> Calculation of lattice sums for \(l \leq 2l_{pot}\)
!> \begin{equation}
!> \sum_{rm} \frac{Y_{lm}\left(q(i)-q(j)-rm\right)}{\left| q(i)-q(j)-rm\right|^{l+1}}
!> \end{equation}
!> in the case of \(i = j\), \(rm = 0\) is omitted. Rhe ewald method is used to 
!> perform the lattice summations the splitting parameter \(\lambda\) is set equal 
!> \(\frac{\sqrt(\pi)}{a_{lat}}\) (\(a_{lat}\) is the lattice constant).
!> If the contribution of the last shell of the direct and the reciprocal lattice 
!> is greater than `1.0e-8` a message is written 
!------------------------------------------------------------------------------------
!> @note V. Popescu May 2004: Dimension of arrays GN,RM changed from `(4,*)` to 
!> `(3,*)`, the 4th one not being used (see also `lattice3d`).
!> @endnote
!------------------------------------------------------------------------------------
module mod_strmat

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculation of lattice sums for \(l \leq 2l_{pot}\)
  !> Author: B. Drittler
  !> Category: electrostatics, geomertry, k-points, KKRhost
  !> Deprecated: False 
  !> Calculation of lattice sums for \(l \leq 2l_{pot}\)
  !> \begin{equation}
  !> \sum_{rm} \frac{Y_{lm}\left(q(i)-q(j)-rm\right)}{\left| q(i)-q(j)-rm\right|^{l+1}}
  !> \end{equation}
  !> in the case of \(i = j\), \(rm = 0\) is omitted. Rhe ewald method is used to 
  !> perform the lattice summations the splitting parameter \(\lambda\) is set equal 
  !> \(\frac{\sqrt(\pi)}{a_{lat}}\) (\(a_{lat}\) is the lattice constant).
  !> If the contribution of the last shell of the direct and the reciprocal lattice 
  !> is greater than `1.0e-8` a message is written 
  !-------------------------------------------------------------------------------
  !> @note V. Popescu May 2004: Dimension of arrays GN,RM changed from `(4,*)` to 
  !> `(3,*)`, the 4th one not being used (see also `lattice3d`).
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine strmat(alat,lpot,naez,ngmax,nrmax,nsg,nsr,nshlg,nshlr,gn,rm,qi0,smat,  &
    vol,iprint,lassld,lmxspd,naezd)

#ifdef CPP_HYBRID
#define CPP_OMPSTUFF
#endif
#ifdef CPP_OMP
#define CPP_OMPSTUFF
#endif
#ifdef CPP_OMPSTUFF
    use :: omp_lib
#endif

    use :: mod_constants
    use :: mod_datatypes, only: dp
    use :: mod_ymy
    use :: mod_gamfc

    implicit none
    ! ..
    ! .. Parameters ..
    real (kind=dp) :: bound
    parameter (bound=1d-8)
    ! ..
    ! .. Scalar arguments ..
    real (kind=dp) :: alat, vol
    integer :: iprint, lpot, naez, ngmax, nrmax, nshlg, nshlr
    integer :: lassld, lmxspd, naezd
    ! ..
    ! .. Array arguments ..
    real (kind=dp) :: gn(3, *), qi0(3, *), rm(3, *), smat(lmxspd, naezd, *)
    integer :: nsg(*), nsr(*)
    ! ..
    ! .. Local scalars ..
    complex (kind=dp) :: bfac
    real (kind=dp) :: alpha, beta, dq1, dq2, dq3, dqdotg, expbsq, fpi, g1, g2, g3, ga
    real (kind=dp) :: lamda, r, r1, r2, r3, rfac, s
    integer :: i, i1, i2, it, l, lm, lmx, lmxsp, lfmt, m, nge, ngs, nre, nrs, nstart
    character (len=80) :: fmt
    ! ..
    ! .. Local arrays ..
    complex (kind=dp) :: stest(lmxspd)
    real (kind=dp) :: g(0:lassld), ylm(lmxspd), qi(3, naezd)
    ! ..................................................................

    lmx = 2*lpot
    lmxsp = (lmx+1)*(lmx+1)
    fpi = 4.0d0*pi

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, '(5X,2A,/)') '< STRMAT > : ', 'calculating lattice sums'
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT


    ! --> choose proper splitting parameter

    lamda = sqrt(pi)/alat

    ! --> loop over atoms per unit cell -- scale basis atoms with alat

    do i2 = 1, naez
      do i1 = 1, 3
        qi(i1, i2) = qi0(i1, i2)*alat
      end do
    end do

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    if (iprint>2) then
      write (1337, 110)
      do i2 = 1, naez
        write (1337, 120) i2, (qi0(i,i2), i=1, 3)
      end do
      write (1337, 130)
    end if
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

    ! **********************************************************************
#ifdef CPP_OMPSTUFF
    ! $omp parallel do default(shared) private(DQ1, DQ2, DQ3, STEST) &
    ! $omp& private(LM, NSTART, IT, NRS, NGS, NRE, NGE, I, R1, R2) &
    ! $omp& private(R3 , R, YLM, ALPHA, G, L, RFAC, M, G1, G2) &
    ! $omp& private(G3, GA, BETA, EXPBSQ, DQDOTG, BFAC, S, I1, I2)
#endif
    do i1 = 1, naez
      do i2 = 1, naez
        ! ======================================================================
        dq1 = qi(1, i1) - qi(1, i2)
        dq2 = qi(2, i1) - qi(2, i2)
        dq3 = qi(3, i1) - qi(3, i2)

        stest(1) = -sqrt(fpi)/vol/(4d0*lamda*lamda)
        do lm = 2, lmxsp
          stest(lm) = 0.0d0
        end do

        ! --> exclude the origine and add correction if i1.eq.i2

        if (i1==i2) then
          stest(1) = stest(1) - lamda/pi
          nstart = 2
        else
          nstart = 1
        end if
        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! --> loop first over n-1 shells of real and reciprocal lattice - then
        ! add the contribution of the last shells to see convergence

        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        do it = 1, 2
          if (it==1) then
            nrs = nstart
            ngs = 2
            nre = nrmax - nsr(nshlr)
            nge = ngmax - nsg(nshlg)
          else
            nrs = nre + 1
            ngs = nge + 1
            nre = nrmax
            nge = ngmax
          end if

          ! --> sum over real lattice

          ! ---------------------------------------------------------------------
          do i = nrs, nre
            r1 = dq1 - rm(1, i)
            r2 = dq2 - rm(2, i)
            r3 = dq3 - rm(3, i)

            call ymy(r1, r2, r3, r, ylm, lmx)
            alpha = lamda*r
            call gamfc(alpha, g, lmx, r)

            ! IF (R>0.3D0 .or. .not. OPT('VIRATOMS') ) THEN         ! added by bauer VIRTAOM
            do l = 0, lmx
              rfac = g(l)/sqrt(pi)
              do m = -l, l
                lm = l*(l+1) + m + 1
                stest(lm) = stest(lm) + ylm(lm)*rfac
              end do
            end do
            ! ELSE
            ! write(*,*) 'omitting ',R
            ! endif !(R>0.3) THEN         ! added by bauer VIRTAOM

          end do
          ! ---------------------------------------------------------------------

          ! --> sum over reciprocal lattice

          ! ---------------------------------------------------------------------
          do i = ngs, nge
            g1 = gn(1, i)
            g2 = gn(2, i)
            g3 = gn(3, i)

            call ymy(g1, g2, g3, ga, ylm, lmx)
            beta = ga/lamda
            dqdotg = dq1*g1 + dq2*g2 + dq3*g3

            if (beta>50.0_dp) then
              bfac = 0.0_dp
            else
              expbsq = exp(beta*beta/4.0d0)
              bfac = fpi*exp(ci*dqdotg)/(ga*ga*expbsq*vol)
            end if

            do l = 0, lmx
              do m = -l, l
                lm = l*(l+1) + m + 1
                stest(lm) = stest(lm) + ylm(lm)*bfac
              end do
              bfac = bfac*ga/ci/real(2*l+1, kind=dp)
            end do
          end do
          ! ---------------------------------------------------------------------
          if (it==1) then
            do lm = 1, lmxsp
              if (abs(aimag(stest(lm)))>bound) then
                write (6, *) ' ERROR: Imaginary contribution', ' to REAL lattice sum'
                stop
              end if
              smat(lm, i1, i2) = real(stest(lm), kind=dp)
              stest(lm) = 0.0d0
            end do
          else

            ! --> test convergence

#ifdef CPP_OMPSTUFF
            ! $omp critical
#endif
            do lm = 1, lmxsp
              s = real(stest(lm), kind=dp)
              smat(lm, i1, i2) = smat(lm, i1, i2) + s
              if (abs(s)>bound) write (6, fmt=100) i1, i2, lm, abs(s)
            end do
#ifdef CPP_OMPSTUFF
            ! $omp end critical
#endif
          end if
          ! ---------------------------------------------------------------------
        end do
        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      end do
    end do
#ifdef CPP_OMPSTUFF
    ! $omp end parallel do
#endif
    ! **********************************************************************

    if (iprint<2) return

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, 140) min(6, naez)
    fmt = ' '
    lfmt = 0
    do i = 1, min(6, naez)
      fmt = fmt(1:lfmt) // '------------'
      lfmt = lfmt + 12
    end do
    write (1337, '(7X,A)') fmt(1:lfmt)
    do i1 = 1, min(6, naez)
      write (1337, 150)(smat(1,i1,i2), i2=1, min(6,naez))
    end do
    write (1337, '(7X,A,/)') fmt
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

100 format (5x, 'WARNING : Convergence of SMAT(', i2, ',', i2, ') ', ' for LMXSP =', i3, ' is ', 1p, d9.2, ' > 1D-8', /, 15x, 'You should use more lattice vectors (RMAX/GMAX)')
110 format (12x, 47('-'), /, 16x, '      Positions of atomic sites', /, 16x, '    in CARTESIAN coordinates (a.u.)', /, 12x, 47('-'), /, 15x, &
      'IQ       x           y           z     ', /, 12x, 47('-'))
120 format (13x, i5, 3f12.6)
130 format (12x, 47('-'), /)
140 format (8x, 'Lattice sum (LMXSP = 1) up to NAEZ =', i2)
150 format (7x, 6(d12.4))
  end subroutine strmat

end module mod_strmat
