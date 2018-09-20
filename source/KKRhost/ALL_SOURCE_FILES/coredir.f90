module mod_coredir

contains

  !-------------------------------------------------------------------------------
  !> Summary: Solver for radial Dirac equation of core wavefunctions
  !> Author: H. Ebert
  !> Date: 1989
  !> Category: KKRhost, core-electrons
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> ROUTINE TO SOLVE RADIAL SPIN-POLARISED DIRAC EQUATIONS
  !> FOR THE CORE WAVE FUNCTIONS                           
  !>                                                       
  !> SIMILAR TO LOUCKS' METHOD TO SOLVE THE COUPLED SET OF 
  !> DIFFERENTIAL EQUATIONS                                
  !>                                  HE JAN. 1989         
  !>                                                       
  !> ADOPTED FOR FINITE NUCLEUS       MB MAR. 1995
  !-------------------------------------------------------------------------------
  subroutine coredir(it, c, e, l, mj, way, vv, bb, rc, drdic, dovrc, nmatch, nzero, gc, fc, d_p, dq, wp, wq, pow, qow, piw, qiw, cgd, cgmd, cgo, nrc, z, nucleus)

    use :: mod_types, only: t_inc
    use :: mod_datatypes, only: dp
    implicit none

    ! PARAMETER definitions
    real (kind=dp), parameter :: eps = 1.0d-12
    integer :: mpsmax, npemax, invmax
    parameter (mpsmax=20, npemax=20, invmax=3)
    real (kind=dp) :: tol
    parameter (tol=1.0d-9)
    integer :: itmax
    parameter (itmax=50)

    ! Dummy arguments
    real (kind=dp) :: c, cgo, mj
    real (kind=dp) :: e
    integer :: it, l, nmatch, nrc, nucleus, nzero, z
    character (len=3) :: way
    real (kind=dp) :: bb(nrc), cgd(2), cgmd(2), dovrc(nrc), d_p(2, 2, nrc), dq(2, 2, nrc), drdic(nrc), fc(2, 2, nrc), gc(2, 2, nrc), piw(2, 2), pow(2, 2), qiw(2, 2), qow(2, 2), &
      rc(nrc), vv(nrc), wp(2, 2, nrc), wq(2, 2, nrc)

    ! Local variables
    real (kind=dp) :: a11, a12, a21, a22, aa11, aa12, aa21, aa22, alpha, bb1, bb2, beta, bova, bpp, bqq, diffa, diffb, dmue, emvpp, emvqq, w1, w2, w3, w4, w5, w6, w7
    real (kind=dp) :: bc(0:npemax), cg1, cg2, cg4, cg5, cg8, csqr, det, dvc, gam(2), gpm, h24, kap(2), pc(2, 2, 0:mpsmax), pnew(2, 2), pold(2, 2), qc(2, 2, 0:mpsmax), qnew(2, 2), &
      qold(2, 2), rpwgpm, rr, tz, vc(0:npemax)
    integer :: i, iv, j, jcorr, k, kap1, kap2, m, mps, n, nn, nsol
    integer :: int, nint
    save :: a11, a12, a21, a22, aa11, aa12, aa21, aa22, alpha, bb1, bb2, bc, beta, bova, bpp, bqq, cg1, cg2, cg4, cg5, cg8, csqr, det, diffa, diffb, dmue, dvc, emvpp, emvqq, gam, &
      gpm, h24, i, iv, j, jcorr, k, kap, kap1, kap2, m, mps, n, nn, nsol, pc, pnew, pold, qc, qnew, qold, rpwgpm, rr, tz, vc, w1, w2, w3, w4, w5, w6, w7

    ! MB
    ! real (kind=dp) CM(INVMAX,INVMAX),CMI(INVMAX,INVMAX)
    ! MB


    h24 = 1.0d0/24.0d0
    dvc = c
    csqr = dvc*dvc

    ! EXPANSION COEFFICIENTS FOR THE POTENTIAL AND B-FIELD
    ! MB
    if (nucleus==0) then
      tz = real(nint(-vv(1)*rc(1)), kind=dp)
      vc(0) = vv(1) - (-tz)/rc(1)
    else
      tz = 2.0d0*real(z, kind=dp)
      vc(0) = vv(1)
    end if
    do i = 1, 2
      do j = 1, 2
        do k = 1, npemax
          pc(i, j, k) = 0.0d0
          qc(i, j, k) = 0.0d0
        end do
      end do
    end do
    ! MB

    bc(0) = bb(1)


    ! CALCULATE G-COEFFICIENTS OF B-FIELD

    kap1 = -l - 1
    kap2 = +l

    cg1 = -mj/(kap1+0.5d0)
    cg5 = -mj/(-kap1+0.5d0)
    cgd(1) = cg1
    cgmd(1) = cg5
    kap(1) = real(kap1, kind=dp)
    ! MB
    if (nucleus==0) then
      gam(1) = sqrt(kap(1)**2-(tz/dvc)**2)
    else
      gam(1) = abs(kap(1))
    end if
    ! MB
    if (abs(mj)>l) then
      cg2 = 0.0d0
      cg4 = 0.0d0
      cg8 = 0.0d0
      nsol = 1
      cgd(2) = 0.0d0
      cgo = 0.0d0
      cgmd(2) = 0.0d0
      gam(2) = 0.0d0
      kap(2) = 0.0d0
    else
      cg2 = -sqrt(1.0d0-(mj/(kap1+0.5d0))**2)
      cg4 = -mj/(kap2+0.5d0)
      cg8 = -mj/(-kap2+0.5d0)
      nsol = 2
      cgd(2) = cg4
      cgo = cg2
      cgmd(2) = cg8
      kap(2) = real(kap2, kind=dp)
      ! MBA
      if (nucleus==0) then
        gam(2) = sqrt(kap(2)**2-(tz/dvc)**2)
      else
        gam(2) = abs(kap(2))
      end if
      ! MBE
    end if


    if (way=='INW') then

      ! #####################################################################
      ! #####################################################################
      ! #####################################################################

      ! INWARD INTEGRATION

      dmue = sqrt(-e-e*e/csqr)
      bova = -dmue/(1.0d0+e/csqr)

      do n = (nzero-3), nzero

        rr = rc(n)

        do j = 1, nsol
          i = 3 - j
          wp(j, j, n) = exp(-dmue*rr)
          d_p(j, j, n) = -dmue*drdic(n)*wp(j, j, n)
          wq(j, j, n) = bova*wp(j, j, n)
          dq(j, j, n) = bova*d_p(j, j, n)

          wp(i, j, n) = 0.0d0
          wq(i, j, n) = 0.0d0
          d_p(i, j, n) = 0.0d0
          dq(i, j, n) = 0.0d0
        end do
      end do

      ! =============================================================== N ====

      ! CALCULATE ALL NEXT POINTS BY PRE/CORR(ADAMS-MOULTON-BASHFORTH)

      do nn = 1, (nzero-3-nmatch)
        n = nzero - 3 - nn

        ! EVALUATE PREDICTOR

        do j = 1, nsol
          do i = 1, nsol
            pnew(i, j) = wp(i, j, n+1) - h24*(55.0d0*d_p(i,j,n+1)-59.0d0*d_p(i,j,n+2)+37.0d0*d_p(i,j,n+3)-9.0d0*d_p(i,j,n+4))
            qnew(i, j) = wq(i, j, n+1) - h24*(55.0d0*dq(i,j,n+1)-59.0d0*dq(i,j,n+2)+37.0d0*dq(i,j,n+3)-9.0d0*dq(i,j,n+4))
          end do
        end do

        emvqq = (e-vv(n)+csqr)*drdic(n)/csqr
        emvpp = -(e-vv(n))*drdic(n)
        bqq = bb(n)*drdic(n)/csqr
        bpp = bb(n)*drdic(n)

        ! EVALUATE CORRECTOR

        do jcorr = 1, itmax

          do j = 1, nsol
            do i = 1, nsol
              pold(i, j) = pnew(i, j)
              qold(i, j) = qnew(i, j)
              d_p(i, j, n) = -kap(i)*pnew(i, j)*dovrc(n) + (emvqq+bqq*cgmd(i))*qnew(i, j)
              dq(i, j, n) = kap(i)*qnew(i, j)*dovrc(n) + (emvpp+bpp*cgd(i))*pnew(i, j) + bpp*cgo*pnew(3-i, j)

              pnew(i, j) = wp(i, j, n+1) - h24*(9.0d0*d_p(i,j,n)+19.0d0*d_p(i,j,n+1)-5.0d0*d_p(i,j,n+2)+d_p(i,j,n+3))
              qnew(i, j) = wq(i, j, n+1) - h24*(9.0d0*dq(i,j,n)+19.0d0*dq(i,j,n+1)-5.0d0*dq(i,j,n+2)+dq(i,j,n+3))
            end do
          end do

          do j = 1, nsol
            do i = 1, nsol
              diffa = pold(i, j) - pnew(i, j)
              if (abs(diffa)>(tol*abs(pnew(i,j)))) go to 100
              if (abs(diffa)>(tol*abs(pnew(i,j)))) go to 100

              diffb = qold(i, j) - qnew(i, j)
              if (abs(diffb)>(tol*abs(qnew(i,j)))) go to 100
              if (abs(diffb)>(tol*abs(qnew(i,j)))) go to 100
            end do
          end do
          go to 110

100     end do
        if (t_inc%i_write>0) write (1337, 140) kap1, n, rc(n), diffa, diffb, it, l, int(2*mj), ' IN'

        ! SORRY NOT CONVERGED IN  ITMAX  ITERATIONS

110     continue
        do j = 1, nsol
          do i = 1, nsol
            wp(i, j, n) = pnew(i, j)
            wq(i, j, n) = qnew(i, j)
            d_p(i, j, n) = -kap(i)*pnew(i, j)*dovrc(n) + (emvqq+bqq*cgmd(i))*qnew(i, j)
            dq(i, j, n) = kap(i)*qnew(i, j)*dovrc(n) + (emvpp+bpp*cgd(i))*pnew(i, j) + bpp*cgo*pnew(3-i, j)
          end do
        end do

      end do
      ! =============================================================== N ====


      ! NOW TRANSFORM TO THE PROPER WAVEFUNCTIONS

      do n = nmatch, nzero
        do j = 1, nsol
          do i = 1, nsol
            gc(i, j, n) = wp(i, j, n)/rc(n)
            fc(i, j, n) = wq(i, j, n)/(rc(n)*dvc)
          end do
        end do
      end do

      do j = 1, nsol
        do i = 1, nsol
          piw(i, j) = wp(i, j, nmatch)
          qiw(i, j) = wq(i, j, nmatch)
        end do
      end do
      go to 150
    end if

    ! #####################################################################
    ! #####################################################################
    ! #####################################################################

    ! OUTWARD INTEGRATION


    ! DETERMINE HIGHER EXPANSION COEFFICIENTS FOR THE WAVE FUNCTIONS

    mps = 20

    aa12 = -tz/csqr
    aa21 = tz
    emvqq = (e-vc(0)+csqr)/csqr
    emvpp = -e + vc(0)
    bqq = bc(0)/csqr
    ! MBA
    if (nucleus==0) then
      ! MBE
      do j = 1, nsol
        i = 3 - j
        pc(j, j, 0) = sqrt(abs(kap(j))-gam(j))
        qc(j, j, 0) = (kap(j)+gam(j))*(csqr/tz)*pc(j, j, 0)
        pc(i, j, 0) = 0.0d0
        qc(i, j, 0) = 0.0d0
      end do

      do j = 1, nsol

        do m = 1, mps
          do i = 1, nsol
            bb1 = (emvqq+bqq*cgmd(i))*qc(i, j, m-1)
            bb2 = (emvpp+bc(0)*cgd(i))*pc(i, j, m-1) + bc(0)*cgo*pc(3-i, j, m-1)
            aa11 = gam(j) + m + kap(i)
            aa22 = gam(j) + m - kap(i)
            det = aa11*aa22 - aa12*aa21
            pc(i, j, m) = (bb1*aa22-aa12*bb2)/det
            qc(i, j, m) = (aa11*bb2-bb1*aa21)/det
          end do
        end do

      end do
      ! MBA
    else
      ! EXPANSION ADAPTED FOR POTENTIALS WITH FINITE NUCLEUS
      ! EXPANSION OF POTENTIAL actually UP TO zeroth ORDER

      ! DO IV=1,INVMAX
      ! DO N=1,INVMAX
      ! CM(N,IV)=RC(N)**(IV-1)
      ! ENDDO
      ! ENDDO

      ! CALL RINVGJ(CMI,CM,INVMAX,INVMAX)
      do iv = 1, invmax
        vc(iv-1) = 0.0d0
        ! DO N=1,INVMAX
        ! VC(IV-1)=VC(IV-1)+CMI(IV,N)*VV(N)
        ! ENDDO
      end do
      do j = 1, nsol
        i = 3 - j
        if (kap(j)>0) then
          ! ARBITRARY STARTING VALUES
          alpha = 0.0d0
          beta = 0.174d0
        else
          beta = 0.0d0
          alpha = 0.174d0
        end if
        pc(j, j, 0) = alpha
        qc(j, j, 0) = beta
        pc(i, j, 0) = 0.0d0
        qc(i, j, 0) = 0.0d0
      end do

      w4 = bc(0)*cgo
      w2 = vc(1)/csqr
      w5 = vc(1)
      w6 = vc(2)/csqr
      w7 = vc(2)
      do j = 1, nsol
        do i = 1, nsol
          w1 = emvqq + bqq*cgmd(i)
          w3 = -emvpp + bc(0)*cgd(i)
          a11 = gam(j) + kap(i) + 1d0
          a12 = gam(j) - kap(i) + 1d0
          if (abs(a11)>eps) pc(i, j, 1) = w1/a11*qc(i, j, 0)
          if (abs(a12)>eps) qc(i, j, 1) = (-w3*pc(i,j,0)+w4*pc(3-i,j,0))/a12

        end do
      end do
      do j = 1, nsol
        do i = 1, nsol
          w1 = emvqq + bqq*cgmd(i)
          w3 = -emvpp + bc(0)*cgd(i)
          a11 = gam(j) + kap(i) + 2d0
          a12 = gam(j) - kap(i) + 2d0
          if (abs(a11)>eps) pc(i, j, 2) = (w1*qc(i,j,1)-w2*qc(i,j,0))/a11
          if (abs(a12)>eps) qc(i, j, 2) = (-w3*pc(i,j,1)+w4*pc(3-i,j,1)+w5*pc(i,j,0))/a12
        end do
      end do
      do j = 1, nsol
        do m = 3, mps
          do i = 1, nsol
            w1 = emvqq + bqq*cgmd(i)
            w3 = -emvpp + bc(0)*cgd(i)
            a21 = gam(j) + kap(i) + real(m, kind=dp)
            a22 = gam(j) - kap(i) + real(m, kind=dp)
            if (abs(a21)>eps) pc(i, j, m) = (w1*qc(i,j,m-1)-w2*qc(i,j,m-2)-w6*qc(i,j,m-3))/a21
            if (abs(a22)>eps) qc(i, j, m) = (-w3*pc(i,j,m-1)+w4*pc(3-i,j,m-1)+w5*pc(i,j,m-2)+w7*pc(i,j,m-3))/a22
          end do
        end do
      end do
    end if
    ! MBE

    ! PERFORM SUMMATION OVER WAVE FUNCTION - EXPANSION COEFFICIENTS
    ! FOR THE FIRST 4 R - MESH - POINTS

    do n = 1, 4
      rr = rc(n)

      do j = 1, nsol
        rpwgpm = rr**gam(j)

        do i = 1, nsol
          wp(i, j, n) = pc(i, j, 0)*rpwgpm
          wq(i, j, n) = qc(i, j, 0)*rpwgpm
          ! print*,WP(i,j,N), WQ(I,J,N),0,i,j
          d_p(i, j, n) = pc(i, j, 0)*rpwgpm*gam(j)*dovrc(n)
          dq(i, j, n) = qc(i, j, 0)*rpwgpm*gam(j)*dovrc(n)
        end do

        do m = 1, mps
          rpwgpm = rpwgpm*rr
          gpm = gam(j) + m

          do i = 1, nsol
            wp(i, j, n) = wp(i, j, n) + pc(i, j, m)*rpwgpm
            wq(i, j, n) = wq(i, j, n) + qc(i, j, m)*rpwgpm
            ! print*,WP(i,j,N),WQ(I,J,N) ,m,i,j
            d_p(i, j, n) = d_p(i, j, n) + pc(i, j, m)*rpwgpm*gpm*dovrc(n)
            dq(i, j, n) = dq(i, j, n) + qc(i, j, m)*rpwgpm*gpm*dovrc(n)
          end do

        end do
        ! if((nsol.eq.2).and.(N.gt.1))stop
      end do
    end do

    ! =============================================================== N ====
    ! CALCULATE ALL NEXT POINTS BY PRE/CORR(ADAMS-MOULTON-BASHFORTH)

    do n = 5, nmatch

      ! EVALUATE PREDICTOR

      do j = 1, nsol
        do i = 1, nsol
          pnew(i, j) = wp(i, j, n-1) + h24*(55.0d0*d_p(i,j,n-1)-59.0d0*d_p(i,j,n-2)+37.0d0*d_p(i,j,n-3)-9.0d0*d_p(i,j,n-4))
          qnew(i, j) = wq(i, j, n-1) + h24*(55.0d0*dq(i,j,n-1)-59.0d0*dq(i,j,n-2)+37.0d0*dq(i,j,n-3)-9.0d0*dq(i,j,n-4))
        end do
      end do

      emvqq = (e-vv(n)+csqr)*drdic(n)/csqr
      emvpp = -(e-vv(n))*drdic(n)
      bqq = bb(n)*drdic(n)/csqr
      bpp = bb(n)*drdic(n)

      ! EVALUATE CORRECTOR


      do jcorr = 1, itmax

        do j = 1, nsol
          do i = 1, nsol
            pold(i, j) = pnew(i, j)
            qold(i, j) = qnew(i, j)
            d_p(i, j, n) = -kap(i)*pnew(i, j)*dovrc(n) + (emvqq+bqq*cgmd(i))*qnew(i, j)
            dq(i, j, n) = kap(i)*qnew(i, j)*dovrc(n) + (emvpp+bpp*cgd(i))*pnew(i, j) + bpp*cgo*pnew(3-i, j)

            pnew(i, j) = wp(i, j, n-1) + h24*(9.0d0*d_p(i,j,n)+19.0d0*d_p(i,j,n-1)-5.0d0*d_p(i,j,n-2)+d_p(i,j,n-3))
            qnew(i, j) = wq(i, j, n-1) + h24*(9.0d0*dq(i,j,n)+19.0d0*dq(i,j,n-1)-5.0d0*dq(i,j,n-2)+dq(i,j,n-3))
          end do
        end do

        do j = 1, nsol
          do i = 1, nsol
            diffa = pold(i, j) - pnew(i, j)
            if (abs(diffa)>(tol*abs(pnew(i,j)))) go to 120
            if (abs(diffa)>(tol*abs(pnew(i,j)))) go to 120

            diffb = qold(i, j) - qnew(i, j)
            if (abs(diffb)>(tol*abs(qnew(i,j)))) go to 120
            if (abs(diffb)>(tol*abs(qnew(i,j)))) go to 120
          end do
        end do
        go to 130

120   end do
      if (t_inc%i_write>0) write (1337, 140) kap1, n, rc(n), diffa, diffb, it, l, int(2*mj), 'OUT'

      ! SORRY NOT CONVERGED IN  ITMAX  ITERATIONS

130   continue
      do j = 1, nsol
        do i = 1, nsol
          wp(i, j, n) = pnew(i, j)
          wq(i, j, n) = qnew(i, j)
          d_p(i, j, n) = -kap(i)*pnew(i, j)*dovrc(n) + (emvqq+bqq*cgmd(i))*qnew(i, j)
          dq(i, j, n) = kap(i)*qnew(i, j)*dovrc(n) + (emvpp+bpp*cgd(i))*pnew(i, j) + bpp*cgo*pnew(3-i, j)
        end do
      end do

    end do
    ! =============================================================== N ====

    ! NOW TRANSFORM TO THE PROPER WAVEFUNCTIONS

    do n = 1, nmatch
      do j = 1, nsol
        do i = 1, nsol
          gc(i, j, n) = wp(i, j, n)/rc(n)


          fc(i, j, n) = wq(i, j, n)/(rc(n)*dvc)
        end do
      end do
    end do

    do j = 1, nsol
      do i = 1, nsol
        pow(i, j) = wp(i, j, nmatch)
        qow(i, j) = wq(i, j, nmatch)
      end do
    end do

    return

140 format (' P/C NOT CONV. IN <DIRAC> ', 2i4, 2x, f10.7, 2x, 2e12.4, 3i2, '/2 ', a3)
150 continue

  end subroutine coredir

end module mod_coredir
