!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_initldau

  private
  public :: initldau

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculates ULDAU
  !> Author: H. Ebert, V. Popescu and Ph. Mavropoulos
  !> Date: March 2003
  !> Category: KKRhost, initialization
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> fivos will add some comments later
  !>
  !> It is already later, but no comments have been added yet.
  !> But wait up, they're coming...    Munich, Feb.2004, Phivos
  !-------------------------------------------------------------------------------
  subroutine initldau(nsra, ntldau, itldau, lopt, ueff, jeff, erefldau, visp, nspin, r, drdi, z, ipan, ircut, phi, uldau)

    use :: mod_datatypes, only: dp
    use :: mod_types, only: t_inc
    use :: global_variables, only: irmd, natypd, mmaxd, ipand, lmaxd, npotd
    use :: mod_phicalc, only: phicalc
    use :: mod_rinit, only: rinit
    use :: mod_simpk, only: simpk
    use :: mod_soutk, only: soutk
    use :: mod_constants, only: pi
    implicit none

    ! Local variables
    integer :: ntldau, nspin, nsra
    real (kind=dp) :: drdi(irmd, natypd), r(irmd, natypd), visp(irmd, npotd), z(natypd)
    real (kind=dp) :: uldau(mmaxd, mmaxd, mmaxd, mmaxd, natypd)

    complex (kind=dp) :: phi(irmd, natypd)
    integer :: itldau(natypd), lopt(natypd)
    integer :: ipan(natypd), ircut(0:ipand, natypd)

    real (kind=dp) :: aa(mmaxd, mmaxd, mmaxd, mmaxd, 0:2*lmaxd), erefldau(natypd), fact(0:100), fclmb(0:2*lmaxd+1), g12, g34, jeff(natypd), rlop, rpw(irmd, 2*lmaxd+1), scl, &
      sg(irmd), sl(irmd), tmpsum, sumfclmb, tg(irmd), tl(irmd), ueff(natypd), wgtfclmb(0:2*lmaxd+1), wig3j, wint(irmd), w2(irmd)
    integer :: i1, im1, im2, im3, im4, ipan1, ir, irc1, it, kk, l1, lf, lfmax, ll, m1, m2, m3, m4
    

    fact(0) = 1.0_dp
    do i1 = 1, 100
      fact(i1) = fact(i1-1)*real(i1, kind=dp) ! factorial
    end do
    if (t_inc%i_write>0) write (1337, '(/,79("="),/,22X,A,/,79("="))') 'LDA+U:  INITIALISE Coulomb matrix U'

    ! -> Calculate test functions Phi. Phi is already normalised to
    ! int phi**2 dr =1, thus it also contains factor r.

    ! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    ! Loop over atoms which need LDA+U ( LOPT >= 0 )

    do it = 1, ntldau
      i1 = itldau(it)
      if (lopt(i1)+1==0) stop ' this atom should be LDA+U'
      call phicalc(i1, lopt(i1), visp, ipan, ircut, r, drdi, z, erefldau(i1), phi(1,i1), nspin, nsra)
    end do

    if (t_inc%i_write>0) write (1337, '(6X,43("-"),/,6X,A,/,6X,43("-"))') 'Slater integrals F^n'

    do it = 1, ntldau
      i1 = itldau(it)
      ipan1 = ipan(i1)
      irc1 = ircut(ipan1, i1)

      lfmax = 2*lopt(i1)
      do ir = 2, irmd
        rpw(ir, 1) = r(ir, i1)
        do l1 = 2, 2*lmaxd + 1
          rpw(ir, l1) = r(ir, i1)*rpw(ir, l1-1)
        end do
      end do

      ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      ! 
      ! 1.  Calculate slater integrals FCLMB
      !     (here using only the large component of sra wavefunction.
      !     Whoever wants can also do the small component.)
      !
      ! 1a. Calculate slater integrals

      rlop = real(lopt(i1), kind=dp)
      sumfclmb = 0.0_dp
      ! ----------------------------------------------------------------------
      do lf = 2, lfmax, 2
        tl(1) = 0.0_dp
        tg(1) = 0.0_dp

        ! Note on integrand:
        ! Integrals are up to IRC1 because we integrate in sphere, 
        ! without thetas.
        ! In case of cell integration, from IRCUT(1)+1 to IRC1 a convolution
        ! with thetas and gaunts is needed: 
        ! Int dr R_l(r)**2 Sum_L' Gaunt_{lm,lm,l'm'}*thetas_{l'm'}(r).
        ! But then the result is m-dependent. Here, the easy way is used!

        do ir = 2, irc1
          wint(ir) = real(conjg(phi(ir,i1))*phi(ir,i1), kind=dp)
          w2(ir) = 2.0_dp*drdi(ir, i1)*wint(ir)
          tl(ir) = w2(ir)*rpw(ir, lf)
          tg(ir) = w2(ir)/rpw(ir, lf+1)
        end do

        call soutk(tl, sl, ipan(i1), ircut(0,i1))
        call soutk(tg, sg, ipan(i1), ircut(0,i1))

        sl(1) = 0.0_dp
        do ir = 2, irc1
          sl(ir) = sl(ir)/rpw(ir, lf+1) + (sg(irc1)-sg(ir))*rpw(ir, lf)
        end do

        sg(1) = 0.0_dp

        ! See Note on integrand above.

        do ir = 2, irc1
          sg(ir) = wint(ir)*sl(ir)
        end do

        call simpk(sg, fclmb(lf), ipan1, ircut(0,i1), drdi(1,i1))

        wig3j = (-1)**nint(rlop)*(1.0_dp/sqrt(2.0_dp*rlop+1.0_dp))*cgcrac(fact, rlop, real(lf,kind=dp), rlop, 0.0_dp, 0.0_dp, 0.0_dp)

        wgtfclmb(lf) = ((2*rlop+1)/(2*rlop))*wig3j**2
        sumfclmb = sumfclmb + wgtfclmb(lf)*fclmb(lf)
      end do

      ! ----------------------------------------------------------------------
      !
      ! 1b.   Normalise slater integrals FCLMB

      scl = jeff(i1)/sumfclmb
      fclmb(0) = ueff(i1)

      do lf = 2, lfmax, 2
        fclmb(lf) = scl*fclmb(lf)
      end do

      ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

      ! ======================================================================
      ! 2.   Calculate coefficient matrix AA = a(m1,m2,m3,m4)

      call rinit(mmaxd*mmaxd*mmaxd*mmaxd*(2*lmaxd+1), aa)
      ll = lopt(i1)
      do lf = 0, lfmax, 2
        do m3 = -ll, ll
          im3 = ll + m3 + 1
          do m2 = -ll, ll
            im2 = ll + m2 + 1
            do m1 = -ll, ll
              im1 = ll + m1 + 1
              m4 = m1 - m2 + m3
              if (-ll<=m4 .and. m4<=ll) then
                im4 = ll + m4 + 1
                tmpsum = 0.0_dp

                do kk = -lf, lf
                  g12 = gauntc(fact, ll, m1, lf, kk, ll, m2)
                  g34 = gauntc(fact, ll, m3, lf, -kk, ll, m4)
                  tmpsum = tmpsum + g12*g34*(-1)**abs(kk)
                end do

                aa(im1, im2, im3, im4, lf) = tmpsum*4.0_dp*pi/(2.0_dp*real(lf,kind=dp)+1.0_dp)
              end if
            end do
          end do
        end do
      end do

      ! ======================================================================
      
      ! UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
      ! 3.  Calculate ULDAU

      call rinit(mmaxd*mmaxd*mmaxd*mmaxd, uldau(1,1,1,1,i1))

      do lf = 0, lfmax, 2
        do im4 = 1, 2*ll + 1
          do im3 = 1, 2*ll + 1
            do im2 = 1, 2*ll + 1
              do im1 = 1, 2*ll + 1
                uldau(im1, im2, im3, im4, i1) = uldau(im1, im2, im3, im4, i1) + aa(im1, im2, im3, im4, lf)*fclmb(lf)
              end do
            end do
          end do
        end do
      end do

      ! UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU

      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      if (t_inc%i_write>0) then
        write (1337, '(/,8X,A,I3,/)') 'ATOM: ', i1
        write (1337, '(12X,A,F8.4,A)') 'LDA+U reference energy :', erefldau(i1), ' Ry'
        write (1337, '(12X,A,2F8.4,A,/)') 'Ueff and Jeff = ', ueff(i1), jeff(i1), ' Ry'
        write (1337, '(12X,"Scaling factor for F^n :",F10.6,/)') scl
        write (1337, '(12X,"  n   F^n calculated   F^n scaled ")')
        
        if (scl==0) then
          ! to prevent division by zero
          ! we can do this here because scl is not use anymore in the remainder
          scl = 1.0
        end if

        do lf = 2, lfmax, 2
          write (1337, '(12X,I3,2(2X,F12.8," Ry"))') lf, fclmb(lf)/scl, fclmb(lf)
        end do
        if (it<ntldau) write (1337, '(8X,58("~"))')
      end if
      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    end do
    !AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

    if (t_inc%i_write>0) then
      write (1337, '(/,6X,60("-"),/,6X,A,/,6X,60("-"))') 'Coulomb matrix U(m1,m1,m3,m3)'
      do it = 1, ntldau
        i1 = itldau(it)
        ll = lopt(i1)
        ll = min(3, ll)
        write (1337, '(/,8X,"ATOM :",I3,/)') i1

        ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        do im1 = 1, 2*ll + 1
          write (1337, 100)(uldau(im1,im1,im3,im3,i1), im3=1, 2*ll+1)
        end do
        if (it<ntldau) write (1337, '(8X,58("~"))')
        ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

      end do
      write (1337, '(/,6X,60("-"),/)')
    end if
100 format (6x, 7f10.6)
  end subroutine initldau


  !-------------------------------------------------------------------------------
  !> Summary: CLEBSCH GORDAN COEFFICIENTS FOR ARBITRARY QUANTUM NUMBERS
  !> Author: 
  !> Category: KKRhost, special-functions
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> CLEBSCH GORDAN COEFFICIENTS FOR ARBITRARY QUANTUM NUMBERS  J1,J2 ...
  !> ACCORDING TO THE FORMULA OF   RACAH
  !> SEE: M.E.ROSE ELEMENTARY THEORY OF ANGULAR MOMENTUM
  !>      EQUATION (3.19)
  !>      EDMONDS EQ. (3.6.11) PAGE 45
  !-------------------------------------------------------------------------------
  function cgcrac(fact, j1, j2, j3, m1, m2, m3)

    use :: mod_datatypes, only: dp
    implicit none

    real (kind=dp) :: j1, j2, j3, m1, m2, m3
    real (kind=dp) :: cgcrac
    real (kind=dp) :: fact(0:100)
    integer :: j, n, n1, n2, n3, n4, n5, nbot, ntop
    real (kind=dp) :: rfact
    real (kind=dp) :: s, tmpsum, vf, x, y

    ! INLINE FUNCTION FACTORIAL FOR REAL ARGUMENT
    rfact(x) = fact(nint(x))

    cgcrac = 0.0_dp
    if (abs(m3-(m1+m2))>1.0d-6) return
    if (abs(j1-j2)>j3) return
    if ((j1+j2)<j3) return
    if (abs(m1)>(j1+1.0d-6)) return
    if (abs(m2)>(j2+1.0d-6)) return
    if (abs(m3)>(j3+1.0d-6)) return

    do j = abs(nint(2*(j1-j2))), nint(2*(j1+j2)), 2
      if (j==nint(2*j3)) go to 100
    end do
    return

100 continue

    x = (2.0_dp*j3+1.0_dp)*rfact(j1+j2-j3)*rfact(j1-j2+j3)*rfact(-j1+j2+j3)*rfact(j1+m1)*rfact(j1-m1)*rfact(j2+m2)*rfact(j2-m2)*rfact(j3+m3)*rfact(j3-m3)

    y = rfact(j1+j2+j3+1)

    vf = sqrt(x/y)

    n1 = nint(j1+j2-j3)
    n2 = nint(j1-m1)
    n3 = nint(j2+m2)
    n4 = nint(j3-j2+m1)
    n5 = nint(j3-j1-m2)
    ntop = min(n1, n2, n3)
    nbot = max(0, -n4, -n5)

    n = nbot + 1
    if (n==(2*(n/2))) then
      s = +1.0_dp
    else
      s = -1.0_dp
    end if
    tmpsum = 0.0_dp

    do n = nbot, ntop
      s = -s
      y = fact(n)*fact(n1-n)*fact(n2-n)*fact(n3-n)*fact(n4+n)*fact(n5+n)
      tmpsum = tmpsum + (s/y)
    end do
    cgcrac = vf*tmpsum

  end function cgcrac


  !-------------------------------------------------------------------------------
  !> Summary: GAUNT COEFFICIENTS for complex spherical harmonics
  !> Author: H. Ebert,
  !> Date: 26/01/95 
  !> Category: KKRhost, special-functions
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> GAUNT COEFFICIENTS for complex spherical harmonics  Y[l,m]
  !>
  !>        G = INT dr^  Y[l1,m1]* Y[l2,m2] Y[l3,m3]
  !>
  !> see: M.E.ROSE ELEMENTARY THEORY OF ANGULAR MOMENTUM  Eq. (4.34)
  !-------------------------------------------------------------------------------
  function gauntc(fact, l1, m1, l2, m2, l3, m3)

    use :: mod_datatypes, only: dp
    use :: mod_constants, only: pi
    implicit none

    integer :: l1, l2, l3, m1, m2, m3
    real (kind=dp) :: fact(0:100)
    real (kind=dp) :: gauntc
    real (kind=dp) :: g

    if ((l1<0) .or. (l2<0) .or. (l3<0)) then
      g = 0.0_dp
    else
      g = (real(2*l2+1,kind=dp)*real(2*l3+1,kind=dp)/(4.0_dp*pi*real(2*l1+1,kind=dp)))**0.5_dp*cgcrac(fact, real(l3,kind=dp), real(l2,kind=dp), real(l1,kind=dp), real(m3,kind=dp), &
        real(m2,kind=dp), real(m1,kind=dp))*cgcrac(fact, real(l3,kind=dp), real(l2,kind=dp), real(l1,kind=dp), 0.0_dp, 0.0_dp, 0.0_dp)
    end if
    gauntc = g

  end function gauntc

end module mod_initldau
