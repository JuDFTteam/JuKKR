!------------------------------------------------------------------------------------
!> Summary: Calculates the overlap integral of test function PHI with regular or irregular wavefunction 
!> Author: Ph. Mavropoulos
!> Calculates the overlap integral of test function PHI with regular or irregular
!> wavefunction (PZ, QZ: spherical part of the large component, 
!> PQNS: nonspherical part of the large component) for LDA+U. 
!> The overlap is given out in array `RESULT`
!>                                                                   
!> In the inner region the non-spherical wavefunctions are approximated as: 
!>                                                                   
!>  * The regular one (ir < irmin = irws-irns) :
!>  \begin{equation}
!>   pns\left(ir,lm1,lm2\right) = pz\left(ir,l1\right) ar\left(lm1,lm2\right)
!>  \end{equation}
!> where \(pz\) is the regular wavefunction of the spherically symmetric
!> part of the potential and \(ar\) the alpha matrix (see sub. regns)
!>                                                                   
!>  * The irregular one (ir < irmin) :
!>  \begin{equation}
!>   qns\left(ir,lm1,lm2\right) = pz\left(ir,l1\right) cr\left(lm1,lm2\right) + qz\left(ir,l1\right) * dr\left(lm1,lm2\right) 
!>  \end{equation}
!> where \(pz\) is the regular and \(qz\) is the irregular wavefunction of the 
!> spherically symmetric part of the potential and \(cr\), \(dr\) the matrices 
!> calculated at the point \(irmin\) (see sub. `irwns`)
!>                                                                   
!> The integrand is convoluted with the shape functions:
!> \begin{equation}
!> Result = \int \bar{phi_l} \sum_{L''L'''} R_L''L' Theta_L''' GAUNT_LL''L'''   
!> \end{equation}
!> Finally, the result should then be transformed to complex spherical harmonics basis! 
!------------------------------------------------------------------------------------
!> @note Here only large component is used, but was PHI normalised according to 
!> large + small component? (qldau)
!> @endnote
!> @warning last index of pns,qns is not the spin index but sra index.
!> @endwarning
!------------------------------------------------------------------------------------
module mod_overlap
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculates the overlap integral of test function PHI with regular or irregular wavefunction 
  !> Author: Ph. Mavropoulos
  !> Category: lda+u, KKRhost
  !> Deprecated: False 
  !> Calculates the overlap integral of test function PHI with regular or irregular
  !> wavefunction (PZ, QZ: spherical part of the large component, 
  !> PQNS: nonspherical part of the large component) for LDA+U. 
  !> The overlap is given out in array `RESULT`
  !>                                                                   
  !> In the inner region the non-spherical wavefunctions are approximated as: 
  !>                                                                   
  !>  * The regular one (ir < irmin = irws-irns) :                     
  !>  \begin{equation}
  !>   pns\left(ir,lm1,lm2\right) = pz\left(ir,l1\right) ar\left(lm1,lm2\right)
  !>  \end{equation}
  !> where \(pz\) is the regular wavefunction of the spherically symmetric
  !> part of the potential and \(ar\) the alpha matrix (see sub. regns)
  !>                                                                   
  !>  * The irregular one (ir < irmin) :
  !>  \begin{equation}
  !>   qns\left(ir,lm1,lm2\right) = pz\left(ir,l1\right) cr\left(lm1,lm2\right) + qz\left(ir,l1\right) * dr\left(lm1,lm2\right) 
  !>  \end{equation}
  !> where \(pz\) is the regular and \(qz\) is the irregular wavefunction of the 
  !> spherically symmetric part of the potential and \(cr\), \(dr\) the matrices 
  !> calculated at the point \(irmin\) (see sub. `irwns`)
  !>                                                                   
  !> The integrand is convoluted with the shape functions:
  !> \begin{equation}
  !>    Result = \int \bar{phi_l} \sum_{L''L'''} R_L''L' Theta_L''' GAUNT_LL''L'''   
  !> \end{equation}
  !>  Finally, the result should then be transformed to complex spherical harmonics basis
  !-------------------------------------------------------------------------------
  !> @note Here only large component is used, but was PHI normalised according to 
  !> large + small component? (qldau)
  !> @endnote
  !> @warning last index of pns,qns is not the spin index but sra index.
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine overlap(result,phi,pz,qz,pqns,acr,dr,lirreg,ipan,ircut,drdi,irmin,lphi,& 
    ipand,lmaxd,lmmaxd,mmaxd,lmpotd,irmind,irmd)

    use :: mod_csimpk
    use :: mod_rinit
    use :: mod_cinit
    implicit none
    ! ..
    ! .. Scalar Arguments ..
    integer, intent(in) :: ipan   !! Number of panels in non-MT-region
    integer, intent(in) :: lphi   !! Is the angular momentum of PHI. l-value for LDA+U
    integer, intent(in) :: irmd   !! Maximum number of radial points
    integer, intent(in) :: irmin  !! Max R for spherical treatment
    integer, intent(in) :: ipand  !! Number of panels in non-spherical part
    integer, intent(in) :: lmaxd  !! Maximum l component in wave function expansion
    integer, intent(in) :: mmaxd
    integer, intent(in) :: lmmaxd !! (KREL+KORBIT+1)*(LMAX+1)**2
    integer, intent(in) :: lmpotd !! (lpot+1)**2
    integer, intent(in) :: irmind !! irmd - irnsd
    logical, intent(in) :: lirreg !! Is true if the irrregular wavefunction is to be used
    integer, dimension(0:ipand), intent(in) :: ircut  !! r points of panel borders
    real (kind=dp), dimension(irmd), intent(in) :: drdi !! Derivative dr/di
    complex (kind=dp), dimension(irmd), intent(in) :: phi
    complex (kind=dp), dimension(irmd, 0:lmaxd), intent(in) :: pz
    complex (kind=dp), dimension(irmd, 0:lmaxd), intent(in) :: qz
    complex (kind=dp), dimension(lmmaxd, lmmaxd), intent(in) :: dr
    complex (kind=dp), dimension(lmmaxd, lmmaxd), intent(in) :: acr
    complex (kind=dp), dimension(lmmaxd, lmmaxd, irmind:irmd, 2), intent(in) :: pqns
    ! .. Output variables
    complex (kind=dp), dimension(mmaxd, mmaxd), intent(out) :: result !! Overlap matrix
    ! ..
    ! .. Locals ..
    integer :: lphisq, mmax, irs1, irc1
    integer :: lm1, lm2, lm3, mm1, mm2, mm3, ir
    real (kind=dp), dimension(lmmaxd, lmmaxd, lmpotd) :: gaunt
    complex (kind=dp), dimension(irmd) :: wint
    ! ..
    mmax = 2*lphi + 1
    lphisq = lphi*lphi
    irs1 = ircut(1)
    irc1 = ircut(ipan)

    call rinit(lmmaxd*lmmaxd*lmpotd, gaunt)

    do lm1 = 1, lmmaxd
      gaunt(lm1, lm1, 1) = 1.e0_dp
    end do
    call cinit(mmaxd*mmaxd, result)

    ! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    ! Loop over mm1,mm2:
    do mm1 = 1, mmax
      lm1 = lphisq + mm1
      do mm2 = 1, mmax
        lm2 = lphisq + mm2
        ! Set up integrand
        call cinit(irmd, wint)
        ! ----------------------------------------------------------------------
        do mm3 = 1, mmax
          lm3 = lphisq + mm3

          ! -> First inner part (up to irmin, no thetas)
          do ir = 2, irmin
            wint(ir) = wint(ir) + phi(ir)*pz(ir, lphi)*acr(lm3, lm2)*gaunt(lm1, lm3, 1)
          end do

          if (lirreg) then
            do ir = 2, irmin
              wint(ir) = wint(ir) + phi(ir)*qz(ir, lphi)*dr(lm3, lm2)*gaunt(lm1, lm3, 1)
            end do
          end if

          ! -> Next middle part (from irmin+1 to irc1, no thetas)
          do ir = irmin + 1, irs1
            wint(ir) = wint(ir) + phi(ir)*pqns(lm3, lm2, ir, 1)*gaunt(lm1, lm3, 1)
          end do

          ! -> Finally last part - from irc1+1 to irs1 - with THETAS and proper
          ! GAUNTS if we integrate in cell:
          ! ccc               DO  LM4 = 1,LMPOTD
          ! ccc                  IF (LMSP(LM4).GT.0) THEN
          ! ccc                     DO IR = IRS1+1,IRC1
          ! ccc                        WINT(IR) = WINT(IR)
          ! ccc     &                           + PHI(IR) * PQNS(LM3,LM2,IR,1)
          ! ccc     &                           * GAUNT(LM1,LM3,LM4)
          ! ccc     &                           * THETAS(IR-IRS1,LM4)
          ! ccc                     ENDDO
          ! ccc                  END IF
          ! ccc               ENDDO

          ! or still without THETAS if we integrate in sphere

          do ir = irs1 + 1, irc1
            wint(ir) = wint(ir) + phi(ir)*pqns(lm3, lm2, ir, 1)*gaunt(lm1, lm3, 1)
          end do
        end do
        ! ----------------------------------------------------------------------

        call csimpk(wint, result(mm1,mm2), ipan, ircut, drdi(1))
      end do
    end do
    ! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
  end subroutine overlap

end module mod_overlap
