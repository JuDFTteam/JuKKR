module mod_amemagvec

  private
  public :: amemagvec

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculates the angular matrix elements
  !> Author: 
  !> Category: KKRhost, dirac, physical-observables
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>                                                    
  !> Calculate the angular matrix elements connected with       
  !>                                                            
  !> 1:       < LAM | sigma(ipol) | LAM' >    spin moment       
  !> 2:       < LAM |     l(ipol) | LAM' >    orbital moment    
  !> 3:       < LAM |     T(ipol) | LAM' >    spin dipole moment
  !> 4:       < LAM |  B_hf(ipol) | LAM' >    hyperfine field   
  !>                                                            
  !> ipol= 1,2,3  ==  (+),(-),(z)                          
  !-------------------------------------------------------------------------------
  subroutine amemagvec(irel, iprint, nkm, amemvec, ikmllim1, ikmllim2, imkmtab, cgc, nlmax, nkmmax, nkmpmax, nmvecmax)

    use :: mod_rmatstr, only: rmatstr
    use :: mod_rinit, only: rinit
    use :: mod_datatypes, only: dp

    implicit none

    ! Dummy arguments

    integer :: iprint, irel, nkm, nkmmax, nkmpmax, nlmax, nmvecmax
    real (kind=dp) :: amemvec(nkmmax, nkmmax, 3, nmvecmax), cgc(nkmpmax, 2)
    integer :: ikmllim1(nkmmax), ikmllim2(nkmmax), imkmtab(nkmmax)

    ! Local variables

    character (len=1) :: chpol(3)
    integer :: i, ikm1, ikm2, imkm, imv, imvec, ipol, j1p05, j2p05, k, k1, k2, kap1, kap2, l, l1, l2, lb1, lb2, m2, msm05, mue1m05, mue2m05, nk, nmvec
    integer :: iabs, nint
    character (len=20) :: str20
    real (kind=dp) :: sum, xj, xjm, xjp, xm, xynorm
    character (len=4) :: txtmvec(4)

    data chpol/'+', '-', 'z'/
    data txtmvec/'spin', 'orb ', 'T_z ', 'B_hf'/

    nk = 2*nlmax - 1
    ! XYNORM = DSQRT(2.0D0)
    xynorm = 2.0e0_dp

    call rinit(nkmmax*nkmmax*3*nmvecmax, amemvec)

    if (irel<=1) return

    ! ----------------------------------------------------------------------
    ! find the bounding indices  IKMLLIM1  and  IKMLLIM2  for IKM-loops
    ! assuming that the matrix elements are diagonal with respect to l
    ! this does not hold for B_hf for which there are l-(l+/-2)-terms
    ! ----------------------------------------------------------------------

    i = 0
    do k = 1, nk
      l = k/2
      xjm = l - 0.5e0_dp
      xjp = l + 0.5e0_dp
      if (mod(k,2)==1) then
        xj = l + 0.5e0_dp
      else
        xj = l - 0.5e0_dp
      end if
      do xm = -xj, +xj
        i = i + 1
        ikmllim1(i) = nint(l*2*(xjm+0.5e0_dp)+1)
        ikmllim2(i) = nint(l*2*(xjp+0.5e0_dp)+2*xjp+1)
      end do
    end do

    ! ----------------------------------------------------------------------

    ikm1 = 0
    do k1 = 1, nk
      l1 = k1/2
      if (mod(k1,2)==0) then
        kap1 = l1
        lb1 = l1 - 1
      else
        kap1 = -l1 - 1
        lb1 = l1 + 1
      end if
      j1p05 = iabs(kap1)

      do mue1m05 = -j1p05, j1p05 - 1
        ikm1 = ikm1 + 1
        imkm = lb1*2*j1p05 + j1p05 + mue1m05 + 1
        imkmtab(ikm1) = imkm

        ikm2 = 0
        do k2 = 1, nk
          l2 = k2/2
          if (mod(k2,2)==0) then
            kap2 = l2
            lb2 = l2 - 1
          else
            kap2 = -l2 - 1
            lb2 = l2 + 1
          end if
          j2p05 = iabs(kap2)

          do mue2m05 = -j2p05, j2p05 - 1
            ikm2 = ikm2 + 1
            ! ----------------------------------------------------------------------
            if ((mue1m05-mue2m05)==+1) then
              amemvec(ikm1, ikm2, 1, 1) = xynorm*cgc(ikm1, 2)*cgc(ikm2, 1)

              sum = 0e0_dp
              do msm05 = -1, 0
                m2 = mue2m05 - msm05
                if (abs(m2)<=l2) sum = sum + cgc(ikm1, msm05+2)*cgc(ikm2, msm05+2)*sqrt(real((l2-m2)*(l2+m2+1),kind=dp))
              end do
              amemvec(ikm1, ikm2, 1, 2) = sum
            end if

            if ((mue1m05-mue2m05)==-1) then
              amemvec(ikm1, ikm2, 2, 1) = xynorm*cgc(ikm1, 1)*cgc(ikm2, 2)

              sum = 0e0_dp
              do msm05 = -1, 0
                m2 = mue2m05 - msm05
                if (abs(m2)<=l2) sum = sum + cgc(ikm1, msm05+2)*cgc(ikm2, msm05+2)*sqrt(real((l2+m2)*(l2-m2+1),kind=dp))

              end do
              amemvec(ikm1, ikm2, 2, 2) = sum
            end if

            if ((mue1m05-mue2m05)==0) then
              amemvec(ikm1, ikm2, 3, 1) = cgc(ikm1, 2)*cgc(ikm2, 2) - cgc(ikm1, 1)*cgc(ikm2, 1)

              sum = 0e0_dp
              do msm05 = -1, 0
                m2 = mue2m05 - msm05
                sum = sum + cgc(ikm1, msm05+2)*cgc(ikm2, msm05+2)*m2
              end do
              amemvec(ikm1, ikm2, 3, 2) = sum
            end if

            ! ----------------------------------------------------------------------
          end do
        end do
      end do
    end do

    nmvec = 2
    ! ----------------------------------------------------------------------
    if (iprint<=90) return

    do imvec = 1, nmvec
      imv = min(imvec, 2)
      do ipol = 1, 3
        str20 = 'A  ' // txtmvec(imv) // '  (' // chpol(ipol) // ')'
        call rmatstr(str20, 12, amemvec(1,1,ipol,imvec), nkm, nkmmax, 3, 3, 1e-8_dp, 6)
      end do
    end do
  end subroutine amemagvec

end module mod_amemagvec
