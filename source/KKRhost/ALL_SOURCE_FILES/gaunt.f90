module mod_gaunt
  
  private
  public :: gaunt

contains

  !-------------------------------------------------------------------------------
  !> Summary: Computes gaunt coefficients
  !> Author: 
  !> Category: KKRhost, special-functions
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> - fills the array cleb with the gaunt coeffients ,i.e.
  !>   the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)
  !>   but only for lm2.le.lm1 and lm3>1
  !> - calculate the pointer array jend  to project the indices
  !>   array cleb with the same lm3,l1,l2 values - because of
  !>   the special ordering of array cleb only the last index
  !>   has to be determined.
  !> (the parameter n has to be chosen that l1+l2+l3 .lt. 2*n)
  !> using gaussian quadrature as given by
  !> M. Abramowitz and I.A. Stegun, Handbook of Mathematical Functions,
  !> NBS Applied Mathematics Series 55 (1968), pages 887 and 916
  !> M. Weinert and E. Wimmer
  !> Northwestern University March 1980
  !>
  !> An index array -icleb- is used to save storage place.
  !> fills the array loflm which is used to determine the
  !> l-value of a given lm-value.
  !> @note
  !> This subroutine has to be called only once!
  !> @endnote
  !>
  !> B. Drittler   November 1987
  !>
  !> modified gaunt coefficients are als calculated defined by
  !> the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)*i**(l2-l1+l3)
  !> 
  !>
  !> Attention : ncleb is an empirical factor - it has to be optimized
  !-------------------------------------------------------------------------------
  subroutine gaunt(lmax, lpot, w, yr, cleb, loflm, icleb, iend, jend, ncleb, lmaxd, lmgf0d, lmpotd)

    use :: mod_datatypes, only: dp
    use :: mod_rcstop
    implicit none

    real (kind=dp), parameter :: eps = 1.0e-12_dp
    complex (kind=dp), parameter :: ci = (0.0e0_dp, 1.0e0_dp)
    integer :: lmpotd, lmgf0d, lmaxd, ncleb
    integer :: iend, lmax, lpot
    real (kind=dp) :: cleb(ncleb, 2), w(*), yr(4*lmaxd, 0:4*lmaxd, 0:4*lmaxd)
    integer :: icleb(ncleb, 4), jend(lmpotd, 0:lmaxd, 0:lmaxd), loflm(*)
    real (kind=dp) :: clecg, factor, fci, s
    integer :: i, j, l, l1, l1p, l2, l2p, l3, lm1, lm2, lm3, lm3p, lmpot, m, m1, m1a, m1s, m2, m2a, m2s, m3, m3a, m3s


    i = 1
    do l = 0, 2*lmax
      do m = -l, l
        loflm(i) = l
        i = i + 1
      end do
    end do

    icleb = 0
    cleb = 0e0_dp
    if (lpot==0) then
      iend = 1
      icleb(1, 1) = (lmax+1)**2
      icleb(1, 3) = 1
    end if

    if (lpot/=0) then

      ! ---> set up of the gaunt coefficients with an index field

      i = 1
      do l3 = 1, lpot
        do m3 = -l3, l3

          do l1 = 0, lmax
            do l2 = 0, l1

              if (mod((l1+l2+l3),2)/=1 .and. (l1+l2-l3)>=0 .and. (l1-l2+l3)>=0 .and. (l2-l1+l3)>=0) then

                fci = real(ci**(l2-l1+l3))
                do m1 = -l1, l1
                  do m2 = -l2, l2

                    ! ---> store only gaunt coeffients for lm2.le.lm1

                    lm1 = l1*(l1+1) + m1 + 1
                    lm2 = l2*(l2+1) + m2 + 1
                    if (lm2<=lm1) then

                      m1s = sign(1, m1)
                      m2s = sign(1, m2)
                      m3s = sign(1, m3)

                      if (m1s*m2s*m3s>=0) then

                        m1a = abs(m1)
                        m2a = abs(m2)
                        m3a = abs(m3)

                        factor = 0.0_dp

                        if (m1a+m2a==m3a) factor = factor + real(3*m3s+sign(1,-m3), kind=dp)/8.0e0_dp
                        if (m1a-m2a==m3a) factor = factor + real(m1s, kind=dp)/4.0e0_dp
                        if (m2a-m1a==m3a) factor = factor + real(m2s, kind=dp)/4.0e0_dp

                        if (abs(factor)>eps) then

                          if (m1s*m2s/=1 .or. m2s*m3s/=1 .or. m1s*m3s/=1) factor = -factor

                          s = 0.0_dp
                          do j = 1, 4*lmaxd
                            s = s + w(j)*yr(j, l1, m1a)*yr(j, l2, m2a)*yr(j, l3, m3a)
                          end do
                          clecg = s*factor
                          if (abs(clecg)>1.e-10_dp) then
                            cleb(i, 1) = clecg
                            cleb(i, 2) = fci*clecg
                            icleb(i, 1) = lm1
                            icleb(i, 2) = lm2
                            icleb(i, 3) = l3*(l3+1) + m3 + 1
                            icleb(i, 4) = lm2*lmgf0d - (lm2*lm2-lm2)/2 + lm1 - lmgf0d
                            i = i + 1
                          end if

                        end if

                      end if

                    end if

                  end do
                end do
              end if

            end do
          end do
        end do
      end do
      iend = i - 1
      if (ncleb<iend) then
        write (6, fmt=100) ncleb, iend
        call rcstop('33      ')

      else

        ! ---> set up of the pointer array jend,use explicitly
        ! the ordering of the gaunt coeffients

        lmpot = (lpot+1)*(lpot+1)
        do l1 = 0, lmax
          do l2 = 0, l1
            do lm3 = 2, lmpot
              jend(lm3, l1, l2) = 0
            end do
          end do
        end do

        lm3 = icleb(1, 3)
        l1 = loflm(icleb(1,1))
        l2 = loflm(icleb(1,2))

        do j = 2, iend
          lm3p = icleb(j, 3)
          l1p = loflm(icleb(j,1))
          l2p = loflm(icleb(j,2))

          if (lm3/=lm3p .or. l1/=l1p .or. l2/=l2p) then
            jend(lm3, l1, l2) = j - 1
            lm3 = lm3p
            l1 = l1p
            l2 = l2p
          end if

        end do
        jend(lm3, l1, l2) = iend


      end if

    end if



100 format (13x, 'error stop in gaunt : dimension of NCLEB = ', i10, ' too small ', /, 13x, 'change NCLEB to ', i6)
  end subroutine gaunt

end module mod_gaunt
