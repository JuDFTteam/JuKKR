subroutine shape_corr(lpot, natyp, gsh, ilm_map, imaxsh, lmsp, ntcell, w, yr, &
  lassld, lmpotd, natypd, ngshd)
  use :: mod_datatypes, only: dp
  ! **********************************************************************
  ! *  Prepares shape corrections using gaussian quadrature as given by  *
  ! *  m. abramowitz and i.a. stegun, handbook of mathematical functions *
  ! *  nbs applied mathematics series 55 (1968), pages 887 and 916       *
  ! *                                                                    *
  ! *  the parameter LASSLD has to be chosen such that                   *
  ! *                        l1+l2+l3 .le. 2*LASSLD                      *
  ! *                                                                    *
  ! **********************************************************************

  implicit none
  real (kind=dp), parameter :: eps=1.0D-12
  ! ..
  ! .. Scalar Arguments ..
  integer :: lassld, lmpotd, natypd, ngshd
  integer :: lpot, natyp
  ! ..
  ! .. Array Arguments ..
  real (kind=dp) :: gsh(*), w(lassld), yr(lassld, 0:lassld, 0:lassld)
  integer :: ilm_map(ngshd, 3), imaxsh(0:lmpotd), lmsp(natypd, *), ntcell(*)
  ! ..
  ! .. Local Scalars ..
  real (kind=dp) :: factor, gaunt, s
  integer :: i, iat, icell, isum, j, l1, l2, l3, lm1, lm2, lm3, m1, m1a, m1s, &
    m2, m2a, m2s, m3, m3a, m3s
  logical :: triangle
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic :: abs, real, sign
  ! ..
  ! .. External Subroutines ..
  external :: rcstop, triangle
  ! ..

  ! -> set up of the gaunt coefficients with an index field
  ! so that  c(lm,lm',lm'') is mapped to c(i)
  i = 1
  ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
  do l1 = 0, lpot
    ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    do m1 = -l1, l1

      lm1 = l1*l1 + l1 + m1 + 1
      imaxsh(lm1-1) = i - 1
      ! llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
      do l3 = 0, lpot*2
        ! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
        do m3 = -l3, l3

          lm3 = l3*l3 + l3 + m3 + 1
          isum = 0

          do iat = 1, natyp
            icell = ntcell(iat)
            ! write(*,*) 'test icell=ntcell(iat) in shape_corr.f',
            ! +                icell,iat
            isum = isum + lmsp(icell, lm3)
          end do

          ! ======================================================================
          if (isum>0) then
            do l2 = 0, lpot
              ! ----------------------------------------------------------------------
              if (triangle(l1,l2,l3)) then
                do m2 = -l2, l2

                  lm2 = l2*l2 + l2 + m2 + 1

                  ! -> use the m-conditions for the gaunt coefficients not to
                  ! be 0

                  m1s = sign(1, m1)
                  m2s = sign(1, m2)
                  m3s = sign(1, m3)
                  ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                  if (m1s*m2s*m3s>=0) then
                    m1a = abs(m1)
                    m2a = abs(m2)
                    m3a = abs(m3)
                    factor = 0.0e0_dp

                    if (m1a+m2a==m3a) factor = factor + &
                      real(3*m3s+sign(1,-m3), kind=dp)/8.0e0_dp

                    if (m1a-m2a==m3a) factor = factor + &
                      real(m1s, kind=dp)/4.0e0_dp

                    if (m2a-m1a==m3a) factor = factor + &
                      real(m2s, kind=dp)/4.0e0_dp
                    ! ......................................................................
                    if (abs(factor)>eps) then

                      if (m1s*m2s/=1 .or. m2s*m3s/=1 .or. m1s*m3s/=1) &
                        factor = -factor

                      s = 0.0e0_dp
                      do j = 1, lassld
                        s = s + w(j)*yr(j, l1, m1a)*yr(j, l2, m2a)*yr(j, l3, &
                          m3a)
                      end do

                      gaunt = s*factor
                      if (abs(gaunt)>1e-10_dp) then
                        gsh(i) = gaunt
                        ilm_map(i, 1) = lm1
                        ilm_map(i, 2) = lm2
                        ilm_map(i, 3) = lm3
                        i = i + 1
                      end if
                    end if
                    ! ......................................................................
                  end if
                  ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                end do
              end if
              ! ----------------------------------------------------------------------
            end do
          end if
          ! ======================================================================
        end do
        ! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
      end do
      ! llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
    end do
    ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
  end do
  ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

  imaxsh(lm1) = i - 1
  write (1337, fmt=100) imaxsh(lm1), ngshd
  if (imaxsh(lm1)>ngshd) call rcstop('SHAPE   ')

100 format (' >>> SHAPE : IMAXSH(', i4, '),NGSHD :', 2i6)

end subroutine shape_corr

function triangle(l1, l2, l3)
  implicit none
  integer :: l1, l2, l3
  logical :: triangle
  intrinsic :: mod
  ! ..
  triangle = (l1>=abs(l3-l2)) .and. (l1<=(l3+l2)) .and. (mod((l1+l2+l3),2)==0)
end function triangle
