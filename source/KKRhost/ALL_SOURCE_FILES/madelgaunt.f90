subroutine madelgaunt(lpot, yrg, wg, cleb, icleb, iend, lassld, nclebd)
  implicit none
!..
!.. Scalar arguments
  integer :: lpot, iend
  integer :: lassld, nclebd
!..
!.. Array arguments
!.. Attention: Dimension NCLEBD appears sometimes as NCLEB1
!..            an empirical factor - it has to be optimized
  double precision :: yrg(lassld, 0:lassld, 0:lassld), wg(lassld)
  double precision :: cleb(nclebd)
  integer :: icleb(nclebd, 3)
!..
!.. Local scalars
  double precision :: clecg, factor, s
  integer :: i, j, l1, l2, l3, m1, m1a, m1s, m2, m2a, m2s, m3, m3a, m3s
!..
!.. Intrinsic functions
  intrinsic :: abs, atan, dble, sign

! --> set up of the gaunt coefficients with an index field
!     recognize that they are needed here only for l3=l1+l2

  if (2*lpot>lassld) then
    write (6, *) 'Dim ERROR in MADELGAUNT -- 2*LPOT > LASSLD', 2*lpot, lassld
    stop
  end if

  i = 1
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
  do l1 = 0, lpot
    do l2 = 0, lpot
      l3 = l1 + l2
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      do m1 = -l1, l1
        do m2 = -l2, l2
          do m3 = -l3, l3
            m1s = sign(1, m1)
            m2s = sign(1, m2)
            m3s = sign(1, m3)
! **********************************************************************
            if (m1s*m2s*m3s>=0) then
              m1a = abs(m1)
              m2a = abs(m2)
              m3a = abs(m3)

              factor = 0.0d0
              if (m1a+m2a==m3a) factor = factor + dble(3*m3s+sign(1,-m3))/ &
                8.0d0
              if (m1a-m2a==m3a) factor = factor + dble(m1s)/4.0d0
              if (m2a-m1a==m3a) factor = factor + dble(m2s)/4.0d0
! ======================================================================
              if (factor/=0.0d0) then
                if (m1s*m2s/=1 .or. m2s*m3s/=1 .or. m1s*m3s/=1) &
                  factor = -factor

                s = 0.0d0
                do j = 1, lassld
                  s = s + wg(j)*yrg(j, l1, m1a)*yrg(j, l2, m2a)*yrg(j, l3, m3a &
                    )
                end do

                clecg = s*factor
! ----------------------------------------------------------------------
                if (abs(clecg)>1.d-10) then
                  cleb(i) = clecg
                  icleb(i, 1) = l1*(l1+1) + m1 + 1
                  icleb(i, 2) = l2*(l2+1) + m2 + 1
                  icleb(i, 3) = l3*(l3+1) + m3 + 1
                  i = i + 1
                  if (i>nclebd) then
                    write (6, fmt='(2I10)') i, nclebd
                    stop ' Dim stop in MADELGAUNT '
                  end if
                end if
! ----------------------------------------------------------------------
              end if
! ======================================================================
            end if
! **********************************************************************
          end do
        end do
      end do
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    end do
  end do
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
  iend = i - 1
end subroutine
