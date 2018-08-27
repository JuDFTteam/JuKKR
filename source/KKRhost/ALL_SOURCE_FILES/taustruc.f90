module mod_taustruct

contains

subroutine taustruct(drot, nsym, symunitary, nkm, nq, nqmax, nkmmax, iprint, &
  irel)
  ! ********************************************************************
  ! *                                                                  *
  ! *   find the structure of the site-diagonal TAU - matrices  TAUQ   *
  ! *                                                                  *
  ! ********************************************************************
  use :: mod_datatypes
  use mod_cmatstr
   use mod_cinit
  implicit none

  ! Dummy arguments
  integer :: iprint, irel, nkm, nkmmax, nq, nqmax, nsym
  complex (kind=dp) :: drot(nkmmax, nkmmax, 48)
  logical :: symunitary(48)

  ! Local variables
  real (kind=dp) :: abst, x
  integer :: i, i0, imweight, iq, isym, iw, iwr, j, k, l, lin, nelmt, &
    nkmq(nqmax), nkmtop, nlin, non0(nqmax)
  complex (kind=dp) :: st(nkmmax, nkmmax), tauk(nkmmax, nkmmax, nqmax)

  do iq = 1, nq
    nkmq(iq) = nkm
  end do

  imweight = 0
  nelmt = 0
  nlin = 0
  iw = 6

  do iq = 1, nq
    non0(iq) = 0
    nkmtop = nkmq(iq)

    if (iprint>0) write (1337, 120) iq
    do i = 1, nkmtop
      do j = 1, nkmtop
        st(i, j) = 0.0e0_dp

        call cinit(nkmmax*nkmmax*nqmax, tauk)

        do isym = 1, nsym
          i0 = iq

          if (symunitary(isym)) then
            do l = 1, nkmtop
              do k = 1, nkmtop
                tauk(k, l, i0) = tauk(k, l, i0) + drot(i, k, isym)*conjg(drot( &
                  j,l,isym))
              end do
            end do
          else
            do l = 1, nkmtop
              do k = 1, nkmtop
                tauk(l, k, i0) = tauk(l, k, i0) + drot(i, k, isym)*conjg(drot( &
                  j,l,isym))
              end do
            end do
          end if
        end do

        lin = 0
        iwr = 0
        do k = 1, nkmq(iq)
          do l = 1, nkmq(iq)
            abst = abs(tauk(k,l,iq))
            st(i, j) = st(i, j) + abst
            if (abst>1e-8_dp) then
              if (aimag(tauk(k,l,iq))>1e-5_dp) then
                if (iprint>0) write (1337, *) ' Im(Weight) > 1D-5 ', i, j, k, &
                  l
                imweight = 1
              end if
              x = real(tauk(k,l,iq))/real(nsym, kind=dp)

              if (iprint>1) then
                if (iwr==0) then
                  iwr = 1
                  write (iw, 100) i, j, iq, x, k + (iq-1)*nkm, l + (iq-1)*nkm
                else
                  write (iw, 110) x, k + (iq-1)*nkm, l + (iq-1)*nkm
                end if
              end if
              lin = lin + 1
            end if

          end do
        end do

        if (lin>0) then
          nlin = nlin + lin
          nelmt = nelmt + 1
          non0(iq) = non0(iq) + 1
        end if

        if (abs(st(i,j))>1e-5_dp) st(i, j) = 2

      end do
    end do

    if (iprint>1) call cmatstr('TAU-MAT', 7, st, nkmtop, nkmmax, irel, irel, &
      0, 1e-8_dp, 6)
  end do

  write (1337, 130) nelmt, (non0(iq), iq=1, nq)
  write (1337, 140) nlin

  if (imweight/=0) write (1337, 150)

  ! -----------------------------------------------------------------------

  return
100 format ('     TAUQ(', i2, ',', i2, ',', i2, ') =  ', '   ', f10.4, ' * <', &
    i3, '|T(K)|', i3, '>')
110 format (23x, ' + ', f10.4, ' * <', i3, '|T(K)|', i3, '>')
120 format (/, /, &
    ' ===========================================================', /, &
    '   structure of  TAU-matrix   INT <i|t(k)|j>     IQ=', i3, /, &
    ' ===========================================================', /)
130 format (/, 5x, 'non-0 TAU-elements          ', i5, '   Q:', 80i4)
140 format (5x, 'terms to sum up             ', i5, /)
150 format (/, 5x, 50('#'), /, 5x, 'WARNING: complex TAU weights found', /, &
    5x, 'this may occur for rotated magnetic moments', /, 5x, &
    'relevant only for tetrahedron BZ-integration', /, 5x, 50('#'), /)
end subroutine taustruct

end module mod_taustruct
