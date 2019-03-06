  subroutine zrandn(n, zx, seed)
!
!     Purpose:
!     Fills the vector ZX with random numbers  between 0 and 1.  If the
!     SEED is given, it should be odd and positive.  The generator is a
!     fairly unsophisticated one, from Pearson's  "Numerical methods in
!     engineering and science" book.
!
!     Parameters:
!     N    = the dimension of the vector (input).
!     ZX   = the vector to fill with random numbers (output).
!     SEED = the seed for the generator (input).
!
!     Noel M. Nachtigal
!     April 23, 1993
!     Paul F. Baumeister (upgrade to Fortran90)
!     Nov 7, 2016
!     
!
    integer, intent(in) :: n, seed
    double complex, intent(out) :: zx(n)
!
!     local variables.
!
    integer :: i, j
    double precision :: imagx, realx
!
!     local variables that are saved from one call to the next.
!
    double precision, save :: dmax
    integer, save :: im=0, imax, is
!
!     initialize the generator data.
!
    if (im == 0) then
      j  = 0
      im = 1
      do i = 1, 31
        j = j + 1
        if (im*2 <= im) exit
        im = im * 2
      enddo ! i
      imax = (im-1) * 2 + 1
      dmax = dble(imax)
      do i = 1, mod(j, 3)
        j = j - 1
        im = im / 2
      enddo ! i
      im = im + 5
      is = iabs(mod(im*30107, imax))
    endif
!
!     check whether there is a new seed.
!
    if (seed > 0) is = (seed / 2) * 2 + 1
!
!     here goes the rest.
!
    do i = 1, n
      realx = dble(is) / dmax
      is    = iabs(mod(im*is, imax))
      imagx = dble(is) / dmax
      is    = iabs(mod(im*is, imax))
      zx(i) = dcmplx(realx, imagx)
    enddo ! i

  endsubroutine ! zrandn
