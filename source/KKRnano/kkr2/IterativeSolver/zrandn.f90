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
!
      integer, intent(in) :: n, seed
      double complex, intent(out) :: zx(n)
!
!     local variables.
!
      integer :: i, j
      double precision imagx, realx
!
!     local variables that are saved from one call to the next.
!
      double precision, save :: dmax
      integer, save :: im=0, imax, is
      
      intrinsic :: dble, dcmplx, iabs, mod
!
!     initialize the generator data.
!
      if (im == 0) then
         j  = 0
         im = 1
         do 10 i = 1, 31
            j = j + 1
            if (im*2 <= im) go to 20
            im = im * 2
 10      continue
 20      imax = (im-1) * 2 + 1
         dmax = dble(imax)
         do 30 i = 1, mod(j,3)
            j = j - 1
            im = im / 2
 30      continue
         im = im + 5
         is = iabs(mod(im*30107,imax))
      end if
!
!     check whether there is a new seed.
!
      if (seed > 0) is = (seed / 2) * 2 + 1
!
!     here goes the rest.
!
      do 40 i = 1, n
         realx = dble(is) / dmax
         is    = iabs(mod(im*is,imax))
         imagx = dble(is) / dmax
         is    = iabs(mod(im*is,imax))
         zx(i) = dcmplx(realx, imagx)
 40   continue

      endsubroutine ! zrandn
