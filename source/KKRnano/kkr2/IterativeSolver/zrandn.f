      SUBROUTINE ZRANDN (N,ZX,SEED)
C
C     Purpose:
C     Fills the vector ZX with random numbers  between 0 and 1.  If the
C     SEED is given, it should be odd and positive.  The generator is a
C     fairly unsophisticated one, from Pearson's  "Numerical methods in
C     engineering and science" book.
C
C     Parameters:
C     N    = the dimension of the vector (input).
C     ZX   = the vector to fill with random numbers (output).
C     SEED = the seed for the generator (input).
C
C     Noel M. Nachtigal
C     April 23, 1993
C
C**********************************************************************
c
      intrinsic dble, dcmplx, iabs, mod
c
      integer n, seed
      double complex zx(n)
c
c     local variables.
c
      integer i, j
      double precision imagx, realx
c
c     local variables that are saved from one call to the next.
c
      double precision dmax
      integer im, imax, is
      save dmax, im, imax, is
      data im/0/
c
c     initialize the generator data.
c
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
c
c     check whether there is a new seed.
c
      if (seed > 0) is = (seed / 2) * 2 + 1
c
c     here goes the rest.
c
      do 40 i = 1, n
         realx = dble(is) / dmax
         is    = iabs(mod(im*is,imax))
         imagx = dble(is) / dmax
         is    = iabs(mod(im*is,imax))
         zx(i) = dcmplx(realx,imagx)
 40   continue
c
      return

      end
