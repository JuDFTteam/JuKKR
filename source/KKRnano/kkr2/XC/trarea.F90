      subroutine trarea(a, b, lmax)
!     from complex to real  (differenciated spherical harmonics)
      integer, intent(in) :: lmax
      double complex, intent(in) :: a(:)
      double precision, intent(out) :: b(:)

!     .. locals ..
      double precision, parameter :: mrtwo=-1.414213562373d0
!       , rhalf=1.d0/rtwo
!       double precision, parameter :: rhalf=0.7071067811865476d0
      double complex, parameter :: ci=(0.d0,1.d0)
      double precision :: srtwo
      integer :: ilm, l, m
!     ..
!
!    calculate real the spherical harmonics derivetived
!
      ilm = 0
      do l = 0, lmax
        ilm = ilm + l + 1
        b(ilm-0) = dreal(a(ilm)) ! m==0
        srtwo = mrtwo
        do m = 1, l
          b(ilm-m) = mrtwo*dimag(a(ilm-m))
          b(ilm+m) = srtwo*dreal(a(ilm+m))
          srtwo = -srtwo
        enddo ! m
        ilm = ilm + l
      enddo ! l
      endsubroutine trarea
