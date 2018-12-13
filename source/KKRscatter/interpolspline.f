      subroutine interpolspline(rmesh,rmeshnew,vpot,vpotnew,
     +                          nrmax,nrmaxnew)
      use mod_splint, only: splint_real
      implicit none
!interface
      integer                       :: nrmax
      integer                       :: nrmaxnew
      double precision              :: rmesh(nrmax)
      double precision              :: rmeshnew(nrmaxnew)
      double precision              :: vpot(nrmax)
      double precision              :: vpotnew(nrmaxnew)
!local
      double precision              :: maxa
      double precision              :: spline(nrmax)
      double precision              :: PARSUM, PARSUMDERIV,R0
      integer                       :: ir
      maxa = 1.d35
      call spline_real(nrmax,rmesh,vpot,nrmax,maxa,maxa,spline)  
!           CALL SPLINE(IRMDJJ,R,VM2Z,NR,maxa,maxa,VM2ZB)  

      DO ir = 1,nrmaxnew
         R0 = rmeshnew(IR)
       call splint_real(rmesh,vpot,spline,nrmax,R0,PARSUM,PARSUMDERIV)
         vpotnew(IR) = PARSUM
      END DO
      end subroutine  interpolspline

