      module Quadrature_mod
      implicit none
      private
      
      public :: simpson ! functions
!     public :: CSIMPK, SIMP3, SIMPK ! subroutines (deprecated)
!     public :: csum, ssum ! functions

      interface simpson
        module procedure csimpkf, simpkf, simp3f
      endinterface
      
      double precision, parameter :: a1 = 4.d0/3.d0, a2 = 2.d0/3.d0
      
      contains

      double precision function simpkf(f,ipan,ircut,drdi) result(fint)
      integer, intent(in) :: ipan
      double precision, intent(in) :: f(:)
      double precision, intent(in) :: drdi(:)
      integer, intent(in) :: ircut(0:ipan)
      
      double precision :: t(size(f))
      integer :: ien, ip, ist, n
      
      fint = 0.d0
      do ip = 1, ipan
!---> loop over kinks
        ist = ircut(ip-1) + 1
        ien = ircut(ip)

        t(ist:ien) = f(ist:ien)*drdi(ist:ien)

        if (mod(ien-ist, 2) == 0) then
          fint = fint + (t(ist) - t(ien))/3.d0
          ist = ist + 1
        else
        
!---> four point lagrange integration for the first step
          fint = fint + (9.d0*t(ist) + 19.d0*t(ist+1) - 5.d0*t(ist+2) + t(ist+3))/24.d0 + (t(ist+1) - t(ien))/3.d0
          ist = ist + 2
        endif
        n = (ien-ist+1)/2

!---> calculate with an extended 3-point-simpson
        fint = fint + a1*ssum(n, t(ist:), 2) + a2*ssum(n, t(ist+1:), 2)
      enddo ! ipan

      endfunction ! simpkf

      
      subroutine simpk(f,fint,ipan,ircut,drdi)
! **********************************************************************
!     this subroutine does an integration up to rcut of an
!     real function f with an extended 3-point-simpson :
!
!                            rcut
!                      fint = { f(r') dr'
!                             0
!
!     modified for functions with kinks - at each kink the
!     integration is restarted .
!
!     attention : input f is destroyed !
!
!-----------------------------------------------------------------------
      double precision, intent(out) :: fint ! result
      integer, intent(in) :: ipan
      double precision, intent(in) :: drdi(*)
      double precision, intent(inout) :: f(*)
      integer, intent(in) :: ircut(0:ipan)

      integer :: ien, ip, ist, n
      
      fint = 0.d0
      do ip = 1, ipan

      !---> loop over kinks
        ist = ircut(ip-1) + 1
        ien = ircut(ip)

        f(ist:ien) = f(ist:ien)*drdi(ist:ien)
        if (mod(ien-ist, 2) == 0) then
          fint = fint + (f(ist)-f(ien))/3.d0
          ist = ist + 1
        else
        
!---> four point Lagrange integration for the first step
          fint = fint + (9.d0*f(ist) + 19.d0*f(ist+1) - 5.d0*f(ist+2) + f(ist+3))/24.d0 + (f(ist+1) - f(ien))/3.d0
          ist = ist + 2
        endif
        n = (ien-ist+1)/2

!---> calculate with an extended 3-point-simpson
        fint = fint + a1*ssum(n,f(ist),2) + a2*ssum(n,f(ist+1),2)
      enddo ! ip

      endsubroutine simpk
      
      

      double complex function csimpkf(cf,ipan,ircut,drdi) result(cfint)
      integer, intent(in) :: ipan
      double complex, intent(in) :: cf(:)
      double precision, intent(in) :: drdi(:)
      integer, intent(in) :: ircut(0:ipan)
      
      double complex :: t(size(cf))
      integer :: ien, ip, ist, n
      
      cfint = 0.d0
      do ip = 1, ipan
!---> loop over kinks
        ist = ircut(ip-1) + 1
        ien = ircut(ip)

        t(ist:ien) = cf(ist:ien)*drdi(ist:ien)

        if (mod(ien-ist, 2) == 0) then
          cfint = cfint + (t(ist) - t(ien))/3.d0
          ist = ist + 1
        else
        
!---> four point lagrange integration for the first step
          cfint = cfint + (9.d0*t(ist) + 19.d0*t(ist+1) - 5.d0*t(ist+2) + t(ist+3))/24.d0 + (t(ist+1) - t(ien))/3.d0
          ist = ist + 2
        endif
        n = (ien-ist+1)/2

!---> calculate with an extended 3-point-simpson
        cfint = cfint + a1*csum(n, t(ist:), 2) + a2*csum(n, t(ist+1:), 2)
      enddo ! ipan

      endfunction ! csimpkf
      
      
      
      
      
      subroutine csimpk(cf,cfint,ipan,ircut,drdi)
!-----------------------------------------------------------------------
!     this subroutine does an integration up to rcut of an
!     complex function cf with an extended 3-point-simpson :
!
!                             rcut
!                      cfint = { cf(r') dr'
!                              0
!
!     modified for functions with kinks - at each kink the
!     integration is restarted .
!
!     attention : input cf is destroyed !
!
!-----------------------------------------------------------------------
      double complex, intent(out) :: cfint
      integer, intent(in) :: ipan
      double complex, intent(inout) :: cf(*)
      double precision, intent(in) :: drdi(*)
      integer, intent(in) :: ircut(0:ipan)

      integer :: ien, ip, ist, n
      
      cfint = 0.d0
      do ip = 1, ipan
!---> loop over kinks
        ist = ircut(ip-1) + 1
        ien = ircut(ip)

        cf(ist:ien) = cf(ist:ien)*drdi(ist:ien)

        if (mod(ien-ist, 2) == 0) then
          cfint = cfint + (cf(ist) - cf(ien))/3.d0
          ist = ist + 1
        else
        
!---> four point lagrange integration for the first step
          cfint = cfint + (9.d0*cf(ist) + 19.d0*cf(ist+1) - 5.d0*cf(ist+2) + cf(ist+3))/24.d0 + (cf(ist+1) - cf(ien))/3.d0
          ist = ist + 2
        endif
        n = (ien-ist+1)/2

!---> calculate with an extended 3-point-simpson
        cfint = cfint + a1*csum(n,cf(ist),2) + a2*csum(n,cf(ist+1),2)
      enddo ! ipan
!
      endsubroutine ! csimpk
      
      
      
      double precision function simp3f(f,istart,iend,drdi) result(fint)
        integer, intent(in) :: iend, istart
        double precision, intent(in) :: drdi(*), f(*)
        call simp3(f,fint,istart,iend,drdi)
      endfunction
      
      
      subroutine simp3(f,fint,istart,iend,drdi)
!-----------------------------------------------------------------------
!     this subroutine does an integration from istart to iendof
!     the real function f with an extended 3-point-simpson :
!
!                          r(istart)
!
!                       fint = { f(r') dr'
!
!                           r(iend)
!
!-----------------------------------------------------------------------
      double precision, intent(out) :: fint ! result
      integer, intent(in) :: iend, istart
      double precision, intent(in) :: drdi(*), f(*)

      integer :: i, ist

!---> initialize fint
      if (mod(iend-istart, 2) == 0) then
        fint = f(istart)*drdi(istart)/3.d0
        ist = istart + 1
      else
        fint = ( f(istart+3)*drdi(istart+3) &
                 - 5.d0*f(istart+2)*drdi(istart+2) &
                 + 19.d0*f(istart+1)*drdi(istart+1) &
                 + 9.d0*f(istart)*drdi(istart) )/24.d0 &
               + f(istart+1)*drdi(istart+1)/3.d0
        ist = istart + 2
      endif
      
!---> calculate with an extended 3-point-simpson
      do i = ist, iend-1, 2
        fint = fint + a1*f(i)*drdi(i) + a2*f(i+1)*drdi(i+1)
      enddo ! i
      fint = fint - f(iend)*drdi(iend)/3.d0

      endsubroutine simp3
      
      
      

      
      
      double precision function ssum(n,v,iv)
! **********************************************************************
!        sum up the first n elements of the double precision
!        array v(*) with a stepwidth of iv
! ----------------------------------------------------------------------
      integer, intent(in) :: iv, n
      double precision, intent(in) :: v(*)

      integer :: i, ibot, itop

      if (iv >= 0) then
        ibot = 1
        itop = 1 + (n-1)*iv
      else
        ibot = 1 - (n-1)*iv
        itop = 1
      endif

      ssum = 0.d0
      do i = ibot, itop, iv
        ssum = ssum + v(i)
      enddo ! i
      endfunction
      

      double complex function csum(n,v,iv)
! **********************************************************************
!        sum up the first n elements of the double complex
!        array v(*) with a stepwidth of iv
! ----------------------------------------------------------------------
      integer, intent(in) :: iv, n
      double complex, intent(in) :: v(*)

      integer :: i, ibot, itop

      if (iv >= 0) then
        ibot = 1
        itop = 1 + (n-1)*iv
      else
        ibot = 1 - (n-1)*iv
        itop = 1
      endif

      csum = dcmplx(0.d0, 0.d0)
      do i = ibot, itop, iv
        csum = csum + v(i)
      enddo ! i
      
      endfunction ! csum
      
      
      endmodule Quadrature_mod