MODULE MOD_SIMPK

CONTAINS

  !-------------------------------------------------------------------------------
  !> Summary: Modified extended 3-point-Simpson integration with restart at kinks
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This subroutine does an integration up to rcut of an
  !> real function f with an extended 3-point-Simpson :
  !>
  !>                        rcut
  !>                  fint = { f(r') dr'
  !>                         0
  !>
  !> modified for functions with kinks - at each kink the
  !> integration is restarted.
  !>
  !> @warning Input f is destroyed! @endwarning
  !-------------------------------------------------------------------------------
  subroutine simpk(f,fint,ipan,ircut,drdi,ipand)
    implicit none
    !     .. scalar arguments ..
    double precision fint
    integer ipan, ipand
    !     .. array arguments ..
    double precision drdi(*),f(*)
    integer ircut(0:ipand)
    !     .. local scalars ..
    double precision a1,a2
    integer i,ien,ip,ist,n
    !     .. intrinsic functions ..
    intrinsic mod


    a1 = 4.0d0/3.0d0
    a2 = 2.0d0/3.0d0
    fint = 0.0d0
    
    do ip = 1,ipan
    
      !---> loop over kinks
      
      ist = ircut(ip-1) + 1
      ien = ircut(ip)
      
      do  i = ist,ien
        f(i) = f(i)*drdi(i)
      end do
      if (mod(ien-ist,2).eq.0) then
        fint = fint + (f(ist)-f(ien))/3.0d0
        ist = ist + 1
        n = (ien-ist+1)/2

      else
        !---> four point lagrange integration for the first step
        fint = fint + (9.0d0*f(ist)+19.0d0*f(ist+1)-5.0d0*f(ist+2)+ f(ist+3))/24.0d0 + (f(ist+1)-f(ien))/3.0d0
        ist = ist + 2
        n = (ien-ist+1)/2
      end if
      
      !---> calculate with an extended 3-point-simpson
      
      fint = fint + a1*simpk_ssum(n,f(ist),2) + a2*simpk_ssum(n,f(ist+1),2)
    end do
  
  end subroutine simpk

  !-------------------------------------------------------------------------------
  !> Summary: Sum up the first N elements of the double precision array V(*) with a stepwidth of IV
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  double precision function simpk_ssum(n,v,iv)
    !     .. scalar arguments ..
    integer iv,n
    !     .. array arguments ..
    double precision v(*)
    !     .. local scalars ..
    double precision vsum
    integer i,ibot,itop


    if (iv.ge.0) then
      ibot = 1
      itop = 1 + (n-1)*iv

    else
      ibot = 1 - (n-1)*iv
      itop = 1
    end if

    vsum = 0.0d0
    do i = ibot,itop,iv
      vsum = vsum + v(i)
    end do
    simpk_ssum = vsum
  end function simpk_ssum


end module mod_simpk


