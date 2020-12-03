
!-----------------------------------------------------------------------------------------!
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_patch_intel

  private
  public :: patch_intel


  !-------------------------------------------------------------------------------
  !> Summary: interface to patch_intel.c which make mkl believe it works on a intel CPU
  !> Author: 
  !> Category: 
  !> Deprecated: False 
  !> taken from fleur code, seems to give better performance on AMD hardware
  !> than unpatched MKL or AMD's BLIS+FLAME libraries.
  !-------------------------------------------------------------------------------
contains

  subroutine patch_intel()
    !we try to patch the intel libraries to overwrite determination of 'INTEL' brand
    !otherwise performance on AMD CPUs is bad.
    INTERFACE
      subroutine mkl_patch() BIND(C, name="intel_mkl_patch")
      END subroutine
    END INTERFACE
    INTERFACE
      subroutine cpu_patch() BIND(C, name="intel_cpu_patch")
      END subroutine
    END INTERFACE

    print *,"INTEL PATCH applied"

    call cpu_patch()
    call mkl_patch()

  end subroutine patch_intel

end module mod_patch_intel