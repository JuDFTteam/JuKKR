!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module Exceptions_mod
!-------------------------------------------------------------------------------
!> Summary: Unified treatment of exceptions
!> Author: Paul F Baumeister
!> Category: KKRnano
!-------------------------------------------------------------------------------
  ! this modules purpose is to make the following functionality public
  use Errors_mod, only: die
  use Warnings_mod, only: launch_warning
  use StringHelpers_mod, only: operator(-), operator(+)
  implicit none
  public !
  
endmodule Exceptions_mod
