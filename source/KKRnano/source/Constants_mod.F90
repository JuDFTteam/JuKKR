!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module Constants_mod
!-------------------------------------------------------------------------------
!> Summary: Hard-coded constants used allover the code
!> Author: Paul F Baumeister
!> Category: KKRnano
!-------------------------------------------------------------------------------
implicit none
  public ! all module vars are constant and hence read-only

  double precision, parameter :: pi = 4.d0*atan(1.d0) ! 3.1415926535897932d0

endmodule ! Constants_mod


