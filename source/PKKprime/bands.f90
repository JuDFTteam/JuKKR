!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


program BAND

  use mod_bandstr, only: calc_bandstr

  implicit none
!test
  call calc_bandstr()

end program BAND
