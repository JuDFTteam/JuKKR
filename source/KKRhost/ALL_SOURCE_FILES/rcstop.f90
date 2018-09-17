module mod_rcstop

contains

  ! ************************************************************************
  subroutine rcstop(c)
    ! ************************************************************************
    ! .. Scalar Arguments ..

    character (len=8), intent (in) :: c

    ! ..
    print *, 'ERROR: STOP AT POSITION ', c
    stop

  end subroutine rcstop

end module mod_rcstop
