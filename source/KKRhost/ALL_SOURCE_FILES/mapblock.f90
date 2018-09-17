module mod_mapblock

contains

  integer function mapblock(ie, ie1, ne, iterstep, nodefirst, nodelast)
    ! **********************************************************************
    ! *                                                                    *
    ! *                                                                    *
    ! *                                                                    *
    ! **********************************************************************
    implicit none

    ! Arguments ..
    integer :: ie, ie1, ne, iterstep
    integer :: nodefirst, nodelast

    ! Locals ..
    integer :: inc, ip, ipp, iproc, je, ke
    integer :: iesort(ne), iproce(ne)
    ! ......................................................................
    ipp = iterstep                 ! dummy use of argument iterstep
    do je = ie1, ne
      iesort(je) = je
      iproce(je) = 0
    end do

    ipp = 0
    do ip = nodefirst, nodelast
      ipp = ipp + 1
    end do
    ! ----------------------------------------------------------------------
    if (ipp>1) then
      iproc = 0
      inc = 1
      do je = ie1, ne - 1
        ke = iesort(je)
        iproc = iproc + inc

        if (iproc==ipp) then
          iproc = 0
          inc = 1
        else if (iproc==-1) then
          iproc = 0
          inc = 1
        end if

        iproce(ke) = iproc
      end do
      mapblock = iproce(ie)
      ! ----------------------------------------------------------------------
    else
      ! ----------------------------------------------------------------------
      mapblock = 0
    end if
    ! ----------------------------------------------------------------------

    return
  end function mapblock

end module mod_mapblock
