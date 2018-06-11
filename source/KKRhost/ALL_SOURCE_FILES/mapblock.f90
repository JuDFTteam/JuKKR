    Integer Function mapblock(ie, ie1, ne, iterstep, nodefirst, nodelast)
! **********************************************************************
! *                                                                    *
! *                                                                    *
! *                                                                    *
! **********************************************************************
      Implicit None

!Arguments ..
      Integer :: ie, ie1, ne, iterstep
      Integer :: nodefirst, nodelast

!Locals ..
      Integer :: inc, ip, ipp, iproc, je, ke
      Integer :: iesort(ne), iproce(ne)
! ......................................................................
      ipp = iterstep !         dummy use of argument iterstep
      Do je = ie1, ne
        iesort(je) = je
        iproce(je) = 0
      End Do

      ipp = 0
      Do ip = nodefirst, nodelast
        ipp = ipp + 1
      End Do
! ----------------------------------------------------------------------
      If (ipp>1) Then
        iproc = 0
        inc = 1
        Do je = ie1, ne - 1
          ke = iesort(je)
          iproc = iproc + inc

          If (iproc==ipp) Then
            iproc = 0
            inc = 1
          Else If (iproc==-1) Then
            iproc = 0
            inc = 1
          End If

          iproce(ke) = iproc
        End Do
        mapblock = iproce(ie)
! ----------------------------------------------------------------------
      Else
! ----------------------------------------------------------------------
        mapblock = 0
      End If
! ----------------------------------------------------------------------

      Return
    End Function
