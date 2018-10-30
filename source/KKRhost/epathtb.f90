!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_epathtb

contains

  !-------------------------------------------------------------------------------
  !> Summary: Driver for emery mesh creation
  !> Author: Ph. Mavropoulos, V. Popescu
  !> Category: KKRhost, undefined
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Generating the energy mesh.
  !>
  !> Calls the routine EMESHT once for the valence contour and once for 
  !> the semicore contour.
  !> In the semicore range, -NPOLSEM is used to create a rectangular
  !> contour.
  !>              Ph. Mavropoulos, V. Popescu Juelich/Munich 2004 
  !-------------------------------------------------------------------------------
  subroutine epathtb(ez, df, efermi, npnt, iesemicore, idosemicore, ebotval, emuval, tkval, npolval, n1val, n2val, n3val, ebotsem, emusem, tksem, npolsem, n1sem, n2sem, n3sem, &
    iemxd)

    use :: mod_types, only: t_inc
    use :: mod_datatypes, only: dp
    use :: mod_emesht, only: emesht
    implicit none
    
    integer :: iemxd
    complex (kind=dp) :: ez(*), df(*), ezsemi(iemxd), dfsemi(iemxd)
    complex (kind=dp) :: ezval(iemxd), dfval(iemxd)
    real (kind=dp) :: ebotsem, emusem, tksem, ebotval, emuval, tkval, efermi
    integer :: npolsem, n1sem, n2sem, n3sem
    integer :: npolval, n1val, n2val, n3val
    integer :: iesemicore, npnt, npntsemi, npntval
    integer :: ie, je
    integer :: idosemicore


    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    if (t_inc%i_write>0) then
      write (1337, *)
      write (1337, '(79("="))')
      write (1337, '(20X,A)') 'EPATHTB: generates a complex E contour'
      write (1337, '(79("="))')
      write (1337, *)
    end if
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

    iesemicore = 0
    if (idosemicore==1) then
      if (t_inc%i_write>0) write (1337, 100) 'semi-core contour'
      call emesht(ezsemi, dfsemi, npntsemi, ebotsem, emusem, efermi, tksem, -npolsem, n1sem, n2sem, n3sem, iemxd)
      iesemicore = npntsemi
      if (t_inc%i_write>0) write (1337, 100) 'valence contour'
    end if
    call emesht(ezval, dfval, npntval, ebotval, emuval, efermi, tkval, npolval, n1val, n2val, n3val, iemxd)

    npnt = iesemicore + npntval

    do ie = 1, iesemicore
      ez(ie) = ezsemi(ie)
      df(ie) = dfsemi(ie)
    end do

    do ie = iesemicore + 1, npnt
      je = ie - iesemicore
      ez(ie) = ezval(je)
      df(ie) = dfval(je)
    end do

100 format (7x, '* ', a, /, 7x, 20('-'), /)
  end subroutine epathtb

end module mod_epathtb
