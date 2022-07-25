!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_gfmask

contains

  !-------------------------------------------------------------------------------
  !> Summary: Prepare mask for GF elements that are needed to be computed
  !> Author: 
  !> Date: 29.02.2000
  !> Category: KKRhost, structural-greensfunction, reference-system
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This subroutine prepares the ICHECK matrix that is used for
  !> calculating the proper off-diagonal GF matrix elements ( e.g.
  !> impurity) in case of no full inversion algorithm
  !>
  !> ICHECK(I,J) points to the block (I,J) of the GF matrix having the
  !> size NPRINCD
  !-------------------------------------------------------------------------------
  subroutine gfmask(icheck, icc, invmod, nsh1, nsh2, naez, nshell, naezd, nprincd)

    use :: mod_runoptions, only: print_ickeck, use_cond_LB

    implicit none

    integer :: naezd, nprincd
    integer :: icc, invmod, nlayer, naez, nshell

    integer :: icheck(naezd/nprincd, naezd/nprincd)
    integer :: nsh1(*), nsh2(*)

    integer :: icouple(naezd, naezd)
    integer :: i, j, k, ii, istep1, ilt1, istep2, ilt2, il2, il1, lfchk
    character (len=80) :: fmtchk
    character (len=35) :: invalg(0:3)


    data invalg/'FULL MATRIX                        ', 'BANDED MATRIX (slab)               ', 'BANDED + CORNERS MATRIX (supercell)', 'godfrin module                     '/


    write (1337, 100)

    ! 21.10.2014 GODFRIN Flaviano
    if ((invmod/=0) .and. (mod(naez,nprincd)/=0) .and. (invmod/=3)) then
      write (6, 110) naez, nprincd
      stop
    end if

    write (1337, 120) invalg(invmod)

    nlayer = naez/nprincd
    ! ----------------------------------------------------------- INVMOD = 1
    ! band-diagonal matrix
    if (invmod==1) then
      do i = 1, nlayer
        do j = 1, nlayer
          if (i==j) then
            icheck(i, j) = 1
          else
            icheck(i, j) = 0
          end if
        end do
      end do
    end if
    ! ----------------------------------------------------------- INVMOD = 2
    ! band-diagonal matrix with corners (slab periodic along z)
    if (invmod==2) then
      do i = 1, nlayer
        do j = 1, nlayer
          if ((i==j) .or. (j==nlayer) .or. (i==nlayer)) then
            icheck(i, j) = 1
          else
            icheck(i, j) = 0
          end if
        end do
      end do
    end if
    ! ================================================= INVMOD = 1, ICC <> 0
    ! band-diagonal matrix, off-diagonal G elements needed

    ! --> prepare the matrix ICOUPLE which has 1 in all nn' blocks
    ! (atomic sites) that are needed

    ! ======================================================================
    if ((icc/=0) .and. (invmod==1)) then
      do i = 1, naez
        do j = 1, naez
          icouple(i, j) = 0

          do ii = 1, nshell
            if (((nsh1(ii)==i) .and. (nsh2(ii)==j)) .or. ((nsh1(ii)==j) .and. (nsh2(ii)==i))) icouple(i, j) = 1
          end do
        end do
      end do
      ! cccC
      ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ! cccC                                               conductivity
      ! calculation
      ! ccc         IF (use_cond_LB) THEN
      ! ccc            DO I=1,NLAYER
      ! ccc               DO J=1,NLAYER
      ! ccc                  ICHECK(I,J)=0
      ! ccc               ENDDO
      ! ccc            ENDDO
      ! ccc            DO I=1,NAEZ
      ! ccc               DO J=1,NAEZ
      ! ccc                  ICOUPLE(I,J) = 0
      ! ccc               END DO
      ! ccc            END DO
      ! ccc            DO II=1,NCONDPAIR
      ! ccc               I = IATCONDL(II)
      ! ccc               J = IATCONDR(II)
      ! ccc               ICOUPLE(I,J) = 1
      ! ccc               ICOUPLE(J,I) = 1
      ! ccc            ENDDO
      ! ccc         END IF                 ! Conductivity calculation
      ! cccC
      ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      ! ----------------------------------------------------------------------
      ! Now given the matrix ICOUPLE prepare the matrix ICHECK which has 1 in
      ! all principal-layer blocks that we need -- this will be used in the
      ! matrix inversion
      ! ----------------------------------------------------------------------
      istep1 = 0
      ilt1 = 1
      ! ----------------------------------------------------------------------
      do il1 = 1, naez
        istep1 = istep1 + 1

        if (istep1>nprincd) then
          ilt1 = ilt1 + 1
          istep1 = 1
        end if

        ilt2 = 1
        istep2 = 0
        ! ......................................................................
        do il2 = 1, naez
          istep2 = istep2 + 1

          if (istep2>nprincd) then
            ilt2 = ilt2 + 1
            istep2 = 1
          end if

          if (icouple(il1,il2)==1) icheck(ilt1, ilt2) = 1
        end do
        ! ......................................................................
      end do
      ! ----------------------------------------------------------------------
      ! in the case of calculation of single blocks it has to put the correct
      ! value to ICHECK in order to calculate all the elements also necessary
      ! to calculate that single block          ?????
      ! ----------------------------------------------------------------------
      do j = 1, nlayer

        ! --> loop over the element ICHECK(I,J) with fixed J and I < J

        if (j/=1) then
          do i = 1, j - 1
            if (icheck(i,j)==1) then
              do k = i + 1, j
                icheck(k, j) = 1
              end do
              do k = j, nlayer
                icheck(k, k) = 1
              end do
            end if
          end do
        end if

        if (.not. use_cond_LB) then

          ! --> loop over the element ICHECK(I,J) with fixed J and I > J

          if (j/=nlayer) then
            do i = nlayer, j + 1, -1
              if (icheck(i,j)==1) then
                do k = i - 1, j, -1
                  icheck(k, j) = 1
                end do
              end if
            end do
          end if
        end if
      end do
      ! ----------------------------------------------------------------------
    end if
    ! ======================================================================

    if (print_ickeck) then

      fmtchk = ' '
      lfchk = 1
      do i = 1, min(35, nlayer)
        fmtchk = fmtchk(1:lfchk) // '--'
        lfchk = lfchk + 2
      end do

      write (1337, '(8X,A,/,8X,A)') 'ICHECK matrix :', fmtchk(1:lfchk)
      do i = 1, nlayer
        write (1337, '(9X,35I2)')(icheck(i,j), j=1, min(35,nlayer))
      end do
      write (1337, '(8X,A,/)') fmtchk(1:lfchk)
    end if

100 format (5x, '< GFMASK > : set KKR matrix inversion algorithm', /)
110 format (6x, 'ERROR: Number of sites (NAEZ) =', i3, ' not an integer multiplier', /, 6x, 'of principal layers (NPRINCD) =', i3, /, 6x, 'Use ONLY  full inversion in this case')
120 format (8x, 'INVERSION algorithm used : ', a35, /)
  end subroutine gfmask

end module mod_gfmask
