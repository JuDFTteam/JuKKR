module mod_impcoefs
  
  private
  public :: impcoefs

contains

  !-------------------------------------------------------------------------------
  !> Summary: Construct `hostimp` mapping and write `impurity.coefs` file
  !> Author: 
  !> Category: KKRhost, input-output, geometry
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Writes out the auxiliary file impurity.coefs which is needed for
  !> impurity calculations                                           
  !> Sets up the array HOSTIMP -- also needed for impurity case      
  !>                                adopted from N. Papanikolaou
  !>
  !> @note
  !> - Writeout not neede anymore?
  !> - Check if array `hostimp` is only set here?
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine impcoefs(natomimp, naez, atomimp, rclsimp, nshell, nsh1, nsh2, ratom, nsymat, isymindex, rotname, hostimp, natypd, lmaxd, nsheld, nsize)
    use :: mod_datatypes, only: dp
    implicit none
    ! ..
    ! .. Scalar arguments
    integer :: lmaxd, naez, natomimp, natypd, nsheld, nsize, nsymat
    ! ..
    ! .. Array arguments
    integer :: atomimp(*), hostimp(0:natypd), isymindex(*), nsh1(*), nsh2(*), nshell(0:nsheld)
    real (kind=dp) :: ratom(3, nsheld), rclsimp(3, *)
    character (len=10) :: rotname(*)
    ! ..
    ! .. Local scalars
    integer :: ai, i, ii, j, nb, ndim, nhost, nrep, ns
    real (kind=dp) :: r1
    ! ..
    ! .. Local arrays
    integer :: imphost(naez), nshout(natomimp)
    logical :: exist(naez)
    ! ..

    ! -->  shells around atom icc are prepared for storing the
    ! cluster-gf in subroutine kkrmat in GMATLL(LMMAXD,LMMAXD,*)

    do i = 1, naez
      exist(i) = .false.
    end do
    do i = 1, natomimp
      exist(atomimp(i)) = .true.
    end do

    nhost = 0
    do i = 1, naez
      imphost(i) = 0
      if (exist(i)) then
        nhost = nhost + 1
        imphost(i) = nhost
        hostimp(nhost) = i
      end if
    end do
    hostimp(0) = nhost
    if (nhost/=naez) write (6, 100)

    nrep = 1
    ndim = 1

    do i = 1, natomimp
      nshout(i) = 1
    end do

    open (58, file='impurity.coefs', form='FORMATTED')
    write (58, 110) nrep, natomimp, lmaxd, natomimp, (nshout(i), i=1, natomimp)
    write (58, 120)
    ! -----------------------------------------------------------------------
    do i = 1, natomimp

      r1 = sqrt(rclsimp(1,i)**2+rclsimp(2,i)**2+rclsimp(3,i)**2)

      if (naez==nhost) then
        ai = atomimp(i)
      else
        ai = imphost(atomimp(i))
      end if

      write (58, 130)(rclsimp(j,i), j=1, 3), ai, i, i, r1, atomimp(i)
    end do
    ! -----------------------------------------------------------------------

    nb = 0
    do ns = 1, nshell(0)
      nb = nb + nshell(ns)
    end do

    write (58, 200) nsize, nb
    write (58, 110) ndim
    write (58, 140) nhost
    write (58, 150)(hostimp(i), i=1, nhost)
    write (58, 160) nsymat
    write (58, 170)(rotname(isymindex(i)), i=1, nsymat)
    write (58, 180)
    write (58, 110) nshell(0)
    write (58, 190)(ns, nsh1(ns), nsh2(ns), (ratom(ii,ns),ii=1,3), nshell(ns), sqrt(ratom(1,ns)**2+ratom(2,ns)**2+ratom(3,ns)**2), ns=1, nshell(0))
    close (58)
    ! ======================================================================
100 format (8x, 'WARNING: Some host atoms are missing in the ', 'impurity cluster', /, 8x, '         Indexing will be changed. Check ', 'impurity.coefs file?', /)
110 format (11i5)
120 format ('     Position of Impurity            Host Imp Shell', '   Dist     Host id in Bulk')
130 format (3f12.8, i4, i4, i5, f10.6, i5)
140 format ('Host order, no of host sites: ', i5)
150 format (12i4)
160 format (i5, '    Symmetries for the Bulk')
170 format (5a10)
180 format ('Shells in the reduced format')
190 format (3i5, 3f12.8, i8, f10.5)
200 format (11i20)
  end subroutine impcoefs

end module mod_impcoefs
