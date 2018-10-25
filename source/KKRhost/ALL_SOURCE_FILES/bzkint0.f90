module mod_bzkint0

  private
  public :: bzkint0

contains

  !-------------------------------------------------------------------------------
  !> Summary: Find k-point mesh with symmetries
  !> Author: 
  !> Category: KKRhost, k-points
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !> 
  !> changes for impurity 20/02/2004 -- v.popescu according to n.papanikolaou
  !-------------------------------------------------------------------------------
  subroutine bzkint0(nshell, naez, natyp, noq, rbasis, kaoez, icc, bravais, recbv, atomimp, rsymat, isymindex, nsymat, ifilimp, natomimp, &
    nsh1, nsh2, rclsimp, ratom, ijtabsym, ijtabsh, ijtabcalc, iofgij, jofgij, nofgij, ish, jsh, rrot, dsymll, para, qmtet, qmphi, symunitary, &
    hostimp, intervx, intervy, intervz, ielast, ez, kmesh, maxmesh, maxmshd, nsymaxd, krel, lmaxd, lmmaxd, kpoibz, naezd, natypd, natomimpd, &
    nsheld, nembd)
    use :: mod_datatypes, only: dp
    use :: mod_gfshells, only: gfshells
    use :: mod_crtstar, only: crtstar
    use :: mod_pointgrp, only: pointgrp
    use :: mod_findgroup, only: findgroup
    use :: mod_bzkmesh, only: bzkmesh
    use :: mod_symtaumat, only: symtaumat
    implicit none
    ! .. Parameters ..
    integer :: nsymaxd, krel, lmaxd, lmmaxd
    integer :: kpoibz, naezd, natypd, natomimpd, nsheld, nembd
    ! ..
    ! .. Scalar Arguments ..
    integer :: icc, naez, natomimp, natyp, nsymat, nofgij
    integer :: intervx, intervy, intervz, maxmesh, maxmshd, ielast
    character (len=40) :: ifilimp
    ! ..
    ! .. Array Arguments ..
    complex (kind=dp) :: dsymll(lmmaxd, lmmaxd, nsymaxd), ez(*)
    real (kind=dp) :: bravais(3, 3), ratom(3, nsheld), rbasis(3, naezd+nembd), rclsimp(3, natomimpd), recbv(3, 3), rrot(48, 3, nsheld), rsymat(64, 3, 3)
    integer :: atomimp(natomimpd), isymindex(nsymaxd), kaoez(natypd, naezd+nembd), noq(naezd), kmesh(*), nsh1(*), nsh2(*), nshell(0:nsheld), ijtabsym(*), ijtabsh(*), ijtabcalc(*), &
      iofgij(*), jofgij(*), ish(nsheld, *), jsh(nsheld, *)

    integer :: hostimp(0:natypd)
    ! ..
    ! .. Local Scalars ..
    integer :: i, ishell, iu, iprint
    logical :: lirr
    ! ..
    ! .. Local Arrays ..
    character (len=10) :: rotname(64)
    ! .. magnetisation angles ..
    real (kind=dp) :: qmtet(naezd), qmphi(naezd)
    ! .. unitary/antiunitary symmetry flag
    logical :: symunitary(nsymaxd), para
    ! ..
    ! .. External Functions ..
    logical :: test, opt
    external :: test, opt

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, '(79("="),/,15X,A)') 'BZKINT0: finding symmetry, setting BZ integration'
    write (1337, '(79("="),/)')
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

    call pointgrp(rsymat, rotname)
    call findgroup(bravais, recbv, rbasis, naez, rsymat, rotname, isymindex, nsymat, para, qmtet, qmphi, symunitary, krel, naezd, nembd, nsymaxd)

    lirr = .true.
    iprint = 0
    if (test('TAUSTRUC')) iprint = 2

    ! --> test: full BZ integration
    if (test('fullBZ  ') .or. opt('NEWSOSOL')) then
      nsymat = 1
      lirr = .false.
      write (1337, '(8X,2A,/)') 'Test option < fullBZ > or Run option < NEWSOSOL >: ', ' overriding NSYMAT, generate full BZ k-mesh'
    end if

    ! --> generate BZ k-mesh
    call bzkmesh(intervx, intervy, intervz, maxmesh, lirr, bravais, recbv, nsymat, rsymat, isymindex, symunitary, ielast, ez, kmesh, iprint, krel, kpoibz, maxmshd)

    call symtaumat(rotname, rsymat, dsymll, nsymat, isymindex, symunitary, naezd, lmmaxd, naez, lmaxd+1, krel, iprint, nsymaxd)

    ! Now DSYMLL hold NSYMAT symmetrization matrices
    ! 20.02.2004
    call gfshells(icc, natomimp, nsh1, nsh2, ijtabsym, ijtabsh, ijtabcalc, iofgij, jofgij, nofgij, ish, jsh, nshell, naez, natyp, noq, rbasis, bravais, ifilimp, ratom, rclsimp, &
      nsymat, isymindex, rsymat, kaoez, atomimp, rotname, hostimp, lmaxd, lmmaxd, naezd, natypd, natomimpd, nembd, nsheld)

    ! -->  creates difference vectors RROT for BZ integration in KKRMAT01
    call crtstar(ratom, nshell(0), rsymat, nsymat, isymindex, rrot)

    ! ----------------------------------------------------------------------
    if (iprint>2) then
      do ishell = 1, nshell(0)
        write (1337, fmt='(I4)') ishell
        write (1337, fmt='((I4,3F10.1))')(iu, (rrot(iu,i,ishell),i=1,3), iu=1, nsymat)
      end do
    end if
    ! ----------------------------------------------------------------------

  end subroutine bzkint0

end module mod_bzkint0
