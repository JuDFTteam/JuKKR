!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Determines the number of different atomic pairs in a cluster by symmetry considerations
!> Author: 
!> Determines the number of different atomic pairs in a cluster by symmetry considerations
!> assigning a "shell" pointer (used to set up the GF matrix), to each representative pair.
!------------------------------------------------------------------------------------
module mod_shellgen2k
  use mod_datatypes, only: dp
  private dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Determines the number of different atomic pairs in a cluster by symmetry considerations
  !> Author:
  !> Category: geometry, input-output, KKRhost
  !> Deprecated: False 
  !> Determines the number of different atomic pairs in a cluster by symmetry considerations
  !> assigning a "shell" pointer (used to set up the GF matrix), to each representative pair.
  !-------------------------------------------------------------------------------
  subroutine shellgen2k(icc,natom,rcls,atom,nofgij,iofgij,jofgij,nrot,rsymat,       &
    isymindex,rotname,nshell,ratom,nsh1,nsh2,ish,jsh,ijtabsym,ijtabsh,ijtabcalc,    &
    iprint,nsheld,natomd)
    ! **********************************************************************
    ! *    Determines the number of different atomic pairs in a cluster by *
    ! * symmetry considerations, assigning a "shell" pointer (used to set  *
    ! * up the GF matrix), to each representative pair.                    *
    ! *                                                                    *
    ! * NATOM       number of atoms in the cluster                         *
    ! * RCLS(3,*)   atom positions in the cluster                          *
    ! * ATOM(*)     corresponding site index in the unit cell              *
    ! * NROT        actual number of symmetry operations                   *
    ! * RSYMAT      symmetry operation matrices                            *
    ! * ISYMINDEX   symmetry operation pointer                             *
    ! * NSHELD      dimension parameter ( max number of different shells)  *
    ! * IJTABCALC   flag to calculate the pair (I,J) - 1/0 for YES/NO      *
    ! *             (e.g. for impurity calc IJTABCALC(I,J) = 1 - delta_ij) *
    ! * NOFGIJ      total number of ij pairs (equals number of non-zero    *
    ! *             IJTABCALC elements                                     *
    ! * IOFGIJ      cluster indices i for pair ij                          *
    ! * JOFGIJ                      j for pair ij                          *
    ! *                                                                    *
    ! * NSHELL(0)   number of different shells (ij pairs)                  *
    ! * NSHELL(NS)  number of equivalent pairs in shell NS                 *
    ! * NSH1(NS),                                                          *
    ! * NSH2(NS)    site indices i,j of shell (representative pair) NS     *
    ! * ISH/JSH     cluster indices i,j of all NSHELL(NS) equivalent pairs *
    ! described by shell NS                                  *
    ! * IJTABSH     the index of the representative shell NS for G_ij      *
    ! * IJTABSYM    the index of the symmetry operation which brings G(NS) *
    ! *             into G_ij                                              *
    ! * RATOM(3,NS) diference vector R_i(NS) - R_j(NS)                     *
    ! *                                                                    *
    ! **********************************************************************
    use mod_constants, only: nsymaxd 
    use mod_runoptions, only: calc_exchange_couplings
    implicit none
    ! ..
    ! .. Scalar arguments
    integer, intent(inout) :: icc
    integer, intent(in) :: natom
    integer, intent(in) :: nofgij
    integer, intent(in) :: nrot
    integer, intent(in) :: iprint
    integer, intent(in) :: nsheld
    integer, intent(in) :: natomd

    ! ..
    ! .. Array arguments
    integer, dimension(natom**2),  intent(in) :: ijtabcalc !! Linear pointer, specifying whether the block (i,j) has to be calculated needs set up for ICC=-1, not used for ICC=1
    integer, dimension(natomd),  intent(in) :: atom
    integer, dimension(nsymaxd), intent(in) :: isymindex
    integer, dimension(nofgij),  intent(in) :: iofgij
    integer, dimension(nofgij),  intent(in) :: jofgij
    integer, dimension(natom**2), intent(out) :: ijtabsym
    integer, dimension(natom**2), intent(out) :: ijtabsh
    integer, dimension(nsheld, 2*nsymaxd), intent(out) :: ish
    integer, dimension(nsheld, 2*nsymaxd), intent(out) :: jsh
    integer, dimension(nsheld), intent(inout) :: nsh1
    integer, dimension(nsheld), intent(inout) :: nsh2
    integer, dimension(0:nsheld), intent(inout) :: nshell
    
    real (kind=dp), dimension(3, natomd), intent(in) :: rcls
    real (kind=dp), dimension(64, 3, 3), intent(in) :: rsymat
    real (kind=dp), dimension(3, nsheld), intent(inout) :: ratom

    character (len=10), dimension(64), intent(in) :: rotname
    ! ..
    ! .. Local scalars
    integer :: ai, aj, i, j, k, ns, nsnew, nsgen, id, isym, ii, ij, igij
    real (kind=dp) :: r1, small
    logical :: lfound
    ! ..
    ! .. Local arrays
    real (kind=dp) :: ri(3), rj(3)
    integer, allocatable :: nsh1i(:), nsh2i(:), nshelli(:)
    real (kind=dp), allocatable :: ratomi(:, :)
    ! ..
    ! .. Data statements
    data small/1.0e-10_dp/
    ! ..

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, 120)
    if (iprint>1) call printijtab(natom, natom**2, ijtabcalc)
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

    if (nofgij<=0) then
      write (6, '(A)') '      ICC set to 0'
      write (6, '(A)') '         maybe you should check your input?'
      icc = 0                      ! Bauer Long 2011-10-11
      return
    end if
    allocate (nsh1i(nsheld), nsh2i(nsheld), nshelli(nsheld), stat=ns)
    if (ns/=0) stop '   < shellgen2k > allocate NSHELLI arrays'
    allocate (ratomi(3,nsheld), stat=ns)
    if (ns/=0) stop '   < shellgen2k > allocate RATOMI array'
    ! ======================================================================
    ! --> initialise number of shells found for this cluster, setup the
    ! working arrays NSH1I,NSH2I,NSHELLI,RATOMI and set the number of
    ! new found shells (NSNEW) to zero

    do i = 1, nshell(0)
      nsh1i(i) = nsh1(i)
      nsh2i(i) = nsh2(i)
      nshelli(i) = nshell(i)
      do j = 1, 3
        ratomi(j, i) = ratom(j, i)
      end do
    end do
    nsnew = 0

    ! **********************************************************************
    ! loop over I,J-pairs in cluster
    do igij = 1, nofgij

      ! --> search for a symmetric equivalent pair of atoms, LFOUND takes
      ! on the value false/true if this equivalent pair is found

      i = iofgij(igij)
      j = jofgij(igij)
      ! take only those shells that have atom in them
      if (i/=0 .or. j/=0) then
        ai = atom(i)
        aj = atom(j)

        lfound = .false.
        nsgen = nshell(0) + nsnew

        ! RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
        do id = 1, nrot
          isym = isymindex(id)
          ! ----------------------------------------------------------------------
          do ii = 1, 3
            ri(ii) = rsymat(isym, ii, 1)*rcls(1, i) + rsymat(isym, ii, 2)*rcls(2, i) + rsymat(isym, ii, 3)*rcls(3, i)

            rj(ii) = rsymat(isym, ii, 1)*rcls(1, j) + rsymat(isym, ii, 2)*rcls(2, j) + rsymat(isym, ii, 3)*rcls(3, j)
          end do
          ! ----------------------------------------------------------------------
          ! --> search for an equivalent pair within the already generated
          ! shells (1..NSHELL(0)+NSNEW)
          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ns = 0
          do while ((.not. lfound) .and. (ns<nsgen))
            ns = ns + 1
            ! ----------------------------------------------------------------------
            ! IF ( ( AI.EQ.NSH1I(NS) .AND. AJ.EQ.NSH2I(NS) ).OR.
            ! &              ( AI.EQ.NSH2I(NS) .AND. AJ.EQ.NSH1I(NS) )  ) THEN
            ! Commented out by Phivos Mavropoulos 31 Oct 2008. The problem is that
            ! if (I,J) and (J,I)
            ! are assigned to the same shell, then G(I,J) should be transposed to
            ! obtain G(J,I).
            ! However, this transposition is not performed in account in kkr1b
            ! (subr. tbxccpljij).
            ! There, only the real-space rotations (DSYMLL) are performed to
            ! generate each pair GF from the
            ! representative pair, but the transposition is forgotten. Thus there
            ! are two ways to resolve this:
            ! Either flag the pairs to be transposed, which is is a little faster
            ! but complicated
            ! to program, or do not consider the (I,J) and (J,I) pairs as
            ! belonging to the same shell,
            ! which is done now:
            if ((ai==nsh1i(ns) .and. aj==nsh2i(ns))) then

              r1 = (ri(1)-rj(1)+ratomi(1,ns))**2 + (ri(2)-rj(2)+ratomi(2,ns))**2 + (ri(3)-rj(3)+ratomi(3,ns))**2

              if (r1<small) then
                lfound = .true.
                nshelli(ns) = nshelli(ns) + 1
                ! this check is only used for jij mode since otherwise ish/jsh arrays are not used anywhere
                if (calc_exchange_couplings) then
                  if (nshelli(ns)>2*nsymaxd) stop 'dimension error in shellgen2k: nshelli > 2*nsymaxd'
                  if (ns<=nshell(0)) write (1337, 130) ai, (rcls(ii,i), ii=1, 3), aj, (rcls(ii,j), ii=1, 3), ns
                  ish(ns, nshelli(ns)) = i
                  jsh(ns, nshelli(ns)) = j
                end if
              end if

            end if
            ! ----------------------------------------------------------------------
          end do                   ! while loop ns = 1..nsgen while .not.lfound
          ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        end do                     ! id=1, nrot
        ! RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR

        ! --> if the rotation and the representative pair (shell) that
        ! identify a pair of atoms was found LFOUND=.TRUE. and
        ! the search for a different pair of atoms starts; otherwise
        ! the pair (I,J) requires a new shell

        if (.not. lfound) then
          nsnew = nsnew + 1
          if (nsnew+nshell(0)>nsheld) then
            write (6, 110) 'global', 'NSHELD', nsnew + nshell(0)
            stop
          end if

          nsh1i(nshell(0)+nsnew) = ai
          nsh2i(nshell(0)+nsnew) = aj
          nshelli(nshell(0)+nsnew) = 1
          if (calc_exchange_couplings) then
            ! ish and jsh only needed fot Jij mode
            ish(nshell(0)+nsnew, 1) = i
            jsh(nshell(0)+nsnew, 1) = j
          end if
          do ii = 1, 3
            ratomi(ii, nshell(0)+nsnew) = rcls(ii, j) - rcls(ii, i)
          end do
        end if                     ! .not. lfound

      end if                       ! (i/=0 .or. j/=0) then

    end do                         ! igij = 1, nofgij
    ! **********************************************************************

    ! --> test number of shells

    if (nsnew+nshell(0)>nsheld) then
      write (6, 110) 'global', 'NSHELD', nsnew + nshell(0)
      stop
    end if

    ! --> update the argument arrays

    do i = 1, nshell(0) + nsnew
      nsh1(i) = nsh1i(i)
      nsh2(i) = nsh2i(i)
      nshell(i) = nshelli(i)
      do j = 1, 3
        ratom(j, i) = ratomi(j, i)
      end do
    end do

    nshell(0) = nshell(0) + nsnew
    deallocate (nsh1i, nsh2i, nshelli, ratomi, stat=ns)
    if (ns/=0) stop '   < shellgen2k > deallocate arrays'

    ! **********************************************************************

    ! --> scan once again the shells to find the corresponding symmetry
    ! index bringing GS(1..NSHELL(0)) to Gij.
    ! Setup the tables IJTABSH  assigning (I,J) --> NS
    ! IJTABSYM assigning (I,J) --> ISYM
    ! G_ij = D^\dagger(ISYM) * G(NS) * D(ISYM)

    ! **********************************************************************
    do i = 1, natom
      ai = (i-1)*natom
      do j = 1, natom
        ij = ai + j
        ijtabsh(ij) = 0
        ijtabsym(ij) = 0
      end do
    end do
    ! **********************************************************************
    do i = 1, natom
      ai = atom(i)
      do j = 1, natom
        aj = atom(j)
        ! =======================================================================
        do ii = 1, nshell(0)
          ! -----------------------------------------------------------------------
          if ( ai==nsh1(ii) .and. aj==nsh2(ii) ) then
            do id = 1, nrot
              isym = isymindex(id)

              do k = 1, 3
                ri(k) = rsymat(isym, k, 1)*ratom(1, ii) + rsymat(isym, k, 2)*ratom(2, ii) + rsymat(isym, k, 3)*ratom(3, ii)
              end do

                r1 = (rcls(1,j)-rcls(1,i)-ri(1))**2 + (rcls(2,j)-rcls(2,i)-ri(2))**2 + (rcls(3,j)-rcls(3,i)-ri(3))**2

                if (r1<small) then
                  ij = (i-1)*natom + j
                  ijtabsh(ij) = ii
                  ijtabsym(ij) = id
                  go to 100
                end if
            end do!id
          end if!ai==nsh1(ii) .and. aj==nsh2(ii)
          ! -----------------------------------------------------------------------
100       continue
        end do!ii
        ! =======================================================================
      end do
    end do
    ! ***********************************************************************
    if (iprint<=0) return

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, 140) 'assigned shells and symmetries'
    do i = 1, natom
      ai = (i-1)*natom
      do j = 1, natom
        ij = ai + j
        if (ijtabcalc(ij)>0) write (1337, 150) i, j, ijtabsh(ij), ijtabsym(ij), rotname(ijtabsym(ij))
      end do
    end do
    write (1337, 160)
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

110 format (6x, 'Dimension ERROR: please increase the ', a, ' parameter', /, 6x, a, ' to a value >=', i5, /)
120 format (9x, '< SHELLGEN2K > : assigning representative pairs', ' (shells) ', /, 26x, 'for the off-diagonal elements Gij', /)
130 format (9x, 'INFO: For the atomic sites   I=', i3, ' :', 3f10.6, /, 9x, 29x, 'J=', i3, ' :', 3f10.6, /, 9x, 6x, 'an already generated equivalent shell (', i3, ') was found', /)
140 format (13x, 30('-'), /, 13x, a, /, 13x, 30('-'), /, 13x, ' I ', 1x, ' J ', ' | ', 'shell', 4x, 'isym', /, 13x, 30('-'))
150 format (13x, i3, 1x, i3, ' | ', 1x, i4, 4x, i2, 2x, a10)
160 format (13x, 30('-'), /)
  end subroutine shellgen2k        ! SUBROUTINE SHELLGEN

  !-------------------------------------------------------------------------------
  !> Summary: Print the pointer indicating if a block of the Greens function has to be calculated
  !> Author: 
  !> Category: geometry, input-output, KKRhost
  !> Deprecated: False 
  !> Print the pointer indicating if a block of the Greens function has to be calculated
  !-------------------------------------------------------------------------------
  subroutine printijtab(natom, nofgij, ijtab)
    implicit none
    ! ..
    integer, intent(in) :: natom !! Number of atoms in the cluster
    integer, intent(in) :: nofgij
    integer, dimension(nofgij), intent(in) :: ijtab !! Linear pointer, specifying whether the block (i,j) has to be calculated needs set up for ICC=-1, not used for ICC=1
    ! ..
    integer :: i, j, ij
    integer :: lgmax
    ! ..
    lgmax = 59
    write (1337, 100, advance='no') '  searched for pairs marked with 1 in the table below'
    do j = 1, min(natom+3, lgmax)
      write (1337, '("-")', advance='no')
    end do
    write (1337, *)
    do i = 1, natom
      write (1337, '(14X,I3," | ")', advance='no') i
      ij = (i-1)*natom
      do j = 1, natom
        write (1337, '(I1)', advance='no') ijtab(ij+j)
      end do
      write (1337, *)
    end do
    write (1337, '(13X,6("-"))', advance='no')
    do j = 1, min(natom+3, lgmax)
      write (1337, '("-")', advance='no')
    end do
    write (1337, '(/)')
    ! ...........................................
100 format (13x, 65('-'), /, 18x, a, /, 13x, 65('-'), /, 13x, '   J |', /, 13x, 'I    | 1..NATCLUS', /, 13x, 6('-'))
  end subroutine printijtab

end module mod_shellgen2k
