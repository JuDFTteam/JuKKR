!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: I/O of the Madelung coefficients
!> Author:
!> For a 2D/3D structures it reads the structure matrices from an unstructured file, and
!> then writes them to a formatted file matrices into a file.
!------------------------------------------------------------------------------------
module mod_madelout
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Reads and writes the Madelung coefficients for a 2D strutcure
  !> Author:
  !> Category: input-output, electrostatics, KKRhost 
  !> Deprecated: False 
  !> Reads the full Madelung coefficients from an unformatted file, and then writes 
  !> the necessary ones to a formatted file.
  !-------------------------------------------------------------------------------
  subroutine madel2out(iprint,naez,lrecamad,lmpotd,nleftoff,nrightoff,nleftall,     &
    nrightall)
    implicit none
    
    integer, intent(in) :: naez       !! Number of atoms in unit cell
    integer, intent(in) :: iprint     !! Printing index controller
    integer, intent(in) :: lmpotd     !! (lpot+1)**2
    integer, intent(in) :: lrecamad   !! Record length for unformatted files
    integer, intent(in) :: nleftoff   !! naez*naez
    integer, intent(in) :: nrightoff  !! nleftoff + naez*nleft*nlbasis
    integer, intent(in) :: nleftall   !! nleft*nlbasis
    integer, intent(in) :: nrightall  !! nright*nrbasis

    ! .. Local variables
    integer :: lfmt, iq1, iq2, lm1, lm2
    real (kind=dp), dimension(6,6) :: smat
    real (kind=dp), dimension(6,200) :: smat1, smat2
    real (kind=dp), dimension(lmpotd,lmpotd) :: avmad
    character (len=80) :: fmt
    integer :: irec

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, 100) 'A(1,1)', min(6, naez), 'avmad.dat'

    fmt = ' '
    lfmt = 0
    do iq1 = 1, min(6, naez)
      fmt = fmt(1:lfmt) // '------------'
      lfmt = lfmt + 12
    end do
    write (1337, '(4X,A,/,8X," Inside the slab ",/,4X,A)')(fmt(1:lfmt), iq1=1, 2)

    open (69, access='direct', recl=lrecamad, file='avmad.unformatted', form='unformatted')
    do iq1 = 1, min(6, naez)
      do iq2 = 1, min(6, naez)
        irec = iq2 + naez*(iq1-1)
        read (69, rec=irec) avmad
        smat(iq1, iq2) = avmad(1, 1)
      end do
      do iq2 = 1, min(200, nleftall)
        irec = iq2 + nleftall*(iq1-1) + nleftoff
        read (69, rec=irec) avmad
        smat1(iq1, iq2) = avmad(1, 1)
      end do
      do iq2 = 1, min(200, nrightall)
        irec = iq2 + nrightall*(iq1-1) + nrightoff
        read (69, rec=irec) avmad
        smat2(iq1, iq2) = avmad(1, 1)
      end do
    end do
    close (69)

    do iq1 = 1, min(6, naez)
      write (1337, 110)(smat(iq1,iq2), iq2=1, min(6,naez))
    end do
    write (1337, '(4X,A,/,8X," Slab - left host",/,4X,A)')(fmt(1:lfmt), iq1=1, 2)
    do iq2 = 1, min(200, nleftall)
      write (1337, 110)(smat1(iq1,iq2), iq1=1, min(6,naez))
    end do
    write (1337, '(4X,A,/,8X," Slab - right host",/,4X,A)')(fmt(1:lfmt), iq1=1, 2)
    do iq2 = 1, min(200, nrightall)
      write (1337, 110)(smat2(iq1,iq2), iq1=1, min(6,naez))
    end do
    write (1337, '(4X,A,/)') fmt(1:lfmt)
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

    if (iprint<3) return

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    open (78, file='avmad.dat', form='formatted')
    write (78, 120) 'A(IQ1,IQ2,LM1,LM2) '

    open (69, access='direct', recl=lrecamad, file='avmad.unformatted', form='unformatted')

    write (78, 150) ' Inside the slab '
    do iq1 = 1, naez
      do iq2 = 1, naez
        write (78, 130) iq1, iq2

        irec = iq2 + naez*(iq1-1)
        read (69, rec=irec) avmad

        do lm1 = 1, lmpotd
          do lm2 = 1, lmpotd
            if (abs(avmad(lm1,lm2))>1e-10_dp) write (78, 140) lm1, lm2, avmad(lm1, lm2)
          end do
        end do
        write (78, '(33("-"))')
      end do
    end do
    write (78, 150) ' Slab - Left Host '
    do iq1 = 1, naez
      do iq2 = 1, nleftall
        write (78, 130) iq1, iq2

        irec = iq2 + nleftall*(iq1-1) + nleftoff
        read (69, rec=irec) avmad

        do lm1 = 1, lmpotd
          do lm2 = 1, lmpotd
            if (abs(avmad(lm1,lm2))>1e-10_dp) write (78, 140) lm1, lm2, avmad(lm1, lm2)
          end do
        end do
        write (78, '(33("-"))')
      end do
    end do
    write (78, 150) ' Slab - Right Host '
    do iq1 = 1, naez
      do iq2 = 1, nrightall
        write (78, 130) iq1, iq2

        irec = iq2 + nrightall*(iq1-1) + nrightoff
        read (69, rec=irec) avmad

        do lm1 = 1, lmpotd
          do lm2 = 1, lmpotd
            if (abs(avmad(lm1,lm2))>1e-10_dp) write (78, 140) lm1, lm2, avmad(lm1, lm2)
          end do
        end do
        write (78, '(33("-"))')
      end do
    end do
    close (69)
    close (78)
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

100 format (5x, 'Madelung potential coefficients ', a, ' up to NAEZ =', i2, /, 5x, '(full matrix (nonzeros) in file ', a, ' if IPRINT.GT.2)')
110 format (4x, 6(d12.4))
120 format (' Madelung potential coefficients ', a, /)
130 format ('IQ1 =', i3, ' IQ2 =', i3, /, 17('-'), /, 5x, '  LM1  LM2  A(LM1,LM2)')
140 format (5x, 2i5, 1p, d18.10)
150 format (70('*'), /, 10x, a, /, 70('*'))
  end subroutine madel2out

  !-------------------------------------------------------------------------------
  !> Summary: Reads and writes the Madelung coefficients for a 3D strutcure
  !> Author:
  !> Category: input-output, electrostatics, KKRhost 
  !> Deprecated: False 
  !> Reads the full Madelung coefficients from an unformatted file, and then writes 
  !> the necessary ones to a formatted file.
  !-------------------------------------------------------------------------------
  subroutine madel3out(iprint, naez, lrecabmad, smat1, smat2, lmpotd)

    implicit none

    integer, intent(in) :: naez       !! Number of atoms in unit cell
    integer, intent(in) :: iprint     !! Printing index controller
    integer, intent(in) :: lmpotd     !! (lpot+1)**2
    integer, intent(in) :: lrecabmad !! Record length for unformatted files
    real (kind=dp), dimension(6,6), intent(in) :: smat1 !! Bulk Madelung coefficients
    real (kind=dp), dimension(6,6), intent(in) :: smat2 !! Bulk Madelung coefficients

    ! .. Local variables
    integer :: irec
    integer :: lfmt, iq1, iq2, lm1, lm2
    character (len=80) :: fmt
    real (kind=dp), dimension(lmpotd) :: bvmad
    real (kind=dp), dimension(lmpotd,lmpotd) :: avmad

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, 100) 'A(1,1)', min(6, naez), 'avmad.dat'

    fmt = ' '
    lfmt = 0
    do iq1 = 1, min(6, naez)
      fmt = fmt(1:lfmt) // '------------'
      lfmt = lfmt + 12
    end do
    write (1337, '(4X,A)') fmt(1:lfmt)
    do iq1 = 1, min(6, naez)
      write (1337, 110)(smat1(iq1,iq2), iq2=1, min(6,naez))
    end do
    write (1337, '(4X,A,/)') fmt(1:lfmt)
    write (1337, 100) 'B(1,1)', min(6, naez), 'bvmad.dat'
    write (1337, '(4X,A)') fmt(1:lfmt)
    do iq1 = 1, min(6, naez)
      write (1337, 110)(smat2(iq1,iq2), iq2=1, min(6,naez))
    end do
    write (1337, '(4X,A,/)') fmt(1:lfmt)
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

    if (iprint<3) return

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    open (78, file='avmad.dat', form='formatted')
    open (79, file='bvmad.dat', form='formatted')
    write (78, 120) 'A(IQ1,IQ2,LM1,LM2) '
    write (79, 120) 'B(IQ1,IQ2,LM)'

    open (69, access='direct', recl=lrecabmad, file='abvmad.unformatted', form='unformatted')
    do iq1 = 1, naez
      do iq2 = 1, naez
        write (78, 130) iq1, iq2
        write (79, 140) iq1, iq2

        irec = iq2 + naez*(iq1-1)
        read (69, rec=irec) avmad, bvmad
        do lm1 = 1, lmpotd
          do lm2 = 1, lmpotd
            if (abs(avmad(lm1,lm2))>1e-10_dp) write (78, 150) lm1, lm2, avmad(lm1, lm2)
          end do
          if (abs(bvmad(lm1))>1e-10_dp) write (79, 160) lm1, bvmad(lm1)
        end do
        write (78, '(33("-"))')
        write (79, '(28("-"))')
      end do
    end do
    close (69)
    close (78)
    close (79)
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

100 format (5x, 'Madelung potential coefficients ', a, ' up to NAEZ =', i2, /, 5x, '(full matrix (nonzeros) in file ', a, ' if IPRINT.GT.2)')
110 format (4x, 6(d12.4))
120 format (' Madelung potential coefficients ', a, /)
130 format ('IQ1 =', i3, ' IQ2 =', i3, /, 17('-'), /, 5x, '  LM1  LM2  A(LM1,LM2)')
140 format ('IQ1 =', i3, ' IQ2 =', i3, /, 17('-'), /, 5x, '   LM  B(LM)')
150 format (5x, 2i5, 1p, d18.10)
160 format (5x, i5, 1p, d18.10)
  end subroutine madel3out

end module mod_madelout
