    Program complexdos
      Use mod_datatypes, Only: dp
      Implicit None
! Principle of DOS here: Two-point contour integration
! for DOS in the middle of the two points. The input DOS
! and energy must be complex. Parameter deltae should be
! of the order of magnitude of eim.


!      <-2*deltae->   _
!           /\        |     DOS=(n(1)+n(2))/2 + (n(1)-n(2))*eim/deltae
!          /  \       |
!        (1)  (2)   2*i*eim=2*i*pi*Kb*Tk
!        /      \     |
!       /        \    |
!------------------------ (Real E axis)
      Integer *4 :: iemaxd, lmaxd
      Parameter (iemaxd=1000, lmaxd=10)
      Integer *4 :: npot, iemax, lmax
      Real (Kind=dp) :: eim, deltae, tk, kb, pi
! Total dos stored in DOS(LMAX+1,IE)
      Complex (Kind=dp) :: dos(0:lmaxd+1, iemaxd), ez(iemaxd)
      Complex (Kind=dp) :: dosnew(0:lmaxd+1)
      Real (Kind=dp) :: temp(2*lmaxd+6)
      Real (Kind=dp) :: enew, ef, ev
      Integer *4 :: ie, ii, ll, iheader
      Character (Len=80) :: text
      Character (Len=45) :: text1
      Character (Len=20) :: text2

      ev = 13.6058_dp

! If only Tk is known, use Bolzmann constant, pi to find eim:
!     Kb=0.6333659D-5
!     pi=3.14159265358979312d0
!     eim=pi*Kb*Tk
      Open (49, File='complex.dos', Form='formatted')
      Open (50, File='new3.dos', Form='formatted')
      Open (51, File='new3_eV_EF.dos', Form='formatted')
      Read (49, *) text2 !dummy readin of header, may be replaced later
      Read (49, *) npot
      Read (49, *) iemax
      Read (49, *) lmax
      If (iemax>iemaxd) Stop 'IEMAX.gt.IEMAXD'
      If (lmax>lmaxd) Stop 'LMAX.gt.LMAXD'


      Do ii = 1, npot
        Write (*, *) 'Reading potential', ii
! Read header:
        Do iheader = 1, 3
          Read (49, 120) text
          Write (50, 120) text
          Write (51, 120) text
          If (iheader==1) Write (*, 120) text
        End Do
        Read (49, Fmt='(A45,F10.6,A20)') text1, ef, text2
        Write (50, Fmt='(A45,F10.6,A20)') text1, ef, text2
        Write (51, Fmt='(A45,F10.6,A20)') text1, ef, text2
        Do iheader = 5, 9
          Read (49, 120) text
          Write (50, 120) text
          Write (51, 120) text
        End Do
! Read dos: (total dos stored at DOS(LMAX+1,IE))
        Do ie = 1, iemax
          Read (49, *) ez(ie), (dos(ll,ie), ll=0, lmax+1)
        End Do
! Compute and write out corrected dos at new (middle) energy points:
        Do ie = 2, iemax - 1
          deltae = real(ez(ie+1)-ez(ie))
          eim = aimag(ez(ie))
          enew = real(ez(ie)) ! Real quantity

          Do ll = 0, lmax + 1
            dosnew(ll) = dos(ll, ie) + 0.5E0_dp*(dos(ll,ie-1)-dos(ll,ie+1))* &
              cmplx(0.E0_dp, eim, kind=dp)/deltae
          End Do

          Write (50, 110) enew, aimag(dosnew(lmax+1)), &
            (aimag(dosnew(ll)), ll=0, lmax)
          Write (51, 110)(enew-ef)*ev, aimag(dosnew(lmax+1))/ev, &
            (aimag(dosnew(ll))/ev, ll=0, lmax)

        End Do ! IE=2,IEMAX-1

        If (ii/=npot) Then
          Read (49, *)
          Write (50, 120) ' '
          Write (51, 120) ' '
        End If

      End Do ! II=1,NPOT

      Close (49)
      Close (50)
      Close (51)

      Stop 'telos'
100   Format (10E12.4)
110   Format (10E14.6)
120   Format (A80)
130   Format (8(2E12.4))
    End Program
