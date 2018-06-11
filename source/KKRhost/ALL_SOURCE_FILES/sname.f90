! ************************************************************************
    Subroutine sname(name, new, band)
      Use mod_datatypes, Only: dp
! ************************************************************************
!.. scalar arguments
      Integer :: band
      Character (Len=40) :: name, new

!.. locals
      Integer :: i, l, lo
      Character (Len=1) :: ch(50), poi
      Character (Len=10) :: s

      Integer :: length
      External :: length
! ------------------------------------------------------------------------
      poi = '.'
      If (band<0) Then
        lo = log(real(-band,kind=dp))/log(10.0E0_dp) + 1
      Else If (band==0) Then
        lo = 0
      Else
        lo = log(real(band,kind=dp))/log(10.0E0_dp)
      End If

!      write(6,*) 'LO ',lo

      Read (name, Fmt='(255a1)')(ch(i), i=1, 40)
      l = length(ch, 40)
!      write(6,*) 'L  ',l

!      write(6,*) 'CH ',(CH(I),I=1,25)

      Write (s, Fmt='(I10)') band
!      write(6,*) 'S  ',s

      Read (s, Fmt='(255A1)')(ch(i), i=l+1, l+10)
!      write(6,*) 'CH ',(CH(I),I=L+1,L+10)

      Write (new, Fmt='(255A1)')(ch(i), i=1, l), poi, (ch(i), i=l+10-lo, l+10)

      Return
    End Subroutine
