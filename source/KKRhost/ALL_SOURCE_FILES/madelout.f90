! **********************************************************************
    Subroutine madel2out(iprint, naez, lrecamad, lmpotd, nleftoff, nrightoff, &
      nleftall, nrightall)
      Use mod_datatypes, Only: dp
      Implicit None
      Integer :: iprint, naez, lmpotd
      Integer :: lrecamad, nleftoff, nrightoff, nleftall, nrightall

      Integer :: lfmt, iq1, iq2, lm1, lm2
      Real (Kind=dp) :: smat(6, 6)
      Real (Kind=dp) :: smat1(6, 200), smat2(6, 200)
      Real (Kind=dp) :: avmad(lmpotd, lmpotd)
      Character (Len=80) :: fmt
      Integer :: irec

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, 100) 'A(1,1)', min(6, naez), 'avmad.dat'

      fmt = ' '
      lfmt = 0
      Do iq1 = 1, min(6, naez)
        fmt = fmt(1:lfmt) // '------------'
        lfmt = lfmt + 12
      End Do
      Write (1337, '(4X,A,/,8X," Inside the slab ",/,4X,A)')(fmt(1:lfmt), &
        iq1=1, 2)

      Open (69, Access='direct', Recl=lrecamad, File='avmad.unformatted', &
        Form='unformatted')
      Do iq1 = 1, min(6, naez)
        Do iq2 = 1, min(6, naez)
          irec = iq2 + naez*(iq1-1)
          Read (69, Rec=irec) avmad
          smat(iq1, iq2) = avmad(1, 1)
        End Do
        Do iq2 = 1, min(200, nleftall)
          irec = iq2 + nleftall*(iq1-1) + nleftoff
          Read (69, Rec=irec) avmad
          smat1(iq1, iq2) = avmad(1, 1)
        End Do
        Do iq2 = 1, min(200, nrightall)
          irec = iq2 + nrightall*(iq1-1) + nrightoff
          Read (69, Rec=irec) avmad
          smat2(iq1, iq2) = avmad(1, 1)
        End Do
      End Do
      Close (69)

      Do iq1 = 1, min(6, naez)
        Write (1337, 110)(smat(iq1,iq2), iq2=1, min(6,naez))
      End Do
      Write (1337, '(4X,A,/,8X," Slab - left host",/,4X,A)')(fmt(1:lfmt), &
        iq1=1, 2)
      Do iq2 = 1, min(200, nleftall)
        Write (1337, 110)(smat1(iq1,iq2), iq1=1, min(6,naez))
      End Do
      Write (1337, '(4X,A,/,8X," Slab - right host",/,4X,A)')(fmt(1:lfmt), &
        iq1=1, 2)
      Do iq2 = 1, min(200, nrightall)
        Write (1337, 110)(smat2(iq1,iq2), iq1=1, min(6,naez))
      End Do
      Write (1337, '(4X,A,/)') fmt(1:lfmt)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

      If (iprint<3) Return

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Open (78, File='avmad.dat', Form='formatted')
      Write (78, 120) 'A(IQ1,IQ2,LM1,LM2) '

      Open (69, Access='direct', Recl=lrecamad, File='avmad.unformatted', &
        Form='unformatted')

      Write (78, 150) ' Inside the slab '
      Do iq1 = 1, naez
        Do iq2 = 1, naez
          Write (78, 130) iq1, iq2

          irec = iq2 + naez*(iq1-1)
          Read (69, Rec=irec) avmad

          Do lm1 = 1, lmpotd
            Do lm2 = 1, lmpotd
              If (abs(avmad(lm1,lm2))>1E-10_dp) Write (78, 140) lm1, lm2, &
                avmad(lm1, lm2)
            End Do
          End Do
          Write (78, '(33(1H-))')
        End Do
      End Do
      Write (78, 150) ' Slab - Left Host '
      Do iq1 = 1, naez
        Do iq2 = 1, nleftall
          Write (78, 130) iq1, iq2

          irec = iq2 + nleftall*(iq1-1) + nleftoff
          Read (69, Rec=irec) avmad

          Do lm1 = 1, lmpotd
            Do lm2 = 1, lmpotd
              If (abs(avmad(lm1,lm2))>1E-10_dp) Write (78, 140) lm1, lm2, &
                avmad(lm1, lm2)
            End Do
          End Do
          Write (78, '(33(1H-))')
        End Do
      End Do
      Write (78, 150) ' Slab - Right Host '
      Do iq1 = 1, naez
        Do iq2 = 1, nrightall
          Write (78, 130) iq1, iq2

          irec = iq2 + nrightall*(iq1-1) + nrightoff
          Read (69, Rec=irec) avmad

          Do lm1 = 1, lmpotd
            Do lm2 = 1, lmpotd
              If (abs(avmad(lm1,lm2))>1E-10_dp) Write (78, 140) lm1, lm2, &
                avmad(lm1, lm2)
            End Do
          End Do
          Write (78, '(33(1H-))')
        End Do
      End Do
      Close (69)
      Close (78)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

100   Format (5X, 'Madelung potential coefficients ', A, ' up to NAEZ =', I2, &
        /, 5X, '(full matrix (nonzeros) in file ', A, ' if IPRINT.GT.2)')
110   Format (4X, 6(D12.4))
120   Format (' Madelung potential coefficients ', A, /)
130   Format ('IQ1 =', I3, ' IQ2 =', I3, /, 17('-'), /, 5X, &
        '  LM1  LM2  A(LM1,LM2)')
140   Format (5X, 2I5, 1P, D18.10)
150   Format (70('*'), /, 10X, A, /, 70('*'))
    End Subroutine
! **********************************************************************

! **********************************************************************

    Subroutine madel3out(iprint, naez, lrecabmad, smat1, smat2, lmpotd)
      Use mod_datatypes, Only: dp

      Implicit None
      Integer :: iprint, naez, lmpotd
      Integer :: lrecabmad
      Integer :: lfmt, iq1, iq2, lm1, lm2
      Real (Kind=dp) :: smat1(6, 6), smat2(6, 6)
      Real (Kind=dp) :: avmad(lmpotd, lmpotd), bvmad(lmpotd)
      Character (Len=80) :: fmt
      Integer :: irec

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, 100) 'A(1,1)', min(6, naez), 'avmad.dat'

      fmt = ' '
      lfmt = 0
      Do iq1 = 1, min(6, naez)
        fmt = fmt(1:lfmt) // '------------'
        lfmt = lfmt + 12
      End Do
      Write (1337, '(4X,A)') fmt(1:lfmt)
      Do iq1 = 1, min(6, naez)
        Write (1337, 110)(smat1(iq1,iq2), iq2=1, min(6,naez))
      End Do
      Write (1337, '(4X,A,/)') fmt(1:lfmt)
      Write (1337, 100) 'B(1,1)', min(6, naez), 'bvmad.dat'
      Write (1337, '(4X,A)') fmt(1:lfmt)
      Do iq1 = 1, min(6, naez)
        Write (1337, 110)(smat2(iq1,iq2), iq2=1, min(6,naez))
      End Do
      Write (1337, '(4X,A,/)') fmt(1:lfmt)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

      If (iprint<3) Return

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Open (78, File='avmad.dat', Form='formatted')
      Open (79, File='bvmad.dat', Form='formatted')
      Write (78, 120) 'A(IQ1,IQ2,LM1,LM2) '
      Write (79, 120) 'B(IQ1,IQ2,LM)'

      Open (69, Access='direct', Recl=lrecabmad, File='abvmad.unformatted', &
        Form='unformatted')
      Do iq1 = 1, naez
        Do iq2 = 1, naez
          Write (78, 130) iq1, iq2
          Write (79, 140) iq1, iq2

          irec = iq2 + naez*(iq1-1)
          Read (69, Rec=irec) avmad, bvmad
          Do lm1 = 1, lmpotd
            Do lm2 = 1, lmpotd
              If (abs(avmad(lm1,lm2))>1E-10_dp) Write (78, 150) lm1, lm2, &
                avmad(lm1, lm2)
            End Do
            If (abs(bvmad(lm1))>1E-10_dp) Write (79, 160) lm1, bvmad(lm1)
          End Do
          Write (78, '(33(1H-))')
          Write (79, '(28(1H-))')
        End Do
      End Do
      Close (69)
      Close (78)
      Close (79)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

100   Format (5X, 'Madelung potential coefficients ', A, ' up to NAEZ =', I2, &
        /, 5X, '(full matrix (nonzeros) in file ', A, ' if IPRINT.GT.2)')
110   Format (4X, 6(D12.4))
120   Format (' Madelung potential coefficients ', A, /)
130   Format ('IQ1 =', I3, ' IQ2 =', I3, /, 17('-'), /, 5X, &
        '  LM1  LM2  A(LM1,LM2)')
140   Format ('IQ1 =', I3, ' IQ2 =', I3, /, 17('-'), /, 5X, '   LM  B(LM)')
150   Format (5X, 2I5, 1P, D18.10)
160   Format (5X, I5, 1P, D18.10)
    End Subroutine
