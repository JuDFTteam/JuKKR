    Subroutine idreals(darry, narry, iprint)
      Use mod_datatypes, Only: dp
      Implicit None

! PARAMETER definitions
      Integer :: nsqr, nmul, divmax
      Parameter (nsqr=7, nmul=5, divmax=15)
      Real (Kind=dp) :: tol
      Parameter (tol=1E-6_dp)

! Dummy arguments
      Integer :: iprint, narry
      Real (Kind=dp) :: darry(narry)

! Local variables
      Real (Kind=dp) :: dabs, dble, dsqrt, dsign
      Integer :: div, i1, i2, idone(narry), imul(nmul), isqr(nsqr)
      Real (Kind=dp) :: dsq, x, xn
      Integer :: iabs, idnint

      Data isqr/2, 3, 5, 6, 7, 8, 10/
      Data imul/3, 7, 11, 13, 17/

! --> mark all numbers as unchecked

      Do i1 = 1, narry
        idone(i1) = 0
      End Do

! --> check darry**2/i integer?, i=1,divmax

      Do div = 1, divmax
        dsq = dble(div)
        Do i2 = 1, narry
          If (idone(i2)==0) Then
            x = darry(i2)*darry(i2)*dsq
            xn = dnint(x)
            If (dabs(x-xn)/dsq<tol .And. xn/=0.E0_dp) Then
              If (iprint>4) Write (1337, 100) dabs(darry(i2)), nint(x), div
              darry(i2) = dsign(1E0_dp, darry(i2))*dsqrt(xn/dsq)
              idone(i2) = 1
            End If
          End If
        End Do
      End Do

! --> check darry/sqrt(n) =?=  i/j
!        n=2,3,5,6,7,8,10      i=1,divmax j=i*n

      Do i1 = 1, nsqr
        Do div = 1, divmax
          dsq = dsqrt(dble(div*div*isqr(i1)))
          Do i2 = 1, narry
            If (idone(i2)==0) Then
              x = darry(i2)*dsq
              xn = dnint(x)
              If (dabs(x-xn)/dsq<tol .And. xn/=0.E0_dp) Then
                If (iprint>4) Write (1337, 110) dabs(darry(i2)), isqr(i1), &
                  iabs(idnint(xn)), iabs(isqr(i1)*div)
                darry(i2) = xn/dsq
                idone(i2) = 1
              End If
            End If
          End Do
        End Do
      End Do

! --> check darry = j/i * n ?
!        n=3,7,11,13,17

      Do i1 = 1, nmul
        Do div = 1, divmax
          dsq = dble(div*imul(i1))
          Do i2 = 1, narry
            If (idone(i2)==0) Then
              x = darry(i2)*dsq
              xn = dnint(x)
              If (dabs(x-xn)/dsq<tol .And. xn/=0.E0_dp) Then
                If (iprint>4) Write (1337, 120) dabs(darry(i2)), imul(i1), &
                  iabs(idnint(xn)), div
                darry(i2) = xn/dsq
                idone(i2) = 1
              End If
            End If
          End Do
        End Do
      End Do
      Return

100   Format (8X, '< IDREALS > : identify ', F12.8, ' as dsqrt(', I3, '/', I3, &
        ')')
110   Format (8X, '< IDREALS > : identify ', F12.8, ' as dsqrt(', I2, ')*', &
        I3, '/', I3)
120   Format (8X, '< IDREALS > : identify ', F12.8, ' as 1/', I2, ' * ', I2, &
        '/', I1)

    End Subroutine
