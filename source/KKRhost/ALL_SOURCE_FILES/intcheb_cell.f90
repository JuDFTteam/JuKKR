    Subroutine intcheb_cell(cden, den, rpan_intervall, ipan_intervall, &
      npan_tot, ncheb, irmdnew)
      Use mod_datatypes, Only: dp
!***********************************************************************
! integrate the complex density of states for LM=1 
! gives the total complex charge which is then
! transformed to the xyz component of the magnetic 
! moment
!***********************************************************************
      Implicit None

      Integer :: ncheb, npan_tot, irmdnew
      Integer :: ipan_intervall(0:npan_tot)
      Real (Kind=dp) :: rpan_intervall(0:npan_tot)
      Complex (Kind=dp) :: cden(irmdnew), den
      Integer :: irstart, irstop, ipan
      Real (Kind=dp) :: widthfac
      Complex (Kind=dp) :: int1

      den = (0.0E0_dp, 0.0E0_dp)

      Do ipan = 1, npan_tot
        irstart = ipan_intervall(ipan-1) + 1
        irstop = ipan_intervall(ipan)
        widthfac = 0.5E0_dp*(rpan_intervall(ipan)-rpan_intervall(ipan-1))
        Call intcheb_complex(ncheb, cden(irstart:irstop), int1)
        den = den + int1*widthfac
      End Do

    End Subroutine

    Subroutine intcheb_complex(ncheb, arr1, result1)
      Use mod_datatypes, Only: dp
      Implicit None
      Integer, Intent (In) :: ncheb
      Complex (Kind=dp), Intent (In) :: arr1(0:ncheb)
      Complex (Kind=dp), Intent (Out) :: result1
      Real (Kind=dp) :: pi
      Real (Kind=dp) :: intweight(0:ncheb)
      Integer :: icheb1, icheb2

      pi = 4E0_dp*atan(1E0_dp)
      intweight = 1.0E0_dp
      Do icheb1 = 0, ncheb
        Do icheb2 = 2, ncheb, 2
          intweight(icheb1) = intweight(icheb1) + (-2.0E0_dp/(icheb2**2- &
            1.0E0_dp))*cos(icheb2*pi*(icheb1+0.5E0_dp)/(ncheb+1))
        End Do
        intweight(icheb1) = intweight(icheb1)*2.0E0_dp/(ncheb+1)
      End Do

      result1 = (0.0E0_dp, 0.0E0_dp)
      Do icheb1 = 0, ncheb
        result1 = result1 + intweight(icheb1)*arr1(icheb1)
      End Do

    End Subroutine
