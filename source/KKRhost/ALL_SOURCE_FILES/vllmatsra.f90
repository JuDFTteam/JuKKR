    Subroutine vllmatsra(vll0, vll, rmesh, lmsize, nrmax, nrmaxd, eryd, lmax, &
      lval_in, cmode)

      Use constants
      Use mod_datatypes, Only: dp
!************************************************************************************
! The perubation matrix for the SRA-equations are set up
!************************************************************************************
      Implicit None

      Integer, Intent (In) :: lmax !< Maximum l component in wave function expansion
      Integer, Intent (In) :: nrmax !< NTOTD*(NCHEBD+1)
      Integer, Intent (In) :: nrmaxd
      Integer, Intent (In) :: lmsize
      Integer, Intent (In) :: lval_in
      Complex (Kind=dp), Intent (In) :: eryd
      Character (Len=*), Intent (In) :: cmode

      Real (Kind=dp), Dimension (nrmaxd), Intent (In) :: rmesh
      Complex (Kind=dp), Dimension (lmsize, lmsize, nrmax), &
        Intent (In) :: vll0
! .. Output variables
      Complex (Kind=dp), Dimension (2*lmsize, 2*lmsize, nrmax), &
        Intent (Out) :: vll
! .. Local variables
      Integer :: ilm, lval, mval, ival, ir
      Integer, Dimension (lmsize) :: loflm
      Complex (Kind=dp) :: mass, mass0

!************************************************************************************
! determine the bounds of the matricies to get the lm-expansion and the max. number
! of radial points
!************************************************************************************

!************************************************************************************
! calculate the index array to determine the L value of an LM index
! in case of spin-orbit coupling 2*(LMAX+1)**2 are used instead of (LMAX+1)**2
! the second half refers to the second spin and has the the same L value
!************************************************************************************
      ilm = 0

      If (lmsize==1) Then
        loflm(1) = lval_in
      Else If ((lmax+1)**2==lmsize) Then
        Do lval = 0, lmax
          Do mval = -lval, lval
            ilm = ilm + 1
            loflm(ilm) = lval
          End Do
        End Do
      Else If (2*(lmax+1)**2==lmsize) Then
        Do ival = 1, 2
          Do lval = 0, lmax
            Do mval = -lval, lval
              ilm = ilm + 1
              loflm(ilm) = lval
            End Do
          End Do
        End Do
      Else
        Stop '[vllmatsra] error'
      End If

      vll = czero

      If (cmode=='Ref=0') Then
        vll(1:lmsize, 1:lmsize, :) = vll0 !/cvlight

        Do ir = 1, nrmax
          Do ival = 1, lmsize
            lval = loflm(ival)
            mass = cone + (eryd-vll0(ival,ival,ir))/cvlight**2
            mass0 = cone + eryd/cvlight**2

!************************************************************************************
! Conventional potential matrix
!************************************************************************************

            vll(lmsize+ival, lmsize+ival, ir) = -vll0(ival, ival, ir)/ &
              cvlight**2 ! TEST 9/22/2011
            vll(ival, ival, ir) = vll(ival, ival, ir) + &
              (1.0E0_dp/mass-1.0E0_dp/mass0)*lval*(lval+1)/rmesh(ir)**2

!************************************************************************************
! The pertubation matrix is changed in the following way
!
!     from  / V11  V12 \   to    / V21  V22 \
!           \ V21  V22 /         \-V11 -V12 /
! because of the convention used for the left solution
!************************************************************************************
          End Do !ival

        End Do !ir
      Else If (cmode=='Ref=Vsph') Then
        vll(lmsize+1:2*lmsize, 1:lmsize, :) = vll0
      End If

    End Subroutine
