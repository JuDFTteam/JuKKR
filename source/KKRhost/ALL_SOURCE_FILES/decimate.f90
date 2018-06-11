!-------------------------------------------------------------------------------
! SUBROUTINE: DECIMATE
!> @brief Decimation method
!> - Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
    Subroutine decimate(gllke, naez, tinvbup, tinvbdown, vacflag, factl, &
      nlbasis, nrbasis)

      Use global_variables
      Use mod_datatypes, Only: dp

      Implicit None

!     .. Parameters ..
!
! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
! *          function, set up in the spin-independent non-relativstic *
! *          (l,m_l)-representation                                   *
! *                                                                   *
! *********************************************************************
!
      Integer, Intent (In) :: naez !< Number of atoms in unit cell
      Integer, Intent (In) :: nlbasis !< Number of basis layers of left host (repeated units)
      Integer, Intent (In) :: nrbasis !< Number of basis layers of right host (repeated units)
      Logical, Dimension (2), Intent (In) :: vacflag
      Complex (Kind=dp), Dimension (lmmaxd, lmmaxd), Intent (In) :: factl
      Complex (Kind=dp), Dimension (lmmaxd, lmmaxd, *), Intent (In) :: tinvbup
      Complex (Kind=dp), Dimension (lmmaxd, lmmaxd, *), &
        Intent (In) :: tinvbdown
      Complex (Kind=dp), Dimension (alm, alm), Intent (Inout) :: gllke
! .. Local Scalars
      Integer :: ldi1, ldi1t, ldi2, ldi2t, lm1, lm2, nlayer, icall
      Integer :: ichck, ihost, ii1, ii2, il1, il2, ip1, ip1t, ip2, ip2t, &
        itermax
      Real (Kind=dp) :: errmax
! .. Local Arrays
      Complex (Kind=dp), Dimension (ndim_slabinv, ndim_slabinv) :: a1, an, b1, bn, c1, cn, x1, &
        xn
! ..
      Data icall/0/
! .. External Subroutines
      External :: bofm, surfgf
! .. Save statement
      Save :: icall, nlayer, itermax, errmax, ichck
! .. External Functions
      Logical :: opt
      External :: opt
! .. Intrinsic Functions
      Intrinsic :: mod
!
      icall = icall + 1
!----------------------------------------------------------------------------
      If (icall==1) Then
        nlayer = naez/nprincd
! Parameters for the "decimation" technique.
        itermax = 300
        errmax = 1.0E-180_dp
        ichck = 1
      End If
!----------------------------------------------------------------------------
      If (.Not. vacflag(1)) Then
!-------------------------------------------------------------------------
! Get the matrix B1
!-------------------------------------------------------------------------
        Call bofm(1, 1, b1, ndim_slabinv, gllke, alm)

! Now Subtract t-mat of left host
        Do ip1 = 1, nprincd
          ihost = mod(ip1-1, nlbasis) + 1
          Do lm1 = 1, lmmaxd
            Do lm2 = 1, lmmaxd
              il1 = lmmaxd*(ip1-1) + lm1
              il2 = lmmaxd*(ip1-1) + lm2
              b1(il1, il2) = (b1(il1,il2)-tinvbup(lm1,lm2,ihost))
            End Do
          End Do
        End Do

        Call bofm(1, 2, c1, ndim_slabinv, gllke, alm)
        Call bofm(2, 1, a1, ndim_slabinv, gllke, alm)

! It performs the 'space decimation' iterative procedure.
        Call surfgf(ndim_slabinv, a1, b1, c1, x1, itermax, errmax, ichck)
! Adds to the matrix GLLKE the elements that couples the
! interface to the two half-spaces.
        Do ip1 = 1, nprincd
          Do ip2 = 1, nprincd
            ii1 = ip1
            ii2 = ip2
            Do lm1 = 1, lmmaxd
              Do lm2 = 1, lmmaxd
                ldi1 = lmmaxd*(ip1-1) + lm1
                il1 = lmmaxd*(ii1-1) + lm1
                ldi2 = lmmaxd*(ip2-1) + lm2
                il2 = lmmaxd*(ii2-1) + lm2
                gllke(il1, il2) = gllke(il1, il2) - x1(ldi1, ldi2)
              End Do
            End Do
          End Do
        End Do
      End If

      If (.Not. vacflag(2)) Then

!  If 'ONEBULK' is activated then it calculates the xn decimated element
!  from the x1 element: this is just in the case of equal bulks on the

        If (.Not. opt('ONEBULK ')) Then

!----------------------------------------------------------------------
! Get the matrix BN
!----------------------------------------------------------------------
          Call bofm(nlayer, nlayer, bn, ndim_slabinv, gllke, alm)

! Now Substract t-mat right host
! Notes : the indexing is easier like that
          Do ip1 = 1, nprincd
            ihost = nrbasis - mod(ip1, nrbasis)
            ihost = mod(ip1-1, nrbasis) + 1
            Do lm1 = 1, lmmaxd
              Do lm2 = 1, lmmaxd
                il1 = lmmaxd*(ip1-1) + lm1
                il2 = lmmaxd*(ip1-1) + lm2
                bn(il1, il2) = (bn(il1,il2)-tinvbdown(lm1,lm2,ihost))
              End Do
            End Do
          End Do

          Call bofm(nlayer, nlayer-1, an, ndim_slabinv, gllke, alm)
          Call bofm(nlayer-1, nlayer, cn, ndim_slabinv, gllke, alm)

! It performs the 'space decimation' iterative procedure.
          Call surfgf(ndim_slabinv, cn, bn, an, xn, itermax, errmax, ichck, lmmaxd)
!
        Else
!
          Do ip1 = 1, nprincd
            Do ip2 = 1, nprincd
              ip1t = (nprincd+1) - ip2
              ip2t = (nprincd+1) - ip1
              Do lm1 = 1, lmmaxd
                Do lm2 = 1, lmmaxd
                  ldi1 = lmmaxd*(ip1-1) + lm1
                  ldi2 = lmmaxd*(ip2-1) + lm2
                  ldi1t = lmmaxd*(ip1t-1) + lm2
                  ldi2t = lmmaxd*(ip2t-1) + lm1
                  xn(ldi1t, ldi2t) = factl(lm1, lm2)*x1(ldi1, ldi2)
                End Do
              End Do
            End Do
          End Do
        End If
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             Added on 1.02.2000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Adds to the matrix GLLKE the elements that couples the
! interface to the two half-spaces.
        Do ip1 = 1, nprincd
          Do ip2 = 1, nprincd
            ii1 = (nlayer-1)*nprincd + ip1
            ii2 = (nlayer-1)*nprincd + ip2
            Do lm1 = 1, lmmaxd
              Do lm2 = 1, lmmaxd
                ldi1 = lmmaxd*(ip1-1) + lm1
                il1 = lmmaxd*(ii1-1) + lm1
                ldi2 = lmmaxd*(ip2-1) + lm2
                il2 = lmmaxd*(ii2-1) + lm2
                gllke(il1, il2) = gllke(il1, il2) - xn(ldi1, ldi2)
              End Do
            End Do
          End Do
        End Do
      End If

      Return

    End Subroutine
