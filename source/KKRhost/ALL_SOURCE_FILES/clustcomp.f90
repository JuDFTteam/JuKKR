    Logical Function clustcomp_tb(rcls, irefpot, atom, iat1, ic1, n1, rcls1, &
      n2, iat2, naclsd)
      Use mod_datatypes, Only: dp
!     This function returns true if cluster ic1 is equal to new cluster
!     RCLS        coordinates of all (already found) clusters
!     IC1         First cluster index
!     N1          Number of atoms in IC1 cluster
!     rcls1       coordinates of new cluster
!     n2          number of atoms in ic2 cluster
!     ATOM:       atom-type at a certain position in a cluster
!     IAT1, IAT2: central atoms of 1st,2nd cluster
!     IREFPOT:    Type of reference potential
      Implicit None
      Integer :: naclsd
      Real (Kind=dp) :: rcls(3, naclsd, *), rcls1(3, naclsd)
      Integer :: atom(naclsd, *), irefpot(*)
      Integer :: ic1, n1, n2, iat1, iat2, n, i
      Real (Kind=dp) :: rd, tol
      Logical :: lreflog

      tol = 1.E-5_dp
      clustcomp_tb = .False.
      lreflog = .True.
      If (n1==n2) Then
        rd = 0.E0_dp
        Do n = 1, n1
! compare ref-potential types
          If (irefpot(abs(atom(n,iat1)))/=irefpot(abs(atom(n,iat2)))) &
            lreflog = .False.
          Do i = 1, 3
            rd = rd + abs(rcls(i,n,ic1)-rcls1(i,n)) ! compare coordinares
          End Do
        End Do
        If (abs(rd)<tol .And. lreflog) clustcomp_tb = .True.
      End If
    End Function
