logical function clustcomp_tb(rcls, irefpot, atom, iat1, ic1, n1, rcls1, n2, &
  iat2, naclsd)
!     This function returns true if cluster ic1 is equal to new cluster
!     RCLS        coordinates of all (already found) clusters
!     IC1         First cluster index
!     N1          Number of atoms in IC1 cluster
!     rcls1       coordinates of new cluster
!     n2          number of atoms in ic2 cluster
!     ATOM:       atom-type at a certain position in a cluster
!     IAT1, IAT2: central atoms of 1st,2nd cluster
!     IREFPOT:    Type of reference potential
  implicit none
  integer :: naclsd
  real *8 :: rcls(3, naclsd, *), rcls1(3, naclsd)
  integer :: atom(naclsd, *), irefpot(*)
  integer :: ic1, n1, n2, iat1, iat2, n, i
  real *8 :: rd, tol
  logical :: lreflog

  tol = 1.d-5
  clustcomp_tb = .false.
  lreflog = .true.
  if (n1==n2) then
    rd = 0.d0
    do n = 1, n1
      if (irefpot(abs(atom(n,iat1)))/= & ! compare ref-potential types
        irefpot(abs(atom(n,iat2)))) lreflog = .false.
      do i = 1, 3
        rd = rd + abs(rcls(i,n,ic1)-rcls1(i,n)) ! compare coordinares
      end do
    end do
    if (abs(rd)<tol .and. lreflog) clustcomp_tb = .true.
  end if
end function
