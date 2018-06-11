    Subroutine addviratoms14(linterface, nvirt, naez, naezd, natypd, nemb, &
      nembd, rbasis, lcartesian, bravais, ncls, nineq, refpot, kaoez, noq, &
      nref, rmtrefat, i25)

      Use mod_datatypes
      use mod_DataTypes
      Implicit None
!interface variables
      Logical :: linterface, lcartesian, labscord
      Integer :: naez, ncls, nineq
      Integer :: naezd, nembd
      Integer :: natypd, nref, nvirt
      Integer :: nemb, naclsd
      Integer :: ivir(1000)
      Parameter (naclsd=1000)

      Integer :: i, i1, j
      Real (Kind=dp) :: rbasis(3, *), rbasisold(3, nemb+naezd), &
        rbasissave(3, nemb+naezd)
      Real (Kind=dp) :: rmtrefat(naezd+nembd)
      Integer :: refpot(*), refpotold(naezd+nemb)
      Integer :: noq(*)
      Integer :: kaoez(natypd, *), kaoezold(1, nemb+naezd)
      Real (Kind=dp) :: diff, rmaxclus, vec1(3), vec2(3, naclsd)
      Integer :: nbr(3), nmax, nmaxz, n1, n2, n3, iq
      External :: getclusnxyz

!local variables
      Character (Len=40) :: i25
      Integer :: nrefold
      Integer :: natomimp
      Real (Kind=dp), Allocatable :: ratomimp(:, :)
!      real (kind=dp),allocatable  :: rbasislist(:,:)
      Integer, Allocatable :: atomimp(:)
      Integer :: iatom, ibasis
      Integer :: ierr
!       real (kind=dp)              :: ratomvtest(3)
      Real (Kind=dp) :: rbasisnew(3)
      Real (Kind=dp), Allocatable :: rclsnew(:, :), rbasisnew1(:, :)

      Real (Kind=dp) :: bravais(3, 3)
      Real (Kind=dp), Allocatable :: bravaisinv(:, :)
      Real (Kind=dp) :: tol
      Integer :: ndim
      Integer :: naeznew

      tol = 1.E-5_dp

      Write (1337, *) 'LINTERFACE', linterface
      Write (1337, *) 'NAEZ', naez
      Write (1337, *) 'NAEZD', naezd
      Write (1337, *) 'NEMB', nemb
      Write (1337, *) 'RBASISOLD'
      Do ibasis = 1, naez + nemb
        Write (1337, *) ibasis, rbasis(:, ibasis)
        refpotold(ibasis) = refpot(ibasis)
        kaoezold(1, ibasis) = kaoez(1, ibasis)
        rbasissave(:, ibasis) = rbasis(:, ibasis)
      End Do



!  -----------------------------------------------------------------
!                      read the scoef file
!  -----------------------------------------------------------------

      Open (Unit=32452345, File=i25, Iostat=ierr)
      Write (1337, *) '*', i25, '*'
      Write (1337, *) '*', ierr, '*'
      If (ierr/=0) Stop '[addvirtual] file not found'
      Read (32452345, *) natomimp
      Write (1337, *) 'natomimp', natomimp
      Allocate (ratomimp(3,natomimp))
      Allocate (atomimp(natomimp))
      Do iatom = 1, natomimp
        Read (32452345, *) ratomimp(:, iatom), atomimp(iatom)
        Write (1337, '(A8,I5,A2,3F25.16)') 'IMPATOM ', iatom, ' :', &
          ratomimp(:, iatom)
      End Do

!  -----------------------------------------------------------------
!                      set bulk/surface
!  -----------------------------------------------------------------
      If (linterface) Then
        ndim = 2
        Write (1337, '(23X,A)') 'ADDVIRTUAL : surface geometry mode'
      Else
        ndim = 3
        Write (1337, '(23X,A)') 'ADDVIRTUAL : bulk geometry mode'
      End If
!  -----------------------------------------------------------------
!                      read bravais vectors
!  -----------------------------------------------------------------

! invert bravais vectors
      Allocate (bravaisinv(ndim,ndim))
      bravaisinv = bravais(1:ndim, 1:ndim)
      Call inverse_d1(bravaisinv(:,:), ndim)


      nrefold = 0
      Do ibasis = 1, naez + nemb
        nrefold = max(nrefold, refpotold(ibasis))
      End Do
      Write (1337, *) 'Number of reference potentials is currently', nrefold

! Change basis vectors to cartesian coordinates in order to calculate distances.
      If (.Not. lcartesian) Then

        If (linterface) Then
          Do i = 1, naez + nemb
            Do j = 1, ndim
              rbasisold(j, i) = (rbasis(1,i)*bravais(j,1)+rbasis(2,i)*bravais( &
                j,2))
            End Do
            rbasisold(3, i) = rbasis(3, i)
          End Do
        Else
          Do i = 1, naez + nemb
            Do j = 1, ndim
              rbasisold(j, i) = (rbasis(1,i)*bravais(j,1)+rbasis(2,i)*bravais( &
                j,2)+rbasis(3,i)*bravais(j,3))
            End Do
          End Do
        End If

      Else

        Do i = 1, naez + nemb
          rbasisold(:, i) = rbasis(:, i)
        End Do

      End If


! If the 1st imp. atom in the list is at (0,0,0) then all coordinates are assumed
! relative to the 1st imp atom, otherwise relative to the lattice coords (absolute coords).
      labscord = .False.
      j = 0
      Do j = 1, 3
        If (abs(ratomimp(j,1))>1E-8_dp) labscord = .True.
      End Do

      Allocate (rclsnew(3,natomimp))
      Allocate (rbasisnew1(3,natomimp))
      Do i = 1, natomimp
        Call dcopy(3, ratomimp(1,i), 1, rclsnew(1,i), 1)
      End Do
      If (.Not. labscord) Then
        iq = atomimp(1)
        Do i = 1, natomimp
          Call daxpy(3, 1E0_dp, rbasisold(1,iq), 1, rclsnew(1,i), 1)
        End Do
      End If
      rmaxclus = 0E0_dp
      Do i = 2, natomimp
        diff = 0E0_dp
        Do j = 1, 3
          diff = diff + (rclsnew(j,i)-rclsnew(j,1))**2
        End Do
        diff = sqrt(diff)
        rmaxclus = max(rmaxclus, diff)
      End Do

      nbr(1:3) = 0
      Call getclusnxyz(rmaxclus, bravais, ndim, diff, nbr)
      nmax = max(nbr(1), nbr(2), nbr(3))
      nmaxz = nmax
      If (ndim==2) nmaxz = 0
      iq = 0
      Do n1 = -nmax, nmax
        Do n2 = -nmax, nmax
          Do n3 = -nmaxz, nmaxz

            vec1(1:3) = real(n1, kind=dp)*bravais(1:3, 1) + &
              real(n2, kind=dp)*bravais(1:3, 2) + real(n3, kind=dp)*bravais(1: &
              3, 3)

            Do i1 = 1, naez
              iq = iq + 1
              diff = 0E0_dp
              vec2(1:3, iq) = vec1(1:3) + rbasisold(1:3, i1)
            End Do

          End Do
        End Do
      End Do

      ibasis = 0
      Do i = 1, natomimp
        Do i1 = 1, iq
          diff = sqrt((rclsnew(1,i)-vec2(1,i1))**2+(rclsnew(2,i)-vec2(2, &
            i1))**2+(rclsnew(3,i)-vec2(3,i1))**2)
          If (diff<=(tol)) Go To 100 ! Position is on lattice, do not set as virtual atom          
        End Do
        Call rtobasis(bravais, rclsnew(:,i), rbasisnew, ndim)
        If (linterface) Then
          Do j = 1, 2
            rbasisnew1(j, i) = rbasisnew(1)*bravaisinv(j, 1) + &
              rbasisnew(2)*bravaisinv(j, 2)
          End Do
          rbasisnew1(3, i) = rbasisnew(3)
        Else
          rbasisnew1(1:3, i) = rbasisnew(1)*bravaisinv(1:3, 1) + &
            rbasisnew(2)*bravaisinv(1:3, 2) + rbasisnew(3)*bravaisinv(1:3, 3)
        End If
        Write (1337, *) 'rnew', rbasisnew1(:, i)
        If (i>1) Then
          Do i1 = 1, i - 1
            diff = sqrt((rbasisnew1(1,i)-rbasisnew1(1,i1))**2+(rbasisnew1(2, &
              i)-rbasisnew1(2,i1))**2+(rbasisnew1(3,i)-rbasisnew1(3,i1))**2)
            If (diff<=1E-05_dp) Go To 100
          End Do
        End If
        ibasis = ibasis + 1
        ivir(ibasis) = i
100   End Do
! IBASIS is the number of virtual atoms
      Write (1337, *) 'ibasis', ibasis, (ivir(j), j=1, ibasis)

      If (ibasis+naez>naezd) Then
        Write (*, *) '[addvirtual] naez increased to ', ibasis
        Write (*, *) '[addvirtual] naezd is', naezd
        Write (*, *) '[addvirtual] naeznew > naezd please change naezd'
        Stop 'addvistual'
      Else
        Write (1337, *) 'NAEZ will soon be increased to : ', ibasis + naez
      End If


      nvirt = ibasis
      nineq = nineq + nvirt
      ncls = ncls + nvirt
      naeznew = nvirt + naez
      If (naeznew>naezd) Then
        Write (*, *) '[addvirtual] naez increased to ', naeznew
        Write (*, *) '[addvirtual] naezd is', naezd
        Write (*, *) '[addvirtual] naeznew > naezd please change naezd'
        Stop 'addvirtual'
      End If

      Do i = naeznew + 1, naeznew + nemb ! Added +1 : Phivos 25.7.2014
        rbasis(:, i) = rbasissave(:, i-ibasis)
        refpot(i) = refpotold(i-ibasis)
        kaoez(1, i) = kaoezold(1, i-ibasis)
      End Do
      Do i = naeznew + nemb, naeznew + 1, -1
        rmtrefat(i) = rmtrefat(i-nvirt) ! Shift values of embedded positions in array
      End Do
      Do i = naeznew - nemb, naeznew
        rmtrefat(i) = 1.E-20_dp
      End Do


      Do i = naez + 1, naeznew
        rbasis(:, i) = rbasisnew1(:, ivir(i-naez))
        refpot(i) = nrefold + 1
        noq(i) = 0
        kaoez(1, i) = -1
      End Do
      Do i = 1, naeznew + nemb
        nref = max(nref, refpot(i))
      End Do

!  -----------------------------------------------------------------
!     write out stuff
!  -----------------------------------------------------------------
      Write (1337, *) &
        'addvirtual: List of new basis atoms including virtual atoms'
      Do j = 1, naeznew + nemb
        Write (1337, *) rbasis(:, j)
      End Do
      Write (1337, *) '-------------------------------------------------'

      Write (1337, *) 'naeznew is now ', naeznew
      Write (1337, *) 'setting naez to naeznew'
      naez = naeznew
      Write (1337, *) 'updating rbasis array with virtual basis sites'

      Do ibasis = 1, naeznew + nemb
        Write (1337, *) 'REFPOT', refpot(ibasis)
        Write (1337, *) 'NOQ', noq(ibasis)
        Write (1337, *) 'KAOEZ', kaoez(1, ibasis)
      End Do

!    stop 'end of ADDVIRTUAL'
!   deallocate(rbasislist,ratomimp,atomimp)
      Deallocate (ratomimp, atomimp)
      Deallocate (bravaisinv)
      Deallocate (rclsnew)
      Deallocate (rbasisnew1)
    End Subroutine

    Logical Function vec_in_list(vec, veclist, bound)
! --------------------------
! checks if the vector vec is in the vector list veclist
! in the range of (1,bound)
! --------------------------
      Use mod_datatypes
      Integer :: bound
      Real (Kind=dp) :: vec(3)
      Real (Kind=dp) :: veclist(3, bound)
      Integer :: ilist
      Real (Kind=dp) :: tempvec(3), diff

      vec_in_list = .False.
      Do ilist = 1, bound
        tempvec = vec - veclist(:, ilist)
        diff = sqrt(tempvec(1)**2+tempvec(2)**2+tempvec(3)**2)
        If (diff<10E-5_dp) vec_in_list = .True.
      End Do
    End Function


    Subroutine rtobasis(bravais, rpos, rbasis, ndim)
! --------------------------
! converts a spacial vector rpos to a basis vector rbasis
! such that rbasis = bravais * n with n in [0,1]^ndim
! --------------------------
      Use mod_datatypes
      Implicit None
      Real (Kind=dp), Intent (In) :: bravais(3, 3)
      Real (Kind=dp), Intent (In) :: rpos(3)
      Integer, Intent (In) :: ndim
      Real (Kind=dp), Intent (Out) :: rbasis(3)

      Real (Kind=dp) :: bravais_inv(ndim, ndim)
      Real (Kind=dp) :: ncoeffreal(ndim)
      Integer :: ncoeffint(ndim)
      Integer :: idim

! --------------------------
! first invert the bravais matrix => bravais_inv
! --------------------------
      bravais_inv = bravais(1:ndim, 1:ndim)
      Call inverse_d1(bravais_inv(:,:), ndim)
! --------------------------
! then do n = Bravais* rpos
! --------------------------
      ncoeffreal = 0.0E0_dp
      Do idim = 1, ndim
        ncoeffreal = ncoeffreal + bravais_inv(1:ndim, idim)*rpos(idim)
      End Do

! ------old method ---------
! take the smaller integer value of n => ncoeffint
! --------------------------
      Do idim = 1, ndim
!       ncoeffint(idim) = floor(ncoeffreal(idim))
      End Do

! ------new method---------
! take the smaller integer value of n + 0.5 => ncoeffint
! --------------------------
      Do idim = 1, ndim
        ncoeffint(idim) = floor(ncoeffreal(idim)+0.5_dp)
      End Do
! --------------------------
! do rbasis = rpos - bravias * ncoeffint
! --------------------------
      rbasis = 0.0E0_dp
      Do idim = 1, ndim
        rbasis = rbasis + bravais(:, idim)*ncoeffint(idim)
      End Do
      rbasis = rpos - rbasis
    End Subroutine

    Subroutine inverse_d1(mat, n)
      Use mod_datatypes
      Implicit None
      Integer, Intent (In) :: n
      Real (Kind=dp), Dimension (n, n), Intent (Inout) :: mat
      Real (Kind=dp), Allocatable :: work(:)
      Integer, Allocatable :: ipiv(:)
      Integer :: info

      Allocate (ipiv(n), work(n), Stat=info)
      If (info/=0) Stop &
        'error allocating work arrays in inverse_d1 of addviratom'
      Call dgetrf(n, n, mat, n, ipiv, info)
      If (info/=0) Stop 'inverse_d1: dpotrf failed.'
      Call dgetri(n, mat, n, ipiv, work, n, info)
      If (info/=0) Stop 'inverse_d1: dpotri failed.'
      Deallocate (ipiv, work, Stat=info)
      If (info/=0) Stop &
        'error allocating work arrays in inverse_d1 of addviratom'
    End Subroutine
