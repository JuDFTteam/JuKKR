module mod_impcheck
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine impcheck(atomimp, natomimp, naez, rclsimp, rbasis, bravais, ndim)
    ! **********************************************************************
    ! * Checking the coordinates and site-index assignments of an impurity *
    ! * cluster read in from a file                                        *
    ! * The size of the Bravais lattice to be generated is determined      *
    ! * dynamically using the radius of the input cluster as reference     *
    ! *                                                                    *
    ! * For an input site not belonging to the Bravais lattice the program *
    ! * stops.                                                             *
    ! * A wrong site-index (unit cell) assignment is corrected and the     *
    ! * execution continues.                                               *
    ! **********************************************************************

    use :: mod_getclusnxyz
    implicit none
    ! ..
    ! .. Scalar arguments
    integer :: naez, natomimp, ndim
    ! .. Array arguments
    integer :: atomimp(*)
    real (kind=dp) :: bravais(3, 3), rbasis(3, *), rclsimp(3, *)
    ! ..
    ! .. Local scalars
    integer :: i, iatok, iposok, iq, j, n1, n2, n3, nmax, nmaxz
    real (kind=dp) :: diff, rmaxclus, rmaxgen
    character (len=6) :: strat, strpos
    ! ..
    ! .. Local arrays
    integer :: ain(natomimp), nbr(3)
    real (kind=dp) :: rclsnew(3, natomimp), vec1(3), vec2(3)
    logical :: latom(natomimp), lpos(natomimp), labscord
    real (kind=dp) :: diffmin(natomimp)
    ! ..

    ! ----------------------------------------------------------------------
    ! initialize diffmin array with high value
    diffmin(:) = 1e+5_dp


    ! LABSCORD - cluster coordinates are absolute atomic positions

    labscord = .false.
    j = 0
    do while (j<3 .and. .not. labscord)
      j = j + 1
      if (abs(rclsimp(j,1))>1e-8_dp) labscord = .true.
    end do

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! --> set cluster coordinates in absolute atomic positions
    ! RCLSNEW(3,*)

    do i = 1, natomimp
      call dcopy(3, rclsimp(1,i), 1, rclsnew(1,i), 1)
    end do
    if (.not. labscord) then
      iq = atomimp(1)
      do i = 1, natomimp
        call daxpy(3, 1e0_dp, rbasis(1,iq), 1, rclsnew(1,i), 1)
      end do
    end if

    ! --> determine the maximum radius of the input cluster
    ! this will be then compared to the maximum generated radius
    ! when testing the positions -- setting NMAX for generating the
    ! lattice

    rmaxclus = 0e0_dp
    do i = 2, natomimp
      diff = 0e0_dp
      do j = 1, 3
        diff = diff + (rclsnew(j,i)-rclsnew(j,1))**2
      end do
      diff = sqrt(diff)
      rmaxclus = max(rmaxclus, diff)
    end do
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    do i = 1, 3
      nbr(i) = 0
    end do
    rmaxclus = 3.0_dp*rmaxclus
    ! RMAXCLUS = 1.5*RMAXCLUS
    call getclusnxyz(rmaxclus, bravais, ndim, diff, nbr)
    nmax = max(nbr(1), nbr(2), nbr(3))
    nmaxz = nmax
    if (ndim==2) nmaxz = 0
    rmaxgen = 0e0_dp
    write (1337, *) nmax, nmaxz

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    write (1337, *) 'rbasis of impurity (in cartesian coordinate)'
    do i = 1, natomimp
      write (1337, *)(rclsnew(j,i), j=1, 3)
      ain(i) = atomimp(i)
      lpos(i) = .false.
      latom(i) = .true.
      ! =======================================================================
      do n1 = -nmax, nmax
        do n2 = -nmax, nmax
          do n3 = -nmaxz, nmaxz

            do j = 1, 3
              vec1(j) = real(n1, kind=dp)*bravais(j, 1) + real(n2, kind=dp)*bravais(j, 2) + real(n3, kind=dp)*bravais(j, 3)
            end do

            ! -----------------------------------------------------------------------
            do iq = 1, naez
              diff = 0e0_dp
              do j = 1, 3
                vec2(j) = vec1(j) + rbasis(j, iq)
                diff = diff + vec2(j)**2
              end do
              rmaxgen = max(rmaxgen, sqrt(diff))

              diff = sqrt((rclsnew(1,i)-vec2(1))**2+(rclsnew(2,i)-vec2(2))**2+(rclsnew(3,i)-vec2(3))**2)

              if (diff<=(1e-5_dp)) then
                atomimp(i) = iq
                if (ain(i)/=iq) latom(i) = .false.
                lpos(i) = .true.
                go to 100
              end if

              if (diff<=diffmin(i)) diffmin(i) = diff

            end do
            ! -----------------------------------------------------------------------
          end do
        end do
      end do
      ! =======================================================================
100 end do
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    write (1337, 110) rmaxclus, rmaxgen
    write (1337, 120) 'Input data for impurity sites - consistency check'


    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    iatok = 0
    iposok = 0
    open (58, file='rimp.dat')
    write (58, *) natomimp, ' NATOMIMP'
    do i = 1, natomimp
      if (lpos(i)) then
        write (strpos, '(A6)') 'OK'
      else
        write (strpos, '(A6)') 'neq BL'
        iposok = iposok + 1
        write (*, *) 'minimal difference for atom', i, '=', diffmin(i)
      end if
      write (strat, '(I3,A3)') atomimp(i), ' <?'
      if (.not. latom(i)) then
        iatok = iatok + 1
      else if (lpos(i)) then
        write (strat, '(I3)') atomimp(i)
      end if
      write (1337, 130) i, (rclsimp(j,i), j=1, 3), ain(i), strpos, strat
      write (58, fmt='(I6,3E16.8,I6)') i, (rclsimp(j,i), j=1, 3), ain(i)
    end do
    close (58)
    write (1337, 140)

    if (iposok/=0) then
      write (6, 150)
      stop
    end if

    do i = 1, natomimp
      if ((atomimp(i)>naez) .or. (atomimp(i)==0)) then
        write (6, 180) i, atomimp(i)
        stop
      end if
    end do

    if (iatok/=0) then
      write (1337, 160)
    else
      write (1337, 170)
    end if
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

110 format (12x, 'input-cluster R  : ', f11.6, /, 12x, 'test-cluster  R  : ', f11.6, /)
120 format (13x, 63('-'), /, 15x, a, /, 13x, 63('-'), /, 13x, ' imp |               READ IN DATA         host  |', ' CHECKED DATA ', /, 13x, &
      'index|       x           y           z    site  |', '  pos.   site ', /, 13x, 63('-'))
130 format (13x, i3, 2x, '|', 3(f12.6), 1x, i4, 1x, '|', a6, 2x, a6)
140 format (13x, 63('-'))
150 format (/, 6x, 'ERROR: At least one of your input sites does not', ' belong to the Bravais ', /, 13x, 'lattice (neq BL). ', 'Please check your input file', /)
160 format (13x, 'WARNING: At least one inconsistent assignment of ', 'site indices', /, 13x, '         was found in your input. The program will', ' override the', /, 13x, &
      '         input data. Crosscheck?  ', /, 13x, 63('-'), /)
170 format (13x, 'Your cluster data is consistent', /, 13x, 63('-'), /)
180 format (/, 6x, 'ERROR: Wrong assignment of impurity site ', i3, ' to the unit-cell site ', i3, /, 13x, 'Please check your input file', /)
  end subroutine impcheck

end module mod_impcheck
