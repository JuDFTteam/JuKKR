!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_findgroup
  
  interface findgroup
    module procedure findgroup_host
    module procedure findgroup_pkkprime
  end interface findgroup

  private
  public :: findgroup

contains

  !-------------------------------------------------------------------------------
  !> Summary: Find rotation matrices checking all 48 point group symmetries
  !> Author: 
  !> Category: KKRhost, geometry
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This subroutine finds the rotation matrices that leave the
  !> real lattice unchanged.
  !> input:  bravais(i,j)    true bravais lattice vectors
  !> i = x,y,z ; j = A, B, C (a.u.)
  !> recbv(i,j)      reciprocal basis vectors
  !> rbasis          coordinates of basis atoms
  !> nbasis          number of basis atoms
  !> rsymat          all 64 rotation matrices.
  !> rotname         names for the rotation matrices
  !> output: nsymat          number of rotations that restore the lattice.
  !> ISYMINDEX       index for the symmeties found
  !>
  !> This sub makes all 64 rotations in the basis vectors and bravais
  !> vectors and checks if the new rotated vectror belongs in the
  !> lattice. The proper rotation must bring all vectors to a lattice
  !> vector. Information about the rotations found is printed in the end.
  !> The array ISYMINDEX holds the numbers of the symmetry operations
  !> that are stored in array RSYMAT
  !> 
  !> in case of relativistic calculation: take account of
  !> direction of the magnetic moment specified by (QMTET,QMPHI)
  !> if the PARA(magnetic) flag is set to .FALSE.
  !-------------------------------------------------------------------------------
  subroutine findgroup_host(bravais, recbv, rbasis, nbasis, rsymat, rotname,      &
    isymindex, nsymat, para, qmtet, qmphi, symunitary, krel, naezd, nembd, nsymaxd)

    use :: mod_datatypes, only: dp
    use :: mod_ddet33, only: ddet33
    use :: mod_latvec, only: latvec
    use :: mod_constants, only : pi
    implicit none

    real (kind=dp), parameter :: eps = 1.0e-12_dp
    integer :: krel
    integer :: naezd, nembd, nsymaxd
    ! ..
    integer :: nbasis, nsymat
    integer :: isymindex(nsymaxd)
    real (kind=dp) :: bravais(3, 3), rbasis(3, naezd+nembd)
    real (kind=dp) :: rsymat(64, 3, 3), recbv(3, 3)
    real (kind=dp) :: qmtet(naezd), qmphi(naezd)
    logical :: symunitary(nsymaxd), para
    ! ..
    ! .. Local variables
    real (kind=dp) :: r(3, 4), rotrbas(3, naezd+nembd)
    real (kind=dp) :: bravais1(3, 3)
    integer :: i, j, isym, nsym, i0, ia
    real (kind=dp) :: mdotmp, mvecq(3, naezd), mvecqp(3, naezd)
    real (kind=dp) :: mrotr(3, 3), symdet, summdotmp
    real (kind=dp) :: stet
    character (len=10) :: rotname(64)
    character (len=10) :: char(64)
    logical :: llatbas, lbulk
    ! external functions
    real (kind=dp), external :: ddot


    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, 100)
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

    nsym = 0
    do isym = 1, nsymaxd
      symunitary(isym) = .true.
    end do
    ! - ---------------------------------
    do i = 1, 3
      do j = 1, 3
        bravais1(j, i) = bravais(j, i)
      end do
    end do
    ! Check for surface mode. If so, set bravais1(3,3) very large, so
    ! that only the in-plane symmetries are found.
    ! not checked, be careful of z--> -z!

    lbulk = .true.
    ! Now check the bravais vectors if they have a z component
    if ((abs(bravais(1,3))<eps) .and. (abs(bravais(2,3))<eps) .and. (abs(bravais(3,3))<eps)) then
      lbulk = .false.
    end if

    do isym = 1, 64

      ! --------------------------------- store rotation matrix

      do i = 1, 3
        do j = 1, 3
          mrotr(i, j) = rsymat(isym, i, j)
        end do
      end do

      summdotmp = 0e0_dp

      symdet = ddet33(mrotr)

      ! rotate bravais lattice vectors

      ! In the case of slab/interface geometry look only for
      ! symmetry opperations that preserve the z axis..

      if (lbulk .or. (abs(rsymat(isym,3,3)-1.0_dp)<eps)) then
        ! do rotation only in case bulk or if slab and z axis is restored..


        do i = 1, 3                ! Loop on bravais vectors
          do j = 1, 3              ! Loop on coordinates
            r(j, i) = rsymat(isym, j, 1)*bravais1(1, i) + rsymat(isym, j, 2)*bravais1(2, i) + rsymat(isym, j, 3)*bravais1(3, i)
          end do
        end do

        ! rotate the basis atoms p and take RSYMAT.p - p then
        ! find if R = (RSYMAT.bravais + RSYMAT.p - p) belongs to the
        ! lattice. This is done by function latvec by checking
        ! if R.q = integer (q reciprocal lattice vector)

        llatbas = .true.
        do ia = 1, nbasis          ! Loop on basis atoms
          do j = 1, 3              ! Loop on coordinates
            rotrbas(j, ia) = rsymat(isym, j, 1)*rbasis(1, ia) + rsymat(isym, j, 2)*rbasis(2, ia) + rsymat(isym, j, 3)*rbasis(3, ia)

            rotrbas(j, ia) = rotrbas(j, ia) - rbasis(j, ia)
            r(j, 4) = rotrbas(j, ia)
          end do

          if (.not. latvec(4,recbv,r)) llatbas = .false.

          if ((krel==1) .and. (.not. para)) then
            stet = sin(qmtet(ia)*pi/180e0_dp)
            mvecq(1, ia) = stet*cos(qmphi(ia)*pi/180e0_dp)
            mvecq(2, ia) = stet*sin(qmphi(ia)*pi/180e0_dp)
            mvecq(3, ia) = cos(qmtet(ia)*pi/180e0_dp)

            call dgemv('N', 3, 3, 1e0_dp, mrotr, 3, mvecq(1,ia), 1, 0e0_dp, mvecqp(1,ia), 1)

            call dscal(3, real(symdet,kind=dp), mvecqp(1,ia), 1)

            mdotmp = ddot(3, mvecq(1,ia), 1, mvecqp(1,ia), 1)
            summdotmp = summdotmp + mdotmp

          end if

        end do                     ! ia=1,nbasis

        if ((krel==1) .and. (.not. para)) then
          if (abs(abs(summdotmp)-nbasis)>0.00001e0_dp) then
            llatbas = .false.
          else
            if (summdotmp>0.00001e0_dp) then
              symunitary(nsym+1) = .true.
            else
              symunitary(nsym+1) = .false.
            end if
          end if
        end if

        ! if llatbas=.true. the rotation does not change the lattice

        if (llatbas) then
          nsym = nsym + 1
          isymindex(nsym) = isym
        end if
      end if                       ! (LBULK .OR. (RSYMAT(ISYM,3,3).EQ.1) )
    end do                         ! isym=1,nmatd
    ! nsym symmetries were found
    ! the ISYMINDEX array has the numbers of the symmetries found


    nsymat = nsym

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, '(8X,60("-"))')
    if (lbulk) then
      write (1337, 110)
    else
      write (1337, 120)
    end if
    write (1337, 130) nsymat
    do i = 1, nsymat
      i0 = isymindex(i)
      char(i) = rotname(i0)
    end do
    nsym = nsymat/5
    do i = 1, nsym + 1
      i0 = (i-1)*5
      isym = min(5, nsymat-i0)
      write (1337, 140)(char(j), j=i0+1, i0+isym)
    end do
    write (1337, 150)

100 format (5x, '< FINDGROUP > : Finding symmetry operations', /)
110 format (8x, '3D symmetries:')
120 format (8x, 'surface symmetries:')
130 format (' found for this lattice: ', i2, /, 8x, 60('-'))
140 format (8x, 5(a10,2x))
150 format (8x, 60('-'), /)

  end subroutine findgroup_host

  !-------------------------------------------------------------------------------
  !> Summary: Find symmetry operations that leave crystal lattice invariant
  !> Author: 
  !> Category: PKKprime, geometry
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> @note copied from host code @endnote
  !>
  !> This subroutine finds the rotation matrices that leave the
  !> real lattice unchanged. 
  !> input:  bravais(i,j)    true bravais lattice vectors
  !>                         i = x,y,z ; j = A, B, C (a.u.)
  !>         recbv(i,j)      reciprocal basis vectors
  !>         rbasis          coordinates of basis atoms
  !>         nbasis          number of basis atoms
  !>         rfctor          alat/4/pi
  !> output: rotmat          all 64 rotation matrices.
  !>         rotname         names for the rotation matrices
  !>         nsymat          number of rotations that restore the lattice.
  !>         isymindex       index for the symmeties found
  !>
  !> This sub makes all 64 rotations in the basis vectors and bravais
  !> vectors and checks if the new rotated vectror belongs in the 
  !> lattice. The proper rotation must bring all vectors to a lattice
  !> vector. Information about the rotations found is printed in the end.
  !> The array isymindex holds the numbers of the symmetry operations
  !> that are stored in array RSYMAT
  !-------------------------------------------------------------------------------
  subroutine findgroup_pkkprime( nbasis,naezd,nembd,bravais,rbasis,  &
    & rfctor,recbv,nbzdim,                &
    & rotmat,rotname,nsymat,isymindex_out )

    use :: mod_datatypes, only: dp
    use mod_mympi, only: nranks, myrank, master
    use mod_pointgrp, only : pointgrp
    implicit none

    integer,          intent(in) :: nbasis, naezd, nembd, nbzdim
    real (kind=dp), intent(in) :: bravais(3,3),rbasis(3,naezd+nembd)
    real (kind=dp), intent(in) :: rfctor,recbv(3,3)

    real (kind=dp),     intent(out) :: rotmat(64,3,3)
    character(len=10),    intent(out) :: rotname(64)
    integer,              intent(out) :: nsymat
    integer, allocatable, intent(out) :: isymindex_out(:)

    integer, parameter :: NSYMAXD=48
    integer :: isymindex(NSYMAXD)

    ! Local variables
    real (kind=dp) :: r(3,4),rotrbas(3,naezd+nembd)
    real (kind=dp) :: bravais1(3,3),tol
    integer :: i,j,isym,nsym,i0,ia
    character(len=10) :: charstr(64)
    logical :: llatbas,lbulk
    tol=1e-5
    !     -------------------------------------------------------------
    nsym = 0
    call pointgrp(rotmat,rotname)
    do i=1,3
      do j=1,3
        bravais1(j,i) = bravais(j,i)
      end do
    end do
    !Check for surface mode. If so, set bravais1(3,3) very large, so
    !that only the in-plane symmetries are found. Not checked, be careful of z--> -z!
    if(nbzdim==2)then
      lbulk=.false.
    else!nbzdim==2
      lbulk=.true.
    end if!nbzdim==2
    !     !Now check the bravais vectors if they have a z component
    !     if ((bravais(1,3).eq.0.d0).and.(bravais(2,3).eq.0.d0).and.&
    !        & (bravais(3,3).eq.0.d0)) then
    !        lbulk=.false.
    !     end if

    !     write(100,*) 'bravais:'
    !     write(100,'(3ES25.16)') bravais
    !     write(100,*) 'lbulk=', lbulk

    do isym=1,64
    !rotate bravais lattice vectors

      !In the case of slab/interface geometry look only for
      !symmetry opperations that preserve the z axis..
      if (lbulk .or. abs(rotmat(isym,3,3)-1.0)<tol ) then
      !do rotation only in case bulk or if slab and z axis is restored..
        do i=1,3            ! Loop on bravais vectors
          do j=1,3         ! Loop on coordinates
            r(j,i) = rotmat(isym,j,1)*bravais1(1,i) +   &
                    & rotmat(isym,j,2)*bravais1(2,i) +  &
                    & rotmat(isym,j,3)*bravais1(3,i)
          enddo
        enddo

        !rotate the basis atoms p and take RSYMAT.p - p then
        !find if R = (RSYMAT.bravais + RSYMAT.p - p) belongs to the
        !lattice. This is done by function latvec by checking
        !if R.q = integer (q reciprocal lattice vector)

        llatbas = .true.
        do ia=1,nbasis      ! Loop on basis atoms
          do j=1,3         ! Loop on coordinates
            rotrbas(j,ia) = rotmat(isym,j,1)*rbasis(1,ia) +   &
                           & rotmat(isym,j,2)*rbasis(2,ia) +  &
                           & rotmat(isym,j,3)*rbasis(3,ia)

            rotrbas(j,ia) = rotrbas(j,ia) - rbasis(j,ia)
            r(j,4) = rotrbas(j,ia) 
          enddo
          if (.not.latvec(4,recbv,r)) llatbas=.false.
        enddo               ! ia=1,nbasis

        !if llatbas=.true. the rotation does not change the lattice 
        if (llatbas) then
          nsym = nsym + 1
          isymindex(nsym) = isym
        end if
      end if                 ! (lbulk .OR. (rotmat(isym,3,3).EQ.1) )
    end do                    ! isym=1,nmatd

    !nsym symmetries were found
    !the isymindex array has the numbers of the symmetries found
    nsymat = nsym
    allocate(isymindex_out(nsymat))
    isymindex_out(:) = isymindex(1:nsymat)

    !write info to the screen
    if(myrank==master)then
      write(6,*) 'Information from FindGroup'
      if (.not.lbulk) write(6,*) 'Surface Symmetries '
      write(6,1020) nsymat
      do i=1,nsymat
        I0 = isymindex(i)
        charstr(i) =  rotname(I0)
      end do
      write(6,1010) (charstr(i),i=1,nsymat)
      write(6,*) '----------- * findgroup ends here * ---------------'
      write(6,*)
    end if!myrank==master
    1010 FORMAT(5(A10,2X))
    1020 FORMAT(' Symmetries found for this lattice: ',I5)

  end subroutine findgroup_pkkprime

  !-------------------------------------------------------------------------------
  !> Summary: Checks if a set of vectors are lattice vectors
  !> Author: 
  !> Category: PKKprime, geometry
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> @note copied from host code @endnote
  !> Inputs:                                                              
  !>   n     :number of vectors                                           
  !>   qlat  :primitive translation vectors in reciprocal space           
  !>   vec   :double-precision vector                                     
  !> Outputs:                                                             
  !>   latvec:.true. if all vectors are lattice vectors
  !-------------------------------------------------------------------------------
  logical function latvec(n,qlat,vec) 
    
    use :: mod_datatypes, only: dp
  
    implicit none 
    ! Passed parameters:                                                    
    integer,          intent(in) :: n 
    real (kind=dp), intent(in) :: qlat(3,3),vec(3,n) 
    ! Local parameters:                                                     
    integer :: i,m 
    real (kind=dp) :: vdiff 
    real (kind=dp), parameter :: tol=1.0_dp-6

    latvec=.false. 
    do i=1,n
      do m=1,3
        vdiff=vec(1,i)*qlat(1,m)+vec(2,i)*qlat(2,m)+vec(3,i)*qlat(3,m) 
        vdiff=abs(vdiff-nint(vdiff)) 
        if (vdiff.gt.tol) return 
      enddo 
    enddo 
    latvec=.true. 

  end function latvec

end module mod_findgroup
