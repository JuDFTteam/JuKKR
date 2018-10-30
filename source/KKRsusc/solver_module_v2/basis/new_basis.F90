  subroutine new_basis(itc,my_rank)
! Driver for basis construction
! All ref wfns for atom ia must have already been provided
  use global

  implicit none

  integer(kind=i4b), intent(in) :: itc
! Mpi  
  integer(kind=i4b), intent(in) :: my_rank
! ----------------------------------------------------------------------
  integer(kind=i4b) :: ia, ib, il, is, nb, nr, nlmsbnew, nsum
  complex(kind=c8b), allocatable :: swap(:,:,:,:,:)

! ********************
! Read basis from file
  if (lbasisfreeze) then
    nlmsbnew = 0
    do ia=1,nasusc
      nsum = 0
      do is=1,nsmax
        do il=0,nlmax
          nb = iwsusc(il,is,ia)
          nsum = nsum + (2*il+1)*nb
          if (nb > 0) call read_wfns(ia,il,is)
          nobasis(il,is,ia) = .false.
        end do
      end do
      nlmsbnew = max(nlmsbnew,nsum)
    end do
! ********************
! Construct new basis
  else
! ********************
    nlmsbnew = 0
    do ia=1,nasusc
      nsum = 0
!     separate basis for each l and spin
      if (ibasis == 1) then
        do il=0,nlmax
          do is=1,nsmax
            call new_basis1(nb,ia,il,is,my_rank)
            if (my_rank == 0) then
              if (itc == 1) write(*,'("New basis1:  ",4i6)') ia, il, is, nb
            end if ! my_rank
            nsum = nsum + nb*(2*il+1)
!           this writes the .wfn file
            if (my_rank == 0) then
              if (lhdio) call save_wfns(ia,il,is)
            end if ! my_rank
          end do
        end do
!     separate basis for each l, same for both spins
      else if (ibasis == 2) then
        do il=0,nlmax
          call new_basis2(nb,ia,il,my_rank)
          if (my_rank == 0 .and. loutbasis) then
            if (itc == 1) write(*,'("New basis2:  ",3i6)') ia, il, nb
          end if ! my_rank
          nsum = nsum + nb*(2*il+1)*nsmax
          do is=1,nsmax
            if (my_rank == 0) then
              if (lhdio) call save_wfns(ia,il,is)
            end if ! my_rank
          end do
        end do
!     same basis for all l and spin
      else if (ibasis == 3) then
        call new_basis3(nb,ia,my_rank)
        if (my_rank == 0) then
          if (itc == 1) write(*,'("New basis3:  ",2i6)') ia, nb
        end if ! my_rank
        nsum = nsum + nb*lmmax*nsmax
        do is=1,nsmax
          do il=0,nlmax
            if (my_rank == 0) then
              if (lhdio) call save_wfns(ia,il,is)
            end if ! my_rank
          end do
        end do
      end if
      nlmsbnew = max(nlmsbnew,nsum)
    end do
! ******
  end if
! ******
  if (ibasismethod == 0) stop 'printing basis functions only'
  if (my_rank == 0 .and. loutbasis) then
    write(*,'("newbasis: nlmsb is changed to",i8,", reallocating arrays")') nlmsbnew
  end if ! my_rank
  nlmsb = nlmsbnew
! Needed for the new SOC version
! Allocate memory for projection coefficients
  call init_gfpartsfit(my_rank)
! If not using I/O
  if (.not.lhdio) then
!   Set up overlaps for the product basis
    call overlaps_gf
!   Set up the angular momentum matrices
    call orbmoment
  end if
! All done!
  end subroutine new_basis
