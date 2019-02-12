!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculates the norm of the wavefunctions with full potential and spin orbit coupling.
!> Author:
!> Calculates the norm of the wavefunctions with full potential and
!> spin orbit coupling. This is needed for the normalization of the
!> coefficients c_Lks.
!------------------------------------------------------------------------------------
!> @note 
!> Added mode (can be 0/1) to determine if operator is written out using
!> host/impurity wavefunctions. In case of mode==1 the index array
!> `t_imp%atomimp` is used to determine which position in the host
!> corresponds to each impurity poisition (i.e. layer index)
!> @endnote
!> @warning The gaunt coeffients are stored in index array (see subroutine gaunt)
!> @endwarning
!------------------------------------------------------------------------------------
module mod_normcoeff_so

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculates the norm of the wavefunctions with full potential and spin orbit coupling.
  !> Author: 
  !> Category: spin-orbit-coupling, KKRhost
  !> Deprecated: False 
  !> Calculates the norm of the wavefunctions with full potential and
  !> spin orbit coupling. This is needed for the normalization of the
  !> coefficients c_Lks.
  !-------------------------------------------------------------------------------
  !> @note
  !> Added mode (can be 0/1) to determine if operator is written out using
  !> host/impurity wavefunctions. In case of mode==1 the index array
  !> `t_imp%atomimp` is used to determine which position in the host
  !> corresponds to each impurity poisition (i.e. layer index)
  !> @endnote
  !> @warning The gaunt coeffients are stored in index array (see subroutine gaunt)
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine normcoeff_so(natom,ircut,lmmax,pns,thetas,ntcell,ifunm,ipan,lmsp,ksra, &
    cleb,icleb,iend,drdi,irws,nspoh,mode)
#ifdef CPP_MPI
    use :: mpi
    use :: mod_types, only: t_mpi_c_grid, t_inc, t_imp
#else
    use :: mod_types, only: t_inc, t_imp
#endif
    use :: mod_mympi, only: myrank, master
    use :: mod_datatypes, only: dp
    use :: global_variables, only: korbit, lmmaxd, natypd, ipand, ncleb, lmpotd, irmd, irid, nfund, nspind, ipand
    use :: mod_calc_rho_ll_ss, only: calc_rho_ll_ss
    use :: mod_constants, only: pi, czero

    implicit none

    integer, intent(in) :: iend
    integer, intent(in) :: mode
    integer, intent(in) :: ksra
    integer, intent(in) :: natom
    integer, intent (in) :: lmmax !! (LMAX+1)^2
    integer, dimension(*), intent(in) :: irws   !! R point at WS radius
    integer, dimension(*), intent(in) :: ntcell !! Index for WS cell
    integer, dimension(natypd), intent(in) :: ipan  !! Number of panels in non-MT-region
    integer, dimension(natypd, *), intent(in)       :: lmsp !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
    integer, dimension(0:ipand, natypd), intent(in) :: ircut !! R points of panel borders
    integer, dimension(ncleb, 4), intent(in)        :: icleb !! Pointer array
    integer, dimension(natypd, lmpotd), intent(in)  :: ifunm
    real (kind=dp), dimension(*), intent(in) :: cleb  !! GAUNT coefficients (GAUNT)
    real (kind=dp), dimension(irmd, natypd), intent(in) :: drdi !! Derivative dr/di
    real (kind=dp), dimension(irid, nfund, *), intent(in) :: thetas
    complex (kind=dp), dimension(nspind*lmmax, nspind*lmmax, irmd, 2, natom), intent(in) :: pns
    ! .. Local Scalars ..
    complex (kind=dp) :: norm
    integer :: nspod
    integer :: ir, lm1, lm2, lm1p, lm2p, i1, i1sp1, i1sp2, lmsp1, lmsp2, i2sp1 
    integer :: i2sp2, insra, nsra, nspoh, isigma, i2
    ! MPI stuff
    integer :: ierr, ihelp, i1_start, i1_end
    ! ..Local Arrays..
    complex (kind=dp), dimension(:, :, :), allocatable :: rll_12
    complex (kind=dp), dimension(:, :, :, :), allocatable :: rhod
    complex (kind=dp), dimension(:, :, :, :, :, :, :), allocatable :: rll
    complex (kind=dp), dimension(:, :, :, :, :, :, :), allocatable :: dens
    ! MPI stuff
    complex (kind=dp), dimension(:, :, :, :, :, :, :), allocatable :: work

    nspod = 1 + korbit

    if (t_inc%i_write>0) write (1337, *) 'KSRA', ksra
    if (ksra>=1) then              ! previously this was .GT. which is wrong for kvrel=1
      nsra = 2
    else
      nsra = 1
    end if
    if (t_inc%i_write>0) then
      write (1337, *) 'NSRA', nsra
      write (1337, *) 'LMMAX', lmmax
      write (1337, *) 'lmmaxd', lmmaxd
    end if

    allocate (rll(irmd,lmmax,lmmax,nspoh,nspoh,nspoh,natom))
    allocate (rll_12(irmd,lmmax,lmmax))
    allocate (dens(lmmax,lmmax,nspind,nspind,nspind,nspind,natom))

    rll = czero
    dens = czero

    ! determine MPI work division for loop over atoms
#ifdef CPP_MPI
    i1_start = t_mpi_c_grid%ioff_pt1(t_mpi_c_grid%myrank_ie) + 1
    i1_end = t_mpi_c_grid%ioff_pt1(t_mpi_c_grid%myrank_ie) + t_mpi_c_grid%ntot_pt1(t_mpi_c_grid%myrank_ie)
#else
    i1_start = 1
    i1_end = natom
#endif

    ! rewrite the wavefunctions in RLL arrays of 1,2*LMMAX
    do i1 = i1_start, i1_end
      if (t_inc%i_write>0) write (1337, *) 'ATOM', i1, i1_start, i1_end

      ! use I2 as index to map for mode==1 each impurity position to the corresponding layer index of the host
      if (mode==1) then
        i2 = t_imp%atomimp(i1)
      else                         ! for mode==0 I2 and I1 are the same
        i2 = i1
      end if

      do insra = 1, nsra
        do ir = 1, irmd
          do i1sp1 = 1, nspoh
            do i1sp2 = 1, nspoh
              do lm1 = 1, lmmax
                lmsp1 = (i1sp1-1)*lmmax + lm1
                do lm2 = 1, lmmax
                  lmsp2 = (i1sp2-1)*lmmax + lm2
                  rll(ir, lm2, lm1, i1sp2, i1sp1, insra, i1) =                      &
                    pns(lmsp2, lmsp1, ir, insra, i1)
                end do             ! LM1=1,LMMAX
              end do               ! LM1=1,LMMAX
            end do                 ! ISP1=1,2
          end do                   ! ISP1=1,2
        end do                     ! IR
      end do                       ! INSRA
      ! set up the array R*_L1L2 R_L3L4
      do i1sp1 = 1, nspoh
        do i1sp2 = 1, nspoh
          do i2sp1 = 1, nspoh
            do i2sp2 = 1, nspoh
              do lm1 = 1, lmmax
                do lm2 = 1, lmmax
                  do insra = 1, nsra
                    rll_12 = czero
                    do lm1p = 1, lmmax
                      do lm2p = 1, lmmax
                        do ir = 1, irmd
                          rll_12(ir, lm1p, lm2p) =                                  &
                            conjg(rll(ir,lm1p,lm1,i1sp1,i1sp2,insra,i1))*           &
                            rll(ir, lm2p, lm2, i2sp1, i2sp2, insra, i1)
                        end do     ! IR
                      end do       ! LM2P
                    end do         ! LM1P
                    call calc_rho_ll_ss(lmmax, rll_12, ircut(0:ipand,i2), ipan(i2), &
                      ntcell(i2), thetas, cleb, icleb, iend, ifunm, lmsp, irws(i2), &
                      drdi(:,i2), norm)
                    dens(lm1, lm2, i1sp1, i1sp2, i2sp1, i2sp2, i1) =                &
                      dens(lm1, lm2, i1sp1, i1sp2, i2sp1, i2sp2, i1)+ norm
                  end do           ! NSRA
                end do             ! LM2
              end do               ! LM1
            end do                 ! I2SP2
          end do                   ! I2SP1
        end do                     ! I1SP2
      end do                       ! I1SP1
    end do                         ! I1
    deallocate (rll)
    deallocate (rll_12)

#ifdef CPP_MPI
    ! finally gather DENS on master in case of MPI run
    allocate (work(lmmax,lmmax,nspind,nspind,nspind,nspind,natom), stat=ierr)
    if (ierr/=0) stop 'Error allocating work for MPI comm of DENS in normcoeff_SO'
    ihelp = lmmax*lmmax*nspind*nspind*nspind*nspind*natom
    call mpi_reduce(dens, work, ihelp, mpi_double_complex, mpi_sum, master, t_mpi_c_grid%mympi_comm_ie, ierr)
    if (ierr/=mpi_success) stop 'Error in MPI comm of DENS in normcoeff_SO'
    dens(:, :, :, :, :, :, :) = work(:, :, :, :, :, :, :)
    deallocate (work, stat=ierr)
    if (ierr/=0) stop 'Error deallocating work for MPI comm of DENS in normcoeff_SO'
#endif

    if (myrank==master) then       ! do last part and writeout only on master
      write (*, *) 'collect terms and writeout'
      ! calculate rho
      allocate (rhod(lmmaxd,lmmaxd,natom,4))
      if (nspoh/=1) then
        do isigma = 1, 4
          do i1 = 1, natom
            if (isigma==1) then
              do i1sp1 = 1, nspod
                do i1sp2 = 1, nspod
                  do lm1 = 1, lmmax
                    do lm2 = 1, lmmax
                      rhod((i1sp2-1)*lmmax+lm2, (i1sp1-1)*lmmax+lm1, i1, isigma) =  &
                        dens(lm2, lm1, 1, i1sp2, 1, i1sp1, i1) +                    &
                        dens(lm2, lm1, 2, i1sp2, 2, i1sp1, i1)
                    end do
                  end do
                end do
              end do
            else if (isigma==2) then
              do i1sp1 = 1, nspod
                do i1sp2 = 1, nspod
                  do lm1 = 1, lmmax
                    do lm2 = 1, lmmax
                      rhod((i1sp2-1)*lmmax+lm2, (i1sp1-1)*lmmax+lm1, i1, isigma) =  &
                        dens(lm2, lm1, 2, i1sp2, 1, i1sp1, i1) +                    &
                        dens(lm2, lm1, 1, i1sp2, 2, i1sp1, i1)
                    end do
                  end do
                end do
              end do
            else if (isigma==3) then
              do i1sp1 = 1, nspod
                do i1sp2 = 1, nspod
                  do lm1 = 1, lmmax
                    do lm2 = 1, lmmax
                      rhod((i1sp2-1)*lmmax+lm2, (i1sp1-1)*lmmax+lm1, i1, isigma) =  &
                        -(0d0, 1d0)*(dens(lm2,lm1,2,i1sp2,1,i1sp1,i1)               &
                        -dens(lm2,lm1,1,i1sp2,2,i1sp1,i1))
                    end do
                  end do
                end do
              end do
            else if (isigma==4) then
              do i1sp1 = 1, nspod
                do i1sp2 = 1, nspod
                  do lm1 = 1, lmmax
                    do lm2 = 1, lmmax
                      rhod((i1sp2-1)*lmmax+lm2, (i1sp1-1)*lmmax+lm1, i1, isigma) =  &
                        -dens(lm2, lm1, 1, i1sp2, 1, i1sp1, i1) +                   &
                        dens(lm2, lm1, 2, i1sp2, 2, i1sp1, i1)
                    end do
                  end do
                end do
              end do
            end if
          end do
        end do

        ! write to the file
        if (mode==0) then
          open (unit=12, file='TBkkr_rhod.txt', form='formatted', action='write')
        else                       ! mode == 1
          open (unit=12, file='TBkkr_rhod_imp.txt', form='formatted', action='write')
        end if
        do isigma = 1, 4
          do i1 = 1, natom
            do lm2 = 1, lmmaxd
              do lm1 = 1, lmmaxd
                write (12, '(2ES25.16)') rhod(lm1, lm2, i1, isigma)
              end do
            end do
          end do
        end do
        close (12)
      end if                       ! NSPOH!=1
      deallocate (dens)
      deallocate (rhod)
    end if                         ! (myrank==master)

  end subroutine normcoeff_so

end module mod_normcoeff_so
