!------------------------------------------------------------------------------------
!> Summary: Calculates the KKR matrix elements for the torque operator
!> Author: Guillaume Géranton
!> Calculates the KKR matrix elements for the torque operator, i.e.,
!> \begin{equation}
!> \int dr\left[R^\mu_{Ls} \right]^\dagger T^\mu \left(r\right) R^\mu_{L's'}
!> \end{equation}
!------------------------------------------------------------------------------------
!> @note Details are in http://arxiv.org/pdf/1602.03417v1.pdf
!> This subroutine was adapted `from NORMCOEFF_SO`.
!> @endnote
!------------------------------------------------------------------------------------
module mod_normcoeff_so_torq

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculates the KKR matrix elements for the torque operator
  !> Author: Guillaume Géranton
  !> Category: physical-observables, KKRhost
  !> Deprecated: False 
  !> Calculates the KKR matrix elements for the torque operator, i.e.,
  !> \begin{equation}
  !> \int dr\left[R^\mu_{Ls} \right]^\dagger T^\mu \left(r\right) R^\mu_{L's'}
  !> \end{equation}
  !-------------------------------------------------------------------------------
  !> @note Details are in http://arxiv.org/pdf/1602.03417v1.pdf
  !> This subroutine was adapted from `NORMCOEFF_SO`.
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine normcoeff_so_torq(natom,ircut,lmmax,pns,ntcell,ifunm,ipan,lmsp,ksra,   &
    cleb,icleb,iend,drdi,irws,visp,nspin,vins,irmin,mode)
#ifdef CPP_MPI
    use :: mpi
#endif
    use :: mod_mympi, only: myrank, master
#ifdef CPP_MPI
    use :: mod_types, only: t_mpi_c_grid, t_inc, t_imp
#else
    use :: mod_types, only: t_inc, t_imp
#endif

    use :: mod_datatypes, only: dp
    use :: global_variables
    use :: mod_calc_torq_ll_ss
    use :: mod_constants, only: czero

    implicit none
    real (kind=dp), parameter :: eps = 1.0d-12
    ! ..
    integer, intent(in) :: iend
    integer, intent(in) :: mode
    integer, intent(in) :: ksra
    integer, intent(in) :: natom
    integer, intent(in) :: nspin  !! Counter for spin directions
    integer, intent (in) :: lmmax !! (LMAX+1)^2
    integer, dimension(*), intent(in) :: irws   !! R point at WS radius
    integer, dimension(*), intent(in) :: irmin  !! Max R for spherical treatment
    integer, dimension(*), intent(in) :: ntcell !! Index for WS cell
    integer, dimension(natypd), intent(in) :: ipan  !! Number of panels in non-MT-region
    integer, dimension(natypd, *), intent(in)       :: lmsp !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
    integer, dimension(0:ipand, natypd), intent(in) :: ircut !! R points of panel borders
    integer, dimension(ncleb, 4), intent(in)        :: icleb !! Pointer array
    integer, dimension(natypd, lmpotd), intent(in)  :: ifunm
    real (kind=dp), dimension(*), intent(in) :: cleb  !! GAUNT coefficients (GAUNT)
    real (kind=dp), dimension(irmd, natypd), intent(in) :: drdi !! Derivative dr/di
    real (kind=dp), dimension(irmd, *), intent(in) :: visp !! spherical part of the potential
    real (kind=dp), dimension(irmind:irmd, lmpotd, *), intent(in) ::  vins !! non-spher. part of the potential
    complex (kind=dp), dimension(nspind*lmmax, nspind*lmmax, irmd, 2, natom), intent(in) :: pns
    ! .. Array Arguments ..
    real (kind=dp) :: theta, phi, theta_tmp, phi_tmp
    real (kind=dp), dimension(3) :: sqa
    ! .. Local Scalars ..
    complex (kind=dp) :: norm
    integer :: i2,lm1,lm2,lm1p,lm2p,ir,i1,i1sp1,i1sp2,lmsp1,lmsp2,isigma,i2sp1,i2sp2,insra,nsra
    logical :: lexist
    ! MPI stuff
    integer :: ierr, ihelp, i1_start, i1_end
    ! ..Local Arrays..
    complex (kind=dp), dimension(:, :, :), allocatable :: rll_12
    complex (kind=dp), dimension(:, :, :, :), allocatable :: torq
    complex (kind=dp), dimension(:, :, :, :, :, :, :), allocatable :: rll
    complex (kind=dp), dimension(:, :, :, :, :, :, :), allocatable :: dens
    ! MPI stuff
    complex (kind=dp), dimension(:, :, :, :, :, :, :), allocatable :: work

    if (t_inc%i_write>0) write (1337, *) 'KSRA', ksra
    if (ksra>=1) then              ! previously this was .GT. which is wrong for kvrel=1
      nsra = 2
    else
      nsra = 1
    end if

    if (t_inc%i_write>0) then
      write (1337, *) 'NSRA', nsra
      write (1337, *) 'LMMAX', lmmax
      write (1337, *) 'LMMAXSO', lmmaxso
    end if

    allocate (rll(irmd,lmmax,lmmax,2,2,2,natom))
    allocate (rll_12(irmd,lmmax,lmmax))
    allocate (dens(lmmax,lmmax,2,2,2,2,natom))
    allocate (torq(lmmaxso,lmmaxso,natom,3))

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

          do i1sp1 = 1, 2
            do i1sp2 = 1, 2
              do lm1 = 1, lmmax
                lmsp1 = (i1sp1-1)*lmmax + lm1
                do lm2 = 1, lmmax
                  lmsp2 = (i1sp2-1)*lmmax + lm2
                  rll(ir, lm2, lm1, i1sp2, i1sp1, insra, i1) = pns(lmsp2, lmsp1, ir, insra, i1)
                end do             ! LM1=1,LMMAX
              end do               ! LM1=1,LMMAX
            end do                 ! ISP1=1,2
          end do                   ! ISP1=1,2

        end do                     ! IR
      end do                       ! INSRA


      ! set up the array R*_L1L2 R_L3L4

      do i1sp1 = 1, 2
        do i1sp2 = 1, 2
          do i2sp1 = 1, 2
            do i2sp2 = 1, 2

              do lm1 = 1, lmmax
                do lm2 = 1, lmmax

                  do insra = 1, nsra

                    rll_12 = czero

                    do lm1p = 1, lmmax
                      do lm2p = 1, lmmax

                        do ir = 1, irmd
                          rll_12(ir, lm1p, lm2p) = conjg(rll(ir,lm1p,lm1,i1sp1,i1sp2,insra,i1))*rll(ir, lm2p, lm2, i2sp1, i2sp2, insra, i1)
                        end do     ! IR

                      end do       ! LM2P
                    end do         ! LM1P

                    call calc_torq_ll_ss(lmmax, rll_12, ircut(0:ipand,i2), ipan(i2), ntcell(i2), cleb, icleb, iend, ifunm, lmsp, irws(i2), drdi(:,i2), norm, visp, nspin, i1, vins, &
                      irmin(i2))

                    dens(lm1, lm2, i1sp1, i1sp2, i2sp1, i2sp2, i1) = dens(lm1, lm2, i1sp1, i1sp2, i2sp1, i2sp2, i1) + norm
                  end do           ! NSRA

                end do             ! LM2
              end do               ! LM1

            end do                 ! I2SP2
          end do                   ! I2SP1
        end do                     ! I1SP2
      end do                       ! I1SP1

    end do                         ! I1

#ifdef CPP_MPI
    ! finally gather DENS on master in case of MPI run
    allocate (work(lmmax,lmmax,2,2,2,2,natom), stat=ierr)
    if (ierr/=0) stop 'Error allocating work for MPI comm of DENS in normcoeff_torq'
    ihelp = lmmax*lmmax*2*2*2*2*natom
    call mpi_reduce(dens, work, ihelp, mpi_double_complex, mpi_sum, master, t_mpi_c_grid%mympi_comm_ie, ierr)
    if (ierr/=mpi_success) stop 'Error in MPI comm of DENS in normcoeff_torq'
    dens(:, :, :, :, :, :, :) = work(:, :, :, :, :, :, :)
    deallocate (work, stat=ierr)
    if (ierr/=0) stop 'Error deallocating work for MPI comm of DENS in normcoeff_torq'
#endif


    if (myrank==master) then       ! do last part and writeout only on master

      write (*, *) 'collect terms and writeout'


      ! reads spin quantization axis from file
      inquire (file='nonco_angle.dat', exist=lexist)
      if (lexist) then
        open (unit=11, file='nonco_angle.dat', form='FORMATTED')

        do i1 = 1, natom
          read (11, *) theta_tmp, phi_tmp
          if (i1>1) then
            if ((abs(theta_tmp-theta)>eps) .or. (abs(phi_tmp-phi)>eps)) stop 'in a non-colinear system this is not implemented yet.'
          end if
          theta = theta_tmp
          phi = phi_tmp
        end do
        close (11)
      else
        theta = 0.0_dp
        phi = 0.0_dp
      end if

      sqa(1) = sin(theta)*cos(phi)
      sqa(2) = sin(theta)*sin(phi)
      sqa(3) = cos(theta)

      torq = czero

      do isigma = 1, 3             ! ISIGMA == 1 --> T_x
        ! ISIGMA == 2 --> T_y
        ! ISIGMA == 3 --> T_z

        write (6, *) 'ISIGMA', isigma
        do i1 = 1, natom

          ! TEMPORARY IMPLEMENTATION OF THE TORQUE OPERATOR FOR B // z
          if (isigma==1) then      ! T_x

            do i1sp1 = 1, 2
              do i1sp2 = 1, 2
                do lm1 = 1, lmmax
                  do lm2 = 1, lmmax
                    torq((i1sp2-1)*lmmax+lm2, (i1sp1-1)*lmmax+lm1, i1, isigma) = (0d0, 1d0)*(dens(lm2,lm1,1,i1sp2,2,i1sp1,i1)-dens(lm2,lm1,2,i1sp2,1,i1sp1,i1))*sqa(3) - &
                      (-dens(lm2,lm1,1,i1sp2,1,i1sp1,i1)+dens(lm2,lm1,2,i1sp2,2,i1sp1,i1))*sqa(2)
                  end do           ! LM2
                end do             ! LM1
              end do               ! I1SP2
            end do                 ! I1SP1

          else if (isigma==2) then ! T_y

            do i1sp1 = 1, 2
              do i1sp2 = 1, 2
                do lm1 = 1, lmmax
                  do lm2 = 1, lmmax
                    torq((i1sp2-1)*lmmax+lm2, (i1sp1-1)*lmmax+lm1, i1, isigma) = (-dens(lm2,lm1,1,i1sp2,1,i1sp1,i1)+dens(lm2,lm1,2,i1sp2,2,i1sp1,i1))*sqa(1) - &
                      (dens(lm2,lm1,2,i1sp2,1,i1sp1,i1)+dens(lm2,lm1,1,i1sp2,2,i1sp1,i1))*sqa(3)
                  end do           ! LM2
                end do             ! LM1
              end do               ! I1SP2
            end do                 ! I1SP1

          else if (isigma==3) then ! T_z

            do i1sp1 = 1, 2
              do i1sp2 = 1, 2
                do lm1 = 1, lmmax
                  do lm2 = 1, lmmax
                    torq((i1sp2-1)*lmmax+lm2, (i1sp1-1)*lmmax+lm1, i1, isigma) = (dens(lm2,lm1,2,i1sp2,1,i1sp1,i1)+dens(lm2,lm1,1,i1sp2,2,i1sp1,i1))*sqa(2) - &
                      (0d0, 1d0)*(-dens(lm2,lm1,2,i1sp2,1,i1sp1,i1)+dens(lm2,lm1,1,i1sp2,2,i1sp1,i1))*sqa(1)
                  end do           ! LM2
                end do             ! LM1
              end do               ! I1SP2
            end do                 ! I1SP1

          end if                   ! (ISIGMA=1,2,3)

        end do                     ! I1
      end do                       ! ISIGMA

      ! writeout
      if (mode==0) then
        open (unit=12, file='TBkkr_torq.txt', form='formatted', action='write')
      else                         ! mode==1
        open (unit=12, file='TBkkr_torq_imp.txt', form='formatted', action='write')
      end if
      do isigma = 1, 3
        do i1 = 1, natom
          do lm2 = 1, lmmaxso
            do lm1 = 1, lmmaxso
              write (12, '(2ES25.16)') torq(lm1, lm2, i1, isigma)
            end do
          end do
        end do
      end do
      close (12)

    end if                         ! (myrank==master)

    deallocate (rll)
    deallocate (dens)
    deallocate (rll_12)
    deallocate (torq)

  end subroutine normcoeff_so_torq


end module mod_normcoeff_so_torq
