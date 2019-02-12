!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculates the KKR matrix elements for the spin flux operator
!> Author: Guillaume Géranton
!> Calculates the KKR matrix elements for the spin flux operator, i.e.,
!> \begin{equation}
!> \int dr\left[R^\mu_{Ls} \right]^\dagger Q_s R^\mu_{L's'}
!> \end{equation}
!------------------------------------------------------------------------------------
!> @note Details are in http://arxiv.org/pdf/1602.03417v1.pdf
!> This subroutine was adapted from `NORMCOEFF_SO`.
!> @endnote
!------------------------------------------------------------------------------------
module mod_normcoeff_so_spinflux

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculates the KKR matrix elements for the spin flux operator
  !> Author: Guillaume Géranton
  !> Category: physical-observables, KKRhost
  !> Deprecated: False 
  !> Calculates the KKR matrix elements for the spin flux operator, i.e.,
  !> \begin{equation}
  !> \int dr\left[R^\mu_{Ls} \right]^\dagger Q_s R^\mu_{L's'}
  !> \end{equation}
  !-------------------------------------------------------------------------------
  !> @note Details are in http://arxiv.org/pdf/1602.03417v1.pdf
  !> This subroutine was adapted from `NORMCOEFF_SO`.
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine normcoeff_so_spinflux(natom, ircut, lmmax0d, pns, ksra, drdi, mode)

#ifdef CPP_MPI
    use :: mpi
    use :: mod_types, only: t_mpi_c_grid, t_inc, t_imp
#else
    use :: mod_types, only: t_inc, t_imp
#endif
    use :: mod_mympi, only: myrank, master
    use :: mod_datatypes, only: dp
    use :: global_variables, only: lmmaxd, ipand, natypd, irmd, nspind
    use :: mod_constants, only: czero

    implicit none

    ! .. Input variables
    integer, intent(in) :: mode
    integer, intent(in) :: ksra
    integer, intent(in) :: natom
    integer, intent (in) :: lmmax0d !! (LMAX+1)^2
    integer, dimension(0:ipand, natypd), intent(in) :: ircut !! R points of panel borders
    real (kind=dp), dimension(irmd, natypd), intent(in) :: drdi !! Derivative dr/di
    complex (kind=dp), dimension(nspind*lmmax0d, nspind*lmmax0d, irmd, 2, natom), intent(in) :: pns
    ! .. Local Scalars ..
    integer :: lm1, lm2, lm1p, ir, i1, i1sp1, i1sp2, lmsp1, lmsp2, isigma, i2sp1, i2sp2, insra, nsra
    integer :: i2
    complex (kind=dp) :: delta1, delta2
    ! MPI stuff
    integer :: ierr, ihelp, i1_start, i1_end
    ! ..Local Arrays..
    complex (kind=dp), dimension(:), allocatable :: rll_12
    complex (kind=dp), dimension(:, :, :, :), allocatable :: spinflux
    complex (kind=dp), dimension(:, :, :, :, :, :, :), allocatable :: dens
    complex (kind=dp), dimension(:, :, :, :, :, :, :), allocatable :: rll
    ! MPI stuff
    complex (kind=dp), dimension(:, :, :, :, :, :, :), allocatable :: work
    ! ..

    if (t_inc%i_write>0) write (1337, *) 'KSRA', ksra
    if (ksra>=1) then              ! previously this was .GT. which is wrong for kvrel=1
      nsra = 2
    else
      nsra = 1
    end if

    if (t_inc%i_write>0) then
      write (1337, *) 'NSRA', nsra
      write (1337, *) 'lmmax0d', lmmax0d
      write (1337, *) 'lmmaxd', lmmaxd
    end if

    allocate (rll(irmd,lmmax0d,lmmax0d,2,2,2,natom))
    allocate (rll_12(lmmax0d))
    allocate (dens(lmmax0d,lmmax0d,2,2,2,2,natom))
    allocate (spinflux(lmmaxd,lmmaxd,natom,3))

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

    ! rewrite the wavefunctions in RLL arrays of 1,2*lmmax0d
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
              do lm1 = 1, lmmax0d
                lmsp1 = (i1sp1-1)*lmmax0d + lm1
                do lm2 = 1, lmmax0d
                  lmsp2 = (i1sp2-1)*lmmax0d + lm2
                  rll(ir, lm2, lm1, i1sp2, i1sp1, insra, i1) = pns(lmsp2, lmsp1, ir, insra, i1)
                end do             ! LM1=1,lmmax0d
              end do               ! LM1=1,lmmax0d
            end do                 ! ISP1=1,2
          end do                   ! ISP1=1,2

        end do                     ! IR
      end do                       ! INSRA


      ! set up the array R*_L1L2 R_L3L4
      do i1sp1 = 1, 2
        do i1sp2 = 1, 2
          do i2sp1 = 1, 2
            do i2sp2 = 1, 2

              do lm1 = 1, lmmax0d
                do lm2 = 1, lmmax0d

                  do insra = 1, nsra

                    rll_12 = czero

                    do lm1p = 1, lmmax0d

                      delta1 = (rll(ircut(1,i2),lm1p,lm2,i2sp1,i2sp2,insra,i1)-rll(ircut(1,i2)-1,lm1p,lm2,i2sp1,i2sp2,insra,i1))/drdi(ircut(1,i2), i2)
                      delta2 = (rll(ircut(1,i2),lm1p,lm1,i1sp1,i1sp2,insra,i1)-rll(ircut(1,i2)-1,lm1p,lm1,i1sp1,i1sp2,insra,i1))/drdi(ircut(1,i2), i2)

                      rll_12(lm1p) = conjg(rll(ircut(1,i2)-1,lm1p,lm1,i1sp1,i1sp2,insra,i1))*delta1 - rll(ircut(1,i2)-1, lm1p, lm2, i2sp1, i2sp2, insra, i1)*conjg(delta2)

                    end do         ! LM1P

                    do lm1p = 1, lmmax0d
                      dens(lm1, lm2, i1sp1, i1sp2, i2sp1, i2sp2, i1) = dens(lm1, lm2, i1sp1, i1sp2, i2sp1, i2sp2, i1) + rll_12(lm1p)
                    end do         ! LM1P

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
    allocate (work(lmmax0d,lmmax0d,2,2,2,2,natom), stat=ierr)
    if (ierr/=0) stop 'Error allocating work for MPI comm of DENS in normcoeff_spinf'
    ihelp = lmmax0d*lmmax0d*2*2*2*2*natom
    call mpi_reduce(dens, work, ihelp, mpi_double_complex, mpi_sum, master, t_mpi_c_grid%mympi_comm_ie, ierr)
    if (ierr/=mpi_success) stop 'Error in MPI comm of DENS in normcoeff_spinflux'
    dens(:, :, :, :, :, :, :) = work(:, :, :, :, :, :, :)
    deallocate (work, stat=ierr)
    if (ierr/=0) stop 'Error deallocating work for MPI comm of DENS in normcoeff_spin'
#endif


    if (myrank==master) then       ! do last part and writeout only on master

      write (*, *) 'collect terms and writeout'

      spinflux = czero

      do isigma = 1, 3             ! ISIGMA == 1 --> Q_x
        ! ISIGMA == 2 --> Q_y
        ! ISIGMA == 3 --> Q_z

        write (6, *) 'ISIGMA', isigma
        do i1 = 1, natom

          if (isigma==1) then      ! Q_x

            do i1sp1 = 1, 2
              do i1sp2 = 1, 2
                do lm1 = 1, lmmax0d
                  do lm2 = 1, lmmax0d
                    spinflux((i1sp2-1)*lmmax0d+lm2, (i1sp1-1)*lmmax0d+lm1, i1, isigma)=&
                    -(0d0, 1d0)*(dens(lm2,lm1,2,i1sp2,1,i1sp1,i1)+dens(lm2,lm1,1,i1sp2,2,i1sp1,i1))/2
                  end do           ! LM2
                end do             ! LM1
              end do               ! I1SP2
            end do                 ! I1SP1

          else if (isigma==2) then ! Q_y

            do i1sp1 = 1, 2
              do i1sp2 = 1, 2
                do lm1 = 1, lmmax0d
                  do lm2 = 1, lmmax0d
                    spinflux((i1sp2-1)*lmmax0d+lm2, (i1sp1-1)*lmmax0d+lm1, i1, isigma)= &
                    -(0d0, 1d0)*(-1)*(0d0, 1d0)*(dens(lm2,lm1,2,i1sp2,1,i1sp1,i1)-dens(lm2,lm1,1,i1sp2,2,i1sp1,i1)) /2
                  end do           ! LM2
                end do             ! LM1
              end do               ! I1SP2
            end do                 ! I1SP1

          else if (isigma==3) then ! Q_z

            do i1sp1 = 1, 2
              do i1sp2 = 1, 2
                do lm1 = 1, lmmax0d
                  do lm2 = 1, lmmax0d
                    spinflux((i1sp2-1)*lmmax0d+lm2, (i1sp1-1)*lmmax0d+lm1, i1, isigma) = (0d0, 1d0)*(dens(lm2,lm1,1,i1sp2,1,i1sp1,i1)-dens(lm2,lm1,2,i1sp2,2,i1sp1,i1))/2
                  end do           ! LM2
                end do             ! LM1
              end do               ! I1SP2
            end do                 ! I1SP1

          end if

        end do                     ! I1
      end do                       ! ISIGMA

      ! write to file
      if (mode==0) then
        open (unit=12, file='TBkkr_spinflux.txt', form='formatted', action='write')
      else                         ! mode==1
        open (unit=12, file='TBkkr_spinflux_imp.txt', form='formatted', action='write')
      end if
      do isigma = 1, 3
        do i1 = 1, natom
          do lm2 = 1, lmmaxd
            do lm1 = 1, lmmaxd
              ! minus sign to get the spin flux into the sphere :
              write (12, '(2ES25.16)') - spinflux(lm1, lm2, i1, isigma)
            end do
          end do
        end do
      end do
      close (12)

    end if                         ! (myrank==master)

    deallocate (rll)
    deallocate (dens)
    deallocate (rll_12)
    deallocate (spinflux)

  end subroutine normcoeff_so_spinflux

end module mod_normcoeff_so_spinflux
