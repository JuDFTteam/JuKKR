  subroutine init_arrays(my_rank)
! Allocate the storage according to the susc input
  use global

  implicit none
  
  integer(kind=i4b) :: is, il, im, ib, i, lm, ia, ilms, jlms, my_rank
  real(kind=r8b)    :: ram

  if (noparameters) stop 'init_arrays: run init_param first!'
! -----------------------------------------------------------------------
!                           Energy mesh          
! -----------------------------------------------------------------------
  allocate(ensusc(nesusc))                            ! energies for susc calculation
  allocate(esusc(nesusc))                             ! energy mesh for susc
  allocate(eksusc(nesusc))                            ! sqrt of e
  allocate(desusc(nesusc))                            ! energy integration weights
  allocate(wsusc(nesusc))                             ! interpolation weights
  esusc = 0.d0; eksusc = 0.d0; desusc = 0.d0
  allocate(escf(nescf))                               ! same things
  allocate(ekscf(nescf))                              ! but to store the energy contour
  allocate(descf(nescf))                              ! used in the SCF calculation
  escf = 0.d0; ekscf = 0.d0; descf = 0.d0
! -----------------------------------------------------------------------
!                      Groundstate quantities 
! -----------------------------------------------------------------------
  allocate(nrpts(nasusc))                             ! number of points in each radial mesh
  allocate(nrpts0(nasusc))                            ! start of each radial mesh
  allocate(nrpts1(nasusc))                            ! end of radial mesh
  nrpts = 0; nrpts0 = 0; nrpts1 = 0
  allocate(rmesh(nrmax,nasusc))                       ! radial mesh for evaluation
  allocate(rsmesh(nrmax,0:nlmax,nasusc))              ! powers of r
  allocate(drmesh(nrmax,nasusc))                      ! weights for radial integration
  allocate(drproj(nrmax,nasusc))                      ! weights for radial projection
  rmesh = 0.d0; drmesh = 0.d0; drproj = 0.d0
  allocate(normesh(nasusc))                           ! was the radial mesh stored?
  normesh      = .true.
  allocate(zat(nasusc))                               ! atomic numbers
  allocate(ri(3,nasusc))                              ! atomic positions
  allocate(magdir(3,nasusc))                          ! magnetization direction
  magdir(1,1:nasusc) = 0.d0
  magdir(2,1:nasusc) = 0.d0
  magdir(3,1:nasusc) = 1.d0
  allocate(newdir(3,nasusc))                          ! magnetization direction
  newdir(:,1:nasusc) = 0.d0
  allocate(iarot(nasusc))                             ! spin direction updates
  iarot(1:nasusc) = 0
  allocate(spinproj(nsmax,nsmax,nasusc))              ! spin projections
  allocate(vr(nrmax,nasusc))                          ! spherical charge potentials
  allocate(br(nrmax,nasusc))                          ! spherical magnetic potentials
  allocate(nrc(nrmax,nasusc))                         ! spherical core charge densities
  allocate(mrc(nrmax,nasusc))                         ! spherical core magnetization densities
  allocate(nrv(nrmax,0:6,nasusc))                     ! spherical valence densities
  zat = 0.d0; vr = 0.d0; br = 0.d0; nrc = 0.d0; mrc = 0.d0
  allocate(old_rho2ns(nrmax,lmmax2,nsmax,nasusc))          ! storage for old groundstate densities
  allocate(gs_qlm(lmmax2,nasusc),gs_mlm(lmmax2,nasusc))    ! storage for groundstate multipoles
  allocate(new_rho2ns(nrmax,lmmax2,nsmax,nasusc),rho_lm(nrmax,lmmax2,0:3,nasusc))          ! storage for new groundstate densities
  allocate(rhomat(nlms,nlms,nasusc))                       ! density matrix integrated over WS sphere
  allocate(kxclm(nrmax,lmmax2,4,4,nasusc))                 ! xc kernel
  allocate(nlmpot(nasusc),i2lmpot(lmmax2,nasusc))          ! lmpot indices (handled by read_rho2ns)
!  allocate(lorb(lmmax,lmmax,3))                       ! storage for angular momentum matrices
  allocate(lorb(lmmax4,lmmax4,3))                       ! storage for angular momentum matrices
  allocate(ebandv(0:nlmax,nsmax,nasusc))              ! storage for recalculated band energy
  allocate(etorque(3,0:nlmax,nasusc))                 ! storage for torque anisotropy energy
  etorque = 0.d0
  allocate(eldau(0:nlmax,nsmax,nasusc))               ! storage for LDA+U energy
  allocate(vshift(0:nlmax,nsmax,nasusc))              ! storage for chemical potential shift
  vshift = 0.d0
! -----------------------------------------------------------------------
!                             Pointers               
! -----------------------------------------------------------------------
! Indices for packing and unpacking arrays
! --> lm (big to accommodate Gaunt numbers)
  allocate(lm2i(-lmmax4:lmmax4,0:lmmax4))
  allocate(i2lm(2,lmmax4))
  i = 0
  do il=0,nlmax4
    do im=-il,il
      i = i + 1
      lm2i(im,il) = i
      i2lm(:,i) = (/im,il/)
    end do
  end do
! --> lms
  allocate(lms2i(lmmax,nsmax))
  allocate(i2lms(2,nlms))
  i = 0
  do lm=1,lmmax
    do is=1,nsmax
      i  = i + 1
      lms2i(lm,is) = i
      i2lms(:,i) = (/lm,is/)
    end do
  end do
  if (my_rank == 0) then
    write(*,'(" nlms   should be ",i8)') i
  end if ! my_rank
! --> alms
  allocate(alms2i(nlms,nasusc))
  allocate(i2alms(2,nalms))
  i = 0
  do ia=1,nasusc
    do ilms=1,nlms
      i = i + 1
      alms2i(ilms,ia) = i
      i2alms(:,i) = (/ilms,ia/)
    end do
  end do
  if (my_rank == 0) then 
    write(*,'(" nalms  should be ",i8)') i
  end if ! my_rank 
! New pointers for KKRSUSC and NEW SOC SOLVER
! --> lms_new
  allocate(lms2i_new(lmmax,nsmax))
  allocate(i2lms_new(2,nlms))
  i = 0
  do is=1,nsmax
    do lm=1,lmmax
      i  = i + 1
      lms2i_new(lm,is) = i
      i2lms_new(:,i) = (/lm,is/)
    end do
  end do
!  write(*,'(" nlms   should be ",i8)') i
! All done!
  end subroutine init_arrays
