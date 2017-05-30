  subroutine groundstate_new()
! Ground state properties without fitting GFs
! The matrix elements of the operators are filled in subroutine observables
  use global

  implicit none

  real(kind=r8b),    parameter :: small = 1.d-6
! auxiliary
  complex(kind=c8b), allocatable :: gfsum(:,:), egfsum(:,:), tgfsum(:,:)
  real(kind=r8b)    :: start, finish
! valence charge, spin and orbital moments
  real(kind=r8b)    :: charge, mspin(3), morb(3), tcharge, tmspin(3), tmorb(3)
  real(kind=r8b)    :: spinlen, spindir(3), orblen, orbdir(3)
! torques
  real(kind=r8b)    :: ttorque(3), dth(3), dph(3), cosp, sinp, cost, sint
! looping
  integer(kind=i4b) :: ia, i, ja, j, nr
! summed currents and non-local orbital moments
  real(kind=r8b)    :: curr_sum(0:3,1:3),morb_non_loc_sum(1:3)
  
  write(*,'(/,"Entering groundstate_new",/)')
  write(*,'("Fermi energy=",f16.8,/)') efscf
  allocate(gfsum(nlmsb,nlmsb),egfsum(nlmsb,nlmsb),tgfsum(nlmsb,nlmsb))
! Prepare fit coefficients for GF
  if (lfit) then
    call cpu_time(start)
    do ja=1,nasusc
      do ia=1,nasusc
        call ratfit_gf2(ia,ja,lgsonsite,lgsstruct)
      end do
    end do
    call cpu_time(finish)
    write(*,'(/," GS fit setup   time=",f10.3," s",/)') finish - start
  end if
  call cpu_time(start)
! ----------------------------------------------------------------------
  tcharge = 0.d0; tmspin = 0.d0; tmorb = 0.d0
!  write(*,'(" ia  |  charge     |  mspin  x    mspin  y    mspin  z   |  morb  x     morb  y     morb  z")')
  write(*,'(" ia  |  charge      mspin       morb       | spin x  spin y  spin z  |  orb x   orb y   orb z")')
! output files for currents (charge=0, (x,y,z) = (1,2,3) )
  if(lcurrent) then
!   Interpolation on a regular grid
    if(lcurrentint) then
      open(unit=220,file="current_int_tot_0.dat")
      open(unit=221,file="current_int_tot_1.dat")
      open(unit=222,file="current_int_tot_2.dat")
      open(unit=223,file="current_int_tot_3.dat")
    end if
!   Net currents of each atom
    open(unit=330,file="current_net_0.dat")
    write(330,'("# positions jx jy jz")')
    open(unit=331,file="current_net_1.dat")
    write(331,'("# positions jx jy jz")')
    open(unit=332,file="current_net_2.dat")
    write(332,'("# positions jx jy jz")')
    open(unit=333,file="current_net_3.dat")
    write(333,'("# positions jx jy jz")') 
    open(unit=440,file="orbital_mag.dat")
    write(440,'("# positions morb morb_non_loc")') 
!   summed currents and summed non-local orbital moment
    curr_sum=0.d0
    morb_non_loc_sum(:)=0.d0
  end if 
  do ia=1,nasusc
! Matrix elements of GF for site-diagonal quantities
    call get_gfsum(ia,gfsum,egfsum,tgfsum,lgsonsite,lgsstruct)
! Charge, magnetization, orbital moments, spin axis
    call integrated_charge(ia,gfsum,egfsum,tgfsum,charge,mspin,morb)
!    write(*,'(i4," |",f12.8," |",3f12.8," |",3f12.8)') ia, charge, mspin, morb
    spinlen = sqrt(dot_product(mspin,mspin))
    if (spinlen < small) then
!     Direction from kkrflex
      spindir =  magdir(:,ia)
    else
      spindir = mspin/spinlen
    end if
!    spindir = mspin/(spinlen + 1.d-12)
    orblen = sqrt(dot_product(morb,morb))
    if (orblen < small) then 
!     Direction from kkrflex
      orbdir =  magdir(:,ia)
    else
      orbdir = morb/orblen
    end if
!    orbdir = morb/(orblen + 1.d-12)
    write(*,'(i4," |",3f12.8," |",3f8.4," |",3f8.4)') ia, charge, spinlen, orblen, spindir, orbdir
    tcharge = tcharge + charge
    tmspin = tmspin + mspin
    tmorb = tmorb + morb
! Charge, magnetization, orbital densities, new rho2ns
    call charge_density(ia,gfsum)
! Current density
    if(lcurrent) then
      nr=nrpts(ia)
      call current_density(ia,nr,gfsum,curr_sum,morb_non_loc_sum,npanat(ia),ircutat(:,ia))
    end if
! Compute DOS
    if (ldos) call density_of_states(ia,lgsonsite,lgsstruct)
!    write(*,*) "before density matrix"
! ----------------------------------------------------------------
! xc-potential and energy
!    call overlaps_susc(ia,gfsum)
!    call get_xc(ia,lmmax4,gfsum,kxclm)
!    call get_rhoden(ia,lmmax2,gfsum,rhoden)
!    call get_xc2(ia,lmmax4,kxclm)
!    call get_xc_basis(ia,lmmax2,gfsum,kxc,kxcbasis(:,:,ia))
! ----------------------------------------------------------------
! Interpolation of charge and magnetization density on regular grid
    if (ia==1) then
      if ((lsusc .AND. lcurrcorr .AND. lcurrcorrint) .OR. (lcurrent .AND. lcurrentint)) call rho_interpolation()
    end if
! ----------------------------------------------------------------------
  end do
  if(lcurrent) then
    close(220); close(221); close(222); close(223)
    close(330); close(331); close(332); close(333)
    close(440)
  end if
  write(*,'(" qe    tot=",f16.8)') tcharge
  write(*,'(" mspin tot=",3f16.8)') tmspin
  write(*,'(" morb  tot=",3f16.8)') tmorb
  if(lcurrent) then
    write(*,'(" morb_non_loc   tot=",3f18.10)') morb_non_loc_sum
    write(*,'(" charge current tot=",3f18.10)') curr_sum(0,:) 
    write(*,'(" spin_x current tot=",3f18.10)') curr_sum(1,:) 
    write(*,'(" spin_y current tot=",3f18.10)') curr_sum(2,:) 
    write(*,'(" spin_z current tot=",3f18.10)') curr_sum(3,:) 
  end if
  write(*,'(" eband tot=",f16.8)') sum(ebandv)
! ----------------------------------------------------------------------
! Output band energies
  do i=1,3
    ttorque(i) = sum(etorque(i,:,:))
  end do
  write(*,'(" etorque tot xyz=",3f16.8)') ttorque
  if (lrot .and. ispinrot == 0) then
    cost = urot(3)
    if (sqrt(1.d0 - cost**2) < 1.d-6) then
      dth = (/1.d0,0.d0,0.d0/)
      dph = (/0.d0,1.d0,0.d0/)
    else
      sint = sqrt((1.d0-cost)*(1.d0+cost))
      cosp = urot(1)/sint
      sinp = urot(2)/sint
      dth = (/cosp*cost,sinp*cost,-sint/) 
      dph = (/-sinp*sint,cosp*sint,0.d0/) 
    end if
    write(*,'(" etorque tot the=",f16.8)') dot_product(ttorque,dth)
    write(*,'(" etorque tot phi=",f16.8)') dot_product(ttorque,dph)
  end if
  if (lldau) write(*,'(" eldau tot=",f16.8)') sum(eldau)
! to file
  open(file='eband.dat',unit=iofile,status='replace')
  write(iofile,'(" ia  |  eband tot  |  eband dn    eband up")')
  do ia=1,nasusc
    write(iofile,'(i4," |",f12.8," |",2f12.8)') ia, sum(ebandv(:,:,ia)), sum(ebandv(:,:,ia),dim=1)
  end do
  write(iofile,'(" ia  |  eband dn l |  eband up l")')
  do ia=1,nasusc
    write(iofile,'(i4," |",4f12.8," |",4f12.8)') ia, ebandv(:,1,ia), ebandv(:,nsmax,ia)
  end do
  write(iofile,'(" ia  | torq xyz  | torq xyz / mspin")')
  do ia=1,nasusc
!    write(iofile,'(i4," |",3f12.8)') ia, sum(etorque(:,:,ia),dim=2) - sum(sum(etorque(:,:,ia),dim=2)*magdir(:,ia))*magdir(:,ia)
    write(iofile,'(i4," |",3f12.8," |",3f12.8)') ia, sum(etorque(:,:,ia),dim=2), sum(etorque(:,:,ia),dim=2)/gs_mlm(1,ia)
  end do
!  write(iofile,'(" ia  | torq x tot  | torq x dn   torq x up")')
!  do ia=1,nasusc
!    write(iofile,'(i4," |",f12.8," |",2f12.8)') ia, sum(etorque(1,:,:,ia)), sum(etorque(1,:,:,ia),dim=1)
!  end do
!  write(iofile,'(" ia  | torq x dn l | torq x up l")')
!  do ia=1,nasusc
!    write(iofile,'(i4," |",4f12.8," |",4f12.8)') ia, etorque(1,:,1,ia), etorque(1,:,nsmax,ia)
!  end do
!  write(iofile,'(" ia  | torq y tot  | torq y dn   torq y up")')
!  do ia=1,nasusc
!    write(iofile,'(i4," |",f12.8," |",2f12.8)') ia, sum(etorque(2,:,:,ia)), sum(etorque(2,:,:,ia),dim=1)
!  end do
!  write(iofile,'(" ia  | torq y dn l | torq y up l")')
!  do ia=1,nasusc
!    write(iofile,'(i4," |",4f12.8," |",4f12.8)') ia, etorque(2,:,1,ia), etorque(2,:,nsmax,ia)
!  end do
!  write(iofile,'(" ia  | torq z tot  | torq z dn   torq z up")')
!  do ia=1,nasusc
!    write(iofile,'(i4," |",f12.8," |",2f12.8)') ia, sum(etorque(3,:,:,ia)), sum(etorque(3,:,:,ia),dim=1)
!  end do
!  write(iofile,'(" ia  | torq z dn l | torq z up l")')
!  do ia=1,nasusc
!    write(iofile,'(i4," |",4f12.8," |",4f12.8)') ia, etorque(3,:,1,ia), etorque(3,:,nsmax,ia)
!  end do
  close(iofile)
  ebandv = ebandv + eldau
! ----------------------------------------------------------------------
! output density matrix for selected atoms and l-channels
  if (ldos .and. ldosdmat) call density_matrix(lgsonsite,lgsstruct)
! ----------------------------------------------------------------------
  call cpu_time(finish)
  deallocate(gfsum,egfsum,tgfsum)
  write(*,'(/," GS quantities  time=",f10.3," s",/)') finish - start
! ----------------------------------------------------------------------
  write(*,'(/,"Leaving groundstate_new",/)')
! All done
  end subroutine groundstate_new
