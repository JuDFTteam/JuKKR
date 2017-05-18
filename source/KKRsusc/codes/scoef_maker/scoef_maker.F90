  program scoef_maker
! generates clusters around chosen sites
! input: BRAVAIS and RBASIS as in the inputcard
! input: site locations and truncation radii

  implicit none

  integer*4 :: nbravais, nbasis, nimp, nmax
  real*8    :: alat
  integer*4, allocatable :: ibasis(:), ijkimp(:,:), nclus(:), iclus(:,:)
  real*8,    allocatable :: rbravais(:,:), rbasis(:,:), rcut(:), zat(:), rclus(:,:,:)
  integer*4 :: i, j
  real*8    :: rxyz(3), r0(3), rij(3), a3d(3,3), b3d(3,3)
  logical   :: lcartesian


! read input
  open(unit=10,file='test_scoef',status='old')
! ----------------------------------------------------------------------
  read(10,*) ! bravais vectors
  read(10,*) nbravais, alat
  write(*,'(" unit of length=",f16.8)') alat
  if (nbravais == 2) then
    write(*,'(" 2D lattice generation")')
  else if (nbravais == 3) then
    write(*,'(" 3D lattice generation")')
  else
    stop 'check nbravais'
  end if
  allocate(rbravais(nbravais,nbravais))
  do i=1,nbravais
    read(10,*) rbravais(1:nbravais,i)
    write(*,'("  ibravais=",i4," rbravais=",3f16.8)') i, rbravais(1:nbravais,i)
  end do
! ----------------------------------------------------------------------
  read(10,*) ! basis vectors
  read(10,*) nbasis, lcartesian
  write(*,'(/," number of basis vectors=",i4,"  cartesian=",l1)') nbasis, lcartesian
  allocate(rbasis(3,nbasis),zat(nbasis))
  write(*,'(" basis vectors divided by unit of length:")')
  do i=1,nbasis
    read(10,*) rbasis(1:3,i), zat(i)
    write(*,'("  ibasis=",i4,"  rbasis=",3f16.8,"  zat=",f6.1)') i, rbasis(1:3,i), zat(i)
  end do
! generate 3D form of bravais vectors
  call bravais3d(nbravais,rbravais,alat,a3d,b3d)
! ----------------------------------------------------------------------
  read(10,*) ! target sites, maximum number of atoms in cluster
  read(10,*) nimp, nmax
  write(*,'(/," number of target sites=",i4,"  maximum number of cluster atoms=",i8)') nimp, nmax
  allocate(ijkimp(nbravais,nimp),rcut(nimp),ibasis(nimp),nclus(nimp),iclus(nmax,nimp),rclus(3,nmax,nimp))
  do i=1,nimp
    read(10,*) ijkimp(1:nbravais,i), rcut(i), ibasis(i)
    write(*,'("  imp=",i4,"  rcut=",f16.8,"  ibasis=",i4," ijkimp=",3i4)') i, rcut(i), ibasis(i), ijkimp(1:nbravais,i)
  end do
  close(10)
! ----------------------------------------------------------------------
! generate truncated cluster around each impurity
  do i=1,nimp
    rxyz = 0.d0
!   basis vector for impurity sublattice
    if (.not.lcartesian) then
      do j=1,3
        rxyz(:) = rxyz(:) + a3d(:,j)*rbasis(j,ibasis(i))
      end do
    else
      rxyz(:) = rbasis(:,ibasis(i))
    end if
!   the sublattice of the first impurity site is taken as the global origin
    if (i == 1) r0(:) = rxyz(:)
!   vector connecting current sublattice to first sublattice
    rij(:) = rxyz(:) - r0(:)
!    fold_vector is causing a bug by shifting the basis of all the impurities to the lattice site closest to the first impurity
!    call fold_vector(nbravais,a3d,rij)
    write(*,'("  imp=",i4," r0=",3f8.4,"  rij=",3f8.4)') i, r0(:), rij(:)
!   now add the lattice vector
    do j=1,nbravais
      rxyz(:) = rxyz(:) + a3d(:,j)*ijkimp(j,i)
    end do
!    write(*,'(/,"  imp=",i8,"  rijkimp, length=",4f16.8)') i, rxyz(1:nbravais), sqrt(dot_product(rxyz(1:nbravais),rxyz(1:nbravais)))
    call truncated_cluster(nbravais,a3d,b3d,lcartesian,nbasis,rbasis,ijkimp(:,i),ibasis(i),rcut(i),nmax,nclus(i),iclus(:,i),rclus(:,:,i))
    write(*,'("  imp=",i4,"  nclus=",i4)') i, nclus(i)
!   add to all cluster positions the vector connecting to the origin
    do j=1,nclus(i)
      rclus(:,j,i) = rclus(:,j,i) + rij(:)
    end do
  end do
! now remove repeated sites and write scoef
  call merge_clusters(nbasis,nimp,nmax,nclus,iclus,rclus,zat)
! All done!

  contains


  subroutine merge_clusters(nbasis,nimp,nmax,nclus,iclus,rclus,zat)
! merges the clusters surrounding each impurity and writes scoef
! origin is the position of the first impurity

  implicit none

! number of sublattices
  integer*4, intent(in) :: nbasis
! number of impurities
  integer*4, intent(in) :: nimp
! maximum number of sites
  integer*4, intent(in) :: nmax
! number of atoms in each cluster surrounding an impurity
  integer*4, intent(in) :: nclus(nimp)
! sublattice for each site in each cluster
  integer*4, intent(in) :: iclus(nmax,nimp)
! site coordinates
  real*8,    intent(in) :: rclus(3,nmax,nimp)
! atomic numbers for atoms in each sublattice
  real*8,    intent(in) :: zat(nbasis)
! ----------------------------------------------------------------------
  real*8, parameter :: tol = 1.d-6
  real*8,    allocatable :: rsites(:,:), zsites(:), dist(:)
  integer*4, allocatable :: isites(:), iimp(:), niimp(:)
  integer*4 :: i, j, k, itot, ntot, ndiff, inew
  real*8    :: diff, rvec(3)

! total number of positions to check
  ntot = sum(nclus(1:nimp))
  allocate(rsites(3,ntot),isites(ntot),zsites(ntot),iimp(ntot),niimp(nimp),dist(ntot))
! ----------------------------------------------------------------------
! first positions are for the impurities
  itot = 0
  do i=1,nimp
    itot = itot + 1
    rsites(:,itot) = rclus(:,1,i)
    isites(itot)   = iclus(1,i)
    zsites(itot)   = zat(iclus(1,i))
    iimp(itot)     = i
    dist(itot)     = 0.d0
    diff = 1.d99
    do j=1,i-1
      diff = min(diff,maxval(abs(rsites(:,j) - rsites(:,i))))
    end do
    if (diff < tol) stop 'merge_clusters: repeated impurity positions'
  end do
  write(*,'(" Distances between each impurity")')
  do j=1,nimp
    do i=1,nimp
      write(*,'(2i4,f16.8)') i, j, sqrt(dot_product(rsites(:,j)-rsites(:,i),rsites(:,j)-rsites(:,i)))
    end do
  end do
! ----------------------------------------------------------------------
! for all other sites
! the positions on each cluster are compared with the current ones in rsites
  itot = nimp
  do i=1,nimp
    ndiff = itot
    do j=2,nclus(i)
      diff = 1.d99
      do k=1,ndiff
        diff = min(diff,maxval(abs(rclus(:,j,i) - rsites(:,k))))
      end do
!     this is a repeated site
      if (diff < tol) cycle
!     this is a distinct site
      itot = itot + 1
      rsites(:,itot) = rclus(:,j,i)
      isites(itot)   = iclus(j,i)
      zsites(itot)   = zat(iclus(j,i))
      iimp(itot)     = i
    end do
  end do
! final tally
  ndiff = itot
! ----------------------------------------------------------------------
! now assign revised impurity indices according to distance
  niimp = 0
  do i=1,ndiff
    diff = 1.d99
    inew = 1
    do j=1,nimp
      dist(i) = sqrt(dot_product(rsites(:,i)-rsites(:,j),rsites(:,i)-rsites(:,j)))
!      if (dist(i) < diff) then
      ! Benedikt: 1) " - 1.d-09"  ---> avoids numerical instabilities when assigning
      !                                to the first impurity and one of the other
      !                                impurities has the same distance
      !           2) " < rcut(j)" ---> no assignment is allowed when the radius around
      !                                the regarded impurity is already exceeded
      if ( (dist(i) < diff - 1.d-09) .and. dist(i) < rcut(j) ) then
        inew = j
        diff = dist(i)
      end if
    end do
    iimp(i) = inew
    niimp(iimp(i)) = niimp(iimp(i)) + 1
    dist(i) = sqrt(dot_product(rsites(:,i)-rsites(:,iimp(i)),rsites(:,i)-rsites(:,iimp(i))))
  end do
  write(*,'(" niimp=",100i6)') niimp(1:nimp)
! ----------------------------------------------------------------------
! sort site positions
  call ugly_sort(ntot,nimp,ndiff,iimp,rsites,isites,zsites,dist)
! ----------------------------------------------------------------------
  write(*,'(/," merge_clusters: ndiff=",i8)') ndiff
! now write scoef
  open(file='scoef_new',unit=10,status='replace')
  write(10,'(i8)') ndiff
  do i=1,ndiff
    write(10,'(3f20.12,i4,f8.1,f20.12,i4)') rsites(:,i), isites(i), zsites(i), dist(i), iimp(i)
  end do
  close(10)
! All done!
  end subroutine merge_clusters



  subroutine ugly_sort(ntot,nimp,ndiff,iimp,rsites,isites,zsites,dist)
! Reorders the list of impurity sites according to some criteria (read comments)

  implicit none

! Dimension of arrays
  integer*4, intent(in)    :: ntot
! Number of target impurities
  integer*4, intent(in)    :: nimp
! Number of distinct sites around all impurities
  integer*4, intent(in)    :: ndiff
! To which impurity is each site assigned
  integer*4, intent(inout) :: iimp(ntot)
! Position of each site
  real*8,    intent(inout) :: rsites(3,ntot)
  integer*4, intent(inout) :: isites(ntot)
! Atomic number of each site
  real*8,    intent(inout) :: zsites(ntot)
! Distance to nearest impurity
  real*8,    intent(inout) :: dist(ntot)
! ----------------------------------------------------------------------
  real*8, parameter :: tol = 1.d-6
  integer*4 :: i, j, inew
  real*8    :: diff, rvec(3)

! reorder according to impurity
  do i=nimp+1,ndiff
    do j=i+1,ndiff
      if (iimp(j) < iimp(i)) then
!       swap positions
        rvec(:) = rsites(:,i)
        rsites(:,i) = rsites(:,j)
        rsites(:,j) = rvec(:)
!       swap layer numbers
        inew = isites(i)
        isites(i) = isites(j)
        isites(j) = inew
!       swap atomic numbers
        diff = zsites(i)
        zsites(i) = zsites(j)
        zsites(j) = diff
!       swap distances
        diff = dist(i)
        dist(i) = dist(j)
        dist(j) = diff
!       swap impurity indices
        inew = iimp(i)
        iimp(i) = iimp(j)
        iimp(j) = inew
      end if
    end do
  end do
! reorder according to distance
  do i=nimp+1,ndiff
    do j=i+1,ndiff
      if (iimp(j) == iimp(i) .and. dist(i) - dist(j) > tol) then
!       swap positions
        rvec(:) = rsites(:,i)
        rsites(:,i) = rsites(:,j)
        rsites(:,j) = rvec(:)
!       swap layer numbers
        inew = isites(i)
        isites(i) = isites(j)
        isites(j) = inew
!       swap atomic numbers
        diff = zsites(i)
        zsites(i) = zsites(j)
        zsites(j) = diff
!       swap distances
        diff = dist(i)
        dist(i) = dist(j)
        dist(j) = diff
!       swap impurity indices
        inew = iimp(i)
        iimp(i) = iimp(j)
        iimp(j) = inew
      end if
    end do
  end do
! reorder according to z-coordinate
  do i=nimp+1,ndiff
    do j=i+1,ndiff
      if (iimp(j) == iimp(i) .and. abs(dist(i)-dist(j)) < tol .and. rsites(3,j) < rsites(3,i)) then
!       swap positions
        rvec(:) = rsites(:,i)
        rsites(:,i) = rsites(:,j)
        rsites(:,j) = rvec(:)
!       swap layer numbers
        inew = isites(i)
        isites(i) = isites(j)
        isites(j) = inew
!       swap atomic numbers
        diff = zsites(i)
        zsites(i) = zsites(j)
        zsites(j) = diff
!       swap distances
        diff = dist(i)
        dist(i) = dist(j)
        dist(j) = diff
!       swap impurity indices
        inew = iimp(i)
        iimp(i) = iimp(j)
        iimp(j) = inew
      end if
    end do
  end do
! reorder according to y-coordinate
  do i=nimp+1,ndiff
    do j=i+1,ndiff
      if (iimp(j) == iimp(i) .and. abs(dist(i)-dist(j)) < tol .and. abs(rsites(3,j)-rsites(3,i)) < tol .and. rsites(2,j) < rsites(2,i)) then
!       swap positions
        rvec(:) = rsites(:,i)
        rsites(:,i) = rsites(:,j)
        rsites(:,j) = rvec(:)
!       swap layer numbers
        inew = isites(i)
        isites(i) = isites(j)
        isites(j) = inew
!       swap atomic numbers
        diff = zsites(i)
        zsites(i) = zsites(j)
        zsites(j) = diff
!       swap distances
        diff = dist(i)
        dist(i) = dist(j)
        dist(j) = diff
!       swap impurity indices
        inew = iimp(i)
        iimp(i) = iimp(j)
        iimp(j) = inew
      end if
    end do
  end do
! reorder according to x-coordinate
  do i=nimp+1,ndiff
    do j=i+1,ndiff
      if (iimp(j) == iimp(i) .and. abs(dist(i)-dist(j)) < tol .and. abs(rsites(3,j)-rsites(3,i)) < tol .and. abs(rsites(2,j)-rsites(2,i)) < tol .and. rsites(1,j) < rsites(1,i)) then
!       swap positions
        rvec(:) = rsites(:,i)
        rsites(:,i) = rsites(:,j)
        rsites(:,j) = rvec(:)
!       swap layer numbers
        inew = isites(i)
        isites(i) = isites(j)
        isites(j) = inew
!       swap atomic numbers
        diff = zsites(i)
        zsites(i) = zsites(j)
        zsites(j) = diff
!       swap distances
        diff = dist(i)
        dist(i) = dist(j)
        dist(j) = diff
!       swap impurity indices
        inew = iimp(i)
        iimp(i) = iimp(j)
        iimp(j) = inew
      end if
    end do
  end do
! All done!
  end subroutine ugly_sort


  subroutine truncated_cluster(nbravais,a3d,b3d,lcartesian,nbasis,rbasis,ijk,ibasis,rcut,nmax,nclus,iclus,rclus)
! generates spherical truncation cluster around chosen site

  implicit none

! number of bravais vectors (2D or 3D)
  integer*4, intent(in)  :: nbravais
! 3D form of bravais vectors for real and reciprocal lattice
  real*8,    intent(in)  :: a3d(3,3), b3d(3,3)
! basis vectors (3D)
  logical,   intent(in)  :: lcartesian
  integer*4, intent(in)  :: nbasis
  real*8,    intent(in)  :: rbasis(3,nbasis)
! lattice position
  integer*4, intent(in)  :: ibasis, ijk(nbravais)
! cutoff radius
  real*8,    intent(in)  :: rcut
! maximum number of sites to generate
  integer*4, intent(in)  :: nmax
! generated truncated cluster
  integer*4, intent(out) :: nclus, iclus(nmax)
  real*8,    intent(out) :: rclus(3,nmax)
! ----------------------------------------------------------------------
  real*8    :: distance, ri(3), rj(3), rij(3), rijkeep(3,nbasis)
  integer*4 :: i, j, k1, k2, k3, ikeep(nbasis), nkeep, nvec(3), ntot


  write(*,'("truncated_cluster:")')
  nclus = 0; rclus(:,:) = 0.d0
! first we have to find which layers/sublattices are reached by rcut
  call coupled_sublattices(nbravais,a3d,lcartesian,nbasis,rbasis,ibasis,rcut,nkeep,ikeep,rijkeep)
! how many multiples of the bravais vectors do we need?
  call set_nvec(nbravais,a3d,rcut,nkeep,nmax,nvec)
! location of the impurity
  ri(:) = (/0.d0,0.d0,0.d0/)
  do i=1,nbravais
    ri(:) = ri(:) + a3d(:,i)*ijk(i)
  end do
! shift the cluster by the lattice position of the impurity
  j = 1
  iclus(j) = ibasis
  rclus(:,j) = ri(:)
  write(*,'(" j, iclus=",2i8,"  rclus, distance=",4f16.8)')  j, iclus(j), rclus(:,j), 0.d0
! for each sublattice
  do i=1,nkeep
!   generate lattice vectors
    do k3=-nvec(3),nvec(3)
      do k2=-nvec(2),nvec(2)
        do k1=-nvec(1),nvec(1)
          rij(:) = rijkeep(:,i) + k1*a3d(:,1) + k2*a3d(:,2) + k3*a3d(:,3)
          distance = sqrt(dot_product(rij(:),rij(:)))
!         the impurity position is listed first
          if (distance < 1.d-8) cycle
          if (distance < rcut) then
            j = j + 1
            iclus(j) = ikeep(i)
!           shift the cluster by the lattice position of the impurity
            rclus(:,j) = rij(:) + ri(:)
            write(*,'(" j, iclus=",2i8,"  rclus, distance=",4f16.8)')  j, iclus(j), rclus(:,j), distance
          end if
        end do
      end do
    end do
  end do
! number of atoms in cluster
  nclus = j
! All done!
  end subroutine truncated_cluster


  subroutine coupled_sublattices(nbravais,a3d,lcartesian,nbasis,rbasis,ibasis,rcut,nkeep,ikeep,rijkeep)
! which sublattices are coupled by rcut

  implicit none

! number of bravais vectors (2D or 3D)
  integer*4, intent(in)  :: nbravais
! 3D form of bravais vectors for real lattice
  real*8,    intent(in)  :: a3d(3,3)
! basis vectors (3D)
  logical,   intent(in)  :: lcartesian
  integer*4, intent(in)  :: nbasis
  real*8,    intent(in)  :: rbasis(3,nbasis)
! lattice position
  integer*4, intent(in)  :: ibasis
! cutoff radius
  real*8,    intent(in)  :: rcut
! how many sublattices couple to selected one
  integer*4, intent(out) :: nkeep, ikeep(nbasis)
! distance between sublattices
  real*8,    intent(out) :: rijkeep(3,nbasis)
! ----------------------------------------------------------------------
  real*8    :: distance, ri(3), rj(3), rij(3)
  integer*4 :: i, j


! location of sublattice containing the origin of the cluster
  if (.not.lcartesian) then
    ri(:) = 0.d0
    do j=1,3
      ri(:) = ri(:) + a3d(:,j)*rbasis(j,ibasis)
    end do
  else
    ri(:) = rbasis(:,ibasis)
  end if
! loop over all sublattices to find which are coupled by rcut
  nkeep = 0
  do i=1,nbasis
    if (.not.lcartesian) then
      rj(:) = 0.d0
      do j=1,3
        rj(:) = rj(:) + a3d(:,j)*rbasis(j,i)
      end do
    else
      rj(:) = rbasis(:,i)
    end if
    rij(:) = rj(:) - ri(:)
!   reduce the sublattice distance vector to the wigner-seitz cell
    call fold_vector(nbravais,a3d,rij)
    distance = sqrt(dot_product(rij(:),rij(:)))
    if (distance < rcut) then
      nkeep = nkeep + 1
      ikeep(nkeep) = i
      rijkeep(:,nkeep) = rij(:)
      write(*,'(" ibasis, i=",2i4," ri, rj=",6f8.4," rij, dist, rcut=",5f8.4)') ibasis, i, ri(:), rj(:), rij(:), distance, rcut
    end if
  end do
! All done!
  end subroutine coupled_sublattices


  subroutine set_nvec(nbravais,a3d,rcut,nkeep,nmax,nvec)
! decides how many multiples of the bravais vectors are needed

  implicit none

! number of bravais vectors (2D or 3D)
  integer*4, intent(in)  :: nbravais
! 3D form of bravais vectors for real lattice
  real*8,    intent(in)  :: a3d(3,3)
! cutoff radius
  real*8,    intent(in)  :: rcut
! number of coupled sublattices
  integer*4, intent(in)  :: nkeep
! maximum number of vectors
  integer*4, intent(in)  :: nmax
! multiples of the bravais vectors to generate
  integer*4, intent(out) :: nvec(3)
! ----------------------------------------------------------------------
  integer*4 :: i, ntot
  real*8    :: length, old(3), new(3)


! First bravais vector:
! which multiple goes beyond rcut
  old(:)  = a3d(:,1)
  length  = sqrt(dot_product(a3d(:,1),a3d(:,1)))
  old(:)  = old(:)/length
  nvec(1) = ceiling(rcut/length)
!  write(*,'(" i, nvec=",2i8)') 1, nvec(1)
! Remaining bravais vectors:
  do i=2,nbravais
!   make it perpendicular to previous bravais vector
    new(:)  = a3d(:,i) - dot_product(a3d(:,i),old(:))*old(:) 
    length  = sqrt(dot_product(new(:),new(:)))
    old(:)  = new(:)/length
    nvec(i) = ceiling(rcut/length)
!    write(*,'(" i, nvec=",2i8)') i, nvec(i)
  end do
  if (nbravais == 2) nvec(3) = 0
! total number of vectors to generate
  ntot = nkeep*product(2*nvec(1:nbravais)+1)
  write(*,'(" ntot=",i4)') ntot
  if (ntot > nmax) stop 'set_nvec: increase nmax'
! All done!
  end subroutine set_nvec


  subroutine bravais3d(nbravais,rbravais,alat,a3d,b3d)
! generate 3D form of bravais vectors for real and reciprocal lattice

  implicit none

! bravais vectors (2D or 3D)
  integer*4, intent(in)    :: nbravais
  real*8,    intent(in)    :: rbravais(nbravais,nbravais)
! unit of length
  real*8,    intent(in)    :: alat
! 3D form of Bravais vectors for real and reciprocal lattice
  real*8,    intent(out)   :: a3d(3,3), b3d(3,3)
! ----------------------------------------------------------------------
  integer*4 :: i, j
  real*8    :: rlen, volume

!  write(*,'(/,"bravais3d:")')
! 2D case
  if (nbravais == 2) then
    a3d(1:2,1:2) = rbravais(1:2,1:2)
    a3d(3,1:2) = 0.d0
    a3d(:,3) = (/0.d0,0.d0,1.d0/)
! 3D case
  else
    a3d(:,:) = rbravais(:,:) 
  end if
!  write(*,'(" Bravais vectors of real lattice:")')
!  write(*,'("  a1=",3f16.8)') a3d(:,1)
!  write(*,'("  a2=",3f16.8)') a3d(:,2)
!  write(*,'("  a3=",3f16.8)') a3d(:,3)
! find bravais vectors of reciprocal lattice
  volume = 0.d0
  volume = volume + a3d(1,1)*(a3d(2,2)*a3d(3,3) - a3d(3,2)*a3d(2,3))
  volume = volume + a3d(2,1)*(a3d(3,2)*a3d(1,3) - a3d(1,2)*a3d(3,3))
  volume = volume + a3d(3,1)*(a3d(1,2)*a3d(2,3) - a3d(2,2)*a3d(1,3))
  volume = abs(volume)
!  write(*,'(" Bravais vectors of reciprocal lattice:")')
  b3d(1,1) = a3d(2,2)*a3d(3,3) - a3d(3,2)*a3d(2,3)
  b3d(2,1) = a3d(3,2)*a3d(1,3) - a3d(1,2)*a3d(3,3)
  b3d(3,1) = a3d(1,2)*a3d(2,3) - a3d(2,2)*a3d(1,3)
  b3d(:,1) = b3d(:,1)/volume
!  write(*,'("  b1=",3f16.8)') b3d(:,1)
  b3d(1,2) = a3d(2,3)*a3d(3,1) - a3d(3,3)*a3d(2,1)
  b3d(2,2) = a3d(3,3)*a3d(1,1) - a3d(1,3)*a3d(3,1)
  b3d(3,2) = a3d(1,3)*a3d(2,1) - a3d(2,3)*a3d(1,1)
  b3d(:,2) = b3d(:,2)/volume
!  write(*,'("  b2=",3f16.8)') b3d(:,2)
  b3d(1,3) = a3d(2,1)*a3d(3,2) - a3d(3,1)*a3d(2,2)
  b3d(2,3) = a3d(3,1)*a3d(1,2) - a3d(1,1)*a3d(3,2)
  b3d(3,3) = a3d(1,1)*a3d(2,2) - a3d(2,1)*a3d(1,2)
  b3d(:,3) = b3d(:,3)/volume
!  write(*,'("  b3=",3f16.8)') b3d(:,3)
!  write(*,'(" a_i . b_j:")')
!  do i=1,3
!    write(*,'(3f16.8)') (dot_product(a3d(:,i),b3d(:,j)),j=1,3)
!  end do
! All done!
  end subroutine bravais3d


  subroutine fold_vector(nbravais,a3d,rxyz)
! adds to input vector the linear combination of bravais vectors that minimizes its length

  implicit none

! number of bravais vectors <=> 2D or 3D lattice
  integer*4, intent(in)    :: nbravais
! 3D form of bravais vectors for real lattice
  real*8,    intent(in)    :: a3d(3,3)
! vector to be folded
  real*8,    intent(inout) :: rxyz(3)
! ----------------------------------------------------------------------
  integer*4 :: k1, k2, k3, nvec(3)
  real*8    :: length, best, rbest(3), rnew(3)

  nvec(:) = 10
  if (nbravais == 2) nvec(3) = 0
  best = sqrt(dot_product(rxyz(:),rxyz(:)))
  rbest(:) = rxyz(:)
  do k1=-nvec(1),nvec(1)
    do k2=-nvec(2),nvec(2)
      do k3=-nvec(3),nvec(3)
        rnew(:) = rxyz(:) - k1*a3d(:,1) - k2*a3d(:,2) - k3*a3d(:,3)
        length = sqrt(dot_product(rnew(:),rnew(:)))
        if (length < best) then
          best = length
          rbest(:) = rnew(:)
        end if
      end do
    end do
  end do
  rxyz(:) = rbest(:)
! All done!
  end subroutine fold_vector


  subroutine lattice_vector(nbravais,a3d,ijk,r0,rxyz)
! generates a 3D lattice vector from multiples of the bravais vectors

  implicit none

! number of bravais vectors <=> 2D or 3D lattice
  integer*4, intent(in)  :: nbravais
! 3D form of bravais vectors for real lattice
  real*8,    intent(in)  :: a3d(3,3)
! multiples of the bravais vectors
  integer*4, intent(in)  :: ijk(3)
! shift
  real*8,    intent(in)  :: r0(3)
! vector to be folded
  real*8,    intent(out) :: rxyz(3)

  rxyz(:) = r0(:) + ijk(1)*a3d(:,1) + ijk(2)*a3d(:,2)
  if (nbravais == 3) rxyz(:) = rxyz(:) + ijk(3)*a3d(:,3)
! All done!
  end subroutine lattice_vector


  end program scoef_maker
