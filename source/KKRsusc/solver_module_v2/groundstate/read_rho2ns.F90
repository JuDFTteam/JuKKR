  subroutine read_rho2ns()
! Read groundstate non-spherical charge and magnetization densities
! There is an extra factor of sqrt(4pi) included
! Sets the initial spin quantization axis before applying corrections
  use global

  implicit none

  integer(kind=i4b) :: ia, nr, ir, lmaxpot, ilm, nlm
  integer(kind=i4b) :: nonzero(lmmax2)
  real(kind=r8b)    :: r, fac
  character*14      :: filename
  character*1       :: dummy

  write(iodb,'(/,"read_rho2ns: setting spin axes!")')
! **************
  do ia=1,nasusc
! **************
!   Read file for rho2ns
    write(filename,'("rho2ns",i4.4,".dat")') ia
    open(file=filename,unit=iofile,status='old')
    read(iofile,*)  ! header
    read(iofile,*) dummy, nr, lmaxpot, nlm
    if (lhdio) write(iodb,'("# ",3i8)') nr, lmaxpot, nlm
    if (nr /= nrpts(ia))   stop 'read_rho2ns: nr /= nrpts'
    if (lmaxpot /= nlmax2) stop 'read_rho2ns: lmaxpot /= nlmax2'
    read(iofile,*) dummy, nonzero(1:nlm)
    if (lhdio) write(iodb,'("# read_rho2ns: ia=",i4)') ia
    if (lhdio) write(iodb,'("# lm =",100i4)') nonzero(1:nlm)
!   Prepare pointers
    nlmpot(ia) = nlm
    i2lmpot(1:nlm,ia) = nonzero(1:nlm)
!   Nonzero multipoles
    gs_qlm(:,ia) = 0.d0; gs_mlm(:,ia) = 0.d0
    read(iofile,*) dummy, gs_qlm(nonzero(1:nlm),ia)
    if (lhdio) write(iodb,'("# qlm=",100f10.6)') gs_qlm(nonzero(1:nlm),ia)
    read(iofile,*) dummy, gs_mlm(nonzero(1:nlm),ia)
    if (lhdio) write(iodb,'("# mlm=",100f10.6)') gs_mlm(nonzero(1:nlm),ia)
!   Nonzero components of rho2ns
    old_rho2ns(:,:,:,ia) = 0.d0
    do ir=1,nr
      read(iofile,*) r, old_rho2ns(ir,nonzero(1:nlm),:,ia)
    end do
    close(iofile)
!   Spin quantization axis
!   - collinear: no need to do anything
!   - noncollinear: the spin direction is needed to set up the spin rotations correctly
!   - New noncollinear magnetism from kkrflex
    if (.not.lsoc_new) magdir(:,ia) = (/0.d0,0.d0,1.d0/)
    fac = 1.d0
    if (lrot .and. gs_mlm(1,ia) < 0.d0) then
      magdir(:,ia) = (/0.d0,0.d0,-1.d0/)
      gs_mlm(:,ia) = -gs_mlm(:,ia)
      old_rho2ns(:,:,:,ia) = -old_rho2ns(:,:,:,ia)
      br(:,ia) = -br(:,ia)
    end if
    write(iodb,'("  ia=",i4," magdir=",3f8.4)') ia, magdir(:,ia)
! ******
  end do
! ******
  write(iodb,*)
! All done!
  end subroutine read_rho2ns
