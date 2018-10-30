  subroutine save_rmesh(ia,nr0,nr1,rr,rs,wr,z,vup,vdn,ncup,ncdn,numpan,numrcut,my_rank)
! Save radial mesh and potentials
  use global

  implicit none

! Start and end of radial mesh; susceptibility atom
  integer(kind=i4b), intent(in) :: nr0, nr1, ia
! Radial mesh, powers of r and radial integration weights
  real(kind=r8b),    intent(in) :: rr(nrmax), wr(nrmax), rs(nrmax,0:nlmax)
! Atomic number, spin up and spin down potentials, core densities
  real(kind=r8b),    intent(in) :: z, vup(nrmax), vdn(nrmax), ncup(nrmax), ncdn(nrmax)
! Number of panels, start and end of each
  integer(kind=i4b), intent(in) :: numpan, numrcut(numpan+1)
! Mpi
  integer(kind=i4b), intent(in) :: my_rank
! -----------------------------------------------------------------
  integer(kind=i4b) :: nr, ir, ip
  real(kind=r8b)    :: magdummy(3)
  complex(kind=c8b) :: work(nrmax), norm
  character*12      :: filename

  nr = nr1 - nr0 + 1
  nrpts(ia)  = nr
  nrpts0(ia) = nr0
  nrpts1(ia) = nr1
  magdummy = (/0.d0,0.d0,1.d0/)
! Check how the radial mesh is defined, what is the first point
  rmesh(1:nr,ia) = rr(nr0:nr1)
  rsmesh(1:nr,0:nlmax,ia) = rs(nr0:nr1,0:nlmax)
! Original radial integration weights
  drmesh(1:nr,ia) = wr(nr0:nr1)
! Modified weights for projection
  drproj(1:nr,ia) = drmesh(1:nr,ia)
! Number of panels > 1 
  npanat(ia) = numpan 
  ircutat(1:numpan+1,ia) = numrcut(1:numpan+1)
! Save atomic number and potentials
  zat(ia) = z
  vr(1:nr,ia) = 0.5d0*(vup(nr0:nr1) + vdn(nr0:nr1))
!  vr(1:nr,ia) = vr(1:nr,ia) - 2.d0*z/rr(nr0:nr1)
  br(1:nr,ia) = 0.5d0*(vup(nr0:nr1) - vdn(nr0:nr1))
! Core charge and magnetization densities
  nrc(1:nr,ia) = ncup(nr0:nr1) + ncdn(nr0:nr1)
  mrc(1:nr,ia) = ncup(nr0:nr1) - ncdn(nr0:nr1)
! Write potential file for susc
  if (my_rank == 0) then
    if (lhdio) then    
      write(filename,'("susc",i4.4,".pot")') ia
      open(file=filename,unit=iofile,status='replace')
      write(iofile,'("# Susc potential: z, nr, lmax; then magdummy; then rmesh, rs, dr, vr, br, nrc, mrc")')
      write(iofile,'("# ",f6.1,2i8)') z, nr, nlmax  ! pot can be plotted
      write(iofile,'("# ",3f12.8)') magdummy  ! magnetization direction
      do ir=1,nr
        write(iofile,'(100es16.8)') rmesh(ir,ia), rsmesh(ir,0:nlmax,ia), drmesh(ir,ia), vr(ir,ia), br(ir,ia), nrc(ir,ia), mrc(ir,ia)
      end do
      close(iofile)
    end if 
  end if ! my_rank
! Output mesh and potentials
  normesh(ia) = .false.
! All done!
  end subroutine save_rmesh
