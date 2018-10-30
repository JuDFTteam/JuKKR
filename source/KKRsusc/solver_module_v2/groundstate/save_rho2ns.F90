  subroutine save_rho2ns(lmaxpot,irmkd,lmpotd,natypd,rho2ns)
! Save non-spherical charge and magnetization densities
  use global

  implicit none

! lmax for potential; dimensions of rho2ns array
  integer(kind=i4b), intent(in) :: lmaxpot, irmkd, lmpotd, natypd
! Radial mesh, powers of r and radial integration weights
  real(kind=r8b),    intent(in) :: rho2ns(irmkd,lmpotd,natypd,nsmax)
! -----------------------------------------------------------------
  real(kind=r8b), parameter :: sqrt4pi = 4.d0*sqrt(atan(1.d0))
  real(kind=r8b), parameter :: tol = 1.d-7
  integer(kind=i4b) :: ia, ih, nr0, nr1, nr, ir, lmpot, ilm, lm
  integer(kind=i4b) :: nonzero(lmpotd)
  real(kind=r8b)    :: qlm(lmpotd), mlm(lmpotd), q, m
  complex(kind=c8b) :: work(nrmax), norm
  character*14      :: filename
  complex(kind=c8b), external :: radint

  lmpot = (lmaxpot + 1)**2
  if (lmpot > lmpotd) stop 'save_rho2ns: lmpot > lmpotd'
  write(*,*) nasusc
! **************
  do ia=1,nasusc
! **************
    if (normesh(ia))  stop 'save_rho2ns: save rmesh first!'
    ih  = iasusc(ia)
    nr  = nrpts(ia)
    nr0 = irmkd - nr + 1 !nrpts0(ia)
    nr1 = irmkd !nrpts1(ia)
!   Check non-zero multipoles
    ilm = 0
    do lm=1,lmpot
      work(1:nr) = sqrt4pi*rho2ns(nr0:nr1,lm,ih,1)
      q = real(radint(nr,work(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia)))
      work(1:nr) = sqrt4pi*rho2ns(nr0:nr1,lm,ih,2)
      m = real(radint(nr,work(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia)))
      if (abs(q) > tol .or. abs(m) > tol) then
        ilm = ilm + 1
        nonzero(ilm) = lm
        qlm(ilm) = q; mlm(ilm) = m
      end if
    end do
!   Write file for rho2ns
    write(filename,'("rho2ns",i4.4,".dat")') ia
    open(file=filename,unit=iofile,status='replace')
    write(iofile,'("# rho2ns: nr, lmaxpot, number of multipoles; then nonzero multipoles; then qlm and mlm; then rmesh and rho2ns")')
    write(iofile,'("# ",3i8)') nr, lmaxpot, ilm
    write(iofile,'("# ",100i4)') nonzero(1:ilm)
    write(iofile,'("# ",100f12.8)') qlm(1:ilm)
    write(iofile,'("# ",100f12.8)') mlm(1:ilm)
    do ir=1,nr
      write(iofile,'(100es16.8)') rmesh(ir,ia), sqrt4pi*rho2ns(nr0+ir-1,nonzero(1:ilm),ih,:)
    end do
    close(iofile)
! ******
  end do
! ******
! All done!
  end subroutine save_rho2ns
