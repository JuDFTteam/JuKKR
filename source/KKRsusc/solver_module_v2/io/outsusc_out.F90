  subroutine outsusc_out()
! write initial info to output file
  use global

  implicit none

  integer(kind=i4b) :: istart, iend, ng, il, ib, ispin, iatom 

  if (noparameters) stop 'open_outsusc: run init_param first!'
  open(file='outsusc.dat',unit=iomain,status='replace')
  write(iomain,'("------------------------------------------------------------------------")')
  write(iomain,'("         Data for dynamical susceptibility calculation")')
  write(iomain,'("------------------------------------------------------------------------")')
  write(iomain,'(" ikkr, ldos, lsoc, lldau, lbfield, lsusc, lrot, lfit=")')
  write(iomain,'(i4,7l4)') ikkr, ldos, lsoc, lldau, lbfield, lsusc, lrot, lfit
  write(iomain,'(" ibasis, ibasismethod, basistol, lbasisfreeze=")')
  write(iomain,'(2i4,es12.4,l4)') ibasis, ibasismethod, basistol, lbasisfreeze
  write(iomain,'(" idos, nedos, e0dos, e1dos=")')
  write(iomain,'(i4,i8,4es12.4)') idos, nedos, dose0, dose1
  write(iomain,'(" lrhomat, ldaumix, ldauitc=")')
  write(iomain,'(l4,es12.4,i8)') lrhomat, ldaumix, ldauitc
  write(iomain,'(" ldynsusc,nomega, omegamin, omegamax, domega, nlmax0=")')
  write(iomain,'(l4,i8,5es12.4,i4)') ldynsusc, nomega, omegamin, omegamax, domega, nlmax0
  write(iomain,'(" lkha, lkxc, lcartesian, lanalytic, lnonanalytic, lenhanced, itermax, lambdamix=")')
  write(iomain,'(6l4,i4,es12.4)') lkha, lkxc, lcartesian, lanalytic, lnonanalytic, lenhanced, itermax, lambdamix
  write(iomain,'(" ispinrot, urot(1:3), dirmix=")')
  write(iomain,'(i4,4es12.4)') ispinrot, urot(1:3), dirmix3
  write(iomain,'(" ifit, numd, dend, eshift, lregf, fudge=")')
  write(iomain,'(3i4,2es12.4,l4,es12.4)') ifit, numd, dend, eshift, lregf, fudge
  write(iomain,'("------------------------------------------------------------------------")')
  write(iomain,'(" nesusc, nasusc, nbmax, nlmax, isra, nrmax, nsmax=")')
  nbmax = maxval(iwsusc)
  write(iomain,'(8i6)') nesusc, nasusc, nbmax, nlmax, isra, nrmax, nsmax
  write(iomain,'(" ngroup, igroup(1:ngroup)=")')
  write(iomain,'(i4,/,1000i4)') ngroup, igroup(1:ngroup)
  iend = 0
  do ng=1,ngroup
    istart = iend + 1
    iend   = iend + igroup(ng)
    write(iomain,'(" atom labels in group")')  
    write(iomain,'(1000i6)') iasusc(istart:iend)
    write(iomain,'(" isoc, socscaling, ildau, ueff(0:lmax), jeff(0:lmax)=")')
    write(iomain,'(i4,es12.4,i4,8es12.4)') isoc(istart), socscaling(istart), ildau(istart), ueff(0:nlmax,istart), jeff(0:nlmax,istart)
    write(iomain,'(" ibfield, blen, bdir(1:3)=")')
    write(iomain,'(i4,4es12.4)') ibfield(istart), blen(istart), bdir(:,istart)
    write(iomain,'(" isusc, ikha, ikxc=")')
    write(iomain,'(4i4)') isusc(istart), ikha(istart), ikxc(istart)
    write(iomain,'(" nspin, iwsusc(0:lmax), ewsusc(1:nbmax,0:lmax,1:nspin)")')
    write(iomain,'(i2,20i6)') issusc(istart), iwsusc(:,:,istart)
  end do
! All done!
  end subroutine outsusc_out
