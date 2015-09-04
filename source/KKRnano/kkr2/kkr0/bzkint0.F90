      subroutine bzkint0(naez, rbasis, bravais, recbv, nsymat, isymindex, &
                         dsymll, intervx, intervy, intervz, ielast, ez, kmesh, maxmesh, maxmshd, lmax, iemxd, krel, kpoibz, ekmd)
!
      implicit none
      integer, parameter :: nsymaxd=48

      integer lmax
      integer iemxd
      integer krel
      integer kpoibz
      integer ekmd
      integer naez,nsymat
      integer intervx,intervy,intervz,maxmesh,maxmshd,ielast
      double complex dsymll((lmax+1)**2,(lmax+1)**2,nsymaxd), ez(iemxd) ! dsymll(lmmaxd,lmmaxd,nsymaxd),ez(iemxd)
      double precision bravais(3,3), rbasis(3,naez), recbv(3,3), rsymat(64,3,3)
      integer isymindex(nsymaxd)
      integer kmesh(iemxd)

      logical, external :: test
      external :: bzkmesh, findgroup, pointgrp, symtaumat
      
      integer iprint
      logical lirr
      character(len=10) :: rotname(64)
      integer :: lmmaxd

      lmmaxd = (lmax+1)**2

      write(6,'(79(1h=),/,15x,a)') 'BZKINT0: finding symmetry, setting BZ integration'
      write(6,'(79(1h=),/)')

      call pointgrp(rsymat, rotname) ! set up the default symmetries and their names
      call findgroup(bravais, recbv, rbasis, naez, rsymat, rotname, isymindex, nsymat, naez)

      lirr = .true.
      iprint = 0    
      if (test('TAUSTRUC')) iprint = 2

! --> test: full bz integration
      if ( test('fullBZ  ') ) then
         nsymat = 1
         lirr = .false.
         write(6,'(8x,2a,/)') 'Test option < fullBZ > : overriding NSYMAT,', ' generate full BZ k-mesh'
      endif

! --> generate bz k-mesh
      call bzkmesh(intervx, intervy, intervz, maxmesh, lirr, bravais, recbv, nsymat, rsymat, isymindex, ielast, ez, kmesh, iprint, maxmshd, iemxd, kpoibz, ekmd)
!
      call symtaumat(rotname, rsymat, dsymll, nsymat, isymindex, naez, lmmaxd, naez, lmax+1, krel, iprint, nsymaxd)
!
! now dsymll hold nsymat symmetrization matrices
! ----------------------------------------------------------------------
      endsubroutine bzkint0
