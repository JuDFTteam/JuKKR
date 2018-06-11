subroutine drvrho_qdos(ldorhoef, rho2ns, r2nef, den, dmuorb, rhotborb, iecurr, &
  eryd, we, ielast, gmatll, vt, bt, r, drdi, r2drdi, zat, jws, ishift, solver, &
  soctl, ctl, qmtet, qmphi, itermvdir, mvevil, mvevilef, lmmaxd, lmaxd, irmd, &
  lmpotd, iemxd, nmvecmax, i1, nqdos) ! qdos ruess
!   ********************************************************************
!   *                                                                  *
!   * driving routine to call relativistic routines                    *
!   *          < SSITE >, < SCFCHRDNS >, < CALCMVEC >                  *
!   * to calculate the charge and spin density in the REL mode         *
!   * v.popescu, munich, may 2004                                      *
!   *                                                                  *
!   ********************************************************************
  use :: mod_types, only: t_tgmat
      Use mod_datatypes, Only: dp
  implicit none

! PARAMETER definitions
  integer :: nrmax
  parameter (nrmax=900)
  integer :: nlamax, nqmax, ntmax, nmmax
  parameter (nlamax=1, nqmax=1, ntmax=1, nmmax=1)
  integer :: nlmax, nkmmax, nmuemax, nkmpmax, nkmax, linmax
  parameter (nlmax=5) ! this should be >= LMAXD + 1
  parameter (nkmmax=2*nlmax**2, nkmax=2*nlmax-1)
  parameter (nkmpmax=nkmmax+2*nlmax, nmuemax=2*nlmax)
  parameter (linmax=2*nlmax*(2*nlmax-1))
  complex (kind=dp) :: cone, czero
  parameter (cone=(1.0d0,0.0d0), czero=(0.0d0,0.0d0))
  real (kind=dp) :: dzero
  parameter (dzero=0.0d0)

! Dummy arguments
  integer :: lmaxd, lmmaxd, irmd, ielast
  integer :: zat(ntmax), jws(nmmax), ishift
  integer :: lmpotd, iemxd, i1
  logical :: ldorhoef
  complex (kind=dp) :: we, eryd
  real (kind=dp) :: rho2ns(irmd, lmpotd, 2), r2nef(irmd, lmpotd, 2)
!  real (kind=dp) VT(NRMAX,NTMAX),BT(NRMAX,NTMAX)
  real (kind=dp) :: vt(nrmax), bt(nrmax)
  real (kind=dp) :: r(nrmax, nmmax), r2drdi(nrmax, nmmax)
  real (kind=dp) :: drdi(nrmax, nmmax), soctl(ntmax, nlmax)
  real (kind=dp) :: ctl(ntmax, nlmax)
  complex (kind=dp) :: gmatll(lmmaxd, lmmaxd, iemxd), den(0:lmaxd+1, 2*ielast)
! l-resolved orbital polarisation 
  complex (kind=dp) :: dmuorb(0:lmaxd, 3)
! orbital density
  real (kind=dp) :: rhotborb(irmd)

! Local variables
  real (kind=dp) :: ameopo(nkmmax, nkmmax, nlamax, 3), at(nrmax, nlamax, 3, ntmax), &
    bcor(ntmax), bcors(ntmax), conc(ntmax), dos(ntmax), dosi(ntmax), efermi, &
    hff(ntmax), hffi(ntmax), mueorb, muespn, nvaltot, omt(ntmax), omti(ntmax), &
    qel(ntmax), rhoorb(nrmax, ntmax), rhochr(nrmax, ntmax), &
    rhospn(nrmax, ntmax)
  real (kind=dp) :: shftef, smt(ntmax), smti(ntmax), pi, sqpi, totdos
  complex (kind=dp) :: bzj(linmax, ntmax), bzz(linmax, ntmax), dosint(nlmax, ntmax), &
    dosl0(nlmax, ntmax), dosm(nmuemax), dzj(linmax, ntmax), &
    dzz(linmax, ntmax), eband, ebandt(ntmax), hffint(nlmax, ntmax), &
    hffl0(nlmax, ntmax), hffm(nmuemax), msst(nkmmax, nkmmax, ntmax), &
    omtint(nlmax, ntmax), omtl0(nlmax, ntmax), omtm(nmuemax), &
    ozj(linmax, ntmax), ozz(linmax, ntmax), p, qzj(linmax, ntmax), &
    qzz(linmax, ntmax), smtint(nlmax, ntmax), smtl0(nlmax, ntmax), &
    smtm(nmuemax), szj(linmax, ntmax), szz(linmax, ntmax)
  complex (kind=dp) :: taut(nkmmax, nkmmax, ntmax), omtls0(nlmax, ntmax, 2)
  complex (kind=dp) :: ozzs(linmax, ntmax, 2), ozjs(linmax, ntmax, 2)
  complex (kind=dp) :: tautlin(linmax, ntmax), tsst(nkmmax, nkmmax, ntmax), &
    tsstlin(linmax, ntmax), tzj(linmax, ntmax), tzz(linmax, ntmax)
  logical :: calcint, getirrsol
  real (kind=dp) :: cgc(nkmpmax, 2)
  real (kind=dp) :: gdia(nkmmax), gmdia(nkmmax), goff(nkmmax), gmoff(nkmmax)
  real (kind=dp) :: fdia(nkmmax), fmdia(nkmmax), foff(nkmmax), fmoff(nkmmax)
  integer :: i, iecurr, ihyper, ikm1lin(linmax), ikm2lin(linmax), il, &
    imt(ntmax), imue, ip, iprint, iq, iqat(nqmax, ntmax), irel, it, iwrirrwf, &
    iwrregwf, j, lin, lopt(ntmax), mmax, nat(ntmax), netab, nkm, nkmq(nqmax), &
    nl, nlinq(nqmax), nlq(nqmax), nt, nucleus
  integer :: nfilcbwf, iol
  integer :: nsollm(nlmax, nmuemax), ltab(nmuemax), lbtab(nmuemax), &
    kaptab(nmuemax), nmuetab(nmuemax)
  character (len=10) :: solver
  character (len=4) :: txtt(ntmax)
  complex (kind=dp) :: w1(lmmaxd, lmmaxd)
  integer :: icall, iec
!qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos
  complex (kind=dp) :: gmat0(lmmaxd, lmmaxd) !qdos ruess
  integer :: nqdos, irec, ipoint !qdos ruess 
!qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos
  intrinsic :: atan, sqrt

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ITERMDIR

  logical :: itermvdir, splitss
  integer :: nmvecmaxd, nmvecmax
  parameter (nmvecmaxd=4)
  real (kind=dp) :: amemvec(nkmmax, nkmmax, 3, nmvecmaxd), fact(0:100)
  integer :: imkmtab(nkmmax), ikmllim1(nkmmax), ikmllim2(nkmmax)
  character (len=1) :: txtl(0:nlmax)
  integer :: igrid(2), iepath, nepath

  real (kind=dp) :: qmtet, qmphi ! ARG. LIST
  real (kind=dp) :: qmphiloc(nqmax), qmtetloc(nqmax) ! DUMMY

  complex (kind=dp) :: bmvevdl0(nlmax, ntmax, 3, nmvecmax), &
    bmvevil1(nlmax, ntmax, 3, nmvecmax), mvevdl0(nlmax, ntmax, 3, nmvecmax), &
    mvevil1(nlmax, ntmax, 3, nmvecmax)
  complex (kind=dp) :: mvevil(0:lmaxd, 3, nmvecmax) ! OUTPUT
  complex (kind=dp) :: mvevilef(0:lmaxd, 3, nmvecmax) ! OUTPUT
!.. dummy arrays
  complex (kind=dp) :: mezj(nkmmax, nkmmax, ntmax, nmvecmax), &
    mezz(nkmmax, nkmmax, ntmax, nmvecmax)

! ITERMDIR
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!..
!.. External Subroutines ..
  external :: amemagvec, calccgc, calcgf, calcmvec, cinit, ikmlin, rinit, &
    scfchrdns, ssite, zcopy, zgemm

  data icall/0/

  save :: icall, ikm1lin, ikm2lin, gdia, gmdia, goff, lopt, nlq, nkmq, iqat, &
    irel, bcor, bcors, qel, nat, conc, txtt, imt, shftef, nvaltot, nkm, &
    ihyper, iprint, it, iq, nl, nt, nucleus, cgc, iwrregwf, iwrirrwf, calcint, &
    getirrsol, nfilcbwf, pi, sqpi, amemvec, imkmtab, ikmllim1, ikmllim2, fact, &
    splitss, txtl, igrid, iepath, nepath

  icall = icall + 1

!=======================================================================
!       initialise relativistic and dummy variables and SAVE them
!=======================================================================
  if (icall==1) then

    if (lmaxd>nlmax-1) then
      write (6, *) ' LMAXD = ', lmaxd, ' > NLMAX-1 = ', nlmax - 1
      stop ' Increase NLMAX in < DRVRHO > '
    end if

    if (irmd>nrmax) then
      write (6, *) ' IRMD = ', irmd, ' > NRMAX = ', nrmax
      write (6, *) ' Increase NRMAX in < sprkkr_rmesh.dim > '
      stop ' In < DRVRHO > '
    end if

    if (nmvecmax>nmvecmaxd) then
      write (6, *) ' NMVECMAX = ', nmvecmax, ' > NMVECMAXD ', nmvecmaxd
      write (6, *) ' Increase NVECMAXD in < DRVRHO > ', &
        'or reduce NMVECMAX in < main1c > '
      stop ' In < DRVRHO > '
    end if

    iprint = 0
    nl = lmaxd + 1

    do i = 1, nmuemax
      ltab(i) = i/2
      if (2*ltab(i)==i) then
        lbtab(i) = ltab(i) - 1
        kaptab(i) = ltab(i)
      else
        lbtab(i) = ltab(i) + 1
        kaptab(i) = -ltab(i) - 1
      end if
      nmuetab(i) = 2*abs(kaptab(i))
    end do

    do il = 1, nlmax
      mmax = 2*il
      do imue = 1, mmax
        if ((imue==1) .or. (imue==mmax)) then
          nsollm(il, imue) = 1
        else
          nsollm(il, imue) = 2
        end if
      end do
    end do

    call ikmlin(iprint, nsollm, ikm1lin, ikm2lin, nlmax, nmuemax, linmax, &
      nlmax)

    call calccgc(ltab, kaptab, nmuetab, cgc, nkmax, nmuemax, nkmpmax)

    call calcgf(nkmax, cgc, gdia, gmdia, goff, gmoff, fdia, fmdia, foff, &
      fmoff, ltab, lbtab, kaptab, nmuetab, nmuemax, nkmmax, nkmpmax)

    do it = 1, ntmax
      bcor(it) = 0d0
      bcors(it) = 0d0
      qel(it) = 0d0
      nat(it) = 1
      conc(it) = 1d0
      txtt(it) = '    '
      imt(it) = 1
      lopt(it) = -1 ! this should change for Brooks' OP
    end do

    do iq = 1, nqmax
      nlq(iq) = nl
      nkmq(iq) = lmmaxd
      nlinq(iq) = 2*nlq(iq)*(2*nlq(iq)-1)
      iqat(iq, 1) = 1
    end do

    irel = 3
    shftef = 0d0
    efermi = 0d0
    nvaltot = 0
    pi = 4d0*atan(1d0)
    sqpi = sqrt(pi)

    nkm = lmmaxd
    ihyper = 0
    it = 1
    nt = 1
    iq = 1
    nucleus = 0

    iwrregwf = 1
    iwrirrwf = 1
    calcint = .true.
    getirrsol = .true.
    nfilcbwf = 87
!     Length in Bytes
    iol = 8*4 + 3 + (16*4*nrmax)
    open (nfilcbwf, status='SCRATCH', form='UNFORMATTED', access='DIRECT', &
      recl=iol)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      ITERMDIR

    if (itermvdir) then
      splitss = .false.
      fact(0) = 1.0d0
      do i = 1, 100
        fact(i) = fact(i-1)*dble(i)
      end do

      call amemagvec(irel, iprint+1, nkm, amemvec, ikmllim1, ikmllim2, &
        imkmtab, cgc, nlmax, nkmmax, nkmpmax, nmvecmax)

      do i = 0, nlmax
        txtl(i) = ' '
      end do
      igrid(1) = 5
      iepath = 1
      nepath = 1
    end if

!      ITERMDIR
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  end if ! ICALL.EQ.1
!=======================================================================

  call ssite(iwrregwf, iwrirrwf, nfilcbwf, calcint, getirrsol, soctl, ctl, &
    eryd, p, ihyper, iprint, ikm1lin, ikm2lin, nlq, nkmq, nlinq, nt, nkm, &
    iqat, tsst, msst, tsstlin, dzz, dzj, szz, szj, ozz, ozj, bzz, bzj, qzz, &
    qzj, tzz, tzj, vt, bt, at, zat, nucleus, r, drdi, r2drdi, jws, imt, &
    ameopo, lopt, solver, cgc, ozzs, ozjs, nlmax, nqmax, linmax, nrmax, nmmax, &
    ntmax, nkmmax, nkmpmax, nlamax)

!-----------------------------------------------------------------------
!     get charge density
!-----------------------------------------------------------------------

  netab = iecurr + 1
  iec = iecurr

! Loop over all qdos points specified in qvec.dat
  do ipoint = 1, nqdos ! qdos ruess
!                                                                    ! qdos ruess
! Read in Green function; remember that for the rel. case, nspin = 1 ! qdos ruess
! (without qdos, IPOINT=NQDOS=1)                                     ! qdos ruess
    irec = ipoint + nqdos*(iecurr-1) + nqdos*ielast*(i1-1) ! qdos ruess
    if (t_tgmat%gmat_to_file) then
      read (69, rec=irec) gmat0 ! qdos ruess
    else
      gmat0(:, :) = t_tgmat%gmat(:, :, irec)
    end if
    gmatll(:, :, iecurr) = gmat0(:, :) ! qdos ruess
!                                                                    ! qdos ruess

!-------- GET TAU MATRIX ------------------------
!         TAUT = t G t + t

! ---> taut = t

    do j = 1, nkm
      call zcopy(nkm, tsst(1,j,it), 1, taut(1,j,it), 1)
    end do

! ---> w1 = G * t

    call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, gmatll(1,1,iecurr), &
      lmmaxd, tsst(1,1,it), nkmmax, czero, w1, lmmaxd)

! ---> taut = t * G * t + t = t * w1 + taut

    call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, tsst(1,1,it), nkmmax, &
      w1, lmmaxd, cone, taut(1,1,it), nkmmax)

! ---> store taut in linear array tautlin

    do lin = 1, nlinq(iq)
      tautlin(lin, it) = taut(ikm1lin(lin), ikm2lin(lin), it)
    end do

    call rinit(nrmax, rhochr(1,it))
    call rinit(nrmax, rhospn(1,it))
    call rinit(nrmax, rhoorb(1,it))
    call cinit(nlmax, omtl0(1,it))
    do lin = 1, 2
      call cinit(nlmax, omtls0(1,it,lin))
    end do

    call scfchrdns(nfilcbwf, r2drdi, jws, imt, shftef, totdos, muespn, mueorb, &
      irel, iprint, nt, nl, nkm, eryd, we, efermi, iec, netab, dos, smt, omt, &
      hff, dosi, smti, omti, hffi, dosm, dosl0, dosint, smtm, smtl0, smtint, &
      omtm, omtl0, omtint, hffm, hffl0, hffint, bcor, bcors, dzz, dzj, szz, &
      szj, ozz, ozj, bzz, bzj, ozzs, ozjs, omtls0, tautlin, nvaltot, txtt, &
      conc, nat, rhochr, rhospn, rhoorb, qel, gdia, gmdia, goff, ntmax, nlmax, &
      nmuemax, linmax, nrmax, nmmax, nkmmax, eband, ebandt)

    do i = 1, ishift
      rho2ns(i, 1, 1) = dzero
      rho2ns(i, 1, 2) = dzero
      rhotborb(i) = dzero
    end do

    do i = 1, jws(it)
      ip = i + ishift
      rho2ns(ip, 1, 1) = rho2ns(ip, 1, 1) - 0.5d0*sqpi*rhochr(i, it)*(r(i,1)** &
        2)
      rho2ns(ip, 1, 2) = rho2ns(ip, 1, 2) - 0.5d0*sqpi*rhospn(i, it)*(r(i,1)** &
        2)
      rhotborb(ip) = rhotborb(ip) - 0.5d0*sqpi*rhoorb(i, it)*(r(i,1)**2)
    end do

    do il = 1, nl
      den(il-1, iecurr+ielast) = -0.5d0*(dosl0(il,it)+smtl0(il,it))*pi

      den(il-1, iecurr) = -0.5d0*(dosl0(il,it)-smtl0(il,it))*pi

      do i = 1, 2
        dmuorb(il-1, i) = -omtls0(il, it, i)*pi
      end do

      dmuorb(il-1, 3) = -omtl0(il, it)*pi
    end do
    den(nl, iecurr+ielast) = czero
    den(nl, iecurr) = czero

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      ITERMDIR

    if (itermvdir) then

      qmphiloc(iq) = qmphi
      qmtetloc(iq) = qmtet

      call cinit(nlmax*ntmax*3*nmvecmax, mvevdl0)
      call cinit(nlmax*ntmax*3*nmvecmax, bmvevdl0)
      call cinit(nlmax*ntmax*3*nmvecmax, mvevil1)
      call cinit(nlmax*ntmax*3*nmvecmax, bmvevil1)

      call calcmvec(nfilcbwf, splitss, iepath, nepath, irel, iprint, nt, nl, &
        mezz, mezj, taut, tsst, iqat, nkmq, nkm, iec, netab, igrid(iepath), &
        we, mvevdl0, mvevil1, bmvevdl0, bmvevil1, r2drdi, jws, imt, amemvec, &
        ikmllim1, ikmllim2, imkmtab, ntmax, nlmax, nmuemax, nqmax, nkmmax, &
        nmmax, nmvecmax, nrmax)

      do i = 1, nmvecmax
        do j = 1, 3
          do il = 1, nl
            mvevil(il-1, j, i) = mvevil(il-1, j, i) - &
              mvevdl0(il, it, j, i)*pi*we
          end do
        end do
      end do
    end if

!      ITERMDIR
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  end do ! IPOINT = 1,NQDOS

  if ((iecurr/=ielast) .or. (.not. ldorhoef)) return

! ======================================================================
!     get the charge at the Fermi energy (IELAST)
!     call SCFCHRDNS with the energy weight CONE --> not overwrite WE

  call rinit(nrmax, rhochr(1,it))
  call rinit(nrmax, rhospn(1,it))

  call scfchrdns(nfilcbwf, r2drdi, jws, imt, shftef, totdos, muespn, mueorb, &
    irel, iprint, nt, nl, nkm, eryd, cone, efermi, iec, netab, dos, smt, omt, &
    hff, dosi, smti, omti, hffi, dosm, dosl0, dosint, smtm, smtl0, smtint, &
    omtm, omtl0, omtint, hffm, hffl0, hffint, bcor, bcors, dzz, dzj, szz, szj, &
    ozz, ozj, bzz, bzj, ozzs, ozjs, omtls0, tautlin, nvaltot, txtt, conc, nat, &
    rhochr, rhospn, rhoorb, qel, gdia, gmdia, goff, ntmax, nlmax, nmuemax, &
    linmax, nrmax, nmmax, nkmmax, eband, ebandt)

  do i = 1, ishift
    r2nef(i, 1, 1) = dzero
    r2nef(i, 1, 2) = dzero
  end do

  do i = 1, jws(it)
    ip = i + ishift
    r2nef(ip, 1, 1) = -0.5d0*sqpi*rhochr(i, it)*(r(i,1)**2)
    r2nef(ip, 1, 2) = -0.5d0*sqpi*rhospn(i, it)*(r(i,1)**2)
  end do

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      ITERMDIR

  if (itermvdir) then
    do i = 1, nmvecmax
      do j = 1, 3
        do il = 1, nl
          mvevilef(il-1, j, i) = -mvevdl0(il, it, j, i)*pi
        end do
      end do
    end do
  end if

!      ITERMDIR
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


! ======================================================================

end subroutine
