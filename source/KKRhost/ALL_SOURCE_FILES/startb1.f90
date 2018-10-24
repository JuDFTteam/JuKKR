!------------------------------------------------------------------------------------
!> Summary: Reads the input potentials
!> Author: B. Drittler
!> Reads the input potentials with units given by:
!>
!> - Rydbergs - units for energy
!> - The lattice constant and all other lengths given in bohr units
!> - The planck constant \(\frac{h}{2\pi}=1\)
!> - The electron charge \(e=\sqrt{2}\)
!> - The electron mass \(m=\frac{1}{2}\)
!> - The speed of light \(c = \frac{2}{\alpha} = 274.0720442\) with the
!> fine structure constant \(\alpha\)
!>
!> In case of shape corrections this routine reads from unit 19 a suitable
!> radial mesh 'xrn',its derivate 'drn' and the shape
!> functions 'thetas'. Thus, the region from the muffin-tin to the
!> circumscribed sphere radii is divided into 'npan'
!> pannels, each one containing 'nm(ipan)' points in order to take care of
!> the discontinuities of the shape-function derivative.
!------------------------------------------------------------------------------------
!> @note 
!> Remember that the input potentials do not include the electro-static contribution 
!> of the nucleus of the cell itself this has to be added explicitly!
!>
!> Modified for bandstructure code
!> @endnote
!------------------------------------------------------------------------------------
module mod_startb1
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Reads the input potentials
  !> Author: B. Drittler
  !> Category: input-output, potential, shape-functions, KKRhost 
  !> Deprecated: False 
  !> Reads the input potentials with units given by:
  !> 
  !> - Rydbergs - units for energy
  !> - The lattice constant and all other lengths given in bohr units
  !> - The planck constant \(\frac{h}{2\pi}=1\)
  !> - The electron charge \(e=\sqrt{2}\)
  !> - The electron mass \(m=\frac{1}{2}\)
  !> - The speed of light \(c = \frac{2}{\alpha} = 274.0720442\) with the
  !> fine structure constant \(\alpha\)
  !>
  !> In case of shape corrections this routine reads from unit 19 a suitable
  !> radial mesh 'xrn',its derivate 'drn' and the shape
  !> functions 'thetas'. Thus, the region from the muffin-tin to the
  !> circumscribed sphere radii is divided into 'npan'
  !> pannels, each one containing 'nm(ipan)' points in order to take care of
  !> the discontinuities of the shape-function derivative.
  !-------------------------------------------------------------------------------
  !> @note 
  !> Remember that the input potentials do not include the electro-static contribution 
  !> of the nucleus of the cell itself this has to be added explicitly!
  !>
  !> Modified for bandstructure code
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine startb1(ifile,ipf,ipfe,ipe,krel,kws,lmax,nbeg,nend,alat,rmtnew,rmt,    &
    ititle,imt,irc,vconst,ins,irns,fpradius,nspin, vins,irmin,kshape,ntcell,ircut,  &
    ipan,thetas,ifunm,nfu,llmsp,lmsp,efermi,vbc,dror,rs,s,vm2z,rws,ecore,lcore,     &
    ncore,drdi,r,zat,a,b,irws,iinfo,lmpot,irmind,irm,lmxspd,ipand,irid,irnsd,natyp, &
    ncelld,nfund,nspotd,ivshift,npotd)

    use :: mod_constants
    use :: mod_potcut
    use :: mod_calrmt
    use :: mod_rinit

    implicit none
    real (kind=dp), parameter :: eps = 1.0e-12_dp
    ! ..
    ! .. Input variables
    integer, intent (in) :: irm    !! Maximum number of radial points
    integer, intent (in) :: kws    !! 0 (MT), 1(ASA)
    integer, intent (in) :: ins    !! 0 (MT), 1(ASA), 2(Full Potential)
    integer, intent (in) :: ipe    !! Not real used, IPFE should be 0
    integer, intent (in) :: ipf    !! Not real used, IPFE should be 0
    integer, intent (in) :: ipfe   !! Not real used, IPFE should be 0
    integer, intent (in) :: irid   !! Shape functions parameters in non-spherical part
    integer, intent (in) :: lmax   !! Maximum l component in wave function expansion
    integer, intent (in) :: nbeg   !! Starting number for reading the potential
    integer, intent (in) :: nend   !! Final number for reading the potential
    integer, intent (in) :: krel   !! Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
    integer, intent (in) :: npotd  !! (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
    integer, intent (in) :: nspin  !! Counter for spin directions
    integer, intent (in) :: ipand  !! Number of panels in non-spherical part
    integer, intent (in) :: irnsd
    integer, intent (in) :: nfund  !! Shape functions parameters in non-spherical part
    integer, intent (in) :: ifile  !! Unit specifier for potential card
    integer, intent (in) :: iinfo
    integer, intent (in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: lmpot  !! (LPOT+1)**2
    integer, intent (in) :: irmind !! IRM-IRNSD
    integer, intent (in) :: lmxspd !! (2*LPOT+1)**2
    integer, intent (in) :: ncelld !! Number of cells (shapes) in non-spherical part
    integer, intent (in) :: nspotd !! Number of potentials for storing non-sph. potentials
    integer, intent (in) :: kshape !! Exact treatment of WS cell
    integer, intent (in) :: ivshift
    real (kind=dp), intent (in) :: vconst !! Potential shift
    integer, dimension (natyp), intent (in) :: ntcell !! Index for WS cell
    real (kind=dp), dimension (natyp), intent (in) :: fpradius !! R point at which full-potential treatment starts
    ! .. In/Out variables
    real (kind=dp), intent (inout) :: alat !! Lattice constant in a.u.
    real (kind=dp), intent (inout) :: efermi !! Fermi energy
    integer, dimension (natyp), intent (inout) :: nfu !! number of shape function components in cell 'icell'
    integer, dimension (natyp), intent (inout) :: imt !! R point at MT radius
    integer, dimension (natyp), intent (inout) :: irc !! R point for potential cutting
    integer, dimension (natyp), intent (inout) :: ipan !! Number of panels in non-MT-region
    integer, dimension (natyp), intent (inout) :: irns !! Position of atoms in the unit cell in units of bravais vectors
    integer, dimension (natyp), intent (inout) :: irws !! R point at WS radius
    integer, dimension (natyp), intent (inout) :: irmin !! Max R for spherical treatment
    integer, dimension (npotd), intent (inout) :: ncore !! Number of core states
    integer, dimension (natyp, lmxspd), intent (inout) :: lmsp !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
    integer, dimension (20, npotd), intent (inout) :: lcore !! Angular momentum of core states
    integer, dimension (natyp, nfund), intent (inout) :: llmsp !! lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
    integer, dimension (0:ipand, natyp), intent (inout) :: ircut !! R points of panel borders
    integer, dimension (natyp, lmxspd), intent (inout) :: ifunm
    integer, dimension (20, npotd), intent (inout) :: ititle !! Titles of the potential card
    real (kind=dp), dimension (natyp), intent (inout) :: a !! Constants for exponential R mesh
    real (kind=dp), dimension (natyp), intent (inout) :: b !! Constants for exponential R mesh
    real (kind=dp), dimension (natyp), intent (inout) :: zat !! Nuclear charge
    real (kind=dp), dimension (2), intent (inout) :: vbc !! Potential constants
    real (kind=dp), dimension (natyp), intent (inout) :: rmt !! Muffin-tin radius of true system
    real (kind=dp), dimension (natyp), intent (inout) :: rws !! Wigner Seitz radius
    real (kind=dp), dimension (natyp), intent (inout) :: rmtnew !! Adapted muffin-tin radius
    real (kind=dp), dimension (0:lmax, natyp), intent (inout) :: s
    real (kind=dp), dimension (irm, natyp), intent (inout) :: r !! Radial mesh ( in units a Bohr)
    real (kind=dp), dimension (irm, natyp), intent (inout) :: drdi !! Derivative dr/di
    real (kind=dp), dimension (irm, natyp), intent (inout) :: dror
    real (kind=dp), dimension (irm, npotd), intent (inout) :: vm2z
    real (kind=dp), dimension (20, npotd), intent (inout) :: ecore !! Core energies
    real (kind=dp), dimension (irm, 0:lmax, natyp), intent (inout) :: rs
    real (kind=dp), dimension (irmind:irm, lmpot, nspotd), intent (inout) :: vins !! Non-spherical part of the potential
    real (kind=dp), dimension (irid, nfund, ncelld), intent (inout) :: thetas !! shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
    ! .. Local Scalars
    integer :: inslpd, lmshapemax
    integer :: j, l, lm, lm1, lmpotp, n, ncell, nfun, nr
    integer :: irminm, irminp, irns1p, irt1p, irws1, isave, ispin, isum
    integer :: i, ia, icell, icore, ifun, ih, imt1, inew, io, ipan1, ir, irc1, iri
    real (kind=dp) :: a1, b1, ea, efnew, s1, z1
    logical :: test
    ! ..
    ! .. Local Arrays
    integer, dimension (ncelld) :: npan
    integer, dimension (ncelld) :: meshn
    integer, dimension (ipand, ncelld) :: nm
    real (kind=dp), dimension (irm) :: u
    real (kind=dp), dimension (ncelld) :: scale
    real (kind=dp), dimension (irid) :: rdummy
    real (kind=dp), dimension (irm) :: vspsme ! dummy for potcut IMPURITY-compatible
    real (kind=dp), dimension (irid, ncelld) :: drn
    real (kind=dp), dimension (irid, ncelld) :: xrn
    ! ..
    ! .. Data statement ..
    integer :: ishape
    data ishape/0/
    ! ----------------------------------------------------------------------------
    ! Output of radial mesh information
    ! ----------------------------------------------------------------------------
    io = 0
    if (iinfo/=0 .and. test('RMESH   ')) io = 1
    ! ----------------------------------------------------------------------------
    ! Set speed of light
    ! ----------------------------------------------------------------------------
    inslpd = (irnsd+1)*lmpot*nspotd
    lmshapemax = (4*lmax+1)**2
    call rinit(inslpd, vins(irmind,1,1))
    ! ----------------------------------------------------------------------------
    ! Read radial mesh information of the shape functions and
    ! shape functions THETAS in the first iteration - if needed
    ! ----------------------------------------------------------------------------
    if ((kshape/=0) .and. (ishape==0)) then
      ishape = 1
      read (19, fmt=100) ncell
      write (1337, fmt=*) '  ncell : ', ncell, ncelld
      ! check consistency with shape numbers from inputcard
      if (maxval(ntcell(1:natyp))>ncell) then
        write (*, *) 'Found ', ncell, 'shapes in shapefun file but need', maxval(ntcell(1:natyp)), 'according to inputcard/default values'
        write (*, *) 'Did you set <SHAPE> correctly in inputcard?'
        stop 'Error consistency shapes from input/shapefun file'
      end if

      if (ncell>ncelld) then
        write (6, *) 'Please, change the parameter ncelld (', ncelld, ') in the inputcard to', ncell
        stop 'STARTB - NCELLD'
      end if

      read (19, fmt=110)(scale(icell), icell=1, ncell)
      do icell = 1, ncell
        read (19, fmt=100) npan(icell), meshn(icell)

        if (npan(icell)+1>ipand) then
          write (6, *) 'Please, change the parameter ipand (', ipand, ') in the inputcard to', npan(icell) + 1
          stop 'STARTB - IPAND'
        end if

        if (meshn(icell)>irid) then
          write (6, *) 'Please, change the parameter irid (', irid, ') in the inputcard to', meshn(icell)
          stop 'STARTB - IRID'
        end if

        read (19, fmt=100)(nm(ipan1,icell), ipan1=2, npan(icell)+1)
        read (19, fmt=110)(xrn(ir,icell), drn(ir,icell), ir=1, meshn(icell))

        read (19, fmt=100) nfu(icell)
        nfun = nfu(icell)
        write (1337, fmt=*) '  nfun  : ', nfun, nfund

        if (nfun>nfund) then
          write (6, *) 'Please, change the parameter nfund (', nfund, ') in the inputcard to', nfun
          stop 'STARTB - NFUND'
        end if

        do lm = 1, lmxspd
          lmsp(icell, lm) = 0
        end do

        do ifun = 1, nfun
          read (19, fmt=100) lm
          if (lm<=lmshapemax) then
            llmsp(icell, ifun) = lm
            lmsp(icell, lm) = 1
            ifunm(icell, lm) = ifun
            read (19, fmt=110)(thetas(n,ifun,icell), n=1, meshn(icell))
          else
            read (19, fmt=110)(rdummy(n), n=1, meshn(icell))
          end if
        end do

      end do
    end if                         ! ((KSHAPE.NE.0) .AND. (IFILE.NE.0))
    ! ----------------------------------------------------------------------------
    ! LMPOT = (LPOT+1)* (LPOT+1)
    do ih = nbeg, nend
      do ispin = 1, nspin
        i = nspin*(ih-1) + ispin

        if (ifile/=0) then
          ircut(0, ih) = 0
          if (ins/=0) then
            ! p.z.            IF (KSHAPE.NE.0) THEN
            icell = ntcell(ih)
            ipan(ih) = 1 + npan(icell)
          else
            ipan(ih) = 1
          end if
          ! -------------------------------------------------------------------
          ! Read title of potential card
          ! -------------------------------------------------------------------
          read (ifile, fmt=120)(ititle(ia,i), ia=1, 20)
          if (iinfo/=0) then
            if (ins==0) then
              write (1337, fmt=180)(ititle(ia,i), ia=1, 20)
            else
              write (1337, fmt=190)(ititle(ia,i), ia=1, 20)
            end if
          end if

          ! -------------------------------------------------------------------
          ! Read muffin-tin radius , lattice constant and new muffin radius
          ! (new mt radius is adapted to the given radial mesh)
          ! -------------------------------------------------------------------
          read (ifile, fmt=*) rmt(ih), alat, rmtnew(ih)
          ! READ (IFILE,FMT=9030) RMT(IH),ALAT,RMTNEW(IH)

          ! -------------------------------------------------------------------
          ! Read nuclear charge , lmax of the core states ,
          ! wigner seitz radius , fermi energy and energy difference
          ! between electrostatic zero and muffin tin zero
          ! -------------------------------------------------------------------
          ! READ (IFILE,FMT=9040) ZAT(IH),RWS(IH),EFNEW,VBC(ISPIN)
          read (ifile, *) z1
          read (ifile, *) rws(ih), efnew, vbc(ispin)

          ! READ (IFILE,*) Z1,RWS(IH),EFNEW,VBC(ISPIN)
          if (zat(ih)<0.e0_dp) zat(ih) = z1
          if (abs(z1-zat(ih))>eps .and. abs(zat(ih))>=0.e0_dp) then
            write (*, *) 'Warning: For atom ', ih, ': ZATOM different in inputcard and in potential.', zat(ih), z1
          end if
          ! -------------------------------------------------------------------
          ! If efermi .eq. 0 use value from in5
          ! -------------------------------------------------------------------
          if (abs(efnew)>eps .and. i==1) efermi = efnew
          ! -------------------------------------------------------------------
          ! Read : number of radial mesh points
          ! (in case of ws input-potential: last mesh point corresponds
          ! to ws-radius, in case of shape-corrected input-potential
          ! last mesh point of the exponential mesh corresponds to
          ! mt-radius/nevertheless this point is always in the array
          ! irws(ih)),number of points for the radial non-muffin-tin
          ! mesh  needed for shape functions, the constants a and b
          ! for the radial exponential mesh : r(i) = b*(exp(a*(i-1))-1)
          ! the no. of different core states and some other stuff
          ! -------------------------------------------------------------------
          ! READ (IFILE,FMT=9050) IRWS(IH),A(IH),B(IH),NCORE(I),INEW
          read (ifile, fmt=*) irws(ih)
          read (ifile, fmt=*) a(ih), b(ih)
          read (ifile, fmt=*) ncore(i), inew

          nr = irws(ih)

          if (nr>irm) then
            write (6, *) 'Increase parameter IRM in the inputcard ', ' to a value .ge. ', nr, ' (= IRWS(', ih, ')).'
            stop 'STARTB1 - IRWS'
          end if
          ! -------------------------------------------------------------------
          ! Read the different core states : l and energy
          ! -------------------------------------------------------------------
          if (ncore(i)>=1) then
            do icore = 1, ncore(i)
              read (ifile, fmt=170) lcore(icore, i), ecore(icore, i)
            end do
          end if

          if (ins<1) then
            ! ----------------------------------------------------------------
            ! Read radial mesh points, its derivative, the spherically averaged
            ! charge density and the input potential without the nuclear pot.
            ! ----------------------------------------------------------------
            if (inew==0) then
              read (ifile, fmt=160)(r(ir,ih), drdi(ir,ih), vm2z(ir,i), ir=1, nr)
            else
              read (ifile, fmt=*)(vm2z(ir,i), ir=1, nr)
            end if
          else                     ! (INS.LT.1)
            ! -------------------------------------------------------------------
            ! Read full potential - the non spherical contribution from irmin
            ! to irt - remember that the lm = 1 contribution is multiplied by
            ! 1/sqrt(4 pi)
            ! -------------------------------------------------------------------
            read (ifile, fmt=200) irt1p, irns1p, lmpotp, isave
            irminp = irt1p - irns1p
            irminm = max(irminp, irmind)
            read (ifile, fmt=210)(vm2z(ir,i), ir=1, nr)
            if (lmpotp>1) then
              lm1 = 2
              do lm = 2, lmpotp
                if (lm1/=1) then
                  if (isave==1) then
                    read (ifile, fmt=200) lm1
                  else
                    lm1 = lm
                  end if

                  if (lm1>1) then
                    read (ifile, fmt=210)(u(ir), ir=irminp, nr)
                    if (lm1<=lmpot) then
                      do ir = irminm, nr
                        vins(ir, lm1, i) = u(ir)
                      end do
                    end if
                  end if
                end if
              end do
            end if
          end if                   ! (INS.LT.1)

          irws1 = irws(ih)
          ! ----------------------------------------------------------------------
          ! Redefine new mt-radius in case of shape corrections
          ! ----------------------------------------------------------------------
          if (ins/=0) then
            ! p.z.      IF (KSHAPE.NE.0) THEN
            rmtnew(ih) = scale(icell)*alat*xrn(1, icell)
            imt1 = nint(log(rmtnew(ih)/b(ih)+1.0e0_dp)/a(ih)) + 1
            ! -------------------------------------------------------------------
            ! For proper core treatment imt must be odd
            ! shift potential by one mesh point if imt is even
            ! -------------------------------------------------------------------
            if (mod(imt1,2)==0) then
              imt1 = imt1 + 1
              do ir = imt1, 2, -1
                vm2z(ir, i) = vm2z(ir-1, i)
              end do
            end if

            imt(ih) = imt1
            b(ih) = rmtnew(ih)/(exp(a(ih)*real(imt1-1,kind=dp))-1.0e0_dp)
          end if                   ! (KSHAPE.NE.0)
          ! ----------------------------------------------------------------------
          ! Generate radial mesh - potential only is stored in potential card
          ! INEW = 1
          ! p. zahn, jan. 99
          ! ----------------------------------------------------------------------
          a1 = a(ih)
          b1 = b(ih)
          r(1, ih) = 0.0e0_dp
          drdi(1, ih) = a1*b1
          do ir = 2, irws1
            ea = exp(a1*real(ir-1,kind=dp))
            r(ir, ih) = b1*(ea-1.0e0_dp)
            drdi(ir, ih) = a1*b1*ea
            dror(ir, ih) = a1/(1.0e0_dp-1.0e0_dp/ea)
          end do
          ! ----------------------------------------------------------------------
          ! Fill cell-type depending mesh points in the non-muffin-tin-region
          ! ----------------------------------------------------------------------
          if (ins/=0) then
            ! p.z.      IF (KSHAPE.NE.0) THEN
            do iri = 1, meshn(icell)
              ir = iri + imt1
              r(ir, ih) = scale(icell)*alat*xrn(iri, icell)
              drdi(ir, ih) = scale(icell)*alat*drn(iri, icell)
              dror(ir, ih) = drdi(ir, ih)/r(ir, ih)
            end do
          end if

          rws(ih) = r(irws1, ih)
          ! ----------------------------------------------------------------------
          ! Kshape.eq.0 : calculate new rmt adapted to exp. mesh
          ! ----------------------------------------------------------------------
          call calrmt(ipf, ipfe, ipe, imt(ih), zat(ih), rmt(ih), rws(ih), rmtnew(ih), alat, drdi(1,ih), a(ih), b(ih), irws1, r(1,ih), io, ins)
          ! p.z. +                  R(1,IH),IO,KSHAPE)

          if (ins>0) then
            ! p.z.            IF (KSHAPE.GT.0) THEN
            ircut(1, ih) = imt(ih)
            isum = imt(ih)
            do ipan1 = 2, ipan(ih)
              isum = isum + nm(ipan1, icell)
              ircut(ipan1, ih) = isum
            end do
            nr = isum
            if (irt1p/=nr) then
              write (*, *) 'STARTB1: Error: IRT1P.NE.NR', irt1p, nr, ' for atom', ih
              stop 'STARTB1: IRT1P.NE.NR'
            end if
          else                     ! (KSHAPE.GT.0)
            nr = irws(ih)
            if (kws>=1) then
              ircut(1, ih) = irws1
            else
              ircut(1, ih) = imt(ih)
            end if
          end if                   ! (KSHAPE.GT.0)

          irc(ih) = ircut(ipan(ih), ih)
          ! ----------------------------------------------------------------------
          ! Fill array irmin in case of full potential
          ! ----------------------------------------------------------------------
          if (ins/=0) then
            if (fpradius(ih)>=0.e0_dp) then
              irmin(ih) = min(floor(log(fpradius(ih)/b(ih)+1.0e0_dp)/a(ih))+1, imt(ih))
              irns(ih) = nr - irmin(ih)
            else if (irns(ih)>=meshn(icell)) then
              irmin(ih) = nr - irns(ih)
            else
              irns(ih) = irns1p
              irmin(ih) = nr - irns(ih)
            end if
          end if
          ! ----------------------------------------------------------------------
          ! Generate arrays for the calculation of the wave functions
          ! ----------------------------------------------------------------------
          z1 = zat(ih)
          do l = 0, lmax
            if (krel>=1) then
              s1 = sqrt(real(l*l+l+1,kind=dp)-4.0e0_dp*z1*z1/(cvlight*cvlight))
              if (abs(z1)<eps) s1 = real(l, kind=dp)
            else
              s1 = real(l, kind=dp)
            end if

            s(l, ih) = s1
            rs(1, l, ih) = 0.0e0_dp
            do ir = 2, nr
              rs(ir, l, ih) = r(ir, ih)**s1
            end do
          end do                   ! L = 0,LMAX
          ! -------------------------------------------------------------------
          ! Cut input potential at rmt if given only at exponential mesh
          ! -------------------------------------------------------------------
          if (kshape==1) then
            imt1 = imt(ih)
            irc1 = ircut(ipan(ih), ih)
            call potcut(imt1, irc1, ins, lmpot, r(1,ih), vm2z(1,i), vspsme, vins(irmind,1,i), zat(ih), irm, irmind)
          end if
          ! -------------------------------------------------------------------
          ! First iteration : shift all potentials (only for test purpose)
          ! in case of test option 'atptshft' shift only potential of atom at
          ! position ivshift
          ! -------------------------------------------------------------------
          if (test('atptshft') .and. (ih==ivshift)) then
            write (1337, *) 'atptshft', ih, ivshift, vconst, nr, irmin(ih)
            do j = 1, irmin(ih)
              vm2z(j, i) = vm2z(j, i) + vconst
            end do
          else if (.not. test('atptshft') .and. abs(vconst)>eps) then
            write (1337, *) 'shifting potential by VCONST=', vconst
            do j = 1, nr
              vm2z(j, i) = vm2z(j, i) + vconst
            end do
          end if
        end if                     ! (ifile.ne.0)

        if (kshape==0 .and. kws==0) then
          ! -------------------------------------------------------------------
          ! In case of a mt calculation cut potential at mt radius
          ! -------------------------------------------------------------------
          imt1 = imt(ih)
          irws1 = irws(ih)
          call potcut(imt1, irws1, ins, lmpot, r(1,ih), vm2z(1,i), vspsme, vins(irmind,1,i), zat(ih), irm, irmind)
        end if                     ! KSHAPE.EQ.0 .AND. KWS.EQ.0
      end do                       ! ISPIN = 1,NSPIN
    end do                         ! IH = NBEG,NEND

    if (ins/=0) then
      i = 0
      do ih = nbeg, nend
        if (irmin(ih)<irmind) then
          write (*, *) 'IRMIN < IRMIND for atom', ih
          write (*, *) irmin(ih), irmind
          write (*, *) 'Increase dimension IRNSD'
          i = 1
        end if
      end do
      if (i/=0) stop 'stop startb1 IRNS IRNSD'
    end if

    return

100 format (16i5)
110 format (4d20.12)
120 format (20a4)
130 format (3f12.8)
140 format (f10.5, /, f10.5, 2f15.10)
150 format (i3, /, 2d15.8, /, 2i2)
160 format (1p, 2d15.6, 1p, d15.8)
170 format (i5, 1p, d20.11)
    ! 9080 format (10x,20a4)
180 format (' < ', 20a4)
190 format (' <#', 20a4)
200 format (10i5)
210 format (1p, 4d20.13)
  end subroutine startb1

end module mod_startb1
