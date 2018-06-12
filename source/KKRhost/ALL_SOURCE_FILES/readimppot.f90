! -------------------------------------------------------------------------------
! SUBROUTINE: READIMPPOT
! > @brief Reads the potential and shapefun of inpurity
! -------------------------------------------------------------------------------
subroutine readimppot(natomimp, ins, ipf, ipfe, ipe, kws, nspin, lpot, &
  ipanimp, thetasimp, ircutimp, irwsimp, khfeld, hfield, vinsimp, vm2zimp, &
  irminimp, rimp, zimp, irmd, irnsd, irid, nfund, ntotd, ipand)
  ! ************************************************************************
  ! read in impurity potential
  ! n.h.long, May 2013
  ! -----------------------------------------------------------------------
  use :: mod_datatypes, only: dp
  implicit none
  ! .. Parameters ..
  integer :: nspin, natomimp, irmd, irnsd, irid, nfund, ntotd, ipand
  ! ..
  ! .. Scalar Arguments ..
  real (kind=dp) :: alat, hfield, vbc(2)
  integer :: ins, ipe, ipf, ipfe, khfeld, kws, lpot
  ! ..
  ! .. Array Arguments ..
  real (kind=dp) :: a(natomimp), b(natomimp), drdi(irmd, natomimp), &
    dror(irmd, natomimp), ecore(20, nspin*natomimp), rimp(irmd, natomimp), &
    rmt(natomimp), rmtnew(natomimp), rws(natomimp), thetasimp(irid, nfund, &
    natomimp), vinsimp((irmd-irnsd):irmd, (lpot+1)**2, natomimp*nspin), &
    vm2zimp(irmd, natomimp*nspin), zimp(natomimp)
  integer :: imt(natomimp), ipanimp(natomimp), ircutimp(0:ipand, natomimp), &
    irminimp(natomimp), irwsimp(natomimp), ititle(20, nspin*natomimp), &
    lcore(20, nspin*natomimp), ncore(nspin*natomimp), nfu(natomimp)
  ! ..
  ! .. Local Arrays ..
  real (kind=dp) :: dummy2(irmd, natomimp*nspin)
  ! ..
  ! .. Local Scalars ..
  real (kind=dp) :: a1, b1, ea, efnew, dummy
  integer :: i, ia, icell, icore, ifun, ih, imt1, inew, io, ipan1, ir, irc1, &
    iri, irminm, irminp, irns1p, irt1p, irws1, isave, ispin, isum, j, lm, lm1, &
    lmpot, lmpotp, n, ncell, nfun, nr
  logical :: test
  ! ..
  ! .. Local Arrays ..
  real (kind=dp) :: drn(irid, natomimp), scale(1), u(irmd), &
    xrn(irid, natomimp)
  integer :: meshn(natomimp), nm(ipand, natomimp), npan(natomimp)
  ! ..
  ! .. External Subroutines ..
  external :: calrmt, potcut, rinit, test
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic :: anint, exp, log, max, mod, real, sqrt
  ! ..
  ! ------------------------------------------------------------------
  write (1337, *) 'in readimppot'
  vinsimp = 0e0_dp
  ! ------------------------------------------------------------------
  ! read data from shapefun_imp file
  if (ins>0) then
    open (unit=20, file='shapefun_imp', form='FORMATTED')
    read (20, *) ncell
    read (20, *) scale(1)
    do icell = 1, ncell
      read (20, fmt=100) npan(icell), meshn(icell)
      read (20, fmt=100)(nm(ipan1,icell), ipan1=2, npan(icell)+1)
      read (20, fmt=110)(xrn(ir,icell), drn(ir,icell), ir=1, meshn(icell))
      read (20, fmt=100) nfu(icell)
      nfun = nfu(icell)

      do ifun = 1, nfun
        read (20, fmt=100) lm
        if (lm<=(2*lpot+1)**2) then
          read (20, fmt=110)(thetasimp(n,ifun,icell), n=1, meshn(icell))
        else
          read (20, fmt=110)(dummy, n=1, meshn(icell))
        end if
      end do

    end do
  end if                           ! INS.EQ.1

  do icell = 1, ncell
    if (ins/=0) then
      ipanimp(icell) = 1 + npan(icell)
    else
      ipanimp(icell) = 1
    end if
  end do
  ! ------------------------------------------------------------------
  ! read in impurity potential

  open (unit=21, file='potential_imp', form='FORMATTED')
  lmpot = (lpot+1)*(lpot+1)
  do ih = 1, ncell
    do ispin = 1, nspin
      i = nspin*(ih-1) + ispin
      ircutimp(0, ih) = 0

      ! ---> read title of potential card
      read (21, fmt=120)(ititle(ia,i), ia=1, 20)

      ! --->read muffin-tin radius , lattice constant and new muffin radius
      ! READ (21,FMT=9030) RMT(IH),ALAT,RMTNEW(IH)
      read (21, fmt=*) rmt(ih), alat, rmtnew(ih)

      ! ---> read nuclear charge , lmax of the core states ,
      ! wigner seitz radius , fermi energy and energy difference
      ! between electrostatic zero and muffin tin zero

      ! READ (21,FMT=9040) ZIMP(IH),RWS(IH),EFNEW,VBC(ISPIN)
      read (21, fmt=*) zimp(ih)
      read (21, fmt=*) rws(ih), efnew, vbc(ispin)

      ! ---> read : number of radial mesh points
      ! (in case of ws input-potential: last mesh point corresponds
      ! to ws-radius, in case of shape-corrected input-potential
      ! last mesh point of the exponential mesh corresponds to
      ! mt-radius/nevertheless this point is always in the array
      ! irws(ih)),number of points for the radial non-muffin-tin
      ! mesh  needed for shape functions, the constants a and b
      ! for the radial exponential mesh : r(i) = b*(exp(a*(i-1))-1)
      ! the no. of different core states and some other stuff

      read (21, fmt=150) irwsimp(ih)
      ! READ (21,FMT=9051) A(IH),B(IH),NCORE(I),INEW
      read (21, fmt=*) a(ih), b(ih)
      read (21, fmt=*) ncore(i), inew
      nr = irwsimp(ih)
      ! ---> read the different core states : l and energy

      if (ncore(i)>=1) read (21, fmt=180)(lcore(icore,i), ecore(icore,i), &
        icore=1, ncore(i))

      if (ins<1) then

        ! --->  read radial mesh points, its derivative, the spherically
        ! averaged
        ! charge density and the input potential without the nuclear pot.

        if (inew==0) then
          read (21, fmt=170)(rimp(ir,ih), drdi(ir,ih), vm2zimp(ir,i), ir=1, &
            nr)
        else
          read (21, fmt=*)(vm2zimp(ir,i), ir=1, nr)
        end if

      else                         ! (INS.LT.1)

        ! --->  read full potential - the non spherical contribution from
        ! irmin
        ! to irt - remember that the lm = 1 contribution is multiplied by
        ! 1/sqrt(4 pi)

        read (21, fmt=190) irt1p, irns1p, lmpotp, isave
        irminp = irt1p - irns1p
        irminm = max(irminp, irmd-irnsd)
        read (21, fmt=200)(vm2zimp(ir,i), ir=1, nr)
        if (lmpotp>1) then
          lm1 = 2
          do lm = 2, lmpotp
            if (lm1/=1) then
              if (isave==1) then
                read (21, fmt=190) lm1
              else
                lm1 = lm
              end if
              if (lm1>1) then
                read (21, fmt=200)(u(ir), ir=irminp, nr)
                if (lm1<=lmpot) then
                  do ir = irminm, nr
                    vinsimp(ir, lm1, i) = u(ir)
                  end do
                end if
              end if
            end if
          end do
        end if
      end if                       ! (INS.LT.1)
      irws1 = irwsimp(ih)

      ! ---> redefine new mt-radius in case of shape corrections

      if (ins/=0) then
        rmtnew(ih) = scale(1)*alat*xrn(1, ih)
        imt1 = anint(log(rmtnew(ih)/b(ih)+1.0e0_dp)/a(ih)) + 1

        ! ---> for proper core treatment imt must be odd
        ! shift potential by one mesh point if imt is even

        if (mod(imt1,2)==0) then
          imt1 = imt1 + 1
          do ir = imt1, 2, -1
            vm2zimp(ir, i) = vm2zimp(ir-1, i)
          end do
        end if

        imt(ih) = imt1
        b(ih) = rmtnew(ih)/(exp(a(ih)*real(imt1-1,kind=dp))-1.0e0_dp)
      end if                       ! (INS.NE.0)

      ! ---> generate radial mesh - potential only is stored in potential card
      ! INEW = 1
      ! p. zahn, jan. 99

      a1 = a(ih)
      b1 = b(ih)
      rimp(1, ih) = 0.0e0_dp
      drdi(1, ih) = a1*b1
      do ir = 2, irws1
        ea = exp(a1*real(ir-1,kind=dp))
        rimp(ir, ih) = b1*(ea-1.0e0_dp)
        drdi(ir, ih) = a1*b1*ea
        dror(ir, ih) = a1/(1.0e0_dp-1.0e0_dp/ea)
      end do

      ! ---> fill cell-type depending mesh points in the non-muffin-tin-region

      if (ins/=0) then
        do iri = 1, meshn(ih)
          ir = iri + imt1
          rimp(ir, ih) = scale(1)*alat*xrn(iri, ih)
          drdi(ir, ih) = scale(1)*alat*drn(iri, ih)
          dror(ir, ih) = drdi(ir, ih)/rimp(ir, ih)
        end do
      end if

      rws(ih) = rimp(irws1, ih)

      ! ---> kshape.eq.0 : calculate new rmt adapted to exp. mesh

      call calrmt(ipf, ipfe, ipe, imt(ih), zimp(ih), rmt(ih), rws(ih), &
        rmtnew(ih), alat, drdi(1,ih), a(ih), b(ih), irws1, rimp(1,ih), io, &
        ins)

      if (ins>0) then
        ircutimp(1, ih) = imt(ih)
        isum = imt(ih)
        do ipan1 = 2, ipanimp(ih)
          isum = isum + nm(ipan1, ih)
          ircutimp(ipan1, ih) = isum
        end do
        nr = isum
      else                         ! INS.EQ.0
        nr = irwsimp(ih)
        if (kws>=1) then
          ircutimp(1, ih) = irws1
        else
          ircutimp(1, ih) = imt(ih)
        end if
      end if                       ! INS.GT.0

      ! ---> fill array irmin in case of full potential
      if (ins/=0) irminimp(ih) = nr - irns1p

      ! ---> cut input potential at rmt if given only at exponential mesh
      if (ins>1) then
        imt1 = imt(ih)
        irc1 = ircutimp(ipanimp(ih), ih)
        call potcut(imt1, irc1, ins, lmpot, rimp(1,ih), vm2zimp(1,i), dummy2, &
          vinsimp(irmd-irnsd,1,i), zimp(ih), irmd, irmd-irnsd)
      end if

      if (ins==0 .and. kws==0) then
        ! ---> in case of a mt calculation cut potential at mt radius
        imt1 = imt(ih)
        irws1 = irwsimp(ih)
        call potcut(imt1, irws1, ins, lmpot, rimp(1,ih), vm2zimp(1,i), dummy2, &
          vinsimp(irmd-irnsd,1,i), zimp(ih), irmd, irmd-irnsd)

      end if                       ! INS.EQ.0 .AND. KWS.EQ.0
      ! --->       maybe apply a magnetic field
      if (khfeld==1) then
        write (1337, *) 'ATOM', ih, 'SPIN', ispin, 'SHIFTED BY', &
          -real(2*ispin-3, kind=dp)*hfield
        do j = 1, ircutimp(ipanimp(ih), ih)
          vm2zimp(j, i) = vm2zimp(j, i) - real(2*ispin-3, kind=dp)*hfield
        end do
      end if

    end do                         ! ISPIN = 1,NSPIN
  end do                           ! IH = 1,NCELL
  close (20)
  close (21)

  return


100 format (16i5)
110 format (4d20.12)
120 format (20a4)
130 format (3f12.8)
140 format (f10.5, /, f10.5, 2f15.10)
150 format (i4)
160 format (2d15.8, /, 2i2)
170 format (1p, 2d15.6, 1p, d15.8)
180 format (i5, 1p, d20.11)
190 format (10i5)
200 format (1p, 4d20.13)
end subroutine readimppot
