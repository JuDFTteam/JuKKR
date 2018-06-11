!-------------------------------------------------------------------------------
! SUBROUTINE: READIMPPOT
!> @brief Reads the potential and shapefun of inpurity
!-------------------------------------------------------------------------------
    Subroutine readimppot(natomimp, ins, ipf, ipfe, ipe, kws, nspin, lpot, &
      ipanimp, thetasimp, ircutimp, irwsimp, khfeld, hfield, vinsimp, vm2zimp, &
      irminimp, rimp, zimp, irmd, irnsd, irid, nfund, ntotd, ipand)
! ************************************************************************
! read in impurity potential
! n.h.long, May 2013
!-----------------------------------------------------------------------
      Use mod_datatypes, Only: dp
      Implicit None
!.. Parameters ..
      Integer :: nspin, natomimp, irmd, irnsd, irid, nfund, ntotd, ipand
!..
!.. Scalar Arguments ..
      Real (Kind=dp) :: alat, hfield, vbc(2)
      Integer :: ins, ipe, ipf, ipfe, khfeld, kws, lpot
!..
!.. Array Arguments ..
      Real (Kind=dp) :: a(natomimp), b(natomimp), drdi(irmd, natomimp), &
        dror(irmd, natomimp), ecore(20, nspin*natomimp), rimp(irmd, natomimp), &
        rmt(natomimp), rmtnew(natomimp), rws(natomimp), &
        thetasimp(irid, nfund, natomimp), vinsimp((irmd-irnsd):irmd, (lpot+1) &
        **2, natomimp*nspin), vm2zimp(irmd, natomimp*nspin), zimp(natomimp)
      Integer :: imt(natomimp), ipanimp(natomimp), ircutimp(0:ipand, natomimp) &
        , irminimp(natomimp), irwsimp(natomimp), ititle(20, nspin*natomimp), &
        lcore(20, nspin*natomimp), ncore(nspin*natomimp), nfu(natomimp)
!..
!.. Local Arrays ..
      Real (Kind=dp) :: dummy2(irmd, natomimp*nspin)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: a1, b1, ea, efnew, dummy
      Integer :: i, ia, icell, icore, ifun, ih, imt1, inew, io, ipan1, ir, &
        irc1, iri, irminm, irminp, irns1p, irt1p, irws1, isave, ispin, isum, &
        j, lm, lm1, lmpot, lmpotp, n, ncell, nfun, nr
      Logical :: test
!..
!.. Local Arrays ..
      Real (Kind=dp) :: drn(irid, natomimp), scale(1), u(irmd), &
        xrn(irid, natomimp)
      Integer :: meshn(natomimp), nm(ipand, natomimp), npan(natomimp)
!..
!.. External Subroutines ..
      External :: calrmt, potcut, rinit, test
!..
!.. Intrinsic Functions ..
      Intrinsic :: anint, exp, log, max, mod, real, sqrt
!..
!------------------------------------------------------------------
      Write (1337, *) 'in readimppot'
      vinsimp = 0E0_dp
!------------------------------------------------------------------
!read data from shapefun_imp file
      If (ins>0) Then
        Open (Unit=20, File='shapefun_imp', Form='FORMATTED')
        Read (20, *) ncell
        Read (20, *) scale(1)
        Do icell = 1, ncell
          Read (20, Fmt=100) npan(icell), meshn(icell)
          Read (20, Fmt=100)(nm(ipan1,icell), ipan1=2, npan(icell)+1)
          Read (20, Fmt=110)(xrn(ir,icell), drn(ir,icell), ir=1, meshn(icell))
          Read (20, Fmt=100) nfu(icell)
          nfun = nfu(icell)

          Do ifun = 1, nfun
            Read (20, Fmt=100) lm
            If (lm<=(2*lpot+1)**2) Then
              Read (20, Fmt=110)(thetasimp(n,ifun,icell), n=1, meshn(icell))
            Else
              Read (20, Fmt=110)(dummy, n=1, meshn(icell))
            End If
          End Do

        End Do
      End If ! INS.EQ.1

      Do icell = 1, ncell
        If (ins/=0) Then
          ipanimp(icell) = 1 + npan(icell)
        Else
          ipanimp(icell) = 1
        End If
      End Do
!------------------------------------------------------------------
!read in impurity potential

      Open (Unit=21, File='potential_imp', Form='FORMATTED')
      lmpot = (lpot+1)*(lpot+1)
      Do ih = 1, ncell
        Do ispin = 1, nspin
          i = nspin*(ih-1) + ispin
          ircutimp(0, ih) = 0

!---> read title of potential card
          Read (21, Fmt=120)(ititle(ia,i), ia=1, 20)

!--->read muffin-tin radius , lattice constant and new muffin radius
!READ (21,FMT=9030) RMT(IH),ALAT,RMTNEW(IH)
          Read (21, Fmt=*) rmt(ih), alat, rmtnew(ih)

!---> read nuclear charge , lmax of the core states ,
!wigner seitz radius , fermi energy and energy difference
!between electrostatic zero and muffin tin zero

!READ (21,FMT=9040) ZIMP(IH),RWS(IH),EFNEW,VBC(ISPIN)
          Read (21, Fmt=*) zimp(ih)
          Read (21, Fmt=*) rws(ih), efnew, vbc(ispin)

!---> read : number of radial mesh points
!    (in case of ws input-potential: last mesh point corresponds
!    to ws-radius, in case of shape-corrected input-potential
!    last mesh point of the exponential mesh corresponds to
!    mt-radius/nevertheless this point is always in the array
!    irws(ih)),number of points for the radial non-muffin-tin
!    mesh  needed for shape functions, the constants a and b
!    for the radial exponential mesh : r(i) = b*(exp(a*(i-1))-1)
!    the no. of different core states and some other stuff

          Read (21, Fmt=150) irwsimp(ih)
!READ (21,FMT=9051) A(IH),B(IH),NCORE(I),INEW
          Read (21, Fmt=*) a(ih), b(ih)
          Read (21, Fmt=*) ncore(i), inew
          nr = irwsimp(ih)
!---> read the different core states : l and energy

          If (ncore(i)>=1) Read (21, Fmt=180)(lcore(icore,i), ecore(icore,i), &
            icore=1, ncore(i))

          If (ins<1) Then

!--->  read radial mesh points, its derivative, the spherically averaged
!      charge density and the input potential without the nuclear pot.

            If (inew==0) Then
              Read (21, Fmt=170)(rimp(ir,ih), drdi(ir,ih), vm2zimp(ir,i), &
                ir=1, nr)
            Else
              Read (21, Fmt=*)(vm2zimp(ir,i), ir=1, nr)
            End If

          Else ! (INS.LT.1)

!--->  read full potential - the non spherical contribution from irmin
!      to irt - remember that the lm = 1 contribution is multiplied by
!      1/sqrt(4 pi)

            Read (21, Fmt=190) irt1p, irns1p, lmpotp, isave
            irminp = irt1p - irns1p
            irminm = max(irminp, irmd-irnsd)
            Read (21, Fmt=200)(vm2zimp(ir,i), ir=1, nr)
            If (lmpotp>1) Then
              lm1 = 2
              Do lm = 2, lmpotp
                If (lm1/=1) Then
                  If (isave==1) Then
                    Read (21, Fmt=190) lm1
                  Else
                    lm1 = lm
                  End If
                  If (lm1>1) Then
                    Read (21, Fmt=200)(u(ir), ir=irminp, nr)
                    If (lm1<=lmpot) Then
                      Do ir = irminm, nr
                        vinsimp(ir, lm1, i) = u(ir)
                      End Do
                    End If
                  End If
                End If
              End Do
            End If
          End If ! (INS.LT.1)
          irws1 = irwsimp(ih)

!---> redefine new mt-radius in case of shape corrections

          If (ins/=0) Then
            rmtnew(ih) = scale(1)*alat*xrn(1, ih)
            imt1 = anint(log(rmtnew(ih)/b(ih)+1.0E0_dp)/a(ih)) + 1

!---> for proper core treatment imt must be odd
!     shift potential by one mesh point if imt is even

            If (mod(imt1,2)==0) Then
              imt1 = imt1 + 1
              Do ir = imt1, 2, -1
                vm2zimp(ir, i) = vm2zimp(ir-1, i)
              End Do
            End If

            imt(ih) = imt1
            b(ih) = rmtnew(ih)/(exp(a(ih)*real(imt1-1,kind=dp))-1.0E0_dp)
          End If ! (INS.NE.0)

!---> generate radial mesh - potential only is stored in potential card
!     INEW = 1
!     p. zahn, jan. 99

          a1 = a(ih)
          b1 = b(ih)
          rimp(1, ih) = 0.0E0_dp
          drdi(1, ih) = a1*b1
          Do ir = 2, irws1
            ea = exp(a1*real(ir-1,kind=dp))
            rimp(ir, ih) = b1*(ea-1.0E0_dp)
            drdi(ir, ih) = a1*b1*ea
            dror(ir, ih) = a1/(1.0E0_dp-1.0E0_dp/ea)
          End Do

!---> fill cell-type depending mesh points in the non-muffin-tin-region

          If (ins/=0) Then
            Do iri = 1, meshn(ih)
              ir = iri + imt1
              rimp(ir, ih) = scale(1)*alat*xrn(iri, ih)
              drdi(ir, ih) = scale(1)*alat*drn(iri, ih)
              dror(ir, ih) = drdi(ir, ih)/rimp(ir, ih)
            End Do
          End If

          rws(ih) = rimp(irws1, ih)

!---> kshape.eq.0 : calculate new rmt adapted to exp. mesh

          Call calrmt(ipf, ipfe, ipe, imt(ih), zimp(ih), rmt(ih), rws(ih), &
            rmtnew(ih), alat, drdi(1,ih), a(ih), b(ih), irws1, rimp(1,ih), io, &
            ins)

          If (ins>0) Then
            ircutimp(1, ih) = imt(ih)
            isum = imt(ih)
            Do ipan1 = 2, ipanimp(ih)
              isum = isum + nm(ipan1, ih)
              ircutimp(ipan1, ih) = isum
            End Do
            nr = isum
          Else ! INS.EQ.0
            nr = irwsimp(ih)
            If (kws>=1) Then
              ircutimp(1, ih) = irws1
            Else
              ircutimp(1, ih) = imt(ih)
            End If
          End If ! INS.GT.0

!---> fill array irmin in case of full potential
          If (ins/=0) irminimp(ih) = nr - irns1p

!---> cut input potential at rmt if given only at exponential mesh
          If (ins>1) Then
            imt1 = imt(ih)
            irc1 = ircutimp(ipanimp(ih), ih)
            Call potcut(imt1, irc1, ins, lmpot, rimp(1,ih), vm2zimp(1,i), &
              dummy2, vinsimp(irmd-irnsd,1,i), zimp(ih), irmd, irmd-irnsd)
          End If

          If (ins==0 .And. kws==0) Then
!---> in case of a mt calculation cut potential at mt radius
            imt1 = imt(ih)
            irws1 = irwsimp(ih)
            Call potcut(imt1, irws1, ins, lmpot, rimp(1,ih), vm2zimp(1,i), &
              dummy2, vinsimp(irmd-irnsd,1,i), zimp(ih), irmd, irmd-irnsd)

          End If ! INS.EQ.0 .AND. KWS.EQ.0
!--->       maybe apply a magnetic field
          If (khfeld==1) Then
            Write (1337, *) 'ATOM', ih, 'SPIN', ispin, 'SHIFTED BY', &
              -real(2*ispin-3, kind=dp)*hfield
            Do j = 1, ircutimp(ipanimp(ih), ih)
              vm2zimp(j, i) = vm2zimp(j, i) - real(2*ispin-3, kind=dp)*hfield
            End Do
          End If

        End Do ! ISPIN = 1,NSPIN
      End Do ! IH = 1,NCELL
      Close (20)
      Close (21)

      Return


100   Format (16I5)
110   Format (4D20.12)
120   Format (20A4)
130   Format (3F12.8)
140   Format (F10.5, /, F10.5, 2F15.10)
150   Format (I4)
160   Format (2D15.8, /, 2I2)
170   Format (1P, 2D15.6, 1P, D15.8)
180   Format (I5, 1P, D20.11)
190   Format (10I5)
200   Format (1P, 4D20.13)
    End Subroutine
