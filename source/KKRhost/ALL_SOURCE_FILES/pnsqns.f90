!-------------------------------------------------------------------------------
! SUBROUTINE: PNSQNS
!> @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
! Added IRMIN 1.7.2014
    Subroutine pnsqns(ar, cr, dr, drdi, ek, icst, pz, qz, fz, sz, pns, qns, &
      nsra, vins, ipan, irmin, ircut, cleb, icleb, iend, loflm, lkonv, &
      idoldau, lopt, lmlo, lmhi, wldau, wldauav, cutoff, mmaxd, lmpot, irmind, &
      lmmaxd, irm, lmax)

      Use global_variables
      Use mod_datatypes, Only: dp

      Implicit None
!
! .. Input variables
      Integer, Intent (In) :: irm !< Maximum number of radial points
      Integer, Intent (In) :: lmax !< Maximum l component in wave function expansion
      Integer, Intent (In) :: icst !< Number of Born approximation
      Integer, Intent (In) :: iend !< Number of nonzero gaunt coefficients
      Integer, Intent (In) :: ipan !< Number of panels in non-MT-region
      Integer, Intent (In) :: lmlo
      Integer, Intent (In) :: lmhi
      Integer, Intent (In) :: lopt !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
      Integer, Intent (In) :: nsra
      Integer, Intent (In) :: lkonv
      Integer, Intent (In) :: irmin !< Max R for spherical treatment
      Integer, Intent (In) :: mmaxd !< 2*LMAX+1
      Integer, Intent (In) :: lmpot !< (LPOT+1)**2
      Integer, Intent (In) :: irmind !< IRM-IRNSD
      Integer, Intent (In) :: lmmaxd !< (KREL+KORBIT+1)(LMAX+1)^2
      Integer, Intent (In) :: idoldau !< flag to perform LDA+U
      Real (Kind=dp), Intent (In) :: wldauav
      Complex (Kind=dp), Intent (In) :: ek
      Integer, Dimension (0:ipand), Intent (In) :: ircut !< R points of panel borders
      Integer, Dimension (*), Intent (In) :: loflm !< l of lm=(l,m) (GAUNT)
      Integer, Dimension (ncleb, 4), Intent (In) :: icleb !< Pointer array
      Real (Kind=dp), Dimension (irm), Intent (In) :: drdi !< Derivative dr/di
      Real (Kind=dp), Dimension (irm), Intent (In) :: cutoff
      Real (Kind=dp), Dimension (ncleb, 2), Intent (In) :: cleb !< GAUNT coefficients (GAUNT)
      Real (Kind=dp), Dimension (irmind:irm, lmpot), Intent (In) :: vins !< Non-spherical part of the potential
      Real (Kind=dp), Dimension (mmaxd, mmaxd), Intent (In) :: wldau !< potential matrix
      Complex (Kind=dp), Dimension (irm, 0:lmax), Intent (In) :: fz
      Complex (Kind=dp), Dimension (irm, 0:lmax), Intent (In) :: qz
      Complex (Kind=dp), Dimension (irm, 0:lmax), Intent (In) :: sz
      Complex (Kind=dp), Dimension (irm, 0:lmax), Intent (In) :: pz
      Complex (Kind=dp), Dimension (lmmaxd, lmmaxd), Intent (In) :: dr
      Complex (Kind=dp), Dimension (lmmaxd, lmmaxd), Intent (In) :: ar
      Complex (Kind=dp), Dimension (lmmaxd, lmmaxd), Intent (In) :: cr
      Complex (Kind=dp), Dimension (lmmaxd, lmmaxd, irmind:irm, 2), &
        Intent (In) :: pns
      Complex (Kind=dp), Dimension (lmmaxd, lmmaxd, irmind:irm, 2), &
        Intent (In) :: qns
! .. Local Scalars
      Integer :: i, lm1, lm2, lmmkonv, m1, m2, ir, irmax
! .. Local Arrays
      Real (Kind=dp), Dimension (lmmaxd, lmmaxd, irmind:irm) :: vnspll
      Complex (Kind=dp), Dimension (lmmaxd) :: efac
      Complex (Kind=dp), Dimension (lmmaxd, lmmaxd) :: tmatll
      Complex (Kind=dp), Dimension (lmmaxd, lmmaxd, irmind:irm) :: dmat
      Complex (Kind=dp), Dimension (lmmaxd, lmmaxd, irmind:irm) :: cmat
      Complex (Kind=dp), Dimension (lmmaxd, irmind:irm, 2) :: pzlm
      Complex (Kind=dp), Dimension (lmmaxd, irmind:irm, 2) :: qzlm
      Complex (Kind=dp), Dimension (lmmaxd, irmind:irm, 2) :: pzekdr
      Complex (Kind=dp), Dimension (lmmaxd, irmind:irm, 2) :: qzekdr
! .. External Subroutines
      External :: irwns, regns, vllns, wftsca
!
      irmax = ircut(ipan) ! Added IRMAX 1.7.2014
      Call vllns(vnspll, vins, cleb, icleb, iend, irm, ncleb, lmpot, irmind, &
        lmmaxd)
      If (lkonv/=lmax) Then
        lmmkonv = (lkonv+1)*(lkonv+1)
        Do lm1 = 1, lmmaxd
          Do lm2 = lmmkonv + 1, lmmaxd
            Do i = irmind, irm
              vnspll(lm2, lm1, i) = 0.0E0_dp
              vnspll(lm1, lm2, i) = 0.0E0_dp
            End Do
          End Do
        End Do
      Else
        lmmkonv = lmmaxd
      End If
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LDA+U
! Add WLDAU to non-spherical porential VINS in case of LDA+U
! Use the average wldau (=wldauav) and the deviation
! of wldau from this. Use the deviation in the Born series
! for the non-spherical wavefunction, while the average is
! used for the spherical wavefunction.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      If (idoldau==1 .And. lopt>=0) Then
        Do ir = irmind, irm
!----------------------------------------------------------------------
! First add wldau to all elements
!----------------------------------------------------------------------
          Do lm2 = lmlo, lmhi
            m2 = lm2 - lmlo + 1
            Do lm1 = lmlo, lmhi
              m1 = lm1 - lmlo + 1
              vnspll(lm1, lm2, ir) = vnspll(lm1, lm2, ir) + &
                wldau(m1, m2)*cutoff(ir)
            End Do
!-------------------------------------------------------------------
! and then subtract average from diag. elements
!-------------------------------------------------------------------
            vnspll(lm2, lm2, ir) = vnspll(lm2, lm2, ir) - wldauav*cutoff(ir)
          End Do
        End Do
      End If
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LDA+U
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------------
! Get wfts of same magnitude by scaling with efac
!----------------------------------------------------------------------------
      Call wftsca(drdi, efac, pz, qz, fz, sz, nsra, pzlm, qzlm, pzekdr, &
        qzekdr, ek, loflm, irmind, irm, irmin, irmax, lmax, lmmaxd) ! Added IRMIN,IRMAX 1.7.2014
!----------------------------------------------------------------------------
! Determine the irregular non sph. wavefunction
!----------------------------------------------------------------------------
! Added IRMIN,IRMAX 1.7.2014
      Call irwns(cr, dr, efac, qns, vnspll, icst, ipan, ircut, nsra, pzlm, &
        qzlm, pzekdr, qzekdr, qns(1,1,irmind,1), cmat, qns(1,1,irmind,2), &
        dmat, irmind, irm, irmin, irmax, ipand, lmmaxd)
!----------------------------------------------------------------------------
! Determine the regular non sph. wavefunction
!----------------------------------------------------------------------------
! Added IRMIN,IRMAX 1.7.2014
      Call regns(ar, tmatll, efac, pns, vnspll, icst, ipan, ircut, pzlm, qzlm, &
        pzekdr, qzekdr, ek, pns(1,1,irmind,1), cmat, pns(1,1,irmind,2), dmat, &
        nsra, irmind, irm, irmin, irmax, ipand, lmmaxd)
!
      Return

    End Subroutine
