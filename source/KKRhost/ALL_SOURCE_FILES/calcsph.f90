    Subroutine calcsph(nsra, irmdnew, nrmaxd, lmax, nspin, z, c, e, lmpotd, &
      lmmaxso, rnew, vins, ncheb, npan_tot, rpan_intervall, jlk_index, hlk, &
      jlk, hlk2, jlk2, gmatprefactor, tmat, alpha, use_sratrick)

      Use constants
      Use profiling
      Use mod_datatypes, Only: dp

      Implicit None
! construct wavefunctions for spherical potentials
      Integer :: nsra, irmdnew, nrmaxd, nspin, lmax, lmpotd, lmmaxso
      Integer :: ncheb, npan_tot
      Integer :: use_sratrick
      Real (Kind=dp) :: z, c
      Complex (Kind=dp) :: e, gmatprefactor
      Real (Kind=dp) :: rnew(nrmaxd), rpan_intervall(0:npan_tot)
      Real (Kind=dp) :: vins(irmdnew, lmpotd, nspin)
      Complex (Kind=dp) :: hlk(1:4*(lmax+1), irmdnew)
      Complex (Kind=dp) :: jlk(1:4*(lmax+1), irmdnew)
      Complex (Kind=dp) :: hlk2(1:4*(lmax+1), irmdnew)
      Complex (Kind=dp) :: jlk2(1:4*(lmax+1), irmdnew)
      Integer :: jlk_index(2*lmmaxso)

! local
      Integer :: lmsize, lmsize2, nvec, nspintemp
      Integer :: ivec, lval, ir, ispin, lspin, lsra, i, l1, m1, lm1
      Integer, Allocatable :: jlk_indextemp(:)
      Complex (Kind=dp), Allocatable :: vll0(:, :, :)
      Complex (Kind=dp), Allocatable :: vll(:, :, :)
      Complex (Kind=dp), Allocatable :: rlltemp(:, :, :), slltemp(:, :, :)
      Complex (Kind=dp), Allocatable :: hlktemp(:, :), jlktemp(:, :)
      Complex (Kind=dp), Allocatable :: hlk2temp(:, :), jlk2temp(:, :)
      Complex (Kind=dp), Allocatable :: hlknew(:, :), jlknew(:, :)
      Complex (Kind=dp), Allocatable :: tmattemp(:, :)
      Complex (Kind=dp), Allocatable :: alphatemp(:, :) ! LLY
      Complex (Kind=dp) :: tmat(2*(lmax+1))
      Complex (Kind=dp) :: alpha(2*(lmax+1)) ! LLY

      lmsize = 1
      If (nsra==2) Then
        lmsize2 = 2
        nvec = 2
      Else
        lmsize2 = 1
        nvec = 1
      End If
      Allocate (rlltemp(lmsize2,lmsize,irmdnew))
      Allocate (slltemp(lmsize2,lmsize,irmdnew))
      Allocate (hlktemp(nvec,irmdnew))
      Allocate (jlktemp(nvec,irmdnew))
      Allocate (hlk2temp(nvec,irmdnew))
      Allocate (jlk2temp(nvec,irmdnew))
      Allocate (jlk_indextemp(lmsize2))
      Allocate (tmattemp(lmsize,lmsize))
      Allocate (alphatemp(lmsize,lmsize)) ! LLY
      Allocate (hlknew(nvec*nspin*(lmax+1),irmdnew))
      Allocate (jlknew(nvec*nspin*(lmax+1),irmdnew))

      Do ivec = 1, nvec
        jlk_indextemp(ivec) = ivec
      End Do
      Allocate (vll0(lmsize,lmsize,irmdnew))
      If (nsra==2) Then
        Allocate (vll(2*lmsize,2*lmsize,irmdnew))
      Else
        Allocate (vll(lmsize,lmsize,irmdnew))
      End If

      alpha = czero
! spin loop
      nspintemp = nspin
!       IF (NSRA.EQ.2) THEN
!         NSPINTEMP=NSPIN
!       ELSEIF (NSRA.EQ.1) THEN
!         NSPINTEMP=1
!       ENDIF

      Do ispin = 1, nspintemp

        lspin = (lmax+1)*(ispin-1)
        lsra = (lmax+1)*nvec
! each value of l, the Lippmann-Schwinger equation is solved using
! the free-potential wavefunctions and potentials corresponding to l-value
        Do lval = 0, lmax

          Do ir = 1, irmdnew
            vll0(lmsize, lmsize, ir) = vins(ir, 1, ispin) - 2E0_dp*z/rnew(ir)
          End Do
          If (nsra==2) Then
            Call vllmatsra(vll0, vll, rnew, lmsize, irmdnew, nrmaxd, e, lmax, &
              lval, 'Ref=0')
          Else
            vll(:, :, :) = vll0(:, :, :)
          End If

          jlktemp(1, :) = jlk(lval+1, :)
          hlktemp(1, :) = hlk(lval+1, :)
          jlk2temp(1, :) = jlk2(lval+1, :)
          hlk2temp(1, :) = hlk2(lval+1, :)
          If (nsra==2) Then
            jlktemp(2, :) = jlk(lmax+lval+2, :)
            hlktemp(2, :) = hlk(lmax+lval+2, :)
            jlk2temp(2, :) = jlk2(lmax+lval+2, :)
            hlk2temp(2, :) = hlk2(lmax+lval+2, :)
          End If
          tmattemp = (0E0_dp, 0E0_dp)
          alphatemp = (0E0_dp, 0E0_dp)
          Call rllsll(rpan_intervall, rnew, vll, rlltemp, slltemp, tmattemp, &
            ncheb, npan_tot, lmsize, lmsize2, nvec, irmdnew, nrmaxd, nvec, &
            jlk_indextemp, hlktemp, jlktemp, hlk2temp, jlk2temp, &
            gmatprefactor, '1', '1', '0', use_sratrick, alphatemp) ! LLY

          Do ir = 1, irmdnew
            hlknew(lspin+lval+1, ir) = slltemp(1, 1, ir)/rnew(ir)
            jlknew(lspin+lval+1, ir) = rlltemp(1, 1, ir)/rnew(ir)
          End Do
          If (nsra==2) Then
            Do ir = 1, irmdnew
              hlknew(lspin+lsra+lval+1, ir) = slltemp(2, 1, ir)/rnew(ir)
              jlknew(lspin+lsra+lval+1, ir) = rlltemp(2, 1, ir)/rnew(ir)
            End Do
          End If
          tmat(lspin+lval+1) = tmattemp(1, 1)
          alpha(lspin+lval+1) = alphatemp(1, 1) ! LLY
        End Do ! LMAX
      End Do ! NSPIN

      lm1 = 1
      Do ivec = 1, nvec
        Do i = 1, 2
          Do l1 = 0, lmax
            Do m1 = -l1, l1
              jlk_index(lm1) = l1 + (ivec-1)*nspintemp*(lmax+1) + &
                (i-1)*(lmax+1) + 1
              lm1 = lm1 + 1
            End Do
          End Do
        End Do
      End Do
      Do ir = 1, irmdnew
        Do l1 = 1, nvec*(lmax+1)*nspintemp
          hlk(l1, ir) = hlknew(l1, ir)
          jlk(l1, ir) = jlknew(l1, ir)
        End Do
      End Do
      If (nsra==2) Then
        Do ir = 1, irmdnew
          Do l1 = 1, (lmax+1)*nspintemp
            hlk2(l1, ir) = -hlknew(l1+lmax+1, ir)
            jlk2(l1, ir) = -jlknew(l1+lmax+1, ir)
          End Do
          Do l1 = nspintemp*(lmax+1) + 1, nvec*(lmax+1)*nspintemp
            hlk2(l1, ir) = hlknew(l1-(lmax+1)*nspintemp, ir)
            jlk2(l1, ir) = jlknew(l1-(lmax+1)*nspintemp, ir)
          End Do
        End Do
      Else
        Do ir = 1, irmdnew
          Do l1 = 1, nvec*(lmax+1)*nspintemp
            hlk2(l1, ir) = -hlknew(l1, ir)
            jlk2(l1, ir) = -jlknew(l1, ir)
          End Do
        End Do
      End If
      Deallocate (rlltemp)
      Deallocate (slltemp)
      Deallocate (hlktemp)
      Deallocate (jlktemp)
      Deallocate (hlk2temp)
      Deallocate (jlk2temp)
      Deallocate (jlk_indextemp)
      Deallocate (tmattemp)
      Deallocate (alphatemp) ! LLY
      Deallocate (hlknew)
      Deallocate (jlknew)
      Deallocate (vll0)
      Deallocate (vll)

    End Subroutine
