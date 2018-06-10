subroutine calcsph(nsra, irmdnew, nrmaxd, lmax, nspin, z, c, e, lmpotd, &
  lmmaxso, rnew, vins, ncheb, npan_tot, rpan_intervall, jlk_index, hlk, jlk, &
  hlk2, jlk2, gmatprefactor, tmat, alpha, use_sratrick)

  use :: constants
  use :: profiling

  implicit none
! construct wavefunctions for spherical potentials
  integer :: nsra, irmdnew, nrmaxd, nspin, lmax, lmpotd, lmmaxso
  integer :: ncheb, npan_tot
  integer :: use_sratrick
  double precision :: z, c
  double complex :: e, gmatprefactor
  double precision :: rnew(nrmaxd), rpan_intervall(0:npan_tot)
  double precision :: vins(irmdnew, lmpotd, nspin)
  double complex :: hlk(1:4*(lmax+1), irmdnew)
  double complex :: jlk(1:4*(lmax+1), irmdnew)
  double complex :: hlk2(1:4*(lmax+1), irmdnew)
  double complex :: jlk2(1:4*(lmax+1), irmdnew)
  integer :: jlk_index(2*lmmaxso)

! local
  integer :: lmsize, lmsize2, nvec, nspintemp
  integer :: ivec, lval, ir, ispin, lspin, lsra, i, l1, m1, lm1
  integer, allocatable :: jlk_indextemp(:)
  double complex, allocatable :: vll0(:, :, :)
  double complex, allocatable :: vll(:, :, :)
  double complex, allocatable :: rlltemp(:, :, :), slltemp(:, :, :)
  double complex, allocatable :: hlktemp(:, :), jlktemp(:, :)
  double complex, allocatable :: hlk2temp(:, :), jlk2temp(:, :)
  double complex, allocatable :: hlknew(:, :), jlknew(:, :)
  double complex, allocatable :: tmattemp(:, :)
  double complex, allocatable :: alphatemp(:, :) ! LLY
  double complex :: tmat(2*(lmax+1))
  double complex :: alpha(2*(lmax+1)) ! LLY

  lmsize = 1
  if (nsra==2) then
    lmsize2 = 2
    nvec = 2
  else
    lmsize2 = 1
    nvec = 1
  end if
  allocate (rlltemp(lmsize2,lmsize,irmdnew))
  allocate (slltemp(lmsize2,lmsize,irmdnew))
  allocate (hlktemp(nvec,irmdnew))
  allocate (jlktemp(nvec,irmdnew))
  allocate (hlk2temp(nvec,irmdnew))
  allocate (jlk2temp(nvec,irmdnew))
  allocate (jlk_indextemp(lmsize2))
  allocate (tmattemp(lmsize,lmsize))
  allocate (alphatemp(lmsize,lmsize)) ! LLY
  allocate (hlknew(nvec*nspin*(lmax+1),irmdnew))
  allocate (jlknew(nvec*nspin*(lmax+1),irmdnew))

  do ivec = 1, nvec
    jlk_indextemp(ivec) = ivec
  end do
  allocate (vll0(lmsize,lmsize,irmdnew))
  if (nsra==2) then
    allocate (vll(2*lmsize,2*lmsize,irmdnew))
  else
    allocate (vll(lmsize,lmsize,irmdnew))
  end if

  alpha = czero
! spin loop
  nspintemp = nspin
!       IF (NSRA.EQ.2) THEN
!         NSPINTEMP=NSPIN
!       ELSEIF (NSRA.EQ.1) THEN
!         NSPINTEMP=1
!       ENDIF

  do ispin = 1, nspintemp

    lspin = (lmax+1)*(ispin-1)
    lsra = (lmax+1)*nvec
! each value of l, the Lippmann-Schwinger equation is solved using
! the free-potential wavefunctions and potentials corresponding to l-value
    do lval = 0, lmax

      do ir = 1, irmdnew
        vll0(lmsize, lmsize, ir) = vins(ir, 1, ispin) - 2d0*z/rnew(ir)
      end do
      if (nsra==2) then
        call vllmatsra(vll0, vll, rnew, lmsize, irmdnew, nrmaxd, e, lmax, &
          lval, 'Ref=0')
      else
        vll(:, :, :) = vll0(:, :, :)
      end if

      jlktemp(1, :) = jlk(lval+1, :)
      hlktemp(1, :) = hlk(lval+1, :)
      jlk2temp(1, :) = jlk2(lval+1, :)
      hlk2temp(1, :) = hlk2(lval+1, :)
      if (nsra==2) then
        jlktemp(2, :) = jlk(lmax+lval+2, :)
        hlktemp(2, :) = hlk(lmax+lval+2, :)
        jlk2temp(2, :) = jlk2(lmax+lval+2, :)
        hlk2temp(2, :) = hlk2(lmax+lval+2, :)
      end if
      tmattemp = (0d0, 0d0)
      alphatemp = (0d0, 0d0)
      call rllsll(rpan_intervall, rnew, vll, rlltemp, slltemp, tmattemp, &
        ncheb, npan_tot, lmsize, lmsize2, nvec, irmdnew, nrmaxd, nvec, &
        jlk_indextemp, hlktemp, jlktemp, hlk2temp, jlk2temp, gmatprefactor, &
        '1', '1', '0', use_sratrick, alphatemp) ! LLY

      do ir = 1, irmdnew
        hlknew(lspin+lval+1, ir) = slltemp(1, 1, ir)/rnew(ir)
        jlknew(lspin+lval+1, ir) = rlltemp(1, 1, ir)/rnew(ir)
      end do
      if (nsra==2) then
        do ir = 1, irmdnew
          hlknew(lspin+lsra+lval+1, ir) = slltemp(2, 1, ir)/rnew(ir)
          jlknew(lspin+lsra+lval+1, ir) = rlltemp(2, 1, ir)/rnew(ir)
        end do
      end if
      tmat(lspin+lval+1) = tmattemp(1, 1)
      alpha(lspin+lval+1) = alphatemp(1, 1) ! LLY
    end do ! LMAX
  end do ! NSPIN

  lm1 = 1
  do ivec = 1, nvec
    do i = 1, 2
      do l1 = 0, lmax
        do m1 = -l1, l1
          jlk_index(lm1) = l1 + (ivec-1)*nspintemp*(lmax+1) + (i-1)*(lmax+1) + &
            1
          lm1 = lm1 + 1
        end do
      end do
    end do
  end do
  do ir = 1, irmdnew
    do l1 = 1, nvec*(lmax+1)*nspintemp
      hlk(l1, ir) = hlknew(l1, ir)
      jlk(l1, ir) = jlknew(l1, ir)
    end do
  end do
  if (nsra==2) then
    do ir = 1, irmdnew
      do l1 = 1, (lmax+1)*nspintemp
        hlk2(l1, ir) = -hlknew(l1+lmax+1, ir)
        jlk2(l1, ir) = -jlknew(l1+lmax+1, ir)
      end do
      do l1 = nspintemp*(lmax+1) + 1, nvec*(lmax+1)*nspintemp
        hlk2(l1, ir) = hlknew(l1-(lmax+1)*nspintemp, ir)
        jlk2(l1, ir) = jlknew(l1-(lmax+1)*nspintemp, ir)
      end do
    end do
  else
    do ir = 1, irmdnew
      do l1 = 1, nvec*(lmax+1)*nspintemp
        hlk2(l1, ir) = -hlknew(l1, ir)
        jlk2(l1, ir) = -jlknew(l1, ir)
      end do
    end do
  end if
  deallocate (rlltemp)
  deallocate (slltemp)
  deallocate (hlktemp)
  deallocate (jlktemp)
  deallocate (hlk2temp)
  deallocate (jlk2temp)
  deallocate (jlk_indextemp)
  deallocate (tmattemp)
  deallocate (alphatemp) ! LLY
  deallocate (hlknew)
  deallocate (jlknew)
  deallocate (vll0)
  deallocate (vll)

end subroutine
