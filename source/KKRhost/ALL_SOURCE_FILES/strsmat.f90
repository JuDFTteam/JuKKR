module mod_strsmat
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine strsmat(lmax, cgc, srrel, nrrel, irrel, nkmmax, nkmpmax)
    ! ********************************************************************
    ! *                                                                  *
    ! *    INITIALIZE TRANSFORMATION MATRIX THAT TAKES MATRICES FROM     *
    ! *    RELATIVISTIC  TO  REAL SPERICAL HARM.  REPRESENTATION         *
    ! *                                                                  *
    ! *    ONLY THE NON-0 ELEMENTS OF THE MATRIX ARE STORED              *
    ! *                                                                  *
    ! * 25/10/95  HE  proper convention of trans. matrix introduced      *
    ! ********************************************************************

    use :: mod_cinit
    implicit none

    ! PARAMETER definitions
    complex (kind=dp) :: ci, c1, c0
    parameter (ci=(0.0e0_dp,1.0e0_dp), c1=(1.0e0_dp,0.0e0_dp), c0=(0.0e0_dp,0.0e0_dp))

    ! Dummy arguments
    integer :: lmax, nkmmax, nkmpmax
    real (kind=dp) :: cgc(nkmpmax, 2)
    integer :: irrel(2, 2, nkmmax), nrrel(2, nkmmax)
    complex (kind=dp) :: srrel(2, 2, nkmmax)

    ! Local variables
    complex (kind=dp) :: crel(nkmmax, nkmmax), rc(nkmmax, nkmmax), rrel(nkmmax, nkmmax)
    integer :: i, ikm, j, jp05, k, l, lam, lm, lnr, lr, m, muem05, muep05, nk, nkm, nlm, ns1, ns2
    real (kind=dp) :: w

    nk = 2*(lmax+1) + 1
    nlm = (lmax+1)**2
    nkm = 2*nlm
    ! ===================================================
    ! INDEXING:
    ! IKM  = L*2*(J+1/2) + J + MUE + 1
    ! LM   = L*(L+1)     +     M   + 1
    ! ===================================================

    ! ----------------------------------------------------------------------
    ! CREL  transforms from  COMPLEX (L,M,S)  to  (KAP,MUE) - representation
    ! |LAM> = sum[LC] |LC> * CREL(LC,LAM)
    ! ----------------------------------------------------------------------
    call cinit(nkmmax*nkmmax, crel)

    lm = 0
    do lnr = 0, lmax
      do m = -lnr, lnr
        lm = lm + 1

        ikm = 0
        do k = 1, nk
          l = k/2
          if (2*l==k) then
            jp05 = l
          else
            jp05 = l + 1
          end if

          do muem05 = -jp05, (jp05-1)
            muep05 = muem05 + 1
            ikm = ikm + 1

            if (l==lnr) then
              if (muep05==m) crel(lm, ikm) = cgc(ikm, 1)
              if (muem05==m) crel(lm+nlm, ikm) = cgc(ikm, 2)
            end if

          end do
        end do

      end do
    end do

    ! ----------------------------------------------------------------------
    ! RC  transforms from  REAL to  COMPLEX (L,M,S) - representation
    ! |LC> = sum[LR] |LR> * RC(LR,LC)
    ! ----------------------------------------------------------------------
    call cinit(nkmmax*nkmmax, rc)

    w = 1.0e0_dp/sqrt(2.0e0_dp)

    do l = 0, lmax
      do m = -l, l
        i = l*(l+1) + m + 1
        j = l*(l+1) - m + 1

        if (m<0) then
          rc(i, i) = -ci*w
          rc(j, i) = w
          rc(i+nlm, i+nlm) = -ci*w
          rc(j+nlm, i+nlm) = w
        end if
        if (m==0) then
          rc(i, i) = c1
          rc(i+nlm, i+nlm) = c1
        end if
        if (m>0) then
          rc(i, i) = w*(-1.0e0_dp)**m
          rc(j, i) = ci*w*(-1.0e0_dp)**m
          rc(i+nlm, i+nlm) = w*(-1.0e0_dp)**m
          rc(j+nlm, i+nlm) = ci*w*(-1.0e0_dp)**m
        end if
      end do
    end do

    ! ----------------------------------------------------------------------
    ! RREL  transforms from   REAL (L,M,S)  to  (KAP,MUE) - representation
    ! |LAM> = sum[LR] |LR> * RREL(LR,LAM)
    ! ----------------------------------------------------------------------
    call zgemm('N', 'N', nkm, nkm, nkm, c1, rc, nkmmax, crel, nkmmax, c0, rrel, nkmmax)

    ! ---------------------------------------------------
    ! store the elements of  RREL
    ! ---------------------------------------------------
    do lam = 1, nkm
      ns1 = 0
      ns2 = 0

      do lr = 1, 2*nlm
        ! IF ( CDABS(RREL(LR,LAM)).GT.1D-6 ) THEN
        if (abs(rrel(lr,lam))>1e-4_dp) then
          if (lr<=nlm) then
            ns1 = ns1 + 1
            if (ns1>2) stop ' IN <STRSMAT>   NS1 > 2'
            srrel(ns1, 1, lam) = rrel(lr, lam)
            irrel(ns1, 1, lam) = lr
          else
            ns2 = ns2 + 1
            if (ns2>2) stop ' IN <STRSMAT>   NS2 > 2'
            srrel(ns2, 2, lam) = rrel(lr, lam)
            irrel(ns2, 2, lam) = lr - nlm
          end if
        end if
      end do

      nrrel(1, lam) = ns1
      nrrel(2, lam) = ns2
    end do

  end subroutine strsmat

end module mod_strsmat
