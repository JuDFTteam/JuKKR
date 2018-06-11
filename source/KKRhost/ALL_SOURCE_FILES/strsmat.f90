    Subroutine strsmat(lmax, cgc, srrel, nrrel, irrel, nkmmax, nkmpmax)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *    INITIALIZE TRANSFORMATION MATRIX THAT TAKES MATRICES FROM     *
!   *    RELATIVISTIC  TO  REAL SPERICAL HARM.  REPRESENTATION         *
!   *                                                                  *
!   *    ONLY THE NON-0 ELEMENTS OF THE MATRIX ARE STORED              *
!   *                                                                  *
!   * 25/10/95  HE  proper convention of trans. matrix introduced      *
!   ********************************************************************

      Implicit None

! PARAMETER definitions
      Complex (Kind=dp) :: ci, c1, c0
      Parameter (ci=(0.0E0_dp,1.0E0_dp), c1=(1.0E0_dp,0.0E0_dp), &
        c0=(0.0E0_dp,0.0E0_dp))

! Dummy arguments
      Integer :: lmax, nkmmax, nkmpmax
      Real (Kind=dp) :: cgc(nkmpmax, 2)
      Integer :: irrel(2, 2, nkmmax), nrrel(2, nkmmax)
      Complex (Kind=dp) :: srrel(2, 2, nkmmax)

! Local variables
      Complex (Kind=dp) :: crel(nkmmax, nkmmax), rc(nkmmax, nkmmax), &
        rrel(nkmmax, nkmmax)
      Integer :: i, ikm, j, jp05, k, l, lam, lm, lnr, lr, m, muem05, muep05, &
        nk, nkm, nlm, ns1, ns2
      Real (Kind=dp) :: w

      nk = 2*(lmax+1) + 1
      nlm = (lmax+1)**2
      nkm = 2*nlm
!     ===================================================
!     INDEXING:
!     IKM  = L*2*(J+1/2) + J + MUE + 1
!     LM   = L*(L+1)     +     M   + 1
!     ===================================================

! ----------------------------------------------------------------------
! CREL  transforms from  COMPLEX (L,M,S)  to  (KAP,MUE) - representation
!                 |LAM> = sum[LC] |LC> * CREL(LC,LAM)
! ----------------------------------------------------------------------
      Call cinit(nkmmax*nkmmax, crel)

      lm = 0
      Do lnr = 0, lmax
        Do m = -lnr, lnr
          lm = lm + 1

          ikm = 0
          Do k = 1, nk
            l = k/2
            If (2*l==k) Then
              jp05 = l
            Else
              jp05 = l + 1
            End If

            Do muem05 = -jp05, (jp05-1)
              muep05 = muem05 + 1
              ikm = ikm + 1

              If (l==lnr) Then
                If (muep05==m) crel(lm, ikm) = cgc(ikm, 1)
                If (muem05==m) crel(lm+nlm, ikm) = cgc(ikm, 2)
              End If

            End Do
          End Do

        End Do
      End Do

! ----------------------------------------------------------------------
!    RC  transforms from  REAL to  COMPLEX (L,M,S) - representation
!                 |LC> = sum[LR] |LR> * RC(LR,LC)
! ----------------------------------------------------------------------
      Call cinit(nkmmax*nkmmax, rc)

      w = 1.0E0_dp/sqrt(2.0E0_dp)

      Do l = 0, lmax
        Do m = -l, l
          i = l*(l+1) + m + 1
          j = l*(l+1) - m + 1

          If (m<0) Then
            rc(i, i) = -ci*w
            rc(j, i) = w
            rc(i+nlm, i+nlm) = -ci*w
            rc(j+nlm, i+nlm) = w
          End If
          If (m==0) Then
            rc(i, i) = c1
            rc(i+nlm, i+nlm) = c1
          End If
          If (m>0) Then
            rc(i, i) = w*(-1.0E0_dp)**m
            rc(j, i) = ci*w*(-1.0E0_dp)**m
            rc(i+nlm, i+nlm) = w*(-1.0E0_dp)**m
            rc(j+nlm, i+nlm) = ci*w*(-1.0E0_dp)**m
          End If
        End Do
      End Do

! ----------------------------------------------------------------------
! RREL  transforms from   REAL (L,M,S)  to  (KAP,MUE) - representation
!                 |LAM> = sum[LR] |LR> * RREL(LR,LAM)
! ----------------------------------------------------------------------
      Call zgemm('N', 'N', nkm, nkm, nkm, c1, rc, nkmmax, crel, nkmmax, c0, &
        rrel, nkmmax)

!     ---------------------------------------------------
!     store the elements of  RREL
!     ---------------------------------------------------
      Do lam = 1, nkm
        ns1 = 0
        ns2 = 0

        Do lr = 1, 2*nlm
!             IF ( CDABS(RREL(LR,LAM)).GT.1D-6 ) THEN
          If (abs(rrel(lr,lam))>1E-4_dp) Then
            If (lr<=nlm) Then
              ns1 = ns1 + 1
              If (ns1>2) Stop ' IN <STRSMAT>   NS1 > 2'
              srrel(ns1, 1, lam) = rrel(lr, lam)
              irrel(ns1, 1, lam) = lr
            Else
              ns2 = ns2 + 1
              If (ns2>2) Stop ' IN <STRSMAT>   NS2 > 2'
              srrel(ns2, 2, lam) = rrel(lr, lam)
              irrel(ns2, 2, lam) = lr - nlm
            End If
          End If
        End Do

        nrrel(1, lam) = ns1
        nrrel(2, lam) = ns2
      End Do

    End Subroutine
