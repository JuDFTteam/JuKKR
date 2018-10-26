!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_bastrmat

  private
  public :: bastrmat

contains

  !-------------------------------------------------------------------------------
  !> Summary: transformation kappa, mu - l,m
  !> Author: Hubert Ebert
  !> date: 13/01/98
  !> Category: KKRhost, dirac, special-functions
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> INITIALIZE TRANSFORMATION MATRIX THAT TAKES MATRICES FROM
  !> RELATIVISTIC  TO  REAL SPERICAL HARM.  REPRESENTATION    
  !>                                                          
  !> this is a special version of <STRSMAT> passing the       
  !> full BASis TRansformation MATrices  RC, CREL and RREL    
  !-------------------------------------------------------------------------------
  subroutine bastrmat(lmax, cgc, rc, crel, rrel, nkmmax, nkmpmax)
    use :: mod_cinit, only: cinit
    use :: mod_datatypes, only: dp
    implicit none

    ! PARAMETER definitions
    complex (kind=dp) :: ci, c1, c0
    parameter (ci=(0.0e0_dp,1.0e0_dp), c1=(1.0e0_dp,0.0e0_dp), c0=(0.0e0_dp,0.0e0_dp))

    ! Dummy arguments
    integer :: lmax, nkmmax, nkmpmax
    real (kind=dp) :: cgc(nkmpmax, 2)
    complex (kind=dp) :: crel(nkmmax, nkmmax), rc(nkmmax, nkmmax), rrel(nkmmax, nkmmax)

    ! Local variables
    integer :: i, ikm, j, jp05, k, l, lm, lnr, m, muem05, muep05, nk, nkm, nlm
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

  end subroutine bastrmat

end module mod_bastrmat
