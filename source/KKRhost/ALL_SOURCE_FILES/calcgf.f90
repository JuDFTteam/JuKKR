module mod_calcgf

  private
  public :: calcgf

contains

  !-------------------------------------------------------------------------------
  !> Summary: G and F coefficients of kappa-mue representation
  !> Author: 
  !> Category: KKRhost, dirac
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> G- AND F-COEFFICIENTS                                    
  !>                                                          
  !> G(K,K',MUE) = CGC(K,MUE,2)*CGC(K',MUE,2) - CGC(K,MUE,1)*CGC(K',MUE,1)
  !>                                                          
  !> GM = G(-K,-K',MUE)                                       
  !>                                                          
  !> F(K,K',MUE) = (MUE-1/2) * CGC(K,MUE,2)*CGC(K',MUE,2) + (MUE+1/2) * CGC(K,MUE,1)*CGC(K',MUE,1)  
  !>                                                          
  !> FM = F(-K,-K',MUE)                                       
  !>                                                          
  !> ..DIA/..OFF ARE THE ELEMENTS FOR  K=K'/K=-K'-1           
  !> IG NUMBERS THE G'S   COLUMN-WISE STARTING WITH COLUMN  1
  !-------------------------------------------------------------------------------
  subroutine calcgf(nk, cgc, gdia, gmdia, goff, gmoff, fdia, fmdia, foff, fmoff, ltab, lbtab, kaptab, nmuetab, nmuemax, nkmmax, nkmpmax)
  
    use :: mod_datatypes, only: dp
    implicit none

    ! Dummy arguments
    integer :: nk, nkmmax, nkmpmax, nmuemax
    real (kind=dp) :: cgc(nkmpmax, 2), fdia(nkmmax), fmdia(nkmmax), fmoff(nkmmax), foff(nkmmax), gdia(nkmmax), gmdia(nkmmax), gmoff(nkmmax), goff(nkmmax)
    integer :: kaptab(nmuemax), lbtab(nmuemax), ltab(nmuemax), nmuetab(nmuemax)

    ! Local variables
    integer :: ig, ikm1, ikm2, imkm1, imkm2, j1p05, j2p05, k1, k2, kap1, kap2, l1, l1bar, l2, l2bar, m1, m2, mldn, mlup, mue1m05, mue2m05


    ! JP05 = J +0.5     MUEP05 = MUE + 0.5
    ! IKM  = L*2*(J+1/2) + J + MUE + 1

    ikm2 = 0
    do k2 = 1, nk
      l2 = ltab(k2)
      kap2 = kaptab(k2)
      l2bar = lbtab(k2)
      j2p05 = abs(kap2)
      mue2m05 = -j2p05 - 1

      do m2 = 1, nmuetab(k2)
        mue2m05 = mue2m05 + 1
        ikm2 = ikm2 + 1
        ig = ikm2
        goff(ig) = 0.0e0_dp
        gmoff(ig) = 0.0e0_dp
        foff(ig) = 0.0e0_dp
        fmoff(ig) = 0.0e0_dp

        imkm2 = l2bar*2*j2p05 + j2p05 + mue2m05 + 1

        ikm1 = 0
        do k1 = 1, nk
          l1 = ltab(k1)
          kap1 = kaptab(k1)
          l1bar = lbtab(k1)
          j1p05 = abs(kap1)
          mue1m05 = -j1p05 - 1

          do m1 = 1, nmuetab(k1)
            mue1m05 = mue1m05 + 1
            ikm1 = ikm1 + 1
            mlup = mue1m05
            mldn = mlup + 1

            if (l1==l2) then
              if (mue1m05==mue2m05) then

                imkm1 = l1bar*2*j1p05 + j1p05 + mue1m05 + 1

                if (kap1==kap2) then

                  gdia(ig) = cgc(ikm1, 2)*cgc(ikm2, 2) - cgc(ikm1, 1)*cgc(ikm2, 1)
                  gmdia(ig) = cgc(imkm1, 2)*cgc(imkm2, 2) - cgc(imkm1, 1)*cgc(imkm2, 1)
                  fdia(ig) = mlup*cgc(ikm1, 2)*cgc(ikm2, 2) + mldn*cgc(ikm1, 1)*cgc(ikm2, 1)
                  fmdia(ig) = mlup*cgc(imkm1, 2)*cgc(imkm2, 2) + mldn*cgc(imkm1, 1)*cgc(imkm2, 1)

                else

                  goff(ig) = cgc(ikm1, 2)*cgc(ikm2, 2) - cgc(ikm1, 1)*cgc(ikm2, 1)
                  gmoff(ig) = cgc(imkm1, 2)*cgc(imkm2, 2) - cgc(imkm1, 1)*cgc(imkm2, 1)
                  foff(ig) = mlup*cgc(ikm1, 2)*cgc(ikm2, 2) + mldn*cgc(ikm1, 1)*cgc(ikm2, 1)
                  fmoff(ig) = mlup*cgc(imkm1, 2)*cgc(imkm2, 2) + mldn*cgc(imkm1, 1)*cgc(imkm2, 1)

                end if

              end if
            end if

          end do ! m1

        end do ! k1

      end do ! m2

    end do ! k2

  end subroutine calcgf

end module mod_calcgf
