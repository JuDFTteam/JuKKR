!------------------------------------------------------------------------------------
!> Summary: Gaunt coefficient type
!> Author: 
!> Category: KKRimp, special-functions
!> Deprecated: False 
!> Contains various informations about the gauntcoefficients $C_{LL'}^{L''}$.
!> @note The gaunt coefficient are constructed for a certain lmmax and only the non-zero gaunt coefficients are stored.
!> iend gaunt coefficients are stored. icleb points to all of the non-zero gaunt coefficients.
!> @endnote
!> @note The $C_{LL}^{(0,0)}$ gaunt coefficient is not stored in this type!
!> @endnote
!------------------------------------------------------------------------------------
module type_gauntcoeff
use nrtype
  type                            :: gauntcoeff_type
     integer,allocatable             :: icleb(:,:)           !! pointer array, first index loops over all gaunt coefficients, second index 1:3 gives access to $L$,$L'$ and $L''$
     integer,allocatable             :: loflm(:)             !! l of lm=(l,m) (gaunt)
     real(kind=dp),allocatable       :: cleb(:,:)            !! gaunt coefficients (gaunt) (1:iend,??), what is the second index?
     integer                         :: iend                 !! number of non-zero gaunt coefficients
     integer,allocatable             :: jend(:,:,:)          !! pointer array for icleb()
     integer                         :: ncleb                !! maximal number of gaunt coefficents (also the non-zero ones) 
     real(kind=dp),allocatable       :: wg(:)                !!
     real(kind=dp),allocatable       :: yrg(:,:,:)           !!
  end type gauntcoeff_type


end module type_gauntcoeff
