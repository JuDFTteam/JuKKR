module type_inpsusc  ! susc |||| This module was added by Benedikt (2014/06)
  ! --------------------------------------------------- !
  ! -- This type should be used when you prepare for -- !
  ! --          running the KKRSUSC program          -- !
  ! --------------------------------------------------- !

  use nrtype

  type                           ::  inpsusc_type

    integer                        :: ia, nr0, nr1, nr, imesh, ilmax, ib, nb              
    real(kind=dp),     allocatable :: rs(:,:), s(:)                                       
    real(kind=dp),     allocatable :: rsdummy(:,:)                                        
                             
    complex(kind=dpc)              :: e, ek                                                
    complex(kind=dpc), allocatable :: dlogdp(:)                                           
    real(kind=dp),     allocatable :: cutoff(:)                                           
    complex(kind=dpc), allocatable :: alpha(:),fz(:,:), pz(:,:), qz(:,:),sz(:,:), tmat(:), tmatll(:,:,:) 

    integer                        :: lmaxd, nrmaxd

!   Solution ASA (Old solver)
    complex(kind=dpc), allocatable :: pzsusc(:,:,:,:)
    complex(kind=dpc), allocatable :: fzsusc(:,:,:,:)
    complex(kind=dpc), allocatable :: qzsusc(:,:,:,:)
    complex(kind=dpc), allocatable :: szsusc(:,:,:,:)
!   Solution SIMULASA + SOC (New solver)
!   Big components
    complex(kind=dpc), allocatable :: pzsusc_soc(:,:,:,:)
    complex(kind=dpc), allocatable :: qzsusc_soc(:,:,:,:)
    complex(kind=dpc), allocatable :: pzsusc_socleft(:,:,:,:)
    complex(kind=dpc), allocatable :: qzsusc_socleft(:,:,:,:)
!   Small componnets (SRA)
    complex(kind=dpc), allocatable :: fzsusc_soc(:,:,:,:)
    complex(kind=dpc), allocatable :: szsusc_soc(:,:,:,:)
    complex(kind=dpc), allocatable :: fzsusc_socleft(:,:,:,:)
    complex(kind=dpc), allocatable :: szsusc_socleft(:,:,:,:)

!   valence charge&spin density for communication between KKRFLEX and module SOLVER
    real(kind=dp), allocatable       :: rho2ns_tmp(:,:,:,:)
    ! valence single-particle energy
    real(kind=dp), allocatable       :: espv_tmp(:,:,:)
    

  end type inpsusc_type

end module type_inpsusc
