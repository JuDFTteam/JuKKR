!------------------------------------------------------------------------------------
!> Summary: This module is used to solve the Dyson equation in free space 
!> Author: People who wrote it
!>  solve the system of linear equations:
!>   math    :  (1 - G_free * T)*G_tb          = G_free
!> ----------------------------------------------------------------------
!>     GREF= 1 - G_free * T
!>    GREFOUT%MAT = G_free
!>     result : GREFOUT%MAT = tb reference Greensfunction G_ref(n,0) 
!>
!> One can write Latex comments like this \(i\hbar\frac{\partial \psi}{\partial t}=-\mathcal{H}\psi\)
!> or add labeled equations using the standard latex way
!> \begin{equation}
!> \mathbf{A} \mathbf{v}= \eta\mathbf{v}
!> \end{equation}
!> **FORd** also accepts markdown style so you can _write with style_ 
!> 
!> **IMPORTANT**
!> The JM-KKR follows the coding conventions noted in this example, one that is
!> not obvious is that **each level of indentation consists of two spaces**. Please keep this:
!> _do it for the children_.
!> So please keep the conventions.
! These are special boxes for ford, notice that this comment does not appear in the html file.
! These boxes contain important information and should be added when necessary. ALWAYS remember to close the box
! BEFORE oppening a new one or they will be nested.
!------------------------------------------------------------------------------------

MODULE mod_GLL95

CONTAINS

   !-------------------------------------------------------------------------------
   !> Summary:  solution of the DYSON equation for a cluster of potentials
   !>           (TMATLL) centered at positions RATOM in free space
   !> Author: Who wrote this subroutine
   !> Category: Green-function, reference system 
   !> Deprecated: True ! This needs to be set to True for deprecated subroutines
   !> A More detailed explanation with the math, concepts, etc necessary to understand the routine
   !-------------------------------------------------------------------------------


SUBROUTINE GLL95(ICLS,ERYD,TMATLL,TBCLUSTER, &
                 ALAT,OUT,GREFOUT,LMMAXD,NIMPATOM,NTOTATOM,LMAXATOM)
! ----------------------------------------------------------------------
!     solution of the DYSON equation for a cluster of potentials
!     (TMATLL) centered at positions RATOM in free space,
! ----------------------------------------------------------------------
  USE type_tbcluster
  USE type_gref
  USE mod_gauntharmonics, only: gauntcoeff
!   USE mod_grefsy
  USE mod_mathtools, only: linearsolve_dc
  USE mod_gfree
  USE mod_config, only: config_testflag
    IMPLICIT NONE

  INTEGER, intent(in)                    ::    ICLS
  COMPLEX(kind=DPC),intent(in)           ::    ERYD
  COMPLEX(kind=DPC),intent(in)           ::    TMATLL(LMMAXD,LMMAXD,NTOTATOM)
  TYPE(TBCLUSTER_TYPE),intent(in)        ::    TBCLUSTER
  REAL(kind=DP),intent(in)               ::    ALAT
  INTEGER,intent(in)                     ::    OUT
  type(GREF_TYPE)                        ::    GREFOUT
  INTEGER,intent(in)                     ::    LMMAXD
  INTEGER,intent(in)                     ::    NIMPATOM
  INTEGER,intent(in)                     ::    NTOTATOM
  INTEGER,intent(in)                     ::    LMAXATOM(NTOTATOM)
  
  INTEGER                                ::    GFDIM(2)
  INTEGER                                ::    LMAXFREE
  REAL(kind=DP)                          ::    RDIFF(3)
  COMPLEX(kind=DPC),ALLOCATABLE          ::    GREF(:,:),GLL(:,:),GTREF(:,:)
  INTEGER                                ::    LM1,LM2,NLM1,NLM2,IATOM1,IATOM2,I,GFNLM,LMMAXTEMP
  INTEGER                                ::    IERROR
  COMPLEX(KIND=DPC),PARAMETER            ::    CZERO = (0.D0,0.D0),&
                                               CONE  = (1.D0,0.D0)


gfnlm=   tbcluster%nclsgf(icls)
gfdim=(/ tbcluster%nclsgf(icls), (lmaxatom(tbcluster%atom(1,icls))+1)**2 /)

allocate( grefout%mat(gfdim(1), gfdim(2)) )

allocate (      gref( gfdim(1), gfdim(1) ), &
                gll( lmmaxd,lmmaxd ), & 
                gtref( gfdim(1),gfdim(2) ), &
                stat=ierror)

if ( ierror.ne.0 ) then
  write(6,'(6x,"error: failed to allocate array(s) :",a,/)') 'gref/gll/gtref'
  stop '           < gll95 > '
end if

! ----------------------------------------------------------------------
!     set up the free greens functions for the whole cluster ICLS
!     G_free(n1,lm1,n2,lm2) and save the result in GREF
! ----------------------------------------------------------------------
do iatom1 = 1,tbcluster%nacls(icls)
  do iatom2 = 1,tbcluster%nacls(icls)
    rdiff = - (tbcluster%rcls(:,iatom1,icls)-tbcluster%rcls(:,iatom2,icls))*alat
    lmaxfree=max(lmaxatom(tbcluster%atom(iatom1,icls)),lmaxatom(tbcluster%atom(iatom2,icls)))
    if (iatom1.ne.iatom2) then
      call gfree(rdiff,eryd,gll, &
                 gauntcoeff(lmaxfree)%cleb(:,2),gauntcoeff(lmaxfree)%icleb,gauntcoeff(lmaxfree)%loflm,gauntcoeff(lmaxfree)%iend,gauntcoeff(lmaxfree)%ncleb, &
                 lmaxfree,(lmaxfree+1)**2,lmaxfree*2+1,(lmaxfree*2+1)**2 )
    else
      gll = czero
    end if
    lm1  = (lmaxatom(tbcluster%atom(iatom1,icls))+1)**2
    lm2  = (lmaxatom(tbcluster%atom(iatom2,icls))+1)**2

    nlm1 = tbcluster%iatom2index(iatom1,icls)
    nlm2 = tbcluster%iatom2index(iatom2,icls)

    gref(nlm1:nlm1+lm1-1,nlm2:nlm2+lm1-1) = gll(1:lm1,1:lm2)
  end do !iatom2
end do !iatom1

! ----------------------------------------------------------------------
!     save the  impurity atom part of the free greens for later use
! ----------------------------------------------------------------------
grefout%mat=GREF(:,1:GFDIM(2))

if (config_testflag('gfree')==1) then
  open(unit=10000,file='test_gfree')
  write(10000,*) '#icls',icls
  do lm1=1,tbcluster%nclsgf(icls)
    write(10000,'(50000f)') gref(lm1,:)
  end do
  close(10000)
end if

! ----------------------------------------------------------------------
!     do the multiplication : - G_free * T (remember the minus sign!)
!     and save the result in the array GREF
! ----------------------------------------------------------------------
do iatom2 = 1,tbcluster%nacls(icls)
  nlm2 = tbcluster%iatom2index(iatom2,icls)
  lmmaxtemp=(lmaxatom(tbcluster%atom(iatom2,icls))+1)**2
  do lm1=1,lmmaxtemp
    gref(:,nlm2+lm1-1)=-gref(:,nlm2+lm1-1)*tmatll(lm1,lm1,tbcluster%atom(iatom2,icls))
  end do !lm1
end do

if (config_testflag('gtfree')==1) then
  open(unit=10000,file='test_gtfree')
  write(10000,*) '#icls',icls
  do lm1=1,tbcluster%nclsgf(icls)
    write(10000,'(50000f)') gref(lm1,:)
  end do
  close(10000)
end if

! ----------------------------------------------------------------------
!     do the multiplication : 1 - G_free * T
!     and save the result in the array GREF
! ----------------------------------------------------------------------
do nlm1 = 1,tbcluster%nclsgf(icls)
        gref(nlm1,nlm1) = cone + gref(nlm1,nlm1) ! gtmat= 1 - g * t
end do

! ----------------------------------------------------------------------
!     solve the system of linear equations:
!   math    :  (1 - G_free * T)*G_tb          = G_free
! ----------------------------------------------------------------------
!     GREF= 1 - G_free * T
!     GREFOUT%MAT = G_free
!     result : GREFOUT%MAT = tb reference Greensfunction G_ref(n,0) 
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

call linearsolve_dc(gref,grefout%mat)

deallocate (gref,gll,gtref)

end subroutine gll95

end module mod_gll95
