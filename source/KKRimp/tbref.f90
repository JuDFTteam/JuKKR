module mod_tbref

contains

!-------------------------------------------------------------------------------
!> Summary: Compute reference system
!> Author: 
!> Category: KKRimp, reference-system
!> Deprecated: True ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
subroutine tbref(eryd,natom, ntotatom,tbcluster,alatc,vref,rmtref,lmaxatom,lmaxd,lmgf0d)
  use nrtype, only: dp
  use type_tbcluster, only: tbcluster_type
  use type_gref, only: gref_type
  use mod_gll95, only: gll95
  use mod_calctref, only: calctref
  use mod_config, only: config_testflag
    implicit none

  complex(kind=dpc),intent(in)              :: eryd
  integer,intent(in)                        :: natom
  integer,intent(in)                        :: ntotatom
  type(tbcluster_type),intent(in)           :: tbcluster
  real(kind=dp),intent(in)                  :: alatc
  real(kind=dp),intent(in)                  :: vref
  real(kind=dp),intent(in)                  :: rmtref(ntotatom)
  integer,intent(in)                        :: lmaxatom(ntotatom)
  integer,intent(in)                        :: lmaxd
  integer,intent(in)                        :: lmgf0d
  type(gref_type),allocatable               :: gref(:)

  complex(kind=dp)                          :: trefll(lmgf0d,lmgf0d,ntotatom)
  integer                                   :: iatom,lm1,ilmgf0d,icls,lmtmp

allocate(gref(tbcluster%ncls))

! ----------------------------------------------------------------------
! calculate the t-matricies for all atoms in the TBCLUSTER
! ----------------------------------------------------------------------
trefll=0.0_dp
do iatom = 1,ntotatom
  call calctref(eryd,vref,rmtref(iatom),lmaxatom(iatom),lmtmp, &
                trefll(:,:,iatom),lmaxd+1,(lmaxd+1)**2)
end do

if (config_testflag('tfree')==1) then
  open(unit=10000,file='test_tfree')
  do iatom = 1,ntotatom
    write(10000,*) '#atom',iatom
    do lm1=1,(lmaxatom(iatom)+1)**2
      write(10000,'(50000E25.14)') trefll(lm1,:,iatom)
    end do !lm
  end do !iatom
  close(10000)
end if


if (config_testflag('gref')==1) then
  open(unit=20000,file='test_gref.dat')
end if

! ----------------------------------------------------------------------
! calculate the reference greens function Gref(n,0)
! ----------------------------------------------------------------------

do icls = 1,tbcluster%ncls
  write(*,*) 'icls',icls!             i1 = 1
  call gll95(icls,eryd,trefll,tbcluster, &
             alatc,0,gref(icls), (lmaxd+1)**2,natom,ntotatom,lmaxatom)
end do !icls

if (config_testflag('gref')==1) then
  do icls = 1,tbcluster%ncls
    do iatom=1, tbcluster%nacls(icls)
!     write(88188,*) 'rcls= ',rcls(:,iatom,icls)
      do ilmgf0d = 1,(lmaxatom(tbcluster%atom(iatom,icls))+1)**2
!         write(*,'(50000E25.14)') gref(icls)%mat( tbcluster%iatom2index(iatom,icls)-1+ &
!                                                   ilmgf0d,:)
        write(20000,'(50000E25.14)') gref(icls)%mat( tbcluster%iatom2index(iatom,icls)-1+ &
                                                ilmgf0d,:)
      end do !ilmgf0d
    end do !iatom
  end do
  close(20000)
end if

end subroutine tbref

end module mod_tbref
