module mod_umatrix

contains

  !-------------------------------------------------------------------------------
  !> Summary: Compute U matrix for U tranfsormation
  !> Author: 
  !> Category: KKRimp, geometry
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> see PhD Bauer p. 100 f.
  !-------------------------------------------------------------------------------
  subroutine umatrix(natom,lmaxhost,umat,dimgmathost,svec,eryd)
    use mod_jmtrx, only: jmtrx
    implicit none
    integer :: natom
    integer :: lmaxhost
    integer :: dimgmathost
    double complex :: umat(dimgmathost,dimgmathost)
    double precision :: svec(3,natom)
    double complex :: eryd
    !local
    integer       :: iatom,ilm,lmmax
    logical       :: lcall
    double complex :: umatrix_small( (lmaxhost+1)**2,(lmaxhost+1)**2 )

    lcall=.false.
    lmmax=(lmaxhost+1)**2
    umat=(0.0d0,0.0d0)
    do iatom=1,natom
      umatrix_small=(0.0d0,0.0d0)
      call  jmtrx(svec(1,iatom),svec(2,iatom),svec(3,iatom),eryd,lmaxhost,umatrix_small,lcall)
      do ilm=1,lmmax
        umat((iatom-1)*lmmax+1:(iatom)*lmmax,(iatom-1)*lmmax+ilm)=umatrix_small(:,ilm)
      end do
    end do

  end subroutine umatrix

end module mod_umatrix