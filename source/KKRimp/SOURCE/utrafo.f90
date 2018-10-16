module mod_utrafo

  private
  public :: utrafo

contains

  !-------------------------------------------------------------------------------
  !> Summary: U-transformation for impurity relaxations
  !> Author: 
  !> Category: KKRimp, geometry
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> see PhD Bauer, p. 100.
  !-------------------------------------------------------------------------------
  subroutine utrafo(ielast, ez, natom, lmaxatom, gmatbulk, rimpshift)
    use mod_preconditioning, only: preconditioning_writegreenfn, preconditioning_readgreenfn
    use type_gmatbulk, only: gmatbulk_type
    use mod_config, only: config_testflag
    use mod_umatrix, only: umatrix
    use nrtype, only: wlength
    implicit none
    !interface
    integer         :: ielast
    double complex  :: ez(ielast)
    integer         :: natom
    integer         :: lmaxatom(natom)
    type(gmatbulk_type) :: gmatbulk
    !local
    double precision         :: rimpshift(3,natom)

    double complex,allocatable :: umat(:,:)
    double complex,allocatable :: gmatvoid(:,:)
    double complex,allocatable :: gmatcut(:,:)
    double complex,allocatable :: gmatnew(:,:)
    integer                    :: my_rank
    integer                    :: nspinhost
    integer                    :: nlmhostnew
    integer                    :: iatom,ie,ispin,ilm
    integer                    :: lmmaxhost,lmaxhost
    integer                    :: recl1,recl2
    integer                    :: lmstart,lmstop,lmstart2,lmstop2
    integer                    :: dimgmathost
    double complex,parameter   :: CONE=(1.0D0,0.0D0),CZERO=(0.0D0,0.0D0)

    ! future parallelization
    my_rank=0

    ! read gmatbulk settings from type definition
    natom=gmatbulk%natom
    lmaxhost=gmatbulk%lmax
    lmmaxhost=gmatbulk%lmmax
    nspinhost=gmatbulk%nspin
    dimgmathost=gmatbulk%hostdim

    write(1337,*) 'Translating void Greensfunction with properties:'
    write(1337,*) 'natom ',natom
    write(1337,*) 'lmmaxhost ',lmmaxhost
    write(1337,*) 'nspinhost ',nspinhost

    !calculate the Greens function dimension
    nlmhostnew=0
    do iatom=1,natom
      if ((lmaxatom(iatom)+1)**2>lmmaxhost) then
        write(*,*) 'lmax value of atom ',iatom,' is greater than the lmax value of the host'
        stop
      end if
      nlmhostnew=nlmhostnew+(lmaxatom(iatom)+1)**2  
    end do

    !open files
    recl1=wlength*4*natom*lmmaxhost*natom*lmmaxhost
    recl2=wlength*4*nlmhostnew*nlmhostnew
    open (88,access='direct',recl=recl1,file='kkrflex_greenvoid',form='unformatted')
    open (89,access='direct',recl=recl2,file='kkrflex_greennew',form='unformatted')
    if(config_testflag('gmat_plain')) open (8989,file='kkrflex_greennew.txt')
    if (my_rank==0) write(89,rec=1) natom,nlmhostnew,nspinhost


    allocate( gmatnew(LMMAXhost*natom,LMMAXhost*natom))
    allocate( gmatcut(nlmhostnew,nlmhostnew) )
    allocate( umat(LMMAXhost*natom,LMMAXhost*natom))

    if (config_testflag('write_umat')) open(unit=7654323, file='test_umat')
    do ie=1,ielast
      ! ***************************************************
      !  calculate U(s=s0)
      ! ***************************************************
      call umatrix(natom,lmaxhost,umat,dimgmathost,rimpshift,ez(ie))

      if (config_testflag('write_umat')) then
        do ilm=1,dimgmathost
          write(7654323,'(50000F)') umat(:,ilm)
        end do
      end if

      do ispin=1,nspinhost
        write(1337,*) 'proc = ',my_rank,' IE = ',ie,' ispin= ',ispin
        call preconditioning_readgreenfn(ie,ispin,ielast,lmmaxhost,natom,gmatvoid,'doubleprecision')

        ! ***************************************************
        !  calculate G(s=s0) = U(s=so)**T G(s=0) * U(s=s0)
        ! ***************************************************
        CALL ZGEMM('N','N',dimgmathost,dimgmathost,dimgmathost,CONE, &
        umat,dimgmathost,gmatvoid,dimgmathost,CZERO,gmatnew,dimgmathost)

        CALL ZGEMM('N','T',dimgmathost,dimgmathost,dimgmathost,CONE, &
        gmatnew,dimgmathost,umat,dimgmathost,CZERO,gmatvoid,dimgmathost)


        ! ***************************************************
        !  cut out the Greens function elements bigger 
        !  then lmaxatom
        ! ***************************************************
        gmatnew=(0.0D0,0.0D0)
        lmstart2=1 !-(lmaxatom(1)+1)**2
        do iatom=1,natom
          lmstart=(iatom-1)*lmmaxhost+1
          lmstop = lmstart+(lmaxatom(iatom)+1)**2-1
          lmstop2  = lmstart2-1 + (lmaxatom(iatom)+1)**2

          gmatnew(:,lmstart2:lmstop2)=gmatvoid(:,lmstart:lmstop)

          lmstart2 = lmstart2   + (lmaxatom(iatom)+1)**2
        end do


        lmstart2=1 !-(lmaxatom(1)+1)**2
        do iatom=1,natom
          lmstart=(iatom-1)*lmmaxhost+1
          lmstop = lmstart+(lmaxatom(iatom)+1)**2-1
          lmstop2  = lmstart2-1 + (lmaxatom(iatom)+1)**2

          gmatcut(lmstart2:lmstop2,:)=gmatnew(lmstart:lmstop,:)
          lmstart2 = lmstart2   + (lmaxatom(iatom)+1)**2
        end do

        call preconditioning_writegreenfn(ie,ispin,ielast,gmatcut,nlmhostnew)

      end do
    end do

    close(88)
    close(89)
    if(config_testflag('gmat_plain')) close(8989)

  end subroutine utrafo

end module mod_utrafo

