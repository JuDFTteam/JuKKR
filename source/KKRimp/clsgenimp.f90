!-------------------------------------------------------------------------------
!> Summary: Generate the cluster around each atom
!> Author: 
!> Deprecated: True ! This needs to be set to True for deprecated subroutines
!-------------------------------------------------------------------------------
module mod_clsgenimp
contains
!-------------------------------------------------------------------------------
!> Summary: Generate the cluster around each atom
!> Author: 
!> Category: KKRimp, geometry, deprecated
!> Deprecated: True ! This needs to be set to True for deprecated subroutines
!-------------------------------------------------------------------------------
      subroutine clsgen(natom,rr,nr,rmt,tbcluster,lmaxatom,lmaxd) !,naclsd,nclsd)
      use configparams, only: naclsd, nclsd,rcut
      use nrtype
      use mod_dsort
      use mod_clustcomp
      use type_tbcluster
      implicit none
!interface variables
  integer,intent(in)             ::    natom                  !! number of impurity atoms
  real(kind=dp),intent(in)       ::    rr(3,nr)               !! real space impurity and surrounding atoms
  integer,intent(in)             ::    nr                     !! number of impurity and surrounding atoms rr
  real(kind=dp),intent(in)       ::    rmt(nr)                !! muffin tin radius of atoms
  type(tbcluster_type)           ::    tbcluster
  integer,intent(in)             ::    lmaxatom(nr)
  integer,intent(in)             ::    lmaxd
!local variables
  integer                        ::    ia
  integer                        ::    i,n1,pos, ib,ii,     &
                                       na,number,n,         &
                                       in,iatom,icu,icluster !,nimpcellmax
  integer                        ::    tmpatom(naclsd), &
                                       isort(naclsd)
  integer                        ::    iatcls(nclsd)
  real(kind=dp)                  ::    epsshl,r2, &
                                       rcls1(3,naclsd), &
                                       rg(3,naclsd),tmp(3),rsort(naclsd)
  real(kind=dp)                  ::    rcut2
  integer                        ::    tb_index2iatom(naclsd*(lmaxd+1)**2)
  integer                        ::    tb_iatom2index(naclsd)
  integer                        ::    tb_nclsgf

  real(kind=dp)                  ::    rcls_sorted(3,naclsd) !(:,:,:natom)       !! real space position of atom in cluster
  integer                        ::    atom_sorted(naclsd,nr) !(naclsd,nr)       !! ??????????
  real(kind=dp)                  ::    cluster_rcls_sorted(3,naclsd,natom)

allocate(tbcluster%rcls(3,naclsd,natom), tbcluster%cls(natom), &
         tbcluster%nacls(nclsd),         tbcluster%atom(naclsd,nr), &
         tbcluster%naclsimp(natom),tbcluster%naclsimpmax(nclsd))

tbcluster%naclsimpmax=0

allocate(tbcluster%iatom2index(nr,nclsd), &
         tbcluster%index2iatom(naclsd*lmaxd,nclsd), &
         tbcluster%nclsgf(nclsd) )

epsshl=1.0d-4

icluster = 1
do n = 1,nclsd
  iatcls(n) = 0
  tbcluster%nacls(n) = 0
end do

tbcluster%rcls=0
rcls1=0
rcut2   = (rcut+epsshl)**2

write(*,*) '--------------------------------------------------------------------'
write(*,*) '-------------- CLSGEN: CLUSTER GENERATION --------------------------'
write(*,*) '--------------------------------------------------------------------'


write(*,*) 'rcut2=',rcut2

            
do iatom = 1,natom       ! loop in all atoms or layers

  tbcluster%cls(iatom) = 0   
  number = 0             ! counter for atoms in cluster
  tb_iatom2index=0
  tb_iatom2index(1)=1
  do n=1,nr    ! loop in all bravais vectors    
    do i=1,3
      tmp(i) = rr(i,n) - rr(i,iatom)
    end do
    r2   =  tmp(1)**2+tmp(2)**2+tmp(3)**2 
    if ( (r2.le.rcut2) )  then
      number = number + 1
      if (number+1.gt.naclsd) then 
        WRITE (6,*)  'ERROR: Dimension NACLSD in inc.cls too small', &
                      number, naclsd
        stop '   < CLSGEN99 >'
      end if
      if (n<=natom) then 
        tbcluster%naclsimp(iatom)=number
      end if
      if (number==1) then
        tb_iatom2index(1)=1
      else
        tb_iatom2index(number)=tb_iatom2index(number-1)+(lmaxatom(tbcluster%atom(number-1,iatom))+1)**2
      end if
      tb_index2iatom(tb_iatom2index(number):tb_iatom2index(number)+(lmaxatom(n)+1)**2-1)=n
      tb_nclsgf=tb_iatom2index(number)+(lmaxatom(n)+1)**2-1
      tbcluster%atom(number,iatom) = n !iclsatom ! store the atom in elem cell
      rcls1(:,number) = tmp(:)
    end if
  end do


!--------------------------------------------------------------------
!     the atom iatom has it's cluster first 
!     sort the atoms of the cluster in increasing order. first by distance
!     then by z then by y and save them temporarily
!--------------------------------------------------------------------
  do ia=1,number
    rsort(ia) = sqrt(rcls1(1,ia)**2+      &
                     rcls1(2,ia)**2+      &
                     rcls1(3,ia)**2)      
    rsort(ia) = 100000000.d0*rsort(ia)+   &
                10000.d0*rcls1(3,ia)+   &
                10.d0*rcls1(2,ia)+       &
                0.1d0*rcls1(1,ia) 
  end do      
  isort=0
  call dsort(rsort,isort,number,pos)!,naclsd)
  do ia =1,number
    ib = isort(ia)
    rcls_sorted(:,ia) = rcls1(:,ib)
    atom_sorted(ia,iatom) = tbcluster%atom(ib,iatom)
  end do

!--------------------------------------------------------------------
!-- now sort the number of atoms just for the impurity atoms
!-- and not for the surrounding atoms
!-- this step is important !
!--------------------------------------------------------------------
  call dsort(rsort,isort,tbcluster%naclsimp(iatom),pos)!,naclsd)
  do ia=1,tbcluster%naclsimp(iatom)       
    do i=1,3
      rg(i,ia)    = rcls1(i,ia)
    end do
    tmpatom(ia) = tbcluster%atom(ia,iatom)
  end do    
  do ia =1,tbcluster%naclsimp(iatom)
    ib = isort(ia)
    do i=1,3
      rcls1(i,ia) = rg(i,ib)
    end do
    tbcluster%atom(ia,iatom) = tmpatom(ib)
  end do
!--------------------------------------------------------------------
!-- find out if this cluster already exists
!--------------------------------------------------------------------
  do icu = 1,icluster-1
    n1 = tbcluster%nacls(icu)
    if( clustcomp(cluster_rcls_sorted,rmt,atom_sorted,iatcls(icu),icu,n1,rcls_sorted, &
        number,iatom,naclsd,nr,lmaxatom)) tbcluster%cls(iatom) = icu
  end do
  if (tbcluster%cls(iatom).eq.0) then
!--------------------------------------------------------------------
!-- if not create a new cluster
!--------------------------------------------------------------------
    if (icluster.gt.nclsd) then
      write(6,*) 'Please, increase the parameter NCLSD in', &
                ' inc.cls to a value greater equal ',icluster,' .'
      STOP 'Dimension error.' 
    end if
    tbcluster%cls(iatom) = icluster
    tbcluster%nacls(icluster) = number
    iatcls(icluster) = iatom
    write(*,'("-----------  CLUSTER ",I4,"  ------------------")') icluster
    do in = 1,number
      do ii=1,3
        cluster_rcls_sorted(ii,in,icluster)=rcls_sorted(ii,in)
        tbcluster%rcls(ii,in,icluster) = rcls1(ii,in)
      end do
      write(6,'(2I5,4F8.4,2I)') iatom,TBCLUSTER%atom(in,iatom), &
                                (rcls1(i,in),i=1,3), &
                                sqrt(rcls1(1,in)**2+rcls1(2,in)**2+rcls1(3,in)**2), &
                                tb_iatom2index(in),tb_nclsgf
    end do   
  
    tbcluster%ncls=icluster
    tbcluster%iatom2index(:,icluster)=tb_iatom2index
    tbcluster%index2iatom(:,icluster)=tb_index2iatom
    tbcluster%nclsgf(icluster)=tb_nclsgf
    tbcluster%naclsimpmax(icluster) = tbcluster%naclsimp(iatom)
    icluster = icluster + 1
  else
!--------------------------------------------------------------------
!-- if the cluster already exist asign the cluster to the atom and
!--   write out some stuff
!--------------------------------------------------------------------
    write(*,*) icu,tbcluster%naclsimpmax(icu),tbcluster%naclsimp(iatom)
    if (tbcluster%naclsimpmax(tbcluster%cls(iatom)) < tbcluster%naclsimp(iatom)) &
        tbcluster%naclsimpmax(tbcluster%cls(iatom)) = tbcluster%naclsimp(iatom)
  end if 

end do !iatom

!--------------------------------------------------------------------
!-- write out cluster information
!--------------------------------------------------------------------
write(*,'("--------------------------------------------------")') 
do icluster=1,tbcluster%ncls
  write(*,'("CLUSTER ",I4," has a max. number of",I4," impurity atoms")') icluster,TBCLUSTER%NACLSIMPMAX(ICLUSTER)
end do
write(*,'(I4," clusters have been created")'), TBCLUSTER%NCLS 
write(*,'("-----------  CLUSTER <---> ATOMS  ----------------")') 
do iatom=1,natom
  write(6,'("   Atom",I4," has cluster",I4," with",I4," sites and ",I4," impurity atoms")') &
        IATOM, TBCLUSTER%CLS(IATOM),TBCLUSTER%NACLS(TBCLUSTER%CLS(IATOM)), &
        TBCLUSTER%NACLSIMP(IATOM)
end do 

end subroutine

end module mod_clsgenimp
