module mod_dysonvirtatom
contains

subroutine dysonvirtatom(natom,ntotatom,lmsizehost,gmathost,tmat,killatom, &
                         isvatom,lmaxatom,gmatimp,nlmhostnew,lattice_relax,NSOC)
use mod_mathtools
use mod_config, only: config_testflag
  implicit none
!interface
integer,intent( in)              :: natom      ! Number of atoms in impurity cluster, excl. "killed atoms"
integer,intent( in)              :: ntotatom   ! Number of host sites where GF is read in, incl. "killed atoms"
integer,intent( in)              :: lmsizehost ! Host lm-max * NSOC (see below for NSOC)
double complex,intent( in)       :: gmathost(ntotatom*lmsizehost, ntotatom*lmsizehost)
double complex,intent( in)       :: tmat(lmsizehost,lmsizehost,ntotatom)
integer,intent( in)              :: killatom(ntotatom)
integer,intent( in)              :: isvatom(ntotatom)
integer,intent( in)              :: lmaxatom(ntotatom) ! lmax to be used for specific atom
integer,intent( in)              :: lattice_relax
integer,intent( in)              :: NSOC ! =2 if host-GF is 2x2 in spin space; =1 otherwise
double complex,intent(out),allocatable  :: gmatimp(:,:)

!local
integer                          :: iatom
integer                          :: nlmhostnew,lmstart,lmstart2,lmstop,lmstop2
integer                          :: ilm,iat_imp1,iat_imp2,iat_host1,iat_host2,isoc1,isoc2
integer                          :: lmmax_imp1,lmmax_imp2
integer                          :: istart_host1,istop_host1,istart_host2,istop_host2 
integer                          :: istart_imp1,istart_imp2,istop_imp1,istop_imp2
integer                          :: ilms_start_imp1,ilms_start_imp2

double complex,allocatable                   :: gmatout(:,:)
double complex,allocatable                   :: tmatbig(:,:)
double complex,allocatable                   :: gmatnew(:,:)

allocate(gmatout(ntotatom*lmsizehost, ntotatom*lmsizehost), &
                 tmatbig(lmsizehost*ntotatom,lmsizehost*ntotatom), &
                 gmatnew(lmsizehost*ntotatom,lmsizehost*ntotatom))


tmatbig=(0.0D0,0.0D0)
do iatom=1,ntotatom
  lmstart=lmsizehost*(iatom-1)+1


    if (isvatom(iatom)==0) then !killatom(iatom)==1
      tmatbig(lmstart:lmstart-1+lmsizehost,lmstart:lmstart-1+lmsizehost) = tmat(:,:,iatom)
    end if

! ##############################################################
! better version as the one above maybe replace in the future
! ##############################################################
if(.false.) then
  if (isvatom(iatom)==0) then !killatom(iatom)==1

    call matmat_dcdc(gmathost(:,lmstart:lmstart-1+lmsizehost), tmat(:,:,iatom),gmatnew(:,lmstart:lmstart-1+lmsizehost))
  else
    gmatnew(:,lmstart:lmstart-1+lmsizehost) = (0.0D0,0.0D0)
  end if
end if

 end do

call matmat_dcdc(gmathost,tmatbig,gmatnew)


! *******   calculate        *********
! gmatnew= gmatnew + 1
! ************************************
do ilm=1,ntotatom*lmsizehost
  gmatnew(ilm,ilm)=gmatnew(ilm,ilm)+(1.0D0,0.0D0)
end do

! *******   calculate        *********
! gmatout = inverse(gmatnew)*g
! ************************************

gmatout=gmathost


if (.not. config_testflag('nopredyson')) then
  call linearsolve_dc(gmatnew,gmatout)
end if !(.not. config_testopt('nodyson')) then





allocate( gmatimp(nlmhostnew,nlmhostnew) )



gmatnew=(0.0D0,0.0D0)
lmstart2=1 !-(lmaxatom(1)+1)**2

if (lattice_relax==0) then

   ! Starting point for atom 1 of impurity-block
   ilms_start_imp1 = 1
   iat_imp1 = 1

   ! Loop 1 over host block
   do iat_host1 = 1,ntotatom

      ! Exclude killed atoms
      if (killatom(iat_host1)/=1) then 
         lmmax_imp1 = (lmaxatom(iat_imp1)+1)**2

         ! Starting point for atom 2 of impurity-block
         ilms_start_imp2 = 1
         iat_imp2 = 1
         ! Loop 2 over host block
         do iat_host2 = 1,ntotatom

            ! Exclude killed atoms
            if (killatom(iat_host2)/=1) then
               lmmax_imp2 = (lmaxatom(iat_imp2)+1)**2

               !-------------------------------------------------------------
               ! Loop over spins (both spin directions in case of NSOC=2)
               do isoc1 = 1,NSOC                   ! NSOC=2 if host GF is 2x2 matrix in spin space 
               do isoc2 = 1,NSOC                   ! NSOC=1 otherwise

                  ! Select block from the host GF matrix
                  ! Remember, lmsizehost includes spin index in case of SOC
                  istart_host1 = (iat_host1 - 1)*lmsizehost + (isoc1 - 1)*lmsizehost/NSOC + 1
                  istop_host1 = istart_host1 + lmmax_imp1 - 1   
                  istart_host2 = (iat_host2 - 1)*lmsizehost + (isoc2 - 1)*lmsizehost/NSOC + 1
                  istop_host2 = istart_host2 +  lmmax_imp2 - 1   

                  ! Select block from the impurity GF matrix
                  ! Remember, lmmax_imp does _not_ include spin index
                  istart_imp1 = ilms_start_imp1 + (isoc1 - 1) * lmmax_imp1
                  istop_imp1 = istart_imp1 + lmmax_imp1 - 1
                  istart_imp2 = ilms_start_imp2 + (isoc2 - 1) * lmmax_imp2
                  istop_imp2 = istart_imp2 + lmmax_imp2 - 1


                  ! Write out the copied blocks for test case
                  !write(*,*) 'dysonviratom:--------------------------------'
                  !write(*,*) 'dysonviratom: isoc1 isoc2',isoc1,isoc2
                  !write(*,*) 'dysonviratom:--------------------------------'
                  !write(*,*) 'dysonviratom: atoms',iat_host1,iat_host2,iat_imp1,iat_imp2
                  !write(*,*) 'dysonviratom: host',istart_host1 ,istop_host1, istart_host2 ,istop_host2
                  !write(*,*) 'dysonviratom:  imp',istart_imp1 ,istop_imp1, istart_imp2 ,istop_imp2

                  ! Copy lmmax x lmmax block
                  gmatimp(istart_imp1 :istop_imp1, istart_imp2 :istop_imp2 ) =        &
                & gmatout(istart_host1:istop_host1,istart_host2:istop_host2)

               enddo ! isoc2 = 1,NSOC 
               enddo ! isoc1 = 1,NSOC 
               !-------------------------------------------------------------


            ! Shift starting point of next block for impurity 2 to next atom
            ilms_start_imp2 = ilms_start_imp2 + lmmax_imp2 * NSOC
            endif ! (killatom(iat_host2)/=1)
            iat_imp2 = iat_imp2 + 1 ! This is outside the if-block because the lmax(impurity) from 
                                    ! the kkrflex_atominfo file is read in even for the killed atoms.

         enddo ! iat_host2 = 1,ntotatom

      ! Shift starting point of next block for impurity 1 to next atom
      ilms_start_imp1 = ilms_start_imp1 + lmmax_imp1 * NSOC
      endif ! (killatom(iat_host1)/=1)
      iat_imp1 = iat_imp1 + 1


   enddo ! iat_host1 = 1,ntotatom


else  ! (lattice_relax==0)

  do iatom=1,ntotatom
    lmstart=(iatom-1)*lmsizehost+1
    lmstop = lmstart+lmsizehost-1
    if (killatom(iatom)/=1) then
      lmstop2  = lmstart2-1 + lmsizehost
      gmatnew(:,lmstart2:lmstop2)=gmatout(:,lmstart:lmstop)
      lmstart2 = lmstart2   + lmsizehost
    end if
  end do

  lmstart2=1 !-(lmaxatom(1)+1)**2
  do iatom=1,ntotatom
    lmstart=(iatom-1)*lmsizehost+1
    lmstop = lmstart+lmsizehost-1
    if (killatom(iatom)/=1) then
      lmstop2  = lmstart2-1 + lmsizehost
      gmatimp(lmstart2:lmstop2,:)=gmatnew(lmstart:lmstop,:)
      lmstart2 = lmstart2   + lmsizehost
    end if
  end do

end if  ! (lattice_relax==0)



end subroutine dysonvirtatom





      subroutine matmat_dcdc(mat1,mat2,matout)
      implicit none
      double complex, intent(in) :: mat1(:,:),mat2(:,:)
      double complex             :: matout(size(mat1,1),size(mat2,2))
      integer                :: n1,n,n2
      n1 = size(mat1,1)
      n  = size(mat1,2)
      n2 = size(mat2,2)
      if(size(mat2,1).ne.n) stop 'matmat_zmzm: dimensions of matrices are inconsistent.'
      call zgemm('N','N',n1,n2,n,(1d0,0d0),mat1,n1,mat2,n,(0d0,0d0),matout,n1)
      end subroutine matmat_dcdc



end module !mod_dysonvirtatom
