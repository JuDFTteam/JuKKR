module mod_vinters2010
contains

!-------------------------------------------------------------------------------
!> Summary: 
!> Author: 
!> Category: KKRimp, 
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
!>-----------------------------------------------------------------------
!>     calculate the intercell-potentials and add these to the poten-
!>     tial v  (in the spin-polarized case for each spin-direction
!>     the intercell-potential is the same . )
!>     it uses the structure dependent matrices amat and bmat which
!>     are calculate once in the subroutine amn .
!>     the charge-moments are calculated in the subroutine vintra2 ,
!>     therefore vintra2 has to be called first .
!>     the intercell-potential is expanded into spherical harmonics .
!>     the lm-term of the intercell-potential v of the representive
!>     atom i is given by
!>     
!>     v(r,lm,i) =  (-r)**l * {amat(i1,i2,lm,l'm')*cmom(i2,l'm')
!>     +bmat(i1,i2,lm)*z(i2)}
!>     
!>     summed over i2 (all shells) and l'm' .    (i1=i-natref)
!>     (see notes by b.drittler)
!>     
!>     in case of shape correction the madelung potential of the host
!>     is taken into account . in all other case the madelung poten-
!>     tial of the host is set to be zero .
!>     as actual values for z and cmom the differences between the
!>     values of the  representive atoms and those of the references
!>     are used .
!>     
!>     attention : the first index of cmom (moment of the charge
!>     density - calculated in vintr2) and of z (nuclear
!>     charge of the atoms in the shell) is in the program
!>     different defined , there one has to use :
!>     cmom(natref+i2,l'm')  and  z(natref+i2)
!>     
!>     b.drittler   june 1987
!>--------------------------------------------------------------------
!> In the case of impurities on surfaces, the intercell potential of the
!> reference system is read in from the fxdr file 'intercell_ref'; which
!> is calculated in the surface program. 
!>                  july 1998
!>-----------------------------------------------------------------------
subroutine vinters2010(natom,nspin,lmaxd,lmaxatom,cell,ratom,cmom,alat,ins,nrmaxd,vpot_out,intercell_ach,zat,lattice_relax,rimpshift)
use nrtype, only: DP
use type_cell, only: cell_type
use mod_shftvout, only: SHFTVOUT
use mod_amn2010, only: amn2010
use mod_gauntharmonics, only: gauntcoeff
use mod_config, only: config_testflag
implicit none
!interface
integer                    :: natom
integer                    :: nspin
integer                    :: lmaxd
integer                    :: lmaxatom(natom)
type(cell_type)            :: cell(natom)
real(kind=DP)              :: ratom(3,natom)
real(kind=DP)              :: cmom((2*lmaxd+1)**2,natom) !(lmpotd,ntotatom)
real(kind=DP)              :: zat(natom) !(lmpotd,ntotatom),achnew(lmpotd,ntotatom)
real(kind=DP)              :: alat
integer                    :: ins
integer                    :: nrmaxd
real*8                     :: vpot_out(nrmaxd,(2*lmaxd+1)**2,nspin,natom) !thetas(iri,nfund,*),
integer                    :: lattice_relax
real(kind=DP)              :: rimpshift(3,natom)
!local
integer                    :: ispin,iatom,jatom,lval,mval,lm,lm2,ir,irs1
integer                    :: lmpotd
real(kind=DP),allocatable  :: intercell_ach(:,:) !(lmpotd,ntotatom),achnew(lmpotd,ntotatom)
real(kind=DP)              :: intercell((2*lmaxd+1)**2) !(lmpotd,ntotatom),achnew(lmpotd,ntotatom)
real(kind=DP)              :: intercell_temp((2*lmaxd+1)**2) !(lmpotd,ntotatom),achnew(lmpotd,ntotatom)
real(kind=DP)              :: ac
real(kind=DP),allocatable  :: amat(:,:,:)
real(kind=DP),allocatable  :: bmat(:,:)

if (.not. allocated(gauntcoeff)) stop '[preconditioning_intercell] gauntcoeff not allocated'
lmpotd=(2*lmaxd+1)**2
allocate( amat(natom,lmpotd,lmpotd), &
          bmat(natom,lmpotd ) ) 
amat=0.0D0
bmat=0.0D0

! allocate( ach( (lmaxd+1)**2 ,natom) )

do iatom=1,natom

   if (lattice_relax==0) then
     intercell=intercell_ach(1:(2*lmaxd+1)**2,iatom)
   else
!      write(*,*) intet
     intercell_temp=intercell_ach(1:(2*lmaxd+1)**2,iatom)
     call SHFTVOUT(intercell_temp,intercell,rimpshift(:,iatom), &
                                 2*LMAXD,gauntcoeff(lmaxd)%WG,gauntcoeff(lmaxd)%YRG,(LMAXD+1)**2,4*LMAXD,2*LMAXD,(4*LMAXD+1)**2,gauntcoeff(lmaxd)%NCLEB)
!                                LMAXD,WG,YRG,LMMAXD,LASSLD,LPOTD,LM3D,NCLEB)
   end if

   call amn2010(iatom,amat,bmat,alat,gauntcoeff(lmaxd),2*lmaxd,natom,ratom)

!          DO jatom= 1,natom
!            WRITE(10000,*) AMAT(jatom,:,:)  !Bauer
!            WRITE(10001,*) BMAT(jatom,:)  !Bauer
!          END DO
!          STOP


!    write(*,*) 'amn2010',iatom,amat,bmat,alat,2*lmaxd,natom,ratom
!    write(*,'(A,5000F)') 'ac2',cmom

!   write(*,*) 'leave amn2010'
!   write(*,*) iatom,amat,bmat,alat,lpothost,ntotatom,ratom
  do lval = 0,2*lmaxatom(iatom) !*lmaxd
    do mval = -lval,lval
      lm = lval*lval + lval + mval + 1
      ac = 0.0d0
      if (.not. config_testflag('hostintercell')) then
        do jatom = 1,natom
          do lm2 = 1,(2*lmaxatom(jatom)+1)**2
            ac = ac + amat(jatom,lm,lm2)*cmom(lm2,jatom)
  !       write(*,*) 'ac5',amat(jatom,lm,lm2),cmom(lm2,jatom)
!                       write(70001,*) amat(jatom,lm,lm2),cmom(lm2,jatom)
          end do !lm2
          ac = ac + bmat(jatom,lm)*zat(jatom)
!                       write(70003,*) bmat(jatom,lm),zat(jatom)
  !       write(*,*) 'ac5',bmat(jatom,lm),zat(jatom)
        end do
  !       write(*,*) 'ac1',intercell_ach(lm,iatom)
  !       write(*,*) 'ac2',ac
      end if !(.not. config_testflag('hostintercell')) then
!       write(60000,*) ac,intercell_ach(lm,iatom)
      ac = ac + intercell(lm)
!       write(*,*) 'ac3',ac
!       write(50000,*) ac,intercell_ach(lm,iatom)
      IF (ins.NE.0) THEN
        IRS1 = cell(iatom)%NRCUT(cell(iatom)%NPAN)
      ELSE
        IRS1 = CELL(IATOM)%NRMAX
      END IF

      DO ISPIN = 1,NSPIN
!             IPOT = NSPIN* (IATYP-1) + ISPIN

          ! code has problems evaluating 0**0 (0 to the power of zero)
          IF (LVAL.EQ.0) THEN 
            VPOT_OUT(1,1,ISPIN,IATOM) = VPOT_OUT(1,1,ISPIN,IATOM) + AC
          END IF
          DO IR = 2,IRS1
            VPOT_OUT(IR,LM,ISPIN,IATOM) = VPOT_OUT(IR,LM,ISPIN,IATOM) + (-CELL(IATOM)%RMESH(IR))**LVAL*AC
          END DO !I

      END DO !ISPIN
    end do    ! m
  end do       ! l
end do !iatom=1,ntotatom

end subroutine vinters2010 !(natotatom)

end module mod_vinters2010
