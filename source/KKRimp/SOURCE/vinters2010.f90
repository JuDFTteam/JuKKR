module mod_vinters2010

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculate and add intercell contribution to potential
  !> Author: B. Drittler
  !> Date: June 1987
  !> Category: KKRimp, potential
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> calculate the intercell-potentials and add these to the potential v (in the spin-polarized case for each spin-direction
  !> the intercell-potential is the same) it uses the structure dependent matrices amat and bmat which
  !> are calculate once in the subroutine amn. The charge-moments are calculated in the subroutine vintra2,
  !> therefore vintra2 has to be called first. The intercell-potential is expanded into spherical harmonics.
  !> The lm-term of the intercell-potential v of the representive atom i is given by
  !> 
  !> v(r,lm,i) =  (-r)**l * {amat(i1,i2,lm,l'm')*cmom(i2,l'm') + bmat(i1,i2,lm)*z(i2)}
  !> summed over i2 (all shells) and l'm' .    (i1=i-natref)
  !> (see notes by B. Drittler)
  !> 
  !> In case of shape correction the madelung potential of the host
  !> is taken into account . in all other case the madelung potential of the host is set to be zero.
  !> as actual values for z and cmom the differences between the values of the  representive atoms and those of the references
  !> are used.
  !> 
  !> @warning The first index of cmom (moment of the charge
  !> density - calculated in vintr2) and of z (nuclear
  !> charge of the atoms in the shell) is in the program
  !> different defined , there one has to use :
  !> cmom(natref+i2,l'm')  and  z(natref+i2) @endwarning
  !>
  !> @notes
  !> In the case of impurities on surfaces, the intercell potential of the
  !> reference system is read in from the fxdr file 'intercell_ref'; which
  !> is calculated in the surface program. - July 1998
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine vinters2010(natom,nspin,lmaxd,lmaxatom,cell,ratom,cmom,alat,ins,nrmaxd,vpot_out,intercell_ach,zat,lattice_relax,rimpshift)
    use nrtype, only: dp
    use type_cell, only: cell_type
    use mod_shftvout, only: shftvout
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
    real(kind=dp)              :: ratom(3,natom)
    real(kind=dp)              :: cmom((2*lmaxd+1)**2,natom) !(lmpotd,ntotatom)
    real(kind=dp)              :: zat(natom) !(lmpotd,ntotatom),achnew(lmpotd,ntotatom)
    real(kind=dp)              :: alat
    integer                    :: ins
    integer                    :: nrmaxd
    real*8                     :: vpot_out(nrmaxd,(2*lmaxd+1)**2,nspin,natom) !thetas(iri,nfund,*),
    integer                    :: lattice_relax
    real(kind=dp)              :: rimpshift(3,natom)
    !local
    integer                    :: ispin,iatom,jatom,lval,mval,lm,lm2,ir,irs1
    integer                    :: lmpotd
    real(kind=dp),allocatable  :: intercell_ach(:,:) !(lmpotd,ntotatom),achnew(lmpotd,ntotatom)
    real(kind=dp)              :: intercell((2*lmaxd+1)**2) !(lmpotd,ntotatom),achnew(lmpotd,ntotatom)
    real(kind=dp)              :: intercell_temp((2*lmaxd+1)**2) !(lmpotd,ntotatom),achnew(lmpotd,ntotatom)
    real(kind=dp)              :: ac
    real(kind=dp),allocatable  :: amat(:,:,:)
    real(kind=dp),allocatable  :: bmat(:,:)

    if (.not. allocated(gauntcoeff)) stop '[preconditioning_intercell] gauntcoeff not allocated'
    lmpotd=(2*lmaxd+1)**2
    allocate( amat(natom,lmpotd,lmpotd), bmat(natom,lmpotd ) ) 
    amat=0.0d0
    bmat=0.0d0

    do iatom=1,natom

      if (lattice_relax==0) then
        intercell=intercell_ach(1:(2*lmaxd+1)**2,iatom)
      else
        intercell_temp=intercell_ach(1:(2*lmaxd+1)**2,iatom)
        call shftvout(intercell_temp,intercell,rimpshift(:,iatom), 2*lmaxd,gauntcoeff(lmaxd)%wg,gauntcoeff(lmaxd)%yrg,(lmaxd+1)**2,4*lmaxd,2*lmaxd,(4*lmaxd+1)**2,gauntcoeff(lmaxd)%ncleb)
      end if

      call amn2010(iatom,amat,bmat,alat,gauntcoeff(lmaxd),2*lmaxd,natom,ratom)

      do lval = 0,2*lmaxatom(iatom) !*lmaxd
        do mval = -lval,lval
          lm = lval*lval + lval + mval + 1
          ac = 0.0d0
          if (.not. config_testflag('hostintercell')) then
            do jatom = 1,natom
              do lm2 = 1,(2*lmaxatom(jatom)+1)**2
                ac = ac + amat(jatom,lm,lm2)*cmom(lm2,jatom)
              end do !lm2
              ac = ac + bmat(jatom,lm)*zat(jatom)
            end do
          end if !(.not. config_testflag('hostintercell')) then

          ac = ac + intercell(lm)

          if (ins.ne.0) then
            irs1 = cell(iatom)%nrcut(cell(iatom)%npan)
          else
            irs1 = cell(iatom)%nrmax
          end if

          do ispin = 1,nspin
            ! code has problems evaluating 0**0 (0 to the power of zero)
            if (lval.eq.0) then 
              vpot_out(1,1,ispin,iatom) = vpot_out(1,1,ispin,iatom) + ac
            end if
            do ir = 2,irs1
              vpot_out(ir,lm,ispin,iatom) = vpot_out(ir,lm,ispin,iatom) + (-cell(iatom)%rmesh(ir))**lval*ac
            end do !i

          end do !ispin
        end do    ! m
      end do       ! l
    end do !iatom=1,ntotatom

  end subroutine vinters2010

end module mod_vinters2010
