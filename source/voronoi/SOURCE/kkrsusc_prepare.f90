  subroutine kkrsusc_prepare(natom,nspin,lmaxd,zatom,nrmaxd,vpot,density,cell,config,itscf,inpsusc,my_rank)    ! sexy |||| added by Manuel (2014/12)

    ! prepares for running the KKRsusc program  ! susc |||| added by Benedikt (2014/06)

    ! ---------------------------------------------------------------------------------
    ! used variable type(s):
    use type_cell
    use type_config
    use type_density
    use type_inpsusc
    use mod_physic_params, only: cvlight

    use mod_wfmesh
    use mod_cradwf
    use mod_config                  ! Juba (2016)

    !
    ! ---------------------------------------------------------------------------------
    ! used module(s):
!    use projection
    use global                      ! sexy |||| added by Manuel (2014/12)
    !

    ! ---------------------------------------------------------------------------------
    implicit none

    !***********************************
    ! type variables
    !***********************************

    type(inpsusc_type),intent(out)   :: inpsusc   ! susc |||| added by Benedikt (2014/06)


    ! ---------------------------------------------------------------------------------
    integer, intent(in)              :: natom                                 ! number of impurity atoms
    integer, intent(in)              :: nspin                                 ! number of spin directions
                                                                              ! 1= paramagnetic
                                                                              ! 2= collinear,non collinear, SOC, ...
    integer, intent(in)              :: lmaxd                                 ! maximum lmax of all atoms 
                                                                              ! meaning max(lmaxatom(iatom),iatom=1,natom)
    integer, intent(in)              :: nrmaxd                                ! cell(1)%nrmaxd
    real(kind=dp),intent(in)         :: zatom(natom)                          ! nucleus charge
    type(density_type),intent(in)    :: density(natom)                        ! density
    type(cell_type), intent(in)      :: cell(natom)                           ! cell properties for each atom
    type(config_type),intent(in)     :: config                                ! type containing all config variables
    real(kind=dp),intent(in)         :: vpot(nrmaxd,(2*lmaxd+1)**2,nspin,natom)
    integer, intent(in)              :: itscf                                 ! current iteration    


    integer                          :: iatom, ispin, il, ia, ib, nb, imesh
    complex(kind=dpc)                :: ene_susc, ek_susc


    ! mpi 
    integer, intent(in)              :: my_rank

    ! ---------------------------------------------------------------------------------
    ! projection module to read in variables from input file ('inpsusc.dat'):
!    call susc_input(lmaxd,natom,nspin,nrmaxd) 
    call new_input(lmaxd,natom,nspin,nrmaxd,itscf,my_rank)  ! newinpsusc   ! sexy |||| added by Manuel (2014/12)

    ! ---------------------------------------------------------------------------------
    ! 

    ! ---------------------------------------------------------------------------------
    ! initialize some values:
    inpsusc%nr0 = 11                                      ! susc
    inpsusc%nr1 = nrmaxd                                  ! susc
    inpsusc%nr  = inpsusc%nr1 - inpsusc%nr0 + 1           ! susc
    inpsusc%nrmaxd = nrmaxd

    ! ---------------------------------------------------------------------------------
    ! write statements:
    if (my_rank == 0) then
      if ( itscf==1 ) then
        write(*,'("  [kkrsusc_prepare] KKRimp program uses following vales:")')
        write(*,'("                        lmaxd  = ",i0)') lmaxd
        write(*,'("                        natom  = ",i0)') natom
        write(*,'("                        nasusc = ",i0)') nasusc
        write(*,'("                        nspin  = ",i0)') nspin
        write(*,'("                        nrmaxd = ",i0)') nrmaxd
        write(*,'("  [kkrsusc_prepare] KKRimp program uses radial mesh points:")')
        write(*,'("                        nr1    = ",i0)') inpsusc%nr1
        write(*,'("                        nr0    = ",i0)') inpsusc%nr0
        write(*,'("                        nr     = ",i0)') inpsusc%nr
      end if
    end if ! my_rank 
!    call outsusc_out


    ! ---------------------------------------------------------------------------------
    ! allocate some arrays:
    allocate ( inpsusc%rs(inpsusc%nr1,0:lmaxd) )          ! susc
    allocate ( inpsusc%s(0:lmaxd) )                       ! susc
    allocate ( inpsusc%rsdummy(nrmaxd,0:lmaxd) )          ! susc
    allocate ( inpsusc%dlogdp(0:lmaxd) )                  ! susc
    allocate ( inpsusc%cutoff(inpsusc%nr1) )              ! susc
    allocate ( inpsusc%alpha(0:lmaxd) )                   ! susc
    allocate ( inpsusc%tmat(0:lmaxd) )                    ! susc
    allocate ( inpsusc%fz(inpsusc%nr1,0:lmaxd) )          ! susc
    allocate ( inpsusc%pz(inpsusc%nr1,0:lmaxd) )          ! susc
    allocate ( inpsusc%qz(inpsusc%nr1,0:lmaxd) )          ! susc
    allocate ( inpsusc%sz(inpsusc%nr1,0:lmaxd) )          ! susc
    allocate ( inpsusc%rho2ns_tmp(nrmaxd,(2*lmaxd+1)**(2),natom,nspin) )  ! sexy |||| added by Manuel (2014/12)
    allocate ( inpsusc%espv_tmp(0:lmaxd,natom,nspin) )                    ! sexy |||| added by Manuel (2014/12)
!   Non spin diagonal t-matrix                            ! Juba
    if (config_testflag('tmatnew')) then
      allocate ( inpsusc%tmatll(2*(lmaxd+1)**(2),2*(lmaxd+1)**(2),natom) )  
    else
      allocate ( inpsusc%tmatll((lmaxd+1)**(2),(lmaxd+1)**(2),natom) ) 
    endif 


!   Read the wave function in restart mode
!    if (.not.lrestart) then     
      ! ---------------------------------------------------------------------------------
      ! loop over atoms to construct the basis set for the wave functions:
      !
      do ia=1,nasusc
        iatom = iasusc(ia)
        !
        ! the value below is hardcoded, this is a dummy variable for the wfmesh call.
        ! It is only used to define the variable Ek, which is not needed at this stage.
        ene_susc = 1.0  
        !
        call wfmesh(ene_susc,ek_susc,cvlight,config%nsra,zatom(iatom),cell(iatom)%rmesh(:), &
                  ! >  <    >         >           >                 >
                  & inpsusc%s,inpsusc%rs, inpsusc%nr1,lmaxd) ! susc
                  !    <           <         >          >
        !
        !
        inpsusc%rsdummy = 0.d0
        inpsusc%rsdummy(1:nrmaxd,:) = inpsusc%rs(:,:) 
        !
        !   -> call save_rmesh

        call save_rmesh(ia,inpsusc%nr0,inpsusc%nr1,cell(iatom)%rmesh(:),inpsusc%rsdummy, & ! susc
           cell(iatom)%drmeshdi(:),zatom(iatom),vpot(:,1,nspin,iatom),vpot(:,1,1,iatom), & ! susc
           density(iatom)%rhoc(:,nspin),density(iatom)%rhoc(:,1),cell(iatom)%npan,cell(iatom)%nrcut,my_rank)                         ! susc 
        !
        !

        ! for each element of the basis set compute the basis wave functions and orthogonalize them:
        !   [[[ in ib loop: get energies -> compute wave functions with cradwf ]]]
        !   [[[ after ib loop: call new_basis1 for il and ispin -> othogonalized basis functions ]]]
        !   [[[ store them in file with save_wfn subroutine ]]]
        if (my_rank == 0) then
          if ( itscf==1 ) then
            if (ia==1) write(*,'("  [kkrsusc_prepare] Setting up the basis for the projected Green functions:")')
            write(*,           '("                      Atom, l-value, spin, basisfunction, energy")')
          end if
        end if ! my_rank 
!       »»»»»»»»»»»»»»»»»»»»»»»»     ! sexy |||| added by Manuel (2014/12)
!       only new basis if needed
        if (.not.lbasisfreeze) then  ! newinpsusc
!       »»»»»»»»»»»»»»»»»»»»»»»»

        do il=0,nlmax   
          do ispin=1,nspin
            ! fetch energy for which the wave function is constructed:
            nb= iwsusc(il,ispin,ia)
            do ib=1,nb
              ene_susc = ewsusc(ib,il,ispin,ia)
              if (my_rank == 0) then 
                if ( itscf==1 ) write(*,     '("                  ",i8,i9,i6,i15,2es16.8)') ia, il, ispin, ib, ene_susc
              end if ! my_rank 
              call wfmesh(ene_susc,ek_susc,cvlight,config%nsra,zatom(iatom),cell(iatom)%rmesh(:), &
                        !    >        <       >         >           >                 >
                        & inpsusc%s,inpsusc%rs, inpsusc%nr1,lmaxd) ! susc
                        !    <           <         >          >
              call cradwf(ene_susc,ek_susc,config%nsra,inpsusc%alpha,cell(iatom)%npan,cell(iatom)%nrcut,cvlight,       &
                        !      >
                        & inpsusc%rs,inpsusc%s,inpsusc%pz,inpsusc%fz, inpsusc%qz,inpsusc%sz,inpsusc%tmat,vpot(:,1,ispin,iatom), &
                        !
                        & cell(iatom)%drmeshdi,cell(iatom)%rmesh,zatom(iatom),.true., 0,0,0.d0,inpsusc%cutoff,lmaxd,lmaxd+1,inpsusc%nr1)
                        !
              call ref_wfn(ia,il,ispin,ib,inpsusc%pz(:,il),cell(iatom)%npan,cell(iatom)%nrcut)

            end do !ib=1,nb

!          ! orthogonalize the wave functions (Gram-Schmidt procedure):
!          nb = iwsusc(il,ispin,ia)
!!           write(*,*) "nb = " , nb
!          if (nb .gt. 0) then
!            call new_basis1(nb,ia,il,ispin,1.d-9)
!!             write(*,'(A,4i6)') "New basis1(atom,l,spin,nb): ", ia, il, ispin, nb ! susc
!            call save_wfns(ia,il,ispin) ! susc
!          end if ! nb>0 ! susc
          end do !ispin ! susc
        end do !il ! susc

!      »»»»»»
       end if
!      »»»»»»

      end do ! ia=1,nasusc
   

!     construct the basis out of the ref wfns      ! sexy |||| added by Manuel (2014/12)
!     if not using I/O allocate memory for projection
      call new_basis(itscf,my_rank)
  
!    end if ! lrestart

!   test: overlaps
    if (my_rank == 0) then 
      if (lhdio) then
        do ia=1,nasusc
          call out_overlap(ia,ia)
        end do
      end if
    end if ! my_rank

!   initial info for outsusc.dat
    if (my_rank == 0) then
      if (lhdio) call outsusc_out
    end if ! my_rank

    ! ---------------------------------------------------------------------------------
    ! deallocate some arrays:
    deallocate ( inpsusc%rsdummy, inpsusc%rs, inpsusc%s )


  end subroutine kkrsusc_prepare
