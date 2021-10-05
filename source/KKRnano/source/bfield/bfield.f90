!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Module storing the run options and the paramters for bfields and constraining fields
!> 
!> Written by Sascha Brinker, ported from KKRhost to KKRnano by the other two.
!> 
!> Author: Sascha Brinker, Eduardo Mendive, Nicolas Essing
!------------------------------------------------------------------------------------
module mod_bfield
  
  ! use :: NonCollinearMagnetism_mod, only : rotatematrix
  ! use :: mod_mympi, only: myrank, master

  implicit none

  ! In KKRhost, there is a type named 'type_bfield',
  ! which contains the input parameters and the information from 'bfield_data'.
  ! Here, one instance of 'bfield_data' is stored for each atom, while
  ! in KKRhost the components of that type were arrays (respectively arrays with one
  ! dimension more than here) containing the information for all atoms.
  
  !-------------------------------------------------------------------------------
  !> Summary: A type storing information on magnetic fields for a single atom
  !-------------------------------------------------------------------------------
  type :: bfield_data
    double precision, dimension(3) :: bfield !! external magnetic field in cartesian coordinates
    double precision               :: bfield_strength !! absolute value of the external magnetic field, dimensions 
    double precision, dimension(3) :: bfield_constr !! constraining field in cartesian coordinates
    double precision               :: theta !! polar angle of the magnetic field
    double precision               :: phi   !! azimuthal angle of the magnetic field
    double precision, dimension(:,:,:), allocatable :: thetallmat !! shapefun in the ll' expansion
    double precision, dimension(3) :: mag_torque !! Magnetic torque 
  end type

contains

  subroutine load_bfields_from_disk(bfields, lbfield, lbfield_constr)
    type(bfield_data), allocatable :: bfields(:)
    logical, intent(in) :: lbfield
    logical, intent(in) :: lbfield_constr
    
    integer :: number_of_atoms

    number_of_atoms = size(bfields)
    
    if (lbfield) then
      call read_bfield(bfields, number_of_atoms)
    end if
    if (lbfield_constr) then
      call read_bconstr(bfields, number_of_atoms)
    end if
    
  end subroutine

  subroutine init_bfield(lbfield, lbfield_constr, bfield, lmax, &
                         npan_tot, npan_log, npan_eq, ncheb, ipan_intervall, thetasnew, &
                         iend, icleb, cleb, ifunm)
    logical, intent(in) :: lbfield        !! Whether to use external nonco bfields
    logical, intent(in) :: lbfield_constr !! Whether to use constraint bfields
    type(bfield_data), intent(inout) :: bfield !! The bfield data
    integer, intent(in) :: lmax     !! Angular momentum cutoff
    integer, intent(in) :: npan_tot !! Chebyshev mesh resolution
    integer, intent(in) :: npan_log !! Chebyshev mesh resolution
    integer, intent(in) :: npan_eq  !! Chebyshev mesh resolution
    integer, intent(in) :: ncheb    !! Chebyshev mesh resolution
    integer, dimension(:), intent(in) :: ipan_intervall !! Indices for important radial points in the mesh
    double precision, dimension(:,:), intent(in)  :: thetasnew !! Interpolated shape function in Chebychev mesh
    integer,                          intent (in) :: iend      !! Number of nonzero gaunt coefficients
    integer,          dimension(:,:), intent (in) :: icleb     !! Mapping from array index to angular momentum indices
    double precision, dimension(:),   intent (in) :: cleb      !! gaunt coefficients
    integer, dimension(:),            intent (in) :: ifunm     !! Switch for use of gaunt coefficients
    
    integer :: ncleb, irmdnew
    integer :: i_stat

    ncleb = size(cleb)
    irmdnew = npan_tot*(ncheb+1)

    !TODO? assert size(theasnew) == irmdnew

    ! Calculate the LL' expansion of the shape function in the new mesh which
    ! is needed to convolute the magnetic field (only done once and stored to
    ! safe computing time)
    ! TODO: the ".or. lbfield_constr" is not present in host. Is it only needed for externals?
    if(lbfield .or. lbfield_constr) then
      allocate (bfield%thetallmat((lmax+1)**2,(lmax+1)**2,irmdnew), stat=i_stat)
      call calc_thetallmat(bfield%thetallmat, lmax, ipan_intervall(npan_log+npan_eq) + 1, &
                            iend, irmdnew, thetasnew, ifunm, icleb, cleb)
    end if

  end subroutine init_bfield

  !TODO save_bconst()

  !-------------------------------------------------------------------------------
  !> Summary: Reads the atom-wise constraining field from bconstr.dat
  !> Author: MdSD
  !> Category: KKRhost, bfield
  !> Deprecated: False
  !> the file has the format:
  !> bz  by  bz (in Ry), mspin (in muB)
  !-------------------------------------------------------------------------------
  subroutine read_bconstr(bfields, number_of_atoms)
    type(bfield_data), dimension(:), intent(inout) :: bfields
    integer                        , intent(in)  :: number_of_atoms
    
    integer :: iatom
    integer :: iostat
   
    open(unit=57493215, file='bconstr_in.dat', iostat=iostat)
    if (iostat == 0) then 
      read(57493215,*)  ! skip header
      do iatom = 1, number_of_atoms
        read(57493215,*,iostat=iostat) bfields(iatom)%bfield_constr(:)
        if (iostat > 0) then
          ! iostat should be zero normally and negative in last line, positive is error
          write(*,*) "Error reading bconstr_in.dat"
          stop
        end if
      end do
      if (iostat >= 0) then
        ! still zero means there are more lines
        write(*,*) "Error reading bconstr_in.dat or file too long"
      end if
      close(57493215)
    else
      do iatom = 1, number_of_atoms
        bfields(iatom)%bfield_constr(:) = 0.
      end do
    end if

    write(*,*) '  ############################################################'
    write(*,*) '  input constraining fields'
    write(*,*) '  ############################################################'
    write(*,*) '  iatom      Bx              By              Bz        (in Ry)'
    do iatom = 1, number_of_atoms
      write(*,'(2X,I4,3(E16.8))') iatom, bfields(iatom)%bfield_constr(:)
    end do
  end subroutine read_bconstr


  !-------------------------------------------------------------------------------
  !> Summary: Reads the atom-wise magnetic field from bfield.dat
  !> Author: Sascha Brinker
  !> Category: KKRhost, bfield
  !> Deprecated: False
  !> the file has the format:
  !> theta  phi  bfield_strength (in Tesla)
  !-------------------------------------------------------------------------------
  subroutine read_bfield(bfields,number_of_atoms)
    type(bfield_data), dimension(:), intent(inout) :: bfields
    integer                        , intent(in)   :: number_of_atoms
    
    integer        :: iatom, iostat
    character(256) :: linebuffer
   
    open(unit=57493215, file='bfield.dat', iostat=iostat)
    if (iostat /= 0) then
      write(*,*) "Error reading bfield.dat."
      write(*,*) "Provide correct file or turn off magnetic fields."
      stop
    end if

    iatom = 1
    do while ((iostat == 0) .and. (iatom <= number_of_atoms))
      read(57493215, *, iostat=iostat) linebuffer
      if (linebuffer(1:1) == '#') continue ! input line commented out, get next line
      read(linebuffer, *, iostat=iostat) bfields(iatom)%theta, bfields(iatom)%phi, bfields(iatom)%bfield_strength
      if (iostat /= 0) then
        write(*,*) "Malformatted line in bfield.dat."
        stop
      end if
      iatom = iatom + 1
    end do
    close(57493215)
    if (iatom <= number_of_atoms + 1) then
      write(*,*) "Error reading bfield.dat or file too short"
      stop
    end if

    write(*,*) '  ###############################################'
    write(*,*) '  external non-collinear magnetic fields'
    write(*,*) '  ###############################################'
    write(*,*) '  iatom      theta       phi         bfield (in Ry)'
    do iatom = 1, number_of_atoms
      bfields(iatom)%theta           = bfields(iatom)%theta / 360.0d0 * 8.d0 * datan(1.d0)
      bfields(iatom)%phi             = bfields(iatom)%phi   / 360.0d0 * 8.d0 * datan(1.d0)
      bfields(iatom)%bfield_strength = bfields(iatom)%bfield_strength ! / 235051.787 ! conversion from Tesla to Ry
      bfields(iatom)%bfield(1)       = bfields(iatom)%bfield_strength*cos(bfields(iatom)%phi)*sin(bfields(iatom)%theta)
      bfields(iatom)%bfield(2)       = bfields(iatom)%bfield_strength*sin(bfields(iatom)%phi)*sin(bfields(iatom)%theta)
      bfields(iatom)%bfield(3)       = bfields(iatom)%bfield_strength*cos(bfields(iatom)%theta)
      write(*,'("  ",i4,3es16.8)') iatom, bfields(iatom)%theta, bfields(iatom)%phi, bfields(iatom)%bfield_strength
   end do
  end subroutine read_bfield
  
  !TODO add_bfield()
  
  !------------------------------------------------------------------------------------
  !> Summary: Shape function LL' expansion
  !> Author: Sascha Brinker
  !> Category: KKRhost, geometry, new-mesh, shapefun
  !> Deprecated: False 
  !> Calculates the LL' expansion of the shape function similarly to vllmat_new.f90
  !> @note The input shapefunction (single L) uses pointer arrays for the lm index.
  !> The output does not need pointers!
  !> @endnote
  !------------------------------------------------------------------------------------
  subroutine calc_thetallmat(thetansll, lmax, imt1, iend, irmdnew, thetasnew, ifunm, icleb, cleb)
    integer                     , intent (in) :: lmax    !! Angular momentum cut-off
    integer                     , intent (in) :: imt1    !! index muffin-tin radius
    integer                     , intent (in) :: iend    !! Number of nonzero gaunt coefficients
    integer                     , intent (in) :: irmdnew !! number of radials point on the Cheby mesh
    integer, dimension(:,:), intent (in) :: icleb   !! Pointer array
    double precision, dimension(:),    intent (in)  :: cleb    !! GAUNT coefficients (GAUNT)
    double precision, dimension(:, :), intent (in)  :: thetasnew !! shapefun on the Cheby mesh
    integer, dimension(1:(2*lmax+1)**2)                , intent (in)  :: ifunm     !! pointer array for shapefun     ! Check index and dimensions!!!!!
    double precision, dimension((lmax+1)**2,(lmax+1)**2,irmdnew), intent (out) :: thetansll !! LL' expansion of the shapefunction 
    !------------------------------------------------------------------------------------
    double precision,parameter                                  :: rfpi=3.5449077018110318
    double precision,dimension(1:irmdnew,1:(2*lmax+1)**2)       :: shapefun_mod  
    double precision                                            :: c0ll
    integer                                                     :: lmmax
    integer                                                     :: lmmax2
    integer                                                     :: ilm
    integer                                                     :: ifun
    integer                                                     :: ir
    integer                                                     :: lm1
    integer                                                     :: lm2
    integer                                                     :: lm3
    integer                                                     :: j
            
    c0ll                                                    = 1d0/rfpi
    lmmax                                                   = (lmax+1)**2
    lmmax2                                                  = (2*lmax+1)**2
    shapefun_mod(:,:)                                       = 0.d0
    shapefun_mod(1:imt1,1)                                  = rfpi ! is multipled by C_LL^0
    shapefun_mod(imt1+1:irmdnew,1)                          = thetasnew(imt1+1:irmdnew,1)
    ! convert from pointer to real ilmdo ilm=2,lmmax2
    do ilm=2,lmmax2
      ifun = ifunm(ilm)
      if(.not. ifun == 0) then !shapefun%lmused(ilm)==1) then
        !if(myrank==master .and. debug) write(17*iatom**2,'(" converted ifun= ",i4," to ilm= ",i4)') ifun, ilm
        shapefun_mod(imt1+1:irmdnew,ilm)           = thetasnew(imt1+1:irmdnew,ifun)
      end if
    end do
    !if(myrank==master .and. debug) write(17*iatom**2,'("  cellnew%shapefun = ")')
    !do ilm= 1,lmmax2
    !  if(sum(abs(thetasnew(:,ilm)))>1e-8) then
    !    if(myrank==master .and. debug) write(17*iatom**2,'(i4,1000es16.8)') ilm, thetasnew(:,ilm)
    !  end if
    !end do
    !if(myrank==master .and. debug) write(17*iatom**2,'("  shapefun_mod = ")')
    !do ilm= 1,lmmax2
    !  if(sum(abs(shapefun_mod(:,ilm)))>1e-8) then
    !    if(myrank==master .and. debug) write(17*iatom**2,'(i4,1000es16.8)') ilm, shapefun_mod(:,ilm)
    !  end if
    !end do
    
    thetansll(:,:,:)=0.d0
    ! diagonal part (not contained in gaunt-coeff)
    do ilm = 1,lmmax
      thetansll(ilm,ilm,:) = shapefun_mod(:,1)*c0ll
    end do
    !write(*,'("ifunm=",1000i4)') ifunm(:)
    do j = 1,iend !gauntcoeff%iend
      lm1 = icleb(j, 1)! gauntcoeff%icleb(j,1) ! lmax
      lm2 = icleb(j, 2)! gauntcoeff%icleb(j,2) ! lmax
      lm3 = icleb(j, 3)! gauntcoeff%icleb(j,3) ! 2*lmax
      !write(17*iatom**2,'("lm1,lm2,lm3 =",3i4,2es16.8)') lm1,lm2,lm3,cleb(j)
      if(lm1<= lmmax .and. lm2 <= lmmax .and. lm3<= lmmax2) then
        ifun = ifunm(lm3)
        if(.not. ifun == 0) then !shapefun%lmused(ilm)==1) then
          !write(*,'("lm1,lm2,lm3,cleb",3i4,10e16.8)') lm1,lm2,lm3,cleb(j)
          do ir = 1,irmdnew!cellnew%nrmaxnew
            thetansll(lm1,lm2,ir) = thetansll(lm1,lm2,ir)+cleb(j)*shapefun_mod(ir,lm3)
            thetansll(lm2,lm1,ir) = thetansll(lm2,lm1,ir)+cleb(j)*shapefun_mod(ir,lm3)
          end do
        end if
      end if
    end do
    !do lm1 = 1,lmmax
    !  do lm2 = 1,lm1-1
    !    do ir = 1, irmdnew
    !      thetansll(ir,lm2,lm1) = thetansll(ir,lm1,lm2)
    !    end do
    !  end do
    !end do
    !if(myrank==master .and. debug) write(17*iatom**2,'("  thetansll = ")')
    !do lm1=1,lmmax
    !  do lm2=1,lmmax
    !    if(sum(abs(thetansll(:,lm1,lm2)))>1e-8) then
    !      if(myrank==master .and. debug) write(17*iatom**2,'(2i4,1000es16.8)') lm1, lm2, thetansll(:,lm1,lm2)
    !    end if
    !  end do
    !end do
              
  end subroutine calc_thetallmat

end module mod_bfield





