!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Module storing the run options and the paramters for bfields and constraining fields
!> Author: Sascha Brinker
!> 
!> 
!> 
!> 
!> 
!------------------------------------------------------------------------------------
module mod_bfield
  
  use :: mod_profiling
  use :: mod_datatypes
  use :: global_variables, only: lmmaxd, ncleb, ntotd, nfund
  use :: mod_types, only: t_inc
  use :: mod_rotatespinframe, only: rotatematrix
  use :: mod_mympi, only: myrank, master

  implicit none

  !-------------------------------------------------------------------------------
  !> Summary: Type used in t_params to store all relevant information for bfields and constraining fields
  !> Author: Sascha Brinker
  !> Category: communication, KKRhost, bfield
  !> Deprecated: False
  !> 
  !-------------------------------------------------------------------------------
  type :: type_bfield
    logical :: lbfield = .False. ! external magnetic field (turned on via runoption <noncobfield>) non-collinear magnetic field
    logical :: lbfield_constr = .False. ! constraining fields (turned on via runoption <noncobfield>) non-collinear magnetic field
    logical :: lbfield_all = .False. ! apply same field to all atoms (True) or individual fields to each atom
    logical :: lbfield_trans = .False. ! apply only transversal field
    logical :: lbfield_mt = .False. ! apply only field in muffin-tin
    logical :: ltorque = .False. ! calculate magnetic torque
    integer :: ibfield = 0 ! spin (0), orbital (1), spin+orbial (2) fields
    integer :: ibfield_constr = 0 ! type of contraint (0 = torque, 1 = magnetic moment)
    integer :: itscf0 = 0    ! start magnetic field at iteration itscf0
    integer :: itscf1 = 10000 ! stop applying magnetic field after iteration itscf1
    real (kind=dp), dimension (:,:), allocatable        :: bfield ! external magnetic field in cartesian coordinates, dimensions (natom,3)
    real (kind=dp), dimension (:), allocatable          :: bfield_strength ! absolute value of the external magnetic field, dimensions (natom)
    real (kind=dp), dimension (:,:), allocatable        :: bfield_constr ! constraining field in cartesian coordinates, dimensions (natom,3)
    real (kind=dp), dimension (:), allocatable          :: theta ! polar angle of the magnetic field
    real (kind=dp), dimension (:), allocatable          :: phi   ! azimuthal angle of the magnetic field
    real (kind=dp), dimension (:,:,:,:), allocatable    :: thetallmat ! shapefun in the ll' expansion
    !------------------------------------------------------------------------------------
    ! Magnetic torque 
    !------------------------------------------------------------------------------------
    real(kind=dp),dimension(:,:), allocatable           :: mag_torque

  end type type_bfield

  type (type_bfield), save :: bfield

contains
  !-------------------------------------------------------------------------------
  !> Summary: Allocate initial magnetic field parameters to be broadcasted via mpi
  !> Author: Sascha Brinker
  !> Category: memory-management, profiling, KKRhost, bfield
  !> Deprecated: False
  !-------------------------------------------------------------------------------
  subroutine init_bfield(bfield,natyp,lbfield,lbfield_constr,lbfield_all,lbfield_trans,lbfield_mt,ltorque,ibfield,ibfield_constr,itscf0, &
    itscf1,npan_log,npan_eq,ncheb,ntotd,nfund,ncelld,lmax,iend,ntcell,ipan_intervall,ifunm,icleb,cleb,thetasnew)

      implicit none

      type(type_bfield), intent(inout) :: bfield
      integer, intent(in) :: natyp 
      integer :: i_stat
      
      logical, intent(in) :: lbfield 
      logical, intent(in) :: lbfield_constr 
      logical, intent(in) :: lbfield_all 
      logical, intent(in) :: lbfield_trans
      logical, intent(in) :: lbfield_mt
      logical, intent(in) :: ltorque
      integer, intent(in) :: ibfield 
      integer, intent(in) :: ibfield_constr 
      integer, intent(in) :: itscf0 
      integer, intent(in) :: itscf1 
      integer, intent(in) :: npan_log
      integer, intent(in) :: npan_eq
      integer, intent(in) :: ncheb
      integer, intent(in) :: ntotd
      integer, intent(in) :: nfund
      integer, intent(in) :: ncelld
      integer, intent(in) :: lmax
      integer, intent (in):: iend         ! Number of nonzero gaunt coefficients
      integer, dimension (natyp), intent (in)                           :: ntcell ! pointer from natyp to ndcell
      integer, dimension (0:ntotd, natyp), intent (in)                  :: ipan_intervall
      integer, dimension (1:(2*lmax+1)**2,natyp), intent (in)           :: ifunm        ! pointer array for shapefun     ! Check index and dimensions!!!!!
      integer, dimension (ncleb, 4), intent (in)                        :: icleb !! Pointer array
      real (kind=dp), dimension (ncleb), intent (in)                    :: cleb !! GAUNT coefficients (GAUNT)
      real (kind=dp), dimension (ntotd*(ncheb+1), nfund, ncelld), intent (in)   :: thetasnew !! interpolated shape function in Chebychev radial mesh
      
      integer   :: icell
      integer   :: i1
      logical, dimension(ncelld)   :: celldone

      !write(*,'("Init bfield")')
      ! init basic parameters
      bfield%lbfield = lbfield
      bfield%lbfield_constr = lbfield_constr
      bfield%lbfield_all = lbfield_all
      bfield%lbfield_trans = lbfield_trans
      bfield%lbfield_mt = lbfield_mt
      bfield%ltorque = ltorque
      bfield%ibfield = ibfield
      bfield%ibfield_constr = ibfield_constr
      bfield%itscf0 = itscf0
      bfield%itscf1 = itscf1
      ! allocate arrays and add to memory screening routine
      allocate (bfield%theta(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(bfield%theta))*kind(bfield%theta), 'bfield%theta', 'init_bfield')
      allocate (bfield%phi(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(bfield%phi))*kind(bfield%phi), 'bfield%phi', 'init_bfield')
      allocate (bfield%bfield(natyp,3), stat=i_stat)
      call memocc(i_stat, product(shape(bfield%bfield))*kind(bfield%bfield), 'bfield%bfield', 'init_bfield')
      allocate (bfield%bfield_strength(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(bfield%bfield_strength))*kind(bfield%bfield_strength), 'bfield%bfield_strength', 'init_bfield')
      allocate (bfield%bfield_constr(natyp,3), stat=i_stat)
      call memocc(i_stat, product(shape(bfield%bfield_constr))*kind(bfield%bfield_constr), 'bfield%bfield_constr', 'init_bfield')
      ! init allocated arrays
      ! bfield%bfield_constr(:,:) = 0.d0
      if(lbfield) call read_bfield(bfield,natyp)
      ! MdSD: constraining fields
      if(lbfield_constr) call read_bconstr(natyp,bfield%bfield_constr)
      

      allocate (bfield%mag_torque(natyp,3), stat=i_stat)
      call memocc(i_stat, product(shape(bfield%mag_torque   ))*kind(bfield%mag_torque   ), 'bfield%mag_torque', 'init_bfield')
      ! Calculate the LL' expansion of hte shape function in the new mesh which
      ! is needed to convolute the magnetic field (only done once and stored to
      ! safe computing time) 
      if(lbfield) then
        allocate (bfield%thetallmat((lmax+1)**2,(lmax+1)**2,ntotd*(ncheb+1),ncelld), stat=i_stat)
        call memocc(i_stat, product(shape(bfield%thetallmat   ))*kind(bfield%thetallmat   ), 'bfield%thetallmat', 'init_bfield')
        celldone(:) = .False.
        do i1=1,natyp
          icell = ntcell(i1)
          if(.not. celldone(icell)) then 
            call calc_thetallmat(bfield%thetallmat(:,:,:,icell), lmax, ipan_intervall(npan_log+npan_eq,i1) + 1, iend, ntotd*(ncheb+1), thetasnew(:,:,icell), ifunm(:,i1), icleb, cleb)
          end if
          celldone(icell) = .True.
        end do
      end if
      !write(*,'("Init bfield done")')
  end subroutine init_bfield


  !-------------------------------------------------------------------------------
  !> Summary: Writes the atom-wise constraining field to bconstr_out.dat
  !> Author: MdSD
  !> Category: KKRhost, bfield
  !> Deprecated: False
  !> the file has the format:
  !> bz  by  bz (in Ry), mspin (in muB)
  !-------------------------------------------------------------------------------
  subroutine save_bconstr(natyp,bconstr_in,bconstr_out)

    implicit none

    integer                          , intent(in)    :: natyp
    real(kind=dp), dimension(4,natyp), intent(in)    :: bconstr_in  !! bx, by, bz, mspin
    real(kind=dp), dimension(natyp,3), intent(inout) :: bconstr_out !! bx, by, bz
    ! local
    integer       :: iatom, nrms
    real(kind=dp) :: rms

    open(unit=57493215, file='bconstr_out.dat', status='replace')
    write(57493215,'("# bconstr_x [Ry], bconstr_y [Ry], bconstr_z [Ry], m_spin [mu_B]")')
    rms = 0.d0; nrms = 0
    do iatom=1,natyp
!     check if this atom is actually being constrained
      if (dot_product(bconstr_in(1:3,iatom), bconstr_in(1:3,iatom)) > 1.d-16) then
        rms = rms + dot_product(bconstr_out(iatom,:)-bconstr_in(1:3,iatom), bconstr_out(iatom,:)-bconstr_in(1:3,iatom))
        nrms = nrms + 1
      end if
      write(57493215,'(4es16.8)') bconstr_in(:,iatom)
      bconstr_out(iatom,:) = bconstr_in(1:3,iatom)
    end do
    close(57493215)
    write(1337,'("Number of constrained atoms=",i8,"  rms for constraining fields=",es16.8)') nrms, sqrt(rms)/max(nrms,1)
  end subroutine save_bconstr


  !-------------------------------------------------------------------------------
  !> Summary: Reads the atom-wise constraining field from bconstr.dat
  !> Author: MdSD
  !> Category: KKRhost, bfield
  !> Deprecated: False
  !> the file has the format:
  !> bz  by  bz (in Ry), mspin (in muB)
  !-------------------------------------------------------------------------------
  subroutine read_bconstr(natyp,bconstr_out)

    implicit none

    integer                          , intent(in)  :: natyp
    real(kind=dp), dimension(natyp,3), intent(out) :: bconstr_out !! bx, by, bz
    ! local
    integer                 :: iatom
    logical                 :: file_exists
   
    inquire(file='bconstr_in.dat',exist=file_exists)
    if (file_exists) then 
      open(unit=57493215, file='bconstr_in.dat')
      read(57493215,*)  ! skip header
      do iatom=1,natyp
        read(57493215,*) bconstr_out(iatom,:)
      end do
      close(57493215)
    else
      bconstr_out(:,:) = 0.0_dp
    end if
    write(1337,*) '  ############################################################'
    write(1337,*) '  input constraining fields'
    write(1337,*) '  ############################################################'
    write(1337,*) '  iatom      Bx              By              Bz        (in Ry)'
    do iatom=1,natyp
      write(1337,'("  ",i4,3es16.8)') iatom, bconstr_out(iatom,:)
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
  subroutine read_bfield(bfield,natyp)
    implicit none

    integer            , intent(in)     :: natyp
    type(type_bfield)  , intent(inout)  :: bfield
    !local
    integer                 :: iatom,ios,ierror
    character(len=200)      :: string1
    integer,save            :: first=1
    
    if (first==1) then
    
      open(unit=57493215, file='bfield.dat', status='old', iostat=ierror)
      if (ierror/=0) then
        write(*,*) '[read_bfield] bfield file does not exist'
        write(*,*) '              setting all bfields to zero'
        if (.not.bfield%lbfield_constr) then
          write(*,*) '              disabling magnetic field lbfield = F'
          bfield%lbfield        = .false.
        end if
        do iatom=1,natyp
          bfield%theta(iatom)            = 0.0D0
          bfield%phi(iatom)              = 0.0D0
          bfield%bfield_strength(iatom)  = 0.0D0
          bfield%bfield(iatom,:)         = 0.0D0
        end do
        return
      end if
      call read_numbofbfields(natyp)
    
    end if
    
    
    
    write(1337,*) '  ###############################################'
    write(1337,*) '  non-collinear magnetic fields'
    write(1337,*) '  ###############################################'
    write(1337,*) '  iatom      theta       phi         bfield (in Ry   )'
    do iatom=1,natyp
       string1=this_readline(57493215,ios)
       if (ios/=0) stop '[read_bfield] Error reading atom info2'
       if (ios==-1) stop '[read_bfield] EOF'
       read(string1,*) bfield%theta(iatom),bfield%phi(iatom),bfield%bfield_strength(iatom)
       write(1337,'("  ",i4,3es16.8)') iatom,bfield%theta(iatom),bfield%phi(iatom),bfield%bfield_strength(iatom)
       bfield%theta(iatom)                  = bfield%theta(iatom)/360.0D0*8.0D0*datan(1.0D0)
       bfield%phi(iatom)                    = bfield%phi(iatom)  /360.0D0*8.0D0*datan(1.0D0)
       bfield%bfield_strength(iatom)        = bfield%bfield_strength(iatom)!/235051.787_dp !conversion from Tesla to Ry
       bfield%bfield(iatom,1)               = bfield%bfield_strength(iatom)*cos(bfield%phi(iatom))*sin(bfield%theta(iatom))
       bfield%bfield(iatom,2)               = bfield%bfield_strength(iatom)*sin(bfield%phi(iatom))*sin(bfield%theta(iatom))
       bfield%bfield(iatom,3)               = bfield%bfield_strength(iatom)*cos(bfield%theta(iatom))
       !write(*,*) iatom,bfield%theta(iatom),bfield%phi(iatom),bfield%bfield_strength(iatom),bfield%bfield(iatom,1),bfield%bfield(iatom,2),bfield%bfield(iatom,3)
    end do
    ! close(33952084)
    first=0
  end subroutine read_bfield
  
  subroutine read_numbofbfields(natom)
    implicit none
    integer,intent(in)      :: natom
    integer                 :: ios,linecount,ierror
    character(len=200)      :: string1
    
    open(unit=55493215, file='bfield.dat', status='old', iostat=ierror)
    ios=0
    linecount = -1
    do while (ios/=-1)
       string1=this_readline(55493215,ios)
       linecount=linecount+1
    end do
    
    if (natom/=linecount) then
      print *,'[read_bfield] number of angles given in file bfield.dat'
      print *,'             is not correct'
      print *,'linecount',linecount
      print *,'natyp',natom
      stop
    end if
    close(55493215)
  end subroutine read_numbofbfields
  
  function this_readline(ifile,ios)
    !--------------------------------------------------------
    !--  reads the next line in file unit ifile            --
    !--------------------------------------------------------
    !--  files starting with a dash (#) are treated as     --
    !--  a comment !!!                                     --
    !--  OUTPUT: next line which is not commented out      --
    !--          error variable IOS (should be zero)       --
    !--------------------------------------------------------
    ! input variables
      implicit none
    integer,intent(in)               :: ifile
    integer,intent(out)              :: ios
    ! local variables
    character(len=200)  ::this_readline
    do
      read(unit=ifile,fmt='(A)', iostat=ios) this_readline
      if (ios/=0 .or. this_readline(1:1)/='#') exit
    end do
  end function this_readline

  
  !------------------------------------------------------------------------------------
  !> Summary: Adds magnetic field to the LL' expansion of the potential
  !> Author: Sascha Brinker
  !> Category: KKRhost, bfield
  !> Deprecated: False 
  !> Calculates the LL' expansion of a magnetic field.
  !> There are two magnetic field, which can be added: a homogeneous magnetic
  !> field read in from bfield.dat and a constraining field, which compensates
  !> the magnetic torque.
  !> @note The magnetic field is added in the following form:
  !> $ H = H_0 - \sigma \cdot \vec{B} $
  !> The input magnetic field is read in from bfield.dat.
  !> It is homogeneous within each cell.
  !> @endnote
  !> @note Runs only with the new solver and either spin-orbit coupling or
  !> noncollinear magnetism.
  !> @endnote
  !------------------------------------------------------------------------------------
  subroutine add_bfield(bfield,iatom,lmax,nspin,irmdnew,imt1,iend,ncheb,theta,phi,ifunm,icleb,cleb,thetasnew,mode,vnspll0,vnspll1,thetansll)

    implicit none

    type(type_bfield)    , intent (in)    :: bfield       ! contains all relevant information about the magnetic fields
    integer              , intent (in)    :: iatom        ! atom index
    integer              , intent (in)    :: lmax         ! angular momentum cut-off
    integer              , intent (in)    :: nspin        ! number of spins
    integer              , intent (in)    :: irmdnew      ! number of radials point on the Cheby mesh
    integer              , intent (in)    :: imt1         ! index muffin-tin radius
    integer              , intent (in)    :: iend         ! Number of nonzero gaunt coefficients
    integer              , intent (in)    :: ncheb        ! Number of Chebychev pannels for the new solver
    real (kind=dp)       , intent (in)    :: theta        ! polar angle of the magnetic moment 
    real (kind=dp)       , intent (in)    :: phi          ! azimuthal angle of the magnetic moment 
    integer              , dimension (1:(2*lmax+1)**2)                     , intent (in) :: ifunm        ! pointer array for shapefun     ! Check index and dimensions!!!!!
    integer, dimension (ncleb, 4), intent (in)                                  :: icleb !! Pointer array
    real (kind=dp), dimension (ncleb), intent (in)                                  :: cleb !! GAUNT coefficients (GAUNT)    ! CHECK THE DIMENSION AND HOW IT IS USED!!!
    real (kind=dp)       , dimension (irmdnew, nfund)     , intent (in) :: thetasnew    ! shapefun on the Cheby mesh
    character (len=*)    , intent (in)                                          :: mode         !! either '1' or 'transpose', depending whether SOC potential is constructed for right or left solution
    complex (kind=dp)    , dimension(lmmaxd, lmmaxd, irmdnew) , intent (in)    :: vnspll0       !! input potential in (l,m,s) basis
    complex (kind=dp)    , dimension(lmmaxd, lmmaxd, irmdnew) , intent (out)   :: vnspll1       !! input potential in (l,m,s) basis
    real(kind=dp) , dimension(1:(lmax+1)**2,1:(lmax+1)**2,1:irmdnew), intent(in)   :: thetansll
    !------------------------------------------------------------------------------------
    ! local variables
    !------------------------------------------------------------------------------------
    integer                                     :: lmmax
    integer                                     :: i
    integer                                     :: j
    integer                                     :: ilm1
    integer                                     :: ilm2
    integer                                     :: ir
    integer                                     :: irend
    integer , dimension(2)                      :: lmstart
    integer , dimension(2)                      :: lmend
    real(kind=dp) , dimension(3)                :: bin(3)
    real(kind=dp),dimension(3)                  :: magdir ! direction of the magnetic moment
    double complex , dimension(2,2)             :: bs ! sigma*b_0
    character(len=1024)                         :: filename
    double complex , dimension(1:irmdnew)       :: temp
    double complex                              :: icompl
    double complex                              :: temp2
    logical                                     :: debug=.False.
    !------------------------------------------------------------------------------------
    ! old variables
    !------------------------------------------------------------------------------------
    !double complex                          :: bxcll(2*(lmax+1)**2,2*(lmax+1)**2,nrmaxnew)
    !double complex                          :: templl(2*(lmax+1)**2,2*(lmax+1)**2,nrmaxnew)
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !!write(*,*) bfield       ! contains all relevant information about the magnetic fields
    !write(*,*) iatom        ! atom index
    !write(*,*) lmax         ! angular momentum cut-off
    !write(*,*) nspin        ! number of spins
    !write(*,*) irmdnew      ! number of radials point on the Cheby mesh
    !write(*,*) imt1         ! index muffin-tin radius
    !write(*,*) iend         ! Number of nonzero gaunt coefficients
    !write(*,*) ncheb        ! Number of Chebychev pannels for the new solver
    !write(*,*) theta        ! polar angle of the magnetic moment 
    !write(*,*) phi          ! azimuthal angle of the magnetic moment 
    !write(*,*) ifunm        ! pointer array for shapefun
    !write(*,*) icleb !! Pointer array
    !write(*,*) cleb !! GAUNT coefficients (GAUNT)
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    
    if (t_inc%i_write>1 .and. mode=='1') then
      write(1337,'("===============================================================================")')
      write(1337,'("                      Magnetic fields for atom ",i4)') iatom
      write(1337,'("===============================================================================")')
    end if
    
    !lmmax=lmmaxd ! could be replace it is still here due to copy/paste
    lmmax=(lmax+1)**2
    icompl=(0d0,1d0) 
    lmstart(1)= 1
    lmstart(2)= lmmax+1
    lmend(1)=lmmax
    lmend(2)=2*lmmax

    !if(myrank==master .and. debug) then
    !  write(filename,'("debug_bfield_",i3.3)') iatom
    !  open(unit=17*iatom**2,file=filename)
    !  write(17*iatom**2,'(" mode =")')
    !  write(17*iatom**2,*) mode 
    !  write(17*iatom**2,'("  lmax= ",i4,"  lmmax= ",i4,"  irmdnew= ",i4,"  imt1= ",i4,"  ncheb= ",i4,"  iend= ",i4,"  iatom= ",i4,"  nspin= ",i4,"  theta= ",f16.8,"  phi= ",f16.8)') lmax,lmmax,irmdnew,imt1,ncheb,iend,iatom,nspin,theta,phi
    !  write(17*iatom**2,'("  bfield is =",3f16.8,"  theta, phi= ",2f16.8)') bfield%bfield(iatom,:),bfield%theta(iatom),bfield%phi(iatom)
    !end if
    
    !!!! calc LL' expansion of the shapefunction
    !!!call calc_thetallmat(thetansll, lmax, imt1, iend, irmdnew, thetasnew, ifunm, icleb, cleb)
    !!!
    !!!if(myrank==master) then
    !!!  write(*,'("  thetansll = ")')
    !!!  do ilm1=1,lmmax
    !!!    do ilm2=1,lmmax
    !!!      if(sum(abs(thetansll(ilm1,ilm2,:)))>1e-8 .or. sum(abs(thetansll_new(ilm1,ilm2,:)))>1e-8) then
    !!!        write(*,'(2i4,1000es16.8)') ilm1, ilm2, thetansll(ilm1,ilm2,:)
    !!!        write(*,'(2i4,1000es16.8)') ilm1, ilm2, thetansll_new(ilm1,ilm2,:)
    !!!      end if
    !!!    end do
    !!!  end do
    !!!end if
    
    ! First add inhomogeneous bfield here (constraining field from xc potential)
    !if( density%magmomentfixed == 6 ) then ! add xc field as constraining field here
    !  write(17*iatom**2,'(" Adding xc constraining field")')
    !  write(17*iatom**2,'(" c_xc= ",10es16.8)') bfield%c_xc(:)
    !  !bxcll(1:lmmax,1:lmmax,:) = (vpotll(1:lmmax,1:lmmax,:)-vpotll(lmmax+1:2*lmmax,lmmax+1:2*lmmax,:))/2.d0 
    !  call vllmat(templl,cellnew%vpotnew(:,:,:),lmax,(lmax+1)**2,lmax,1, &
    !                cellnew%nrmaxnew,gauntcoeff,zatom,cellnew%rmeshnew,2*lmmax,1,nspin,ispin,'NS')
    !  bxcll(1:lmmax,1:lmmax,:) = (templl(1:lmmax,1:lmmax,:)-templl(lmmax+1:2*lmmax,lmmax+1:2*lmmax,:))/2.d0 
    !  write(17*iatom**2,'("  bxcll = ")')
    !  do ilm1=1,lmmax
    !    do ilm2=1,lmmax
    !      if(sum(abs(bxcll(ilm1,ilm2,:)))>1e-8) then
    !        write(17*iatom**2,'(2i4,1000es16.8)') ilm1, ilm2, bxcll(ilm1,ilm2,:)
    !      end if
    !    end do
    !  end do
    !  bs(:,:)=0.d0
    !  ! down/down block
    !  bs(1,1)=-bfield%c_xc(3)
    !  ! up/up block
    !  bs(2,2)= bfield%c_xc(3)
    !  ! down/up block
    !  bs(1,2)=(bfield%c_xc(1)+icompl*bfield%c_xc(2))
    !  ! up/down block
    !  bs(2,1)=(bfield%c_xc(1)-icompl*bfield%c_xc(2))
    !  if (ncoll==1) then                                        
    !    call rotatematrix(bs,theta, phi,1,'glob->loc')  
    !  end if                                                  
    !  if (mode=='transpose') then                            ! used for the left solution / Bxc is symmetric in LL' by construction
    !    temp2   = bs(1,2)
    !    bs(1,2) = bs(2,1)
    !    bs(2,1) = temp2
    !  elseif (mode=='1') then
    !  else
    !    stop'[bfield] mode not known'
    !  end if
    !  do i = 1,2
    !    do j = 1,2
    !      vpotll(lmstart(i):lmend(i),lmstart(j):lmend(j),:)=vpotll(lmstart(i):lmend(i),lmstart(j):lmend(j),:)-bs(i,j)*bxcll(:,:,:)
    !    end do
    !  end do
    !end if
    
    ! Add different contributions to the magnetic field
    bin(:) = 0.d0
    if ( t_inc%i_iteration>=bfield%itscf0 .and. t_inc%i_iteration<=bfield%itscf1 ) then ! preconvergence without adding the constraining field
      !if(myrank==master .and. debug) write(*,'("Adding bfield for iatom="i4)') iatom
      if(bfield%lbfield) then
        bin(:) = bfield%bfield(iatom,:)
      end if
      if( bfield%lbfield_constr ) then
        if( bfield%ibfield_constr == 0 .or. bfield%ibfield_constr == 1 ) then ! add homogenous constraining field here
            bin(:) = bin(:) + bfield%bfield_constr(iatom,:)
        end if
      end if
    end if

    if(bfield%lbfield_trans) then
      magdir(1) = sin(theta)*cos(phi)
      magdir(2) = sin(theta)*sin(phi)
      magdir(3) = cos(theta)
      bin(:) = bin(:) - magdir*dot_product(magdir(:),bin(:))
    end if
    if (t_inc%i_write>1 .and. mode=='1') write(1337,'("iatom, bin=",i4,3es16.8)') iatom, bin
    
    bs(:,:)=0.d0
    ! down/down block
    bs(1,1)=-bin(3)
    ! up/up block
    bs(2,2)= bin(3)
    ! down/up block
    bs(1,2)=(bin(1)+icompl*bin(2))
    ! up/down block
    bs(2,1)=(bin(1)-icompl*bin(2))
    if(myrank==master .and. debug) write(17*iatom**2,'(" bs(glob) =",8es16.8)') bs(1,1),bs(1,2),bs(2,1),bs(2,2)
    ! Just always rotate it even for collinear calculations 
    call rotatematrix(bs,theta, phi,1,1)!'glob->loc')  
    if(myrank==master .and. debug) write(17*iatom**2,'(" bs(loc) =",8es16.8)') bs(1,1),bs(1,2),bs(2,1),bs(2,2)
    if (mode=='transpose') then                            ! used for the left solution
      temp2   = bs(1,2)
      bs(1,2) = bs(2,1)
      bs(2,1) = temp2
    elseif (mode=='1') then
    else
      stop '[bfield] mode not known'
    end if
    
    !if(myrank==master .and. debug) then
    !  write(17*iatom**2,'("  vnspll_old = ")')
    !  do ilm1=1,lmmaxd
    !    do ilm2=1,lmmaxd
    !      if(sum(abs(real(vnspll0(ilm1,ilm2,:))))>1e-8) then
    !        write(17*iatom**2,'(2i4,1000es16.8)') ilm1, ilm2, vnspll0(ilm1,ilm2,:)
    !      end if
    !    end do
    !  end do
    !end if

    if(bfield%lbfield_mt) then
      irend = imt1
    else
      irend = irmdnew
    end if
    do i = 1,2
      do j = 1,2
        do ir=1,irend
          vnspll1(lmstart(i):lmend(i),lmstart(j):lmend(j),ir)=vnspll0(lmstart(i):lmend(i),lmstart(j):lmend(j),ir)-bs(i,j)*thetansll(1:lmmax,1:lmmax,ir)
        end do
        if(bfield%lbfield_mt) then !add non-spherical part of normal potential
          do ir=irend+1,irmdnew
            vnspll1(lmstart(i):lmend(i),lmstart(j):lmend(j),ir)=vnspll0(lmstart(i):lmend(i),lmstart(j):lmend(j),ir)
          end do
        end if
      end do
    end do
    
    !if(myrank==master .and. debug) then
    !  write(17*iatom**2,'("  vnspll_new = ")')
    !  do ilm1=1,lmmaxd
    !    do ilm2=1,lmmaxd
    !      if(sum(abs(real(vnspll1(ilm1,ilm2,:))))>1e-8) then
    !        write(17*iatom**2,'(2i4,1000es16.8)') ilm1, ilm2, vnspll1(ilm1,ilm2,:)
    !      end if
    !    end do
    !  end do
    !end if
 
    !close(17*iatom**2)
  end subroutine add_bfield
  
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

    !use :: global_variables, only: ntotd, nfund

    implicit none

    integer                     , intent (in)   :: lmax         ! Angular momentum cut-off
    integer                     , intent (in)   :: imt1         ! index muffin-tin radius
    integer                     , intent (in)   :: iend         ! Number of nonzero gaunt coefficients
    integer                     , intent (in)   :: irmdnew      ! number of radials point on the Cheby mesh
    integer, dimension (ncleb, 4), intent (in)  :: icleb !! Pointer array
    real (kind=dp), dimension (ncleb), intent (in)  :: cleb !! GAUNT coefficients (GAUNT)
    real (kind=dp)              , dimension (irmdnew, nfund)                    , intent (in)           :: thetasnew    ! shapefun on the Cheby mesh
    integer                     , dimension (1:(2*lmax+1)**2)                   , intent (in)           :: ifunm        ! pointer array for shapefun     ! Check index and dimensions!!!!!
    real(kind=dp)               , dimension((lmax+1)**2,(lmax+1)**2,irmdnew)    , intent (out)          :: thetansll    ! LL' expansion of the shapefunction 
    !------------------------------------------------------------------------------------
    real(kind=dp),parameter                                     :: rfpi=3.5449077018110318
    real(kind=dp),dimension(1:irmdnew,1:(2*lmax+1)**2)          :: shapefun_mod  
    real(kind=dp)                                               :: c0ll
    integer                                                     :: lmmax
    integer                                                     :: lmmax2
    integer                                                     :: ilm
    integer                                                     :: ifun
    integer                                                     :: ir
    integer                                                     :: lm1
    integer                                                     :: lm2
    integer                                                     :: lm3
    integer                                                     :: j
    logical                                                     :: debug=.False.
            
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





