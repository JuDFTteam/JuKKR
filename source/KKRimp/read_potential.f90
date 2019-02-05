module mod_read_potential
contains
subroutine read_potential(filename_pot,filename_shape,natom,lmaxd,zatom,lmaxatom,cell,shapefun,corestate,vpot,alat,nspin,ins)
   use nrtype
!    use configparams, only: ins, nspin
   use type_cell
   use type_corestate
   use type_shapefun
   use mod_config, only: config_testflag
   implicit none
!interfae variables
   character(len=*),intent(in)      ::  filename_pot                ! filename potential
   character(len=*),intent(in)      ::  filename_shape              ! filename shape functions
   integer,intent(in)               ::  natom
   integer,intent(in)               ::  lmaxd                      ! maximum lmax from atominfo
   real(kind=dp),intent(in)         ::  zatom(natom)               ! core charge from atom info file
   integer,intent(in)               ::  lmaxatom(natom)                ! lmax per atom
   type(cell_type),allocatable      ::  cell(:)
   type(shapefun_type),allocatable  :: shapefun(:)
   type(corestate_type),allocatable ::  corestate(:)                ! derived type containing all shape information 
   real(kind=dp),allocatable        ::  vpot(:,:,:,:)              ! input potential array
   real (kind=dp)                   ::  alat,alat_temp             ! lattice constant for scaling
   integer                          :: nspin
   integer                          :: ins
!local variabless
  integer                          ::  insguess                      ! error integer

! io files
   integer                          ::  ifile_pot,ifile_shape       ! file units
   integer                          ::  ifile_temp                  ! 
   integer                          ::  ierror                      ! error integer

! temporary variables
   character(len=100)               ::  cline
   real(kind=dp)                    ::  expa
   integer                          ::  itemp,isum,ipan1,nr,lmval
   integer                          ::  inonzerolm
   real(kind=dp)                    ::  vpottemp

! index variables
   integer                          ::  iatom                       ! running index, number of atoms
   integer                          ::  ispin, nspinguess           ! running index, spin index
   integer                          ::  icorestate                  ! running index core states
   integer                          ::  ir,iri, npot                ! radial running index

! geometry
   integer                          ::  nrmax2                      ! 2. value of nrmax for comparison
   integer                          ::  nrmaxd                      ! maximum number of rad. points for allocationg
   integer                          ::  nrcore1
   real(kind=dp)                    ::  a1,b1                       ! temporary a,b

   integer                          ::  nrmin_nsd                      ! 

   integer                          ::  natomshape
   integer                          ::  ipan, ifun
   integer                          ::  npand
   real(kind=dp),allocatable        ::  xrn(:,:),drn(:,:)
! type                                ::  rmesh_type

! 
!    integer,allocatable              ::  shape_index2lm(:,:)            ! shape function : array index -> lm value
!    integer,allocatable              ::  shape_lmused(:,:)          ! shape function : 1 -> lm-shape/=0,
!                                                                    !                  0 -> lm-shape=0 
!    integer,allocatable              ::  shape_lm2index(:,:)            ! shape function : lm value -> array index
! 
!    real(kind=dp),allocatable        ::  thetas(:,:,:)              ! shape function
! end type rmesh_type

!    type(shape_type)               ::  shape                      ! derived type containing all shape information 

! potential
   real(kind=dp)                    ::  zatompot                   ! core charge from potential file for comparison
   integer                          ::  lmmaxpot                   ! maximum (lm)-max from potential file

! core atoms
! type                                ::  corestate_type
    integer                          ::  ncorestated                ! maximum number of core states (for allocating)
!    integer,allocatable              ::  ncorestate(:)              ! number of core states for atom iatom
!    real(kind=dp),allocatable        ::  lcorestate(:,:,:)          ! angular momentum of core state
!    real(kind=dp),allocatable        ::  ecorestate(:,:,:)          ! energy of core state
! end type corestate_type
!    type(corestate_type),allocatable ::  corestate(:)                ! derived type containing all shape information 


! variables not used
   REAL(KIND=DP)                    ::  EFERMI,VMTZERO             ! Fermi Energy and muffin tin zero 
                                                                   ! (both not used)
   integer                           :: ilmcount,ival
IFILE_POT = 1000000
IFILE_SHAPE = 324238230 
IFILE_TEMP = 234569638
NRMIN_NSD=0

! ------------------------  read SHAPEFUN ---------------------

write(1337,*) '-------------------------------------------'
write(1337,*) '-----  POTENTIAL SUBROUTINE            ----'
write(1337,*) '-------------------------------------------'
write(1337,*) 'natom is ',natom
write(1337,*) 'ins is ',ins
write(1337,*) '-------------------------------------------'

ALLOCATE (CELL(NATOM),STAT=IERROR)
IF (IERROR/=0) stop '[read_potential] Error while allocating arrays'

  ALLOCATE (SHAPEFUN(NATOM),STAT=IERROR)
  IF (IERROR/=0) stop '[read_potential] Error while allocating arrays'

! if (ins==1) then 
  ALLOCATE (CORESTATE(NATOM),STAT=IERROR)
  IF (IERROR/=0) stop '[read_potential] Error while allocating arrays'
! end if



if (ins==1) then 
   CALL PRE_READ_SHAPE(FILENAME_SHAPE,NATOM,NPAND,SHAPEFUN(1)%NRSHAPED,SHAPEFUN(1)%NLMSHAPED)
   SHAPEFUN(:)%NRSHAPED  =  SHAPEFUN(1)%NRSHAPED
   SHAPEFUN(:)%NLMSHAPED =  SHAPEFUN(1)%NLMSHAPED
else
   NPAND=1;SHAPEFUN(:)%NRSHAPED=1;SHAPEFUN(:)%NLMSHAPED=1
   do iatom=1,natom
     allocate(SHAPEFUN(IATOM)%THETAS(1,1))
   end do
end if
!                         >          >     <      <         <
write(1337,*) ''
write(1337,*) '-------------------------------------------'
write(1337,*) '-----  pre read shape functions        ----'
write(1337,*) '-------------------------------------------'
write(1337,'(A,I0,A,I0,A,I0)') 'NPAND= ',NPAND,' NRSHAPED=',SHAPEFUN(1)%NRSHAPED,' NLMSHAPED=',SHAPEFUN(1)%NLMSHAPED

if (ins==0) then
  write(1337,*) '-----  skip shape functions            ----'
  DO IATOM=1,NATOM
    CELL(IATOM)%NPAN=1
  END DO !IATOM
elseif (ins==1) then 
  OPEN(UNIT=IFILE_SHAPE, FILE=FILENAME_SHAPE, STATUS='old', IOSTAT=IERROR)
  IF (IERROR/=0) THEN
    WRITE(*,*) '[read_potential] shape file does not exist'
    STOP
  END IF
  
  READ(UNIT=IFILE_SHAPE,FMT=*) NATOMSHAPE  ! number of shape functions
  READ(UNIT=IFILE_SHAPE,FMT=*) CLINE       ! scaling variables not used
! write(*,*) 't'
! write(*,*) natomshape,natom
! stop
  IF (NATOMSHAPE/=NATOM) THEN
    WRITE(*,*) 'NATOM=',NATOM, 'NATOMSHAPE=',NATOMSHAPE
    STOP '[read_potential] number of shape potentials differs from NATOM'
  END IF
  
  
  
  DO IATOM=1,NATOM
  allocate(CELL(IATOM)%NMESHPAN(NPAND))
  CELL(IATOM)%NMESHPAN=0
  cell(iatom)%npand=npand
  ! allocate(NPAN(NATOM),NRSHAPE(NATOM),NMESHPAN(NPAND,NATOM),NLMSHAPE(NATOM))
  END DO !IATOM
  
  
  allocate(XRN(SHAPEFUN(1)%NRSHAPED,NATOM),DRN(SHAPEFUN(1)%NRSHAPED,NATOM))
  DO IATOM=1,NATOM
    allocate(SHAPEFUN(IATOM)%INDEX2LM(SHAPEFUN(IATOM)%NLMSHAPED),SHAPEFUN(IATOM)%LMUSED((4*LMAXD+1)**2), &
             SHAPEFUN(IATOM)%LM2INDEX((4*LMAXD+1)**2))
    allocate(SHAPEFUN(IATOM)%THETAS(SHAPEFUN(IATOM)%NRSHAPED,SHAPEFUN(IATOM)%NLMSHAPED))
    SHAPEFUN(IATOM)%INDEX2LM(:)=0
    SHAPEFUN(IATOM)%LM2INDEX(:)=0
    SHAPEFUN(IATOM)%LMUSED(:)=0
  END DO!IATOM
  
  DO IATOM=1,NATOM
    read(unit=ifile_shape,fmt=*) CELL(IATOM)%NPAN,SHAPEFUN(IATOM)%NRSHAPE                ! number of panels
    CELL(IATOM)%NPAN =  CELL(IATOM)%NPAN+1                                            ! extend the number of panels :
                                                                          ! old: panel1,panel2,...panelN
                                                                          ! new: corepanel,panel2,panel3,...panelN+1
    read(unit=ifile_shape,fmt=*) (CELL(IATOM)%NMESHPAN(IPAN),IPAN=2,CELL(IATOM)%NPAN) ! IPAN=1 is later on assumed to be
    CELL(IATOM)%NMESHPAN(1)=-1                                                   ! the 1st panel is assumed to be inside RCORE
    read(unit=ifile_shape,FMT='(4d20.12)') (XRN(IR,IATOM), &               ! shape function mesh
                                            DRN(IR,IATOM), &               ! and its derivative
                                            IR=1, SHAPEFUN(IATOM)%NRSHAPE)
    read(unit=ifile_shape,fmt=*) SHAPEFUN(IATOM)%NLMSHAPE                           ! max number of (l,m)
    ilmcount=0
    do ifun = 1,SHAPEFUN(IATOM)%NLMSHAPE
      read(unit=ifile_shape,fmt=*) LMVAL
      if (LMVAL<=(4*LMAXATOM(IATOM)+1)**2) then
        ilmcount=ilmcount+1
        SHAPEFUN(IATOM)%INDEX2LM(IFUN)=LMVAL
        SHAPEFUN(IATOM)%LMUSED(LMVAL)=1
        SHAPEFUN(IATOM)%LM2INDEX(LMVAL)= IFUN
      end if
      read(unit=ifile_shape,FMT='(4d20.12)') (SHAPEFUN(IATOM)%THETAS(IR,IFUN),IR=1, SHAPEFUN(IATOM)%NRSHAPE)
    end do !ifun
    if (ilmcount<SHAPEFUN(IATOM)%NLMSHAPE) then
      write(1337,*) 'too many lm-components are defined in the shape file'
      write(1337,*) 'I will cut the lm-components of atom',iatom
      write(1337,*) 'according to its lmax of ',lmaxatom(iatom)
      write(1337,*) ilmcount,'compoentens are used instead of ',SHAPEFUN(IATOM)%NLMSHAPE
      SHAPEFUN(IATOM)%NLMSHAPE = ilmcount
    else if (ilmcount>SHAPEFUN(IATOM)%NLMSHAPE) then
     write(*,*) 'shape file for atom',iatom,'has not sufficient lm-'
     write(*,*) 'components for its lmax of',lmaxatom(iatom)
     stop
    end if
  END DO !IATOM
  close(ifile_shape)
  
  write(1337,*) ''
  write(1337,*) '-------------------------------------------'
  write(1337,*) '-----      shape functions             ----'
  write(1337,*) '-------------------------------------------'
  write(1337,*) '         IATOM      NPAN  NRSHAPE    NMESHPAN(1..NPAN) '
  DO IATOM=1,NATOM
    write(1337,'(I8,I8,I8,100I8)') IATOM, CELL(IATOM)%NPAN,SHAPEFUN(IATOM)%NRSHAPE, CELL(IATOM)%NMESHPAN(:CELL(IATOM)%NPAN)
  END DO
else !ins
  stop 'INS/=1/0'
end if

! --------------------------------------------------------------



! ------------------------  read POTENTIAL ---------------------


! first read in some properties of the potential file
CALL PRE_READ_POTENTIAL(FILENAME_POT,INSGUESS,NCORESTATED,NPOT,NRMAXD)
 corestate(:)%NCORESTATED=NCORESTATED
do iatom = 1,natom
 cell(iatom)%nrmaxd=nrmaxd
end do
write(1337,*) ''
write(1337,*) '-------------------------------------------'
write(1337,*) '-----      pre read potential          ----'
write(1337,*) '-------------------------------------------'
write(1337,'(A,I0,A,I0,A,I0,A,I0,A,I0)') 'INS(guess)=',INSGUESS,' NCORESTATED=',NCORESTATED,' NPOT=',NPOT,' NRMAXD=',NRMAXD


write(1337,*) '-------------------------------------------'
if (natom==npot) then
   write(1337,*) 'According potential file nspin=1 is used'
   nspinguess=1
elseif (natom*2==npot) then
   write(1337,*) 'According potential file nspin=2 is used'
   nspinguess=2
else
   write(*,*) 'natom =',natom,'npot =',npot
   stop '[read_potential] number of atoms and potentials inconsistent'
end if
write(1337,*) '-------------------------------------------'
if (nspinguess/=nspin) then 
  write(*,*) nspinguess, nspin
  stop '[read_potential] conflict between potential file and NSPIN value'
end if

! open the potential file
OPEN(UNIT=ifile_pot, FILE=filename_pot, status='old', iostat=ierror)
if (ierror/=0) then
   write(*,*) '[read_potential] potential file does not exist'
   STOP
end if


! ALLOCATE (CELL(NATOM)  !RMT(NATOM),RCORE(NATOM),NRCORE(NATOM), &
          !RMAX(NATOM),NRMAX(NATOM)
!           ,  &
          !LOGPARAMS(2,NATOM),NRNS(NATOM),NRMIN_NS(NATOM),  & 
!           ,STAT=IERROR)
do iatom=1,natom
ALLOCATE (corestate(iatom)%LCORE(NCORESTATED,NSPIN), &
          corestate(iatom)%ECORE(NCORESTATED,NSPIN),STAT=IERROR)
          IF (IERROR/=0) STOP'[read_potential] Error while allocationg arrays LCORESTATE,ECORESTATE'
end do !iatom


! IF     (INS==0) THEN
!    ALLOCATE (VPOT(NRMAXD,1,NSPIN,NATOM),STAT=IERROR)
!    IF (IERROR/=0) STOP '[read_potential] Error allocating VPOT'
!    VPOT(:,:,:,:)=0.0D0
! ELSEIF (INS==1) THEN
   ALLOCATE (VPOT(NRMAXD,(2*LMAXD+1)**2,NSPIN, NATOM),STAT=IERROR)
   VPOT(:,:,:,:)=0.0D0
   IF (IERROR/=0) STOP '[read_potential] Error allocating VPOT'
! ELSE
!    STOP '[read_potential] INS/=0, 1'
! END IF

DO IATOM= 1, NATOM
  ALLOCATE (CELL(IATOM)%RMESH(NRMAXD), CELL(IATOM)%DRMESHDI(NRMAXD), CELL(IATOM)%DRMESHOR(NRMAXD), &
            CELL(IATOM)%NRCUT(0:NPAND),STAT=IERROR)
            IF (IERROR/=0) STOP'[read_potential] Error while allocationg arrays LCORESTATE,ECORESTATE'
  CELL(IATOM)%NRCUT=0
END DO !IATOM

DO IATOM= 1, NATOM
   DO ISPIN=1,NSPIN

        READ (IFILE_POT,FMT='(A28)') CELL(IATOM)%VPOT_NAME(ISPIN)
!       IF (ISPIN==2) THEN
!         READ (IFILE_POT,FMT='(A28)') VPOT_NAME_TEMP
!         IF (CELL(IATOM)%VPOT_NAME(1:16,1)/=VPOT_NAME_TEMP(1:16)) THEN
!           STOP '[read_potential] Title in the potentialcard for spin up/down do not match'
!         END IF
!       END IF
!       WRITE(*,*) 'ghdfgfdg', CELL(IATOM)%VPOT_NAME

!       READ(UNIT=ifile_pot,FMT=*) CLINE                               ! skip line
      READ(UNIT=ifile_pot,FMT=*) CELL(IATOM)%RMT, ALAT_TEMP, CELL(IATOM)%RCORE     ! muffin tin and core radius
!        write(*,*) alat, alat_temp
      IF (ABS(ALAT-ALAT_TEMP)>10D-14) then
         write(*,*) alat, alat_temp
         STOP '[read_potential] error ALAT value in potential is does not match'
      end if
      READ(UNIT=ifile_pot,FMT=*) ZATOMPOT                            ! atomic charge
      IF (ZATOM(IATOM)/=ZATOMPOT) THEN
         WRITE(*,*) 'ZATOM(',IATOM,')  is   ',ZATOM(IATOM), ' ZATOMPOT  is',ZATOMPOT
         WRITE(*,*) 'The core charge for atom',iatom,'is not consistent with the atominfo file'
         STOP
      END IF
      READ(UNIT=ifile_pot,FMT=*) CELL(IATOM)%RMAX,EFERMI,VMTZERO          ! wigner-seitz radius, efermi(not used)
                                                                    ! muffin tin zero(not used)
      READ(UNIT=ifile_pot,FMT=*) CELL(IATOM)%NRMAX                         ! number of spherical mesh points
      IF (CELL(IATOM)%NRMAX > NRMAXD) stop '[read_potential] NRMAXD is to small '

      READ(UNIT=ifile_pot,FMT=*) CELL(IATOM)%LOGPARAMS(1), &
                                CELL(IATOM)%LOGPARAMS(2)               ! logarithmic mesh parameters
      READ(UNIT=ifile_pot,FMT=*) corestate(IATOM)%NCORE, ITEMP    ! number of core states, file format (1=new)


      DO ICORESTATE=1,corestate(iatom)%NCORE
         READ(UNIT=ifile_pot,FMT=*) corestate(iatom)%LCORE(ICORESTATE,ISPIN), &
                                    corestate(iatom)%ECORE(ICORESTATE,ISPIN)
      END DO


      IF (INS==0 .AND. INSGUESS==0) THEN
         READ (UNIT=ifile_pot,FMT=*) (VPOT(IR,1,ISPIN, IATOM),IR=1,CELL(IATOM)%NRMAX)
      ELSEIF (INS==0 .AND. INSGUESS==1) THEN
      write(*,*) 'INS=0 but potential file is spherical'
      stop
      ELSEIF (INS==1 .AND. INSGUESS==1) THEN
         READ(UNIT=ifile_pot,FMT=*) NRMAX2,CELL(IATOM)%NRNS, &  ! # mesh points, # non sph. mesh points
                                   LMMAXPOT, INONZEROLM             ! maximum (l,m)-components, 
                                                                 ! INONZEROLM=1 => read only non zero values 

         IF (NRMAX2/=CELL(IATOM)%NRMAX) THEN
           stop '[read_potential] the two values for IRMT in potential file are not consistent'
         END IF

         IF (INONZEROLM/=1) THEN 
            STOP '[read_potential] potential file has old format given all lm-component. Use new format'
         END IF

         IF (LMMAXPOT<(2*LMAXATOM(IATOM)+1)**2) then
!             write(*,*) LMMAXPOT-(2*LMAXATOM(IATOM)+1)**2
            write(1337,*) '[read_potential] WARNING: the lmax value in potential file is too low'
         END IF

         CELL(IATOM)%NRMIN_NS = NRMAX2-CELL(IATOM)%NRNS
         IF (CELL(IATOM)%NRMIN_NS>NRMIN_NSD) NRMIN_NSD=CELL(IATOM)%NRMIN_NS


!         READ (UNIT=ifile_pot,FMT=*) (VPOT(IR,1,ISPIN, IATOM),IR=1,NRMAX2)
         READ (UNIT=ifile_pot,FMT='(1p,4d20.13)') (VPOT(IR,1,ISPIN, IATOM),IR=1,NRMAX2)

         IF (LMMAXPOT>1) THEN
!             ival=0
            DO ival=2,LMMAXPOT
!                write(*,*) ival
!                ival=ival+1
!                IF (ival==LMMAXPOT) EXIT !hack py D. Bauer 
               READ (unit=ifile_pot,fmt=*) LMVAL
               IF (LMVAL==1) EXIT
               IF (LMVAL<=(2*LMAXATOM(IATOM)+1)**2) THEN
                  IF (LMVAL<1) stop '[read_potential] LMVAL in potential file < 1'
                  READ (UNIT=ifile_pot,FMT='(4d20.12)') (VPOT(IR,LMVAL,ISPIN, IATOM),IR=CELL(IATOM)%NRMIN_NS,NRMAX2)
               ELSE
                  READ (UNIT=ifile_pot,FMT='(4d20.12)') (VPOTTEMP,IR=CELL(IATOM)%NRMIN_NS,NRMAX2)
               END IF
            END DO
         END IF
      ELSE
          STOP '[read_potential] INS/=1/0 or INSGUESS/=INS'
      END IF


   IF (INS==1) THEN
      !IF (abs(ALAT*XRN(1,IATOM)-CELL(IATOM)%RCORE)>10e-10)  THEN
      IF (abs(ALAT*XRN(1,IATOM)-CELL(IATOM)%RCORE)>10e-7)  THEN
         write(*,*) IATOM,ALAT*XRN(1,IATOM),CELL(IATOM)%RCORE
         stop '[read_potential] RMT inconsistency betweeen shape and potential file'
      END IF

      ! phivos
     CELL(IATOM)%RCORE=ALAT*XRN(1,IATOM)

      NRCORE1 = NINT(LOG(CELL(IATOM)%RCORE/CELL(IATOM)%LOGPARAMS(2)+1.0_DP)/CELL(IATOM)%LOGPARAMS(1)) + 1

      !---> for proper core treatment imt must be odd
      !     shift potential by one mesh point if imt is even
      !

      IF (MOD(NRCORE1,2)==0) THEN
         STOP 'ICORE is not odd'
         !       ICORE1 = ICORE1 + 1
         !       DO IR = ICORE1,2,-1
         !          VPOT(IR,1,ISPIN,IATOM) = VPOT(IR-1,1,ISPIN,IATOM)
         !       END DO
      END IF

      CELL(IATOM)%NRCORE = NRCORE1
!        phivos: not sure if i should include this:
       CELL(IATOM)%LOGPARAMS(2) = CELL(IATOM)%RCORE / & 
                                        (EXP(CELL(IATOM)%LOGPARAMS(1)*DBLE(NRCORE1-1))-1.0D0)
   ELSEIF (INS==0) THEN
      CELL(IATOM)%NRCORE = CELL(IATOM)%NRMAX ! susc ! line added to fix bug in ASA calculation, Julen and Benedikt 2014/08
   ELSE
      STOP 'INS/=1,0'
   END IF !INS==1

   !
   !---> generate radial mesh - potential only is stored in potential card
   !
   A1 = CELL(IATOM)%LOGPARAMS(1)
   B1 = CELL(IATOM)%LOGPARAMS(2)
   CELL(IATOM)%RMESH(1) = 0.0D0
   CELL(IATOM)%DRMESHDI(1) = A1*B1
   DO IR = 2,CELL(IATOM)%NRMAX
     EXPA = EXP(A1*DBLE(IR-1))
     CELL(IATOM)%RMESH(IR) = B1* (EXPA-1.0D0)
     CELL(IATOM)%DRMESHDI(IR) = A1*B1*EXPA
     CELL(IATOM)%DRMESHOR(IR) = A1/ (1.0D0-1.0D0/EXPA)
   END DO !IR
  
   !   c
   !   c---> fill cell-type depending mesh points in the non-muffin-tin-region
   !   c
   IF (INS==1) THEN
     DO IRI = 1,SHAPEFUN(IATOM)%NRSHAPE
         IR = IRI + CELL(IATOM)%NRCORE
         CELL(IATOM)%RMESH(IR)    = ALAT*XRN(IRI,IATOM)
         CELL(IATOM)%DRMESHDI(IR) = ALAT*DRN(IRI,IATOM)
         CELL(IATOM)%DRMESHOR(IR) = CELL(IATOM)%DRMESHDI(IR)/CELL(IATOM)%RMESH(IR)
     END DO
   END IF
   

   IF (INS==1) THEN
     CELL(IATOM)%NRCUT(1) = CELL(IATOM)%NRCORE
     ISUM = CELL(IATOM)%NRCORE
     DO IPAN1 = 2,CELL(IATOM)%NPAN
         ISUM = ISUM + CELL(IATOM)%NMESHPAN(IPAN1)
         CELL(IATOM)%NRCUT(IPAN1) = ISUM
     END DO
     NR = ISUM
   ELSE
     NR = CELL(IATOM)%NRMAX
     CELL(IATOM)%NRCUT(1) = CELL(IATOM)%NRMAX
   END IF

   END DO !ISPIN
END DO !IATOM


write(1337,*) '-------------------------------------------'
write(1337,*) '-----  NRCUT                            ----'
write(1337,*) '-------------------------------------------'
write(1337,*) 'IATOM     NRCUT'
DO IATOM=1,NATOM
   write(1337,'(1000I5)') IATOM, CELL(IATOM)%NRCUT(:)
END DO!IATOM


write(1337,*) '-------------------------------------------'
write(1337,*) '-----  Atom stuff                      ----'
write(1337,*) '-------------------------------------------'
write(1337,'(2A)') ' IATOM   RCORE   NRCORE   ',&
                   'RMT   RMAX   NRMAX   RMIN_NS   NRMIN_NS'
DO IATOM=1,NATOM
  IF (INS==1) THEN
   write(1337,'(I5,F8.3,I6,F8.3,F8.3,I5,F8.3,I)') IATOM,CELL(IATOM)%RCORE,CELL(IATOM)%NRCORE, &
                                   CELL(IATOM)%RMT,CELL(IATOM)%RMAX,CELL(IATOM)%NRMAX,CELL(IATOM)%RMESH(CELL(IATOM)%NRMIN_NS),CELL(IATOM)%NRMIN_NS
  ELSE
   write(1337,'(I5,F8.3,I6,F8.3,F8.3,I5,F8.3,I)') IATOM,CELL(IATOM)%RCORE,CELL(IATOM)%NRCORE, &
                                   CELL(IATOM)%RMT,CELL(IATOM)%RMAX,CELL(IATOM)%NRMAX
  END IF
END DO !NATOM

write(1337,*) ''
write(1337,*) '-------------------------------------------'
write(1337,*) '-----  Parameters for the radial mesh  ----'
write(1337,*) '-------------------------------------------'
write(1337,*) '               A                       B           '
DO IATOM=1,NATOM
write(1337,*) IATOM,CELL(IATOM)%LOGPARAMS(:)
END DO 

IF(CONFIG_TESTFLAG('write_rmesh')) THEN
  OPEN(UNIT=IFILE_TEMP, FILE='test_rmesh', iostat=ierror)
  do IR=1,NRMAXD
    !    do IATOM = 1, NATOM
    write(IFILE_TEMP,'(I0,1000g24.16)') IR,(CELL(IATOM)%rmesh(IR),IATOM=1,NATOM)
    !    end do !IR
  end do !IATOM
  close(IFILE_TEMP)
END IF

IF(CONFIG_TESTFLAG('readpot_write')) THEN
  OPEN(UNIT=IFILE_TEMP, FILE='test_potential', iostat=ierror)
  ! write(*,*) 'POT'
  if (ins==1) then
    do IATOM = 1, NATOM
      do ISPIN = 1, NSPIN
          do LMVAL = 1, (2*LMAXATOM(IATOM)+1)**2
            write(unit=IFILE_TEMP,fmt='(3(A,I0))',advance='no') 'IATOM=',IATOM,'ISPIN=',ISPIN,'LMVAL=',LMVAL
            do IR = 1, CELL(IATOM)%NRMAX
                write(unit=IFILE_TEMP,fmt='(5000g24.16)',advance='no') VPOT(IR,LMVAL,ISPIN,IATOM)
            end do
            write(unit=IFILE_TEMP,fmt=*) ''
          end do
      end do
    end do
  else
    do IATOM = 1, NATOM
      do ISPIN = 1, NSPIN
        write(unit=IFILE_TEMP,fmt='(3I0)',advance='no') IATOM,ISPIN,LMVAL
        do IR = 1, CELL(IATOM)%NRMAX
            write(unit=IFILE_TEMP,fmt='(5000g24.16)',advance='no') VPOT(IR,1,ISPIN,IATOM)
        end do
        write(unit=IFILE_TEMP,fmt=*) ''
      end do
    end do
  end if
  close(IFILE_TEMP)
END IF !(CONFIG_TESTFLAG('readpot_write')) THEN


END SUBROUTINE read_potential


subroutine pre_read_shape(filename_shape,NATOM,NPAND,NRSHAPED,NLMSHAPED)
    use nrtype
    implicit none
    character(len=80)    :: cline
    character(len=*)     :: filename_shape
    integer              :: NATOM,NATOMSHAPE,NPAN,NPAND,NRSHAPE,NRSHAPED,NLMSHAPED, NLMSHAPE,ios,ifile_shape,ierror
    integer              :: iatom,ir,ifun,itemp,ipan1
    real(kind=dp)        :: temp1,temp2

ifile_shape=1000000
npand=0;nrshaped=0;nlmshaped=0

open(unit=ifile_shape, file=filename_shape, status='old', iostat=ierror)
if (ierror/=0) then
   write(*,*) '[read_potential] shape file does not exist'
   stop
end if

read(unit=ifile_shape,fmt=*) natomshape

if(natomshape/=natom) then 
   write(*,*) 'natom=',natom, 'natomshape=',natomshape
   stop '[read_potential] number of shapefun differs from number of atoms'
end if

read(unit=ifile_shape,fmt=*) cline       ! scaling variables not used
do iatom=1,natom
   read(unit=ifile_shape,fmt=*) npan,nrshape
   npan=npan+1
   if (npan>npand) npand=npan
   if (nrshape>nrshaped) nrshaped=nrshape
   read (ifile_shape,fmt=*) (itemp,ipan1=2,npan)  !in  2/3/2012
!    read (19,fmt='(16i5)') (itemp,ipan1=2,npan)  !in  2/3/2012
!  read(unit=ifile_shape,fmt=*) cline           !out 2/3/2012
   read(unit=ifile_shape,fmt='(4d20.12)') (temp1,temp2,ir=1, nrshape)
   read(unit=ifile_shape,fmt=*) nlmshape
   if (nlmshape>nlmshaped) nlmshaped=nlmshape
   do ifun = 1,nlmshape
      read(unit=ifile_shape,fmt=*) cline
      read(unit=ifile_shape,fmt='(4d20.12)',iostat=ios) (temp1,ir=1, nrshape)
   end do 
end do !iatom


close(ifile_shape)
! write(*,*) filename_shape,natom,npand,nrshaped,nlmshaped
end subroutine pre_read_shape

subroutine pre_read_potential(filename_pot,ins,ncorestated,npot,NRMAXD)
use nrtype
implicit none

   character(len=*)      ::  filename_pot
   integer               ::  ins,lmmaxd,ncorestated,npot,nRMAXd,nrmax1
   integer               ::  ios,ifile_pot,ierror
   character(len=200)    ::  string1
   integer               ::  iline,iline2,itemp,i4
   real(kind=dp)         ::  temp4(4)

ins    = -1
lmmaxd = 0
ncorestated = 0
npot=0
nRMAXd=0

ifile_pot=1000000
OPEN(UNIT=ifile_pot, FILE=filename_pot, status='old', iostat=ierror)
if (ierror/=0) then
   write(*,*) '[read_potential] potential file does not exist'
   STOP
end if


ios=0
do while (ios/=-1)
   read(unit=ifile_pot,fmt='(A)', iostat=ios) string1

   if (ios==-1) cycle

   if (string1(6:15)=='POTENTIAL') then 
      npot = npot+1
!       write(*,*) string1
      do iline=1,3
         read(unit=ifile_pot,fmt='(A)', iostat=ios) string1
      end do

      read(unit=ifile_pot,fmt=*, iostat=ios) nrmax1

      do iline=1,2
         read(unit=ifile_pot,fmt='(A)', iostat=ios) string1
      end do


      read(string1,*) itemp
      if (itemp>ncorestated) ncorestated = itemp
      do iline2=1,itemp+1
         read(unit=ifile_pot,fmt='(A)', iostat=ios) string1
      end do !iline2
      read(string1,*) temp4
      if (ins==-1) then 
         do i4=1,4
            ins=1
            if (abs(int(temp4(i4))-temp4(i4))>10D-10) then
               INS=0
            end if
         end do
      else 
         do i4=1,4
            if (abs(int(temp4(i4))-temp4(i4))>10D-10 .and. ins==1)  then
               stop '[pre_read_potential] INS=1/0 are mixed'
            end if
         end do
      end if

!      if nrmax1/=temp4(1) stop '[pre_read_potential] the two NRMAX in potential file are inconsistant'

!      write(*,*) 'nrmaxd',nrmaxd,nrmax1
      if (nRMAXd<temp4(1) ) nRMAXd=nrmax1 !temp4(1)
!       if (lmmaxd<temp4(3) .and. ins/=0) lmmaxd=temp4(3)
   end if!keyword2='POTENTIAL'

end do !ios/=-1
! stop
close(1000000)
! write(*,*) 'PRE ',filename_pot,ins,ncorestated,npot,NRMAXD
end subroutine pre_read_potential

end module mod_read_potential
