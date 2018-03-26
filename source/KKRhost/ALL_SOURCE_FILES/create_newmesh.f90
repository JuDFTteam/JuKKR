!-------------------------------------------------------------------------------
! MODULE: mod_create_newmesh
!> @brief Module for the creation of the integration grid for the new solver
!-------------------------------------------------------------------------------
module mod_create_newmesh

   use Constants
   use Profiling

   implicit none

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: CREATE_NEWMESH
   !> @brief Creation of the integration grid for the new solver
   !> @note changed interface to get rid of inc.p and to be able to use i
   !> create_newmesh in tmatimp routine for GREENIMP option
   !> this is the list of  array dimensions previously importted from inc.p
   !----------------------------------------------------------------------------
   subroutine CREATE_NEWMESH(NATYP,LMAX,LPOT,IRM,IRNSD,IPAND,IRID,NTOTD,   &
      NFUND,NCHEB,IRMDNEW,NSPIN,R,IRMIN,IPAN,IRCUT,R_LOG,NPAN_LOG,NPAN_EQ, &
      NPAN_LOGNEW,NPAN_EQNEW,NPAN_TOT,RNEW,RPAN_INTERVALL,IPAN_INTERVALL,  &
      NCELLD,NTCELL,THETAS,THETASNEW) !< optional arguments

      implicit none

      ! .. Input variables
      integer, intent(in) :: IRM       !< Maximum number of radial points
      integer, intent(in) :: IRID      !< Shape functions parameters in non-spherical part
      integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
      integer, intent(in) :: LPOT      !< Maximum l component in potential expansion
      integer, intent(in) :: NSPIN     !< Counter for spin directions
      integer, intent(in) :: NATYP     !< Number of kinds of atoms in unit cell
      integer, intent(in) :: IRNSD
      integer, intent(in) :: IPAND     !< Number of panels in non-spherical part
      integer, intent(in) :: NTOTD
      integer, intent(in) :: NFUND     !< Shape functions parameters in non-spherical part
      integer, intent(in) :: NCHEB     !< Number of Chebychev pannels for the new solver
      integer, intent(in) :: NCELLD    !< Number of cells (shapes) in non-spherical part
      integer, intent(in) :: IRMDNEW
      integer, intent(in) :: NPAN_EQ   !< Variables for the pannels for the new solver
      integer, intent(in) :: NPAN_LOG  !< Variables for the pannels for the new solver
      double precision, intent(in) :: R_LOG
      integer, dimension(NATYP), intent(in)     :: IPAN   !< Number of panels in non-MT-region
      integer, dimension(NATYP), intent(in)     :: IRMIN  !< Max R for spherical treatment
      integer, dimension(0:IPAND,NATYP), intent(in) :: IRCUT    !< R points of panel borders
      double precision, dimension(IRM,NATYP), intent(in) :: R   !< Radial mesh ( in units a Bohr)
      ! .. Input/Output variables
      integer, dimension(NATYP), intent(inout)  :: NPAN_TOT
      integer, dimension(NATYP), intent(inout)  :: NPAN_EQNEW
      integer, dimension(NATYP), intent(inout)  :: NPAN_LOGNEW
      integer, dimension(0:NTOTD,NATYP), intent(inout) :: IPAN_INTERVALL
      double precision, dimension(IRMDNEW,NATYP), intent(inout):: RNEW
      double precision, dimension(0:NTOTD,NATYP), intent(inout) :: RPAN_INTERVALL
      ! .. Optional arguments, do interpolation when given
      integer, dimension(NATYP), intent(in), optional :: NTCELL   !< Index for WS cell
      double precision, dimension(IRID,NFUND,NCELLD), intent(in), optional :: THETAS   !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
      double precision, dimension(NTOTD*(NCHEB+1),NFUND,NCELLD), intent(inout), optional :: THETASNEW

      ! .. Local variables
      integer :: NPAN_INST,i_stat
      integer :: I1,IPOT,IPOTM,IR2,IP,ICELL
      integer :: IMIN,IMAX,IMINNEW,IMAXNEW,LM1
      integer :: ISHIFT,ILOGPANSHIFT,ILINPANSHIFT,NPAN_LOGTEMP
      double precision :: FAC
      double precision :: RMIN,RMAX,RVAL
      ! .. Allocatable arrays
      double precision, dimension(:,:,:), allocatable :: THETASIN
      ! .. Assignment of values to parameters
      parameter (FAC=2d0)


      ! checks for optional arguments
      if(present(NTCELL).and..not.present(THETAS).or..not.present(THETASNEW))) then
         write(*,*) 'Error in create_newmesh:'
         write(*,*) 'List of optional arguments not complete'
         stop
      end if

      ! allocations
      if( present(NTCELL)) then
         allocate(THETASIN(IRID,NFUND,NCELLD),stat=i_stat)
         call memocc(i_stat,product(shape(THETASIN))*kind(THETASIN),'THETASIN','CREATE_NEWMESH')
         THETASIN=0.0D0
      endif
      if( present(NTCELL) ) THETASNEW=0.0D0

      IPOTM=0

      do I1 = 1,NATYP

         IPOT=NSPIN*(I1-1)+1
         NPAN_INST= IPAN(I1)-1
         NPAN_TOT(I1)= NPAN_LOG+NPAN_EQ+NPAN_INST

         ! log panel
         RMIN=R(2,I1)
         RMAX=R_LOG
         RVAL=0d0
         ISHIFT=0
         if (R_LOG.GT.R(IRMIN(I1),I1)) then
            ILOGPANSHIFT=1
            ILINPANSHIFT=0
         else
            ILOGPANSHIFT=0
            ILINPANSHIFT=1
         endif

         if (ILINPANSHIFT.EQ.1) then
            stop 'non-spherical part of the potential needs to be inside the log panel'
         endif

         do IP=0,NPAN_LOG-ILOGPANSHIFT
            RVAL=(FAC**IP-1d0)/(FAC**(NPAN_LOG-ILOGPANSHIFT)-1d0)
            RPAN_INTERVALL(IP+ISHIFT,I1)= RMIN+RVAL*(RMAX-RMIN)
            IPAN_INTERVALL(IP+ISHIFT,I1)= (IP+ISHIFT)*(NCHEB+1)
            if (ISHIFT.EQ.0.AND.RPAN_INTERVALL(IP,I1).GT.R(IRMIN(I1),I1)) then
               ISHIFT=1
               NPAN_LOGTEMP=IP
               RPAN_INTERVALL(IP+1,I1)=RPAN_INTERVALL(IP,I1)
               IPAN_INTERVALL(IP+1,I1)=(IP+ISHIFT)*(NCHEB+1)
               RPAN_INTERVALL(IP,I1)=R(IRMIN(I1),I1)
               IPAN_INTERVALL(IP,I1)=IP*(NCHEB+1)
            endif
         enddo ! NPAN_LOG

         ! equivalent panel
         ISHIFT=0
         RMIN=R_LOG
         RMAX=R(IRCUT(1,I1),I1)
         do IP=0,NPAN_EQ-ILINPANSHIFT
            RPAN_INTERVALL(IP+ISHIFT+NPAN_LOG,I1)=RMIN+IP*(RMAX-RMIN)/&
               (NPAN_EQ-ILINPANSHIFT)
            IPAN_INTERVALL(IP+ISHIFT+NPAN_LOG,I1)=(NPAN_LOG+IP+ISHIFT)*(NCHEB+1)
         enddo ! NPAN_EQ

         ! intersection zone
         do IP=1,NPAN_INST
            RPAN_INTERVALL(NPAN_LOG+NPAN_EQ+IP,I1)=R(IRCUT(IP+1,I1),I1)
            IPAN_INTERVALL(NPAN_LOG+NPAN_EQ+IP,I1)=(NPAN_LOG+NPAN_EQ+IP)*(NCHEB+1)
         enddo ! NPAN_INST

         NPAN_EQNEW(I1)=NPAN_EQ+NPAN_LOG-NPAN_LOGTEMP
         NPAN_LOGNEW(I1)=NPAN_LOGTEMP

         call CHEBMESH(NPAN_TOT(I1),NCHEB,RPAN_INTERVALL(0:,I1),RNEW(1,I1))

         ! do interpolation only when optional arguments are given
         if( present(NTCELL) ) then
            ! interpolate shape function THETAS to new shape function THETASNEW
            ! save THETAS to THETASIN
            ICELL = NTCELL(I1)
            do LM1=1,NFUND
               THETASIN(:,LM1,ICELL)=THETAS(:,LM1,ICELL)
               IR2=0
               do IP=NPAN_LOGNEW(I1)+NPAN_EQNEW(I1)+1,NPAN_TOT(I1)
                  IR2=IR2+1
                  IMIN=IRCUT(IR2,I1)+1
                  IMAX=IRCUT(IR2+1,I1)
                  IMINNEW=IPAN_INTERVALL(IP-1,I1)+1
                  IMAXNEW=IPAN_INTERVALL(IP,I1)
                  call INTERPOLSPLINE(R(IMIN:IMAX,I1),RNEW(IMINNEW:IMAXNEW,I1),&
                     THETASIN(IMIN-IRCUT(1,I1):IMAX-IRCUT(1,I1),LM1,ICELL),   &
                     THETASNEW(IMINNEW:IMAXNEW,LM1,ICELL),IMAX-IMIN+1,        &
                     IMAXNEW-IMINNEW+1)
               enddo
            enddo
         end if ! present(ntcell)
      enddo ! i1

      if( present(NTCELL) ) then
         i_all=-product(shape(THETASIN))*kind(THETASIN)
         deallocate(THETASIN,stat=i_stat)
         call memocc(i_stat,i_all,'THETASIN','CREATE_NEWMESH')
      endif

   end subroutine CREATE_NEWMESH

   !----------------------------------------------------------------------------
   ! SUBROUTINE: CHEBMESH
   !----------------------------------------------------------------------------
   subroutine CHEBMESH(NPAN,NCHEB,RI,RO)

      implicit none

      integer, intent(in) :: NPAN
      integer, intent(in) :: NCHEB  !< Number of Chebychev pannels for the new solver
      double precision, dimension(0:NPAN), intent(in) :: RI
      double precision, dimension(NPAN*(NCHEB+1)), intent(inout) :: RO
      integer :: I,K,IK
      double precision :: TAU

      do I=1,NPAN
         do K=0,NCHEB
            IK=I*NCHEB+I-K
            TAU=DCOS(((2*K+1)*PI)/(2*(NCHEB+1)))
            TAU=0.5d0*((RI(I)-RI(I-1))*TAU+RI(I)+RI(I-1))
            RO(IK)=TAU
         enddo
      enddo
   end subroutine CHEBMESH

end module mod_create_newmesh
