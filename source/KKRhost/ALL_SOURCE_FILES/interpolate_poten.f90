!-------------------------------------------------------------------------------
! SUBROUTINE: INTERPOLATE_POTEN
!> @brief Routine for the interpolation of the potential in the integration grid.
!> @note Jonathan Chico: Include array dimensions in interface explicitly to get rid of
!> inc.p import and to be able to use routine for different number of atoms
!-------------------------------------------------------------------------------
subroutine INTERPOLATE_POTEN(LPOT,IRM,IRNSD,NATYP,IPAND,LMPOT,    &
      NSPOTD,NTOTD,NCHEBD,IRMDNEW,NSPIN,R,IRMIND,IRWS,IRCUT,VINS, &
      VISP,NPAN_LOG,NPAN_EQ,NPAN_TOT,RNEW,IPAN_INTERVALL,VINSNEW)

   implicit none

   integer, intent(in) :: IRM    !< Maximum number of radial points
   integer, intent(in) :: LPOT   !< Maximum l component in potential expansion
   integer, intent(in) :: LMPOT  !< (LPOT+1)**2
   integer, intent(in) :: IRNSD
   integer, intent(in) :: NTOTD
   integer, intent(in) :: NCHEB  !< Number of Chebychev pannels for the new solver
   integer, intent(in) :: NATYP  !< Number of kinds of atoms in unit cell
   integer, intent(in) :: IPAND  !< Number of panels in non-spherical part
   integer, intent(in) :: LMPOT  !< (LPOT+1)**2
   integer, intent(in) :: NSPIN  !< Counter for spin directions
   integer, intent(in) :: NSPOTD
   integer, intent(in) :: IRMIND !< IRM-IRNSD
   integer, intent(in) :: IRMDNEW
   integer, dimension(NATYP), intent(in) :: IRWS   !< R point at WS radius
   integer, dimension(NATYP), intent(in) :: IRMIN  !< Max R for spherical treatment
   integer, dimension(NATYP), intent(in) :: NPAN_EQ   !< Variables for the pannels for the new solver
   integer, dimension(NATYP), intent(in) :: NPAN_LOG  !< Variables for the pannels for the new solver
   integer, dimension(NATYP), intent(in) :: NPAN_TOT
   integer, dimension(0:IPAND,NATYP), intent(in) :: IRCUT   !< R points of panel borders
   integer, dimension(0:NTOTD,NATYP), intent(in) :: IPAN_INTERVALL
   double precision, dimension(IRM,NATYP), intent(in) :: R
   double precision, dimension(IRM,NSPOTD), intent(in) :: VISP
   double precision, dimension((IRM-IRNSD):IRM,(LPOT+1)**2,NSPOTD), intent(in)  :: VINS
   ! .. Input/Output variables
   double precision, dimension(IRMDNEW,(LPOT+1)**2,NSPOTD), intent(inout) :: VINSNEW

   ! .. Local variables
   integer I1,IPOT,IPOTM,IMIN,IMAX,IP,IR,LM1,ISPIN,IMINNEW,IMAXNEW,IR2
   double precision, dimension(IRMDNEW,NATYP)         :: RNEW
   double precision, dimension(IRM,(LPOT+1)**2,NSPIN) :: VINSIN

   VINSNEW=0d0
   IPOTM=0

   ! interpolate potential to new mesh
   do I1 = 1,NATYP

      IPOT=NSPIN*(I1-1)+1

      ! save input potential to VINSIN
      VINSIN=0d0
      do IR=1,IRWS(I1)
         if (IR.LT.IRMIN(I1)) then
            VINSIN(IR,1,1)=VISP(IR,IPOT)
            VINSIN(IR,1,NSPIN)=VISP(IR,IPOT+NSPIN-1)
         else
            do LM1=1,LMPOT
               if (LM1.EQ.1) then
                  VINSIN(IR,LM1,1)=VISP(IR,IPOT)
                  VINSIN(IR,LM1,NSPIN)=VISP(IR,IPOT+NSPIN-1)
               else
                  VINSIN(IR,LM1,1)=VINS(IR,LM1,IPOT)
                  VINSIN(IR,LM1,NSPIN)=VINS(IR,LM1,IPOT+NSPIN-1)
               endif
            enddo
         endif
      enddo

      do ISPIN=1,NSPIN
         IPOTM=IPOTM+1
         do LM1=1,LMPOT

            IMIN=1
            IMAX=IRMIN(I1)
            do IP=1,NPAN_LOG(I1)
               IMINNEW=IPAN_INTERVALL(IP-1,I1)+1
               IMAXNEW=IPAN_INTERVALL(IP,I1)
               call INTERPOLSPLINE(R(IMIN:IMAX,I1),RNEW(IMINNEW:IMAXNEW,I1),&
                  VINSIN(IMIN:IMAX,LM1,ISPIN),&
                  VINSNEW(IMINNEW:IMAXNEW,LM1,IPOTM),&
                  IMAX-IMIN+1,IMAXNEW-IMINNEW+1)
            enddo

            IMIN=IRMIN(I1)
            IMAX=IRCUT(1,I1)
            do IP=NPAN_LOG(I1)+1,NPAN_LOG(I1)+NPAN_EQ(I1)
               IMINNEW=IPAN_INTERVALL(IP-1,I1)+1
               IMAXNEW=IPAN_INTERVALL(IP,I1)
               call INTERPOLSPLINE(R(IMIN:IMAX,I1),RNEW(IMINNEW:IMAXNEW,I1),&
                  VINSIN(IMIN:IMAX,LM1,ISPIN),&
                  VINSNEW(IMINNEW:IMAXNEW,LM1,IPOTM),&
                  IMAX-IMIN+1,IMAXNEW-IMINNEW+1)
            enddo

            IR2=0
            do IP=NPAN_LOG(I1)+NPAN_EQ(I1)+1,NPAN_TOT(I1)
               IR2=IR2+1
               IMIN=IRCUT(IR2,I1)+1
               IMAX=IRCUT(IR2+1,I1)
               IMINNEW=IPAN_INTERVALL(IP-1,I1)+1
               IMAXNEW=IPAN_INTERVALL(IP,I1)
               call INTERPOLSPLINE(R(IMIN:IMAX,I1),RNEW(IMINNEW:IMAXNEW,I1),&
                  VINSIN(IMIN:IMAX,LM1,ISPIN),&
                  VINSNEW(IMINNEW:IMAXNEW,LM1,IPOTM),&
                  IMAX-IMIN+1,IMAXNEW-IMINNEW+1)
            enddo
         enddo ! lm1
      enddo ! ispin
   enddo ! i1
end subroutine INTERPOLATE_POTEN
