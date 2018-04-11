!-------------------------------------------------------------------------------
!> @brief add core and valence density expanded in spherical harmonics
!>         ( convention see subroutine rholm )
!> @details In the paramagnetic case (nspin=1) the core valence charge times
!> r**2 is added to the valence charge density times r**2
!> then only rho2ns(irmd,lmxtsq,natypd,1) is used .
!> In the spin-polarized case (nspin=2) the spin-splitted core
!> charge density times r**2 is converted into core charge
!> density times r**2 and core spin density times r**2 .
!> then these parts are added to corresponding parts of
!> the valence densities times r**2 , that are rho2ns(...,1)
!> which contains the charge density  and rho2ns(...,2) which
!> contains in that case the spin density .
!> (see notes by b.drittler)
!>
!> Attention : the core density is spherically averaged and multiplied by 4 pi.
!> therefore the core density is only added to l=0 part .
!
!> @author B. Drittler
!> @date   Nov. 1989
!
!> @note -V. Popescu March 2002: Total orbital moment within the WS sphere is also calculated
!> in the relativistic case; orbital density is normalised in the
!> same way as the charge density.
!> @note -Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine RHOTOTB(IPF,NATYP,NAEZ,NSPIN,RHO2NS,RHOC,RHOORB, &
   Z,DRDI,IRWS,                                             &
   IRCUT,LPOT,NFU,LLMSP,THETAS,NTCELL,KSHAPE,IPAN,          &
   CHRGNT,ITC,NSHELL,NOQ,CONC,KAOEZ,CATOM)

   implicit none

   !     .. Parameters ..
   ! .. Input variables
   integer, intent(in) :: ITC
   integer, intent(in) :: IPF
   integer, intent(in) :: IRM    !< Maximum number of radial points
   integer, intent(in) :: NEMB   !< Number of 'embedding' positions
   integer, intent(in) :: LPOT   !< Maximum l component in potential expansion
   integer, intent(in) :: NAEZ   !< Number of atoms in unit cell
   integer, intent(in) :: NATYP  !< Number of kinds of atoms in unit cell
   integer, intent(in) :: NSPIN  !< Counter for spin directions
   integer, intent(in) :: LMPOT  !< (LPOT+1)**2
   integer, intent(in) :: KSHAPE !< Exact treatment of WS cell
   !     .. Array Arguments ..
   integer, dimension(NAEZ), intent(in)      :: NOQ !< Number of diff. atom types located
   integer, dimension(*), intent(in)         :: NFU  !< number of shape function components in cell 'icell'
   integer, dimension(*), intent(in)         :: IPAN !< Number of panels in non-MT-region
   integer, dimension(*), intent(in)         :: IRWS !< R point at WS radius
   integer, dimension(0:NSHELD), intent(in)  :: NSHELL   !< Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
   integer, dimension(*), intent(in)         :: NTCELL !< Index for WS cell

   integer, dimension(0:IPAND,*), intent(in)          :: IRCUT !< R points of panel borders
   integer, dimension(NATYP,*), intent(in)            :: LLMSP  !< lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
   integer, dimension(NATYP,NAEZ+NEMB), intent(in)    :: KAOEZ !< Kind of atom at site in elem. cell
   double precision, dimension(*), intent(in)      :: Z
   double precision, dimension(NATYP), intent(in)  :: CONC  !< Concentration of a given atom
   double precision, dimension(IRM,*), intent(in)  :: DRDI  !< Derivative dr/di
   double precision, dimension(IRM,*), intent(in)  :: RHOC  !< core charge density
   double precision, dimension(IRM*KREL+(1-KREL),NATYP), intent(in) :: RHOORB   !< Orbital density
   double precision, dimension(IRID,NFUND,*), intent(in) :: THETAS   !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion

   ! .. In/Out variables
   double precision, dimension(IRM,LMPOT,NATYP,*), intent(inout) :: RHO2NS
   ! .. Output variables
   double precision, intent(out) :: CHRGNT
   double precision, dimension(NATYP,2*KREL+(1-KREL)*NSPIND), intent(out) :: CATOM
   ! .. Local variables
   integer :: I,I1,IATYP,ICELL,IFUN,IPAN1,IPOTD,IPOTU,IRC1,IRS1,ISPIN,LM,IQEZ,IOEZ
   double precision :: DIFF,FACTOR,RFPI,SUM,TOTSMOM,TOTOMOM,SUMO
   double precision, dimension(NATYP)                       :: OMOM   !< Orbital moment
   double precision, dimension(IRM)                         :: RHO
   double precision, dimension(NAEZ,2*KREL+(1-KREL)*NSPIND) :: CSITE
   double precision, dimension(KREL*NAEZ+(1-KREL))          :: MUOSITE
   !
   logical :: OPT
   ! .. External Subroutines
   external :: SIMP3,SIMPK,OPT
   ! .. Intrinsic Functions
   intrinsic :: ATAN,SQRT
   ! .. Save statement
   SAVE
   !     ..
   RFPI = SQRT(16.0D0*ATAN(1.0D0))

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Loop over atomic sites
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do IQEZ = 1, NAEZ
      !-------------------------------------------------------------------------
      ! Loop over atoms located on IQEZ
      !-------------------------------------------------------------------------
      do ISPIN=1,NSPIN
         CSITE(IQEZ,ISPIN) = 0.0D0
      end do

      if (KREL.EQ.1) MUOSITE(IQEZ) = 0.0D0

      do IOEZ =1, NOQ(IQEZ)

         IATYP = KAOEZ(IOEZ,IQEZ)
         !---------------------------------------------------------------------
         ! Determine the right potential numbers for rhoc
         !---------------------------------------------------------------------
         if (NSPIN.EQ.2) then
            IPOTD = 2*IATYP - 1
            IPOTU = 2*IATYP
            FACTOR = 1.0D0
         else
            IPOTD = IATYP
            IPOTU = IATYP
            FACTOR = 0.5D0
         end if

         if (KSHAPE.NE.0) then
            IPAN1 = IPAN(IATYP)
            IRS1 = IRCUT(1,IATYP)
            IRC1 = IRCUT(IPAN1,IATYP)
         else
            IRS1 = IRWS(IATYP)
            IRC1 = IRS1
         end if
         !-----------------------------------------------------------------------
         do I = 2,IRS1
            !-------------------------------------------------------------------
            ! Convert core density
            !-------------------------------------------------------------------
            SUM = (RHOC(I,IPOTD)+RHOC(I,IPOTU))*FACTOR/RFPI
            DIFF = (RHOC(I,IPOTU)-RHOC(I,IPOTD))/RFPI
            !-------------------------------------------------------------------
            ! Add this to the lm=1 component of rho2ns
            !-------------------------------------------------------------------
            RHO2NS(I,1,IATYP,1) = RHO2NS(I,1,IATYP,1) + SUM
            RHO2NS(I,1,IATYP,NSPIN) = RHO2NS(I,1,IATYP,NSPIN) + DIFF
         end do
         !----------------------------------------------------------------------
         ! Calculate  charge and moment of the atom
         !----------------------------------------------------------------------
         do ISPIN = 1,NSPIN
            !
            if (KSHAPE.EQ.0) then
               !----------------------------------------------------------------
               ! Integrate over wigner seitz sphere - no shape correction
               !----------------------------------------------------------------
               call SIMP3(RHO2NS(1,1,IATYP,ISPIN),SUM,1,IRS1,DRDI(1,IATYP))
               !----------------------------------------------------------------
               ! The result has to be multiplied by sqrt(4 pi)
               ! (4 pi for integration over angle and 1/sqrt(4 pi) for
               ! the spherical harmonic y(l=0))
               !----------------------------------------------------------------
               SUM = SUM*RFPI
            else ! (KSHAPE.EQ.0)
               !----------------------------------------------------------------
               ! convolute charge density with shape function to get the
               ! charge in the exact cell - if kshape .gt. 0
               !----------------------------------------------------------------
               ICELL = NTCELL(IATYP)

               do I = 1,IRS1
                  RHO(I) = RHO2NS(I,1,IATYP,ISPIN)*RFPI
               end do
               !
               do I = IRS1 + 1,IRC1
                  RHO(I) = 0.0D0
               end do
               !
               do IFUN = 1,NFU(ICELL)
                  LM = LLMSP(ICELL,IFUN)
                  if (LM.LE.LMPOT) then
                     do I = IRS1 + 1,IRC1
                        RHO(I) = RHO(I) + RHO2NS(I,LM,IATYP,ISPIN)*  &
                        THETAS(I-IRS1,IFUN,ICELL)
                     end do
                  end if
               end do
               !----------------------------------------------------------------
               ! Integrate over circumscribed sphere
               !----------------------------------------------------------------
               call SIMPK(RHO,SUM,IPAN1,IRCUT(0,IATYP),DRDI(1,IATYP))
            end if                    ! (KSHAPE.EQ.0)

            CATOM(IATYP,ISPIN) = SUM
            CSITE(IQEZ,ISPIN) = CSITE(IQEZ,ISPIN)+ CATOM(IATYP,ISPIN) *CONC(IATYP)

            if (ISPIN.NE.1) then
               !----------------------------------------------------------------
               ! Calculate orbital moment (ASA) and add it to the total
               !----------------------------------------------------------------
               if ((KREL.EQ.1) .AND. (KSHAPE.EQ.0)) then
                  call SIMP3(RHOORB(1,IATYP),SUMO,1,IRS1,DRDI(1,IATYP))
                  SUMO = SUMO*RFPI
                  OMOM(IATYP) = SUMO
                  MUOSITE(IQEZ) = MUOSITE(IQEZ)+ OMOM(IATYP) * CONC(IATYP)
               end if

               if (KSHAPE.NE.0) then
                  write (IPF,FMT=9010) SUM
               else
                  write (IPF,FMT=9050) SUM
                  if (KREL.EQ.1) then
                     write (IPF,FMT=9051) OMOM(IATYP)
                     write (IPF,FMT=9052) SUM+OMOM(IATYP)
                  end if
               end if
            else                      ! (ISPIN.NE.1)
               if (KSHAPE.NE.0) then
                  write (IPF,FMT=9000) IATYP,SUM
               else
                  write (IPF,FMT=9040) IATYP,SUM
               end if
            end if                    ! (ISPIN.NE.1)
         end do                      ! ISPIN = 1,NSPIN
         !----------------------------------------------------------------------
         if (IOEZ.NE.NOQ(IQEZ)) write(IPF,'(2X,77(1H-))')
      end do
      !-------------------------------------------------------------------------
      ! IOEZ = 1, NOQ(IQEZ)
      !-------------------------------------------------------------------------
      if (NOQ(IQEZ).GT.1) then
         write(IPF,'(2X,77(1H=))')
         write(IPF,FMT=9071) IQEZ,CSITE(IQEZ,1)
         if (NSPIN.EQ.2) then
            write(IPF,FMT=9072) CSITE(IQEZ,NSPIN)
            if (KREL.EQ.1) then
               write(IPF,FMT=9073) MUOSITE(IQEZ)
               write(IPF,FMT=9074) CSITE(IQEZ,NSPIN)+MUOSITE(IQEZ)
            end if
         end if
         if (IQEZ.NE.NAEZ) write(IPF,'(2X,77(1H=))')
      else
         if (IQEZ.NE.NAEZ) write(IPF,'(2X,77(1H=))')
      end if
   end do
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! IQEZ = 1, NAEZ
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(IPF,*)

   CHRGNT = 0.0D0
   do I1 = 1,NATYP
      CHRGNT = CHRGNT+ DBLE(NSHELL(I1))*(CATOM(I1,1) - Z(I1))*CONC(I1)
   end do

   write(IPF,'(79(1H+))')
   write (IPF,FMT=9020) ITC,CHRGNT
   write ( 6 ,FMT=9020) ITC,CHRGNT

   if (NSPIN.EQ.2) then
      TOTSMOM = 0.0D0
      if (KREL.EQ.1) TOTOMOM = 0.0D0
      do I1 = 1,NATYP
         TOTSMOM = TOTSMOM+ DBLE(NSHELL(I1))*CATOM(I1,NSPIN)*CONC(I1)
         if (KREL.EQ.1) TOTOMOM = TOTOMOM+ DBLE(NSHELL(I1))*OMOM(I1)*CONC(I1)
      end do

      if (KREL.EQ.0) then
         write (IPF,FMT=9030) TOTSMOM
         write ( 6 ,FMT=9030) TOTSMOM
      else
         write (IPF,FMT=9030) TOTSMOM+TOTOMOM
         write (IPF,FMT=9031) TOTSMOM
         write (IPF,FMT=9032) TOTOMOM
         write ( 6 ,FMT=9030) TOTSMOM+TOTOMOM
         write ( 6 ,FMT=9031) TOTSMOM
         write ( 6 ,FMT=9032) TOTOMOM
      end IFUN
   end if
   write (IPF,*)

   return

   9000 format ('  Atom ',I4,' charge in wigner seitz cell =',f10.6)
   9010 format (7X,'spin moment in wigner seitz cell =',f10.6)
   9040 format ('  Atom ',I4,' charge in wigner seitz sphere =',f10.6)
   9050 format (7X,'spin moment in wigner seitz sphere =',f10.6)
   9051 format (7X,'orb. moment in wigner seitz sphere =',f10.6)
   9052 format (7X,'total magnetic moment in WS sphere =',f10.6)
   9020 format ('      ITERATION',I4,  &
   ' charge neutrality in unit cell = ',f12.6)
   9030 format ('                   ', &
   ' TOTAL mag. moment in unit cell = ',f12.6)
   9031 format ('                   ', &
   '           spin magnetic moment = ',f12.6)
   9032 format ('                   ', &
   '        orbital magnetic moment = ',f12.6)
   9071 format ('      Site ',i3,' total charge =',f10.6)
   9072 format ('         ',' total spin moment =',f10.6)
   9073 format ('         ',' total orb. moment =',f10.6)
   9074 format ('      total magnetic moment =',f10.6)

end subroutine RHOTOTB
