!-------------------------------------------------------------------------------
! SUBROUTINE: VMADELBLK
!> @brief Calculate the madelung potentials and add these to the potential \f$V\f$
!> (in he spin-polarized case for each spin-direction this is the same)
!> @details It uses the structure dependent matrices AVMAD and BVMAD which
!> are calculated once in the subroutine MADELUNG3D() and saved in
!> the DA-file abvmad.unformatted (May 2004)
!> The charge-moments are calculated in the subroutine vintras,
!> therefore vintras has to be called first.
!> The madelung-potential is expanded into spherical harmonics.
!> The lm-term of the potential \f$V\f$ of the atom \f$i\f$ is given by
!> \f$ V(r,lm,i) = \sum_{i2}^{N} \sum_{l'm'} (-r)^l * \left\{avmad(i,i2,lm,l'm')*cmom(i2,l'm') +bvmad(i,i2,lm)*z(i2)\right\}\f$
!> where \f$ N\f$ is the number of atoms
!> @author B. Drittler
!> @date Nov. 1989
!> @note
!> - V. Popescu Feb. 2002: Adopted for the case of more atoms on the same site, summation is done over the occupants of that site, the charge is weighted with the appropriate concentration of the occupant
!> - Impurity-program adopted feb. 2004 (according to N. Papanikalou)
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine VMADELBLK(CMOM,CMINST,LMAX,NSPIN,NAEZ,V,ZAT,R,IRWS,IRCUT,IPAN,KSHAPE, &
   NOQ,KAOEZ,CONC,CATOM,ICC,HOSTIMP,VINTERS,IRM,NEMB,LMPOT,NPOTD,LMMAXD,NATYP)

   use Constants
   use global_variables

   implicit none

   ! .. Input variables
   integer, intent(in) :: ICC       !< Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
   integer, intent(in) :: IRM       !< Maximum number of radial points
   integer, intent(in) :: NAEZ      !< Number of atoms in unit cell
   integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
   integer, intent(in) :: NEMB      !< Number of 'embedding' positions
   integer, intent(in) :: NATYP     !< Number of kinds of atoms in unit cell
   integer, intent(in) :: NSPIN     !< Counter for spin directions
   integer, intent(in) :: LMPOT     !< (LPOT+1)**2
   integer, intent(in) :: NPOTD     !< (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
   integer, intent(in) :: KSHAPE    !< Exact treatment of WS cell
   integer, intent(in) :: LMMAXD    !< (KREL+KORBIT+1)(LMAX+1)^2
   ! .. Array Arguments
   integer, dimension(NAEZ), intent(in)            :: NOQ      !< Number of diff. atom types located
   integer, dimension(NATYP), intent(in)           :: IRWS     !< Position of atoms in the unit cell in units of bravais vectors
   integer, dimension(NATYP), intent(in)           :: IPAN     !< Number of panels in non-MT-region
   integer, dimension(0:NATYP), intent(in)         :: HOSTIMP
   integer, dimension(0:IPAND,NATYP), intent(in)   :: IRCUT    !< R points of panel borders
   integer, dimension(NATYP,NAEZ+NEMB), intent(in) :: KAOEZ    !< Kind of atom at site in elem. cell
   double precision, dimension(NATYP), intent(in)        :: ZAT      !< Nuclear charge
   double precision, dimension(NATYP), intent(in)        :: CONC     !< Concentration of a given atom
   double precision, dimension(NATYP), intent(in)        :: CATOM
   double precision, dimension(IRM,NATYP), intent(in)    :: R        !< Radial mesh ( in units a Bohr)
   double precision, dimension(LMPOT,NATYP), intent(in)  :: CMOM     !< LM moment of total charge
   double precision, dimension(LMPOT,NATYP), intent(in)  :: CMINST   !< charge moment of interstitial
   ! .. Input/Ouput variables
   double precision, dimension(IRM,LMPOT,NPOTD), intent(inout) :: V
   ! .. Output variables
   double precision, dimension(LMPOT,NAEZ), intent(out) :: VINTERS
   ! .. Local Scalars
   integer :: LRECABMAD,IREC
   integer :: I,L,LM,LM2,LMMAX,M,IO1,IO2,IPOT,IQ1,IQ2
   integer :: IRS1,ISPIN,IT1,IT2,NOQVAL
   double precision :: AC
   ! .. Local Arrays
   double precision, dimension(LMPOT) :: BVMAD  !< Structure dependent matrix
   double precision, dimension(LMPOT,LMPOT) :: AVMAD  !< Structure dependent matrix
   logical :: OPT
   ! .. Intrinsic Functions ..
   intrinsic :: SQRT
   !----------------------------------------------------------------------------
   write(1337,FMT=99001)
   write(1337,FMT=99002)
   !
   LRECABMAD = WLENGTH*2*LMPOT*LMPOT + WLENGTH*2*LMPOT
   open (69,ACCESS='direct',RECL=LRECABMAD,FILE='abvmad.unformatted',FORM='unformatted')
   !
   LMMAX = (LMAX+1)*(LMAX+1)
   !
   if (ICC.NE.0) then
      do IQ1=1,NAEZ
         do LM=1,LMPOT
            VINTERS(LM,IQ1) = 0D0
         end do
      end do
   end if
   !----------------------------------------------------------------------------
   ! Loop over all types in unit cell
   !----------------------------------------------------------------------------
   do IQ1 = 1,NAEZ               ! added bauer 2/7/2012
      NOQVAL=NOQ(IQ1)            ! added bauer 2/7/2012
      if (NOQVAL<1) NOQVAL=1     ! added bauer 2/7/2012
      do IO1 = 1,NOQVAL          ! added bauer 2/7/2012
         IT1 = KAOEZ(IO1,IQ1)    ! added bauer 2/7/2012

         !----------------------------------------------------------------------
         ! Take a site occupied by atom IT1
         !----------------------------------------------------------------------
         if (IT1/=-1) then                ! added bauer 2/7/2012
            if ( KSHAPE.NE.0 ) then
               IRS1 = IRCUT(IPAN(IT1),IT1)
            else
               IRS1 = IRWS(IT1)
            end if
         end if                           ! added bauer 2/7/2012
         !----------------------------------------------------------------------
         do L = 0,LMAX
            !-------------------------------------------------------------------
            do M = -L,L
               LM = L*L + L + M + 1
               AC = 0.0D0
               !----------------------------------------------------------------
               if ( NAEZ.EQ.1 ) then
                  IREC = IQ1 + NAEZ*(IQ1-1)
                  read(69,REC=IREC) AVMAD,BVMAD
                  !-------------------------------------------------------------
                  ! Loop over all occupants of site IQ2=IQ1
                  !-------------------------------------------------------------
                  do IO2 = 1,NOQ(IQ1)
                     IT2 = KAOEZ(IO2,IQ1)
                     !----------------------------------------------------------
                     ! lm = 1 component disappears if there is only one host atom
                     ! take moments of sphere
                     !----------------------------------------------------------
                     do LM2 = 2,LMMAX
                        AC = AC + AVMAD(LM,LM2)*CMOM(LM2,IT2)*CONC(IT2)
                     end do
                     !----------------------------------------------------------
                     ! Add contribution of interstial in case of shapes
                     !----------------------------------------------------------
                     if ( KSHAPE.NE.0 ) then
                        do LM2 = 2,LMMAX
                           AC = AC + AVMAD(LM,LM2)*CMINST(LM2,IT2)*CONC(IT2)
                        end do
                     end if
                  end do
                  !-------------------------------------------------------------
               else
                  !-------------------------------------------------------------
                  ! Loop over all sites
                  !-------------------------------------------------------------
                  do IQ2 = 1,NAEZ
                     IREC = IQ2 + NAEZ*(IQ1-1)
                     read(69,REC=IREC) AVMAD,BVMAD
                     !----------------------------------------------------------
                     ! Loop over all occupants of site IQ2
                     !----------------------------------------------------------
                     do IO2 = 1,NOQ(IQ2)
                        !
                        IT2 = KAOEZ(IO2,IQ2)
                        AC = AC + BVMAD(LM)*ZAT(IT2)*CONC(IT2)
                        !-------------------------------------------------------
                        ! Take moments of sphere
                        !-------------------------------------------------------
                        do LM2 = 1,LMMAX
                           AC = AC + AVMAD(LM,LM2)*CMOM(LM2,IT2)*CONC(IT2)
                        end do
                        !-------------------------------------------------------
                        ! Add contribution of interstial in case of shapes
                        !-------------------------------------------------------
                        if ( KSHAPE.NE.0 ) then
                           do LM2 = 1,LMMAX
                              AC = AC + AVMAD(LM,LM2)*CMINST(LM2,IT2)*CONC(IT2)
                           end do
                        end if
                     end do    ! IO2 = 1, NOQ(IQ2)
                     !----------------------------------------------------------
                  end do       ! IQ2 = 1, NAEZ
                  !-------------------------------------------------------------
               end if          ! NAEZ.GT.1
               !----------------------------------------------------------------
               if ( LM.eq.1 ) then
                  write (1337,FMT=99003) IT1,(CATOM(IT1)-ZAT(IT1)),(AC/SQRT(4.D0*PI))
               endif
               !----------------------------------------------------------------
               ! Add to v the intercell-potential
               !----------------------------------------------------------------
               !----------------------------------------------------------------
               ! SPIN
               !----------------------------------------------------------------
               do ISPIN = 1,NSPIN
                  !-------------------------------------------------------------
                  ! Determine the right potential number
                  !-------------------------------------------------------------
                  IPOT = NSPIN*(IT1-1) + ISPIN
                  !-------------------------------------------------------------
                  ! In the case of l=0 : r(1)**l is not defined
                  !-------------------------------------------------------------
                  if (IT1/=-1) then                ! added bauer 2/7/2012
                     if ( L.EQ.0 ) V(1,1,IPOT) = V(1,1,IPOT) + AC
                     do I = 2,IRS1
                        V(I,LM,IPOT) = V(I,LM,IPOT) + (-R(I,IT1))**L*AC
                     end do
                  end if
               end do                              ! added bauer 2/7/2012
               !----------------------------------------------------------------
               ! SPIN
               !----------------------------------------------------------------
               if (ICC.NE.0 .or. OPT('KKRFLEX ')) then
                  LM = L*L + L + M + 1
                  write(1337,*) 'ac',iq1,lm,ac
                  VINTERS(LM,IQ1) = AC
               end if
               !
            end do
            !-------------------------------------------------------------------
         end do
         !----------------------------------------------------------------------
      end do
   end do
   !----------------------------------------------------------------------------
   close(69)
   !----------------------------------------------------------------------------
   write(1337,*) 'ICC in VMADELBLK',ICC
   write(1337,'(25X,30(1H-),/)')
   write(1337,'(79(1H=))')
   !
   if ( (ICC==0) .and. (.not.OPT('KKRFLEX ')) ) return
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Now Prepare output for Impurity calculation
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   open (91,FILE='intercell_ref',STATUS='unknown',FORM='formatted')
   write(1337,*)
   write(1337,*) '                     ','Writing intercell potential for impurity'
   write(1337,'(/,20X,55(1H-))')
   write(1337,99004) HOSTIMP(0),LMMAX
   write(1337,'(20X,55(1H-),/,35X,"  i host lm  Vint")')
   do I=1,HOSTIMP(0)
      write(1337,*)
      LM = 1
      write(1337,'(35X,I4,I4,I3,1X,F10.6)') I, HOSTIMP(I),LM,VINTERS(LM,HOSTIMP(I))
      do LM=2,9
         write (1337,'(43X,I3,1X,F10.6)') LM,VINTERS(LM,HOSTIMP(I))
      end do
      write(1337,'(20X,55(1H-))')
   end do
   write(1337,'(79(1H=),/)')
   !
   write(91,99005) HOSTIMP(0),LMMAX
   do I=1,HOSTIMP(0)
      write(91,99006) (VINTERS(LM,HOSTIMP(I)),LM=1,LMMAX)
   end do
   close(91)

   return
   !
   99001 format (79(1H=),/,18X,' MADELUNG POTENTIALS ','(spherically averaged) ')
   99002 format (/,25X,' ATOM ','  Delta_Q  ','     VMAD',/,25X,30(1H-))
   99003 format (25X,I4,2X,F10.6,1X,F12.6)
   99004 format (22X,I4,' host atoms, LMPOT = ',I2,' output up to LM = 9')
   99005 format (3I6)
   99006 format (4D20.10)
end subroutine VMADELBLK
