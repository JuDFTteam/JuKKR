!-------------------------------------------------------------------------------
! SUBROUTINE: EMESHT
!> @brief This subroutine provides the energy mesh in array EZ and the
!> appropriate integration weights in array DF.
!> @details Poles of the Fermi function C (Matsubara frequencies) and a contour in
!> the complex energy are used as described in (????).
!> The contour consists of three straight lines with NPNT1, NPNT2, and NPNT3
!> integration points and is determined by the input arguments: EBOT, EMU, TK, and NPOL.
!> The three lines are defined by:
!> 1. The line from EBOT to \f$ EBOT+2*NPOL*\pi*i*k*TK \f$ with NPNT1 integration points (Gauss-Legendre rule)
!> 2. The line from \f$ EBOT+2*NPOL*\pi*i*k*TK\f$ to \f$ EMU+(2*NPOL*\pi*i-30)*k*TK\f$ with NPNT2 integration points (Gauss-Legendre rule)
!> 3. The line from \f$ EMU+(2*NPOL*\pi*i-30)*k*TK\f$ to \f$ \infty \f$
!>
!> The total number of integration points is given by: \f$ NPNT=NPNT1+NPNT2+NPNT3+NPOL\f$
!> The integration points and weights on three lines are chosen according to Gauss integration rules. Only in third interval
!> the Fermi function matters since \f$ e^x < 10^{-10} \f$ for \f$ x < -25\f$.
!> There are two special cases determined by NPOL = 0 and NPOL < 0.
!> - NPOL = 0 leads to density-of-states calculations with constant integration weights and equally distributed points
!> between \f$ EBOT - \pi*i*k*TK\f$ and \f$ EMU - \pi*i*k*TK\f$.
!> The total number of integration points is given by: NPNT=NPNT2
!> - NPOL < 0 is meant for calculations where the Fermi-Dirac function is replaced by a step function with step at EMU. When
!> this option is used no poles of the Fermi-Dirac function are used and the contour consists of the three straight lines:
!> 1. The line from \f$EBOT\f$ to \f$ EBOT-2*NPOL*\pi*i*k*TK\f$ with NPNT1 integration points (Gauss-Legendre rule)
!> 2. The line from \f$EBOT-2*NPOL*\pi*i*k*TK\f$ to \f$EMU-2*NPOL*\pi*i*k*TK\f$ with NPNT2 integration points (Gauss-Legendre rule)
!> 3. The line from \f$ EMU-2*NPOL*\pi*i*k*TK\f$ to \f$EMU\f$ with NPNT3 integration points (Gauss-Legendre rule)
!>
!> The total number of integration points is given by: \f$NPNT=NPNT1+NPNT2+NPNT3\f$
!> @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine EMESHT(EZ,DF,NPNT,EBOT,EMU,EFERMI,TK,NPOL,NPNT1,NPNT2,NPNT3,IEMXD)

   use mod_types, only: t_inc
   use Constants

   implicit none
   ! ..
   ! .. Input variables
   integer, intent(in) :: NPOL   !< Number of Matsubara Poles (EMESHT)
   integer, intent(in) :: NPNT1  !< number of E points (EMESHT) for the contour integration
   integer, intent(in) :: NPNT2  !< number of E points (EMESHT) for the contour integration
   integer, intent(in) :: NPNT3  !< number of E points (EMESHT) for the contour integration
   integer, intent(in) :: IEMXD  !< Dimension for energy-dependent arrays
   double precision, intent(in) :: TK     !< Temperature
   double precision, intent(in) :: EMU    !< Top of the contour
   double precision, intent(in) :: EBOT   !< Bottom of the contour
   double precision, intent(in) :: EFERMI !< Fermi energy
   ! .. Input/Output variables
   integer, intent(inout) :: NPNT
   double complex, dimension(IEMXD), intent(inout) :: DF
   double COMPLEX, dimension(IEMXD), intent(inout) :: EZ
   ! .. Local Scalars ..
   integer :: I
   double complex :: DE
   double precision :: ER,ETK
   ! .. Local Arrays ..
   double precision, dimension(128) :: WI,XI
   ! .. External Subroutines ..
   logical :: OPT
   external :: GAUFD,GAULEG,OPT
   ! ..
   ! .. Intrinsic Functions ..
   intrinsic :: DCMPLX
   ! ..
   !----------------------------------------------------------------------------
   ! OUTPUT
   !----------------------------------------------------------------------------
   if(t_inc%i_write>0) then
      write (1337,'(5X,A,F12.6," (Ry)",8X,A,F12.6," (Ry)")') &
         'E min = ',EBOT,'Fermi energy = ',EFERMI
      write (1337,'(5X,A,F12.6," (Ry)",8X,A,F12.6," (K )",/,5X,62(1H-))')  &
         'E max = ',EMU,'Temperature  = ',TK
   endif
   !----------------------------------------------------------------------------
   ! OUTPUT
   !----------------------------------------------------------------------------
   ETK = PI*KB*TK
   !----------------------------------------------------------------------------
   if (NPOL.EQ.0) then
      DE = (EMU-EBOT)
      if (NPNT2.GT.1) then
         DE = DE/(NPNT2-1)
      else
         DE=DCMPLX(1.0D0,0.0D0)
      end if
      NPNT = 0
      do I = 1,NPNT2
         NPNT = NPNT + 1
         if ( NPNT.GT.IEMXD ) then
            write(6,'(/,5X,2A,I4)') 'Dimension ERROR: Increase IEMXD in the inputcard to ',  &
            'at least ',NPNT
            stop '     < EMESHT >'
         end if
         ER = EBOT + (I-1)*DE
         EZ(NPNT) = DCMPLX(ER,ETK)
         DF(NPNT) = DE
      enddo ! I
      if(t_inc%i_write>0) write (1337,FMT=9000) NPNT,ETK,ETK*RYD
   !----------------------------------------------------------------------------
   ! NPOL > 0
   !----------------------------------------------------------------------------
   else if (NPOL.GT.0) then
      call GAULEG(XI,WI,NPNT1)
      DE = NPOL*DCMPLX(0.0D0,ETK)
      NPNT = 0
      do I = 1,NPNT1
         NPNT = NPNT + 1
         if ( NPNT.GT.IEMXD ) then
            write(6,'(/,5X,2A,I4)')'Dimension ERROR: Increase IEMXD in the inputcard to ',  &
            'at least ',NPNT
            stop '     < EMESHT >'
         end if
      EZ(NPNT) = XI(I)*DE + DE + EBOT
      DF(NPNT) = WI(I)*DE
   enddo ! I -> NPNT1
      call GAULEG(XI,WI,NPNT2)
      DE = (EMU-30*KB*TK-EBOT)*0.5D0
      do I = 1,NPNT2
         NPNT = NPNT + 1
         if ( NPNT.GT.IEMXD ) then
            write(6,'(/,5X,2A,I4)')'Dimension ERROR: Increase IEMXD in the inputcard to ',   &
            'at least ',NPNT
            stop '     < EMESHT >'
         end if
         EZ(NPNT) = XI(I)*DE + DE + EBOT + 2*NPOL*DCMPLX(0.0D0,ETK)
         DF(NPNT) = WI(I)*DE
      enddo ! I -> NPTN2
      call GAUFD(XI,WI,NPNT3)
      DE = 30*KB*TK
      do I = 1,NPNT3
         NPNT = NPNT + 1
         if ( NPNT.GT.IEMXD ) then
            write(6,'(/,5X,2A,I4)')'Dimension ERROR: Increase IEMXD in the inputcard to ',  &
            'at least ',NPNT
            stop '     < EMESHT >'
         end if
      EZ(NPNT) = XI(I)*DE + EMU + 2*NPOL*DCMPLX(0.0D0,ETK)
      DF(NPNT) = WI(I)*DE
   enddo ! I - >NPTN3
      do 50 I = NPOL,1,-1
         NPNT = NPNT + 1
         if ( NPNT.gt.IEMXD ) then
            write(6,'(/,5X,2A,I4)')'Dimension ERROR: Increase IEMXD in the inputcard to ',   &
            'at least ',NPNT
            stop '     < EMESHT >'
         end if
         EZ(NPNT) = EMU + (2*I-1)*DCMPLX(0.0D0,ETK)
         DF(NPNT) = -2*DCMPLX(0.0D0,ETK)
      50 continue
      if(t_inc%i_write>0) write(1337,9090) NPNT,NPOL,NPNT1,NPNT2,NPNT3
      !-------------------------------------------------------------------------
      ! NPOL < 0
      !-------------------------------------------------------------------------
   else
      if (NPNT1.gt.0) call GAULEG(XI,WI,NPNT1)
      DE = -NPOL*DCMPLX(0.0D0,ETK)
      NPNT = 0
      do I = 1,NPNT1
         if ( NPNT.GT.IEMXD ) then
            write(6,'(/,5X,2A,I4)')'Dimension ERROR: Increase IEMXD in the inputcard to ',   &
            'at least ',NPNT
            stop '     < EMESHT >'
         end if
         NPNT = NPNT + 1
         EZ(NPNT) = XI(I)*DE + DE + EBOT
         DF(NPNT) = WI(I)*DE
      enddo ! I -> NPNT1
      call GAULEG(XI,WI,NPNT2)
      DE = (EMU-EBOT)*0.5D0
      do I = 1,NPNT2
         NPNT = NPNT + 1
         if ( NPNT.GT.IEMXD ) then
            write(6,'(/,5X,2A,I4)')'Dimension ERROR: Increase IEMXD in the inputcard to ',   &
            'at least ',NPNT
            stop '     < EMESHT >'
         end if
         EZ(NPNT) = XI(I)*DE + DE + EBOT - 2*NPOL*DCMPLX(0.0D0,ETK)
         if (OPT('GF-EF   ')) EZ(NPNT) = EMU + NPOL*DCMPLX(0.0D0,ETK)
         DF(NPNT) = WI(I)*DE
      enddo ! I -> NPNT2
      if (NPNT3.GT.0) call GAULEG(XI,WI,NPNT3)
      DE = -NPOL*DCMPLX(0.0D0,ETK)
      do I = NPNT3,1,-1
         NPNT = NPNT + 1
         if ( NPNT.GT.IEMXD ) then
            write(6,'(/,5X,2A,I4)')'Dimension ERROR: Increase IEMXD in the inputcard to ',   &
            'at least ',NPNT
            stop '     < EMESHT >'
         end if
         EZ(NPNT) = XI(I)*DE + DE + EMU
         DF(NPNT) = -WI(I)*DE
      enddo ! I -> NPNT3
      if(t_inc%i_write>0) write(1337,9091) NPNT,-NPOL,NPNT1,NPNT2,NPNT3
   end if
   !----------------------------------------------------------------------------
   if(t_inc%i_write>0) WRITE(1337,*)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Correction Factor for the weight in the integration according to Phivos Idea
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do I=1,NPNT
      DF(I)=DF(I)
   enddo
   !*(8.5D0/8.48686D0)*(8.75D0/8.74083D0)
   ! GaCrN*(8.5D0/8.49286D0)
   !*(8.5D0/8.48969D0)
   !*(8.5D0/8.48823D0)
   !*(8.75D0/8.73983D0)
   !*(8.5D0/8.48686D0)
   !*(8.75D0/8.75659D0)
   !*(6.5D0/6.55253D0)*(7.5D0/7.47798D0)
   !     *(8.75D0/8.75659D0)
   !*(8.5D0/8.54963D0)
   !*(6.5D0/6.41299D0)
   !*(8.5D0/8.47767D0)
   !*(6.5D0/6.45787D0)
   !*(4.D0/4.01579D0)
   !*(8.8D0/8.80272D0)*(8.8D0/8.78691D0)
   !*(4.0D0/4.0419213D0)*(17.5D0/17.508D0)
   !     &  *(8.0D0/7.9885D0)*(8.75D0/8.74682D0)*(8.0D0/7.9246D0)
   !     &  *(8.25D0/8.24085)
   ! 90        write(*,*)'DF=',I,DF(I)
   !*************************************************************
   !**********************************************************

   return
   9000 format (5X,'Density-of-States calculation',/,    &
   5X,'Number of energy points :',I4,4X,'broadening =',  &
   3P,F9.3,' ( mRy )',/,48X,' =',3P,F9.3,' ( meV )')
   9090 format (5X,'GF integration rectangular contour ( ImE > 0 )',/,  &
   5X,'Number of energy points :',I4,13X,'poles =',I2,/,                &
   23X,'contour: N1 =',I2,', N2 =',I4,', N3 =',I2)
   9091 format (5X,'GF integration rectangular contour ( ImE < 0 )',/,  &
   5X,'Number of energy points :',I4,13X,'poles =',I2,/,                &
   23X,'contour: N1 =',I2,', N2 =',I4,', N3 =',I2)
end subroutine EMESHT
