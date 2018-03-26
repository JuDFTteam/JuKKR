!*==wmatldau.f    processed by SPAG 6.05Rc at 16:58 on 22 Dec 2004
!-------------------------------------------------------------------------------
! SUBROUTINE: WMATLDAU
!> @brief Calculation of Coulomb interaction potential in LDA+U non-relativistic case
!> otherwise matrices DENMAT and VLDAU must have double dimension
!> @details Uses the Coulomb matrix U (array ULDAU), the density matrix \f$n\f$
!> (array DENMAT) and the occupation numbers dentot (total) and \f$n_s\f$  (array DENTOTS) (per spin).
!> The expression evaluated (array VLDAU) is
!> \f$V_{m1,s,m2,s'} = \delta_{ss'} \sum_{s'',m3,m4} U_{m1,m2,m3,m4} n_{m3,s'',m4,s''} - \sum_{m3,m4} U_{m1,m4,m3,m2} n_{m3,s',m4,s} - [Ueff (dentot-1/2) - Jeff (n_s - 1/2)] delta_{ss'} \delta_{m1,m2} \f$
!> @author Ph. Mavropoulos, H. Ebert (Munich)
!> @date 2002-2004
!> @note Modifications by N. Long Xmas Juelich 2015
!-------------------------------------------------------------------------------
subroutine WMATLDAU(NTLDAU,ITLDAU,NSPIN,DENMATC,LOPT, &
      UEFF,JEFF,ULDAU,WLDAU,EU,EDC,MMAXD,NPOTD,NATYP,NSPIND)
   !**********************************************************************
   !*                                                                    *
   !* Calculation of Coulomb interaction potential in LDA+U              *
   !* non-relativistic case -- otherwise matrices DENMAT and VLDAU must  *
   !*                          have double dimension                     *
   !*                                                                    *
   !* Uses the Coulomb matrix U (array ULDAU), the density matrix n      *
   !* (array DENMAT) and the occupation numbers dentot (total) and n_s   *
   !* (array DENTOTS) (per spin).                                        *
   !*                                                                    *
   !* The expression evaluated (array VLDAU) is                          *
   !*                                                                    *
   !*       V_{m1,s,m2,s'} =                                             *
   !* delta_{ss'} Sum_{s'',m3,m4} U_{m1,m2,m3,m4} n_{m3,s'',m4,s''}      *
   !* - Sum_{m3,m4} U_{m1,m4,m3,m2} n_{m3,s',m4,s}                       *
   !* - [Ueff (dentot-1/2) - Jeff (n_s - 1/2)] delta_{ss'} delta_{m1,m2} *
   !*                                                                    *
   !*                  ph. mavropoulos, h.ebert munich/juelich 2002-2004 *
   !*                  n.long, Xmas Juelich 2015                         *
   !**********************************************************************
   use Constants

   implicit none
   ! .. Input variables
   integer, intent(in) :: NSPIN  !< Counter for spin directions
   integer, intent(in) :: MMAXD  !< 2*LMAX+1
   integer, intent(in) :: NPOTD  !< (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
   integer, intent(in) :: NATYP  !< Number of kinds of atoms in unit cell
   integer, intent(in) :: NTLDAU !< number of atoms on which LDA+U is applied
   integer, intent(in) :: NSPIND !< KREL+(1-KREL)*(NSPIN+1)
   integer, dimension(NATYP), intent(in) :: LOPT   !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
   integer, dimension(NATYP), intent(in) :: ITLDAU !< integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cell
   double precision, dimension(NATYP), intent(in) :: UEFF   !< input U parameter for each atom
   double precision, dimension(NATYP), intent(in) :: JEFF   !< input J parameter for each atom
   ! .. Input/Output variables
   double precision, dimension(NATYP), intent(inout) :: EDC !< Double-counting correction
   double precision, dimension(NATYP), intent(inout) :: EU  !< Total energy corrections
   double precision, dimension(MMAXD,MMAXD,NSPIND,NATYP), intent(inout) :: WLDAU !< potential matrix
   double precision, dimension(MMAXD,MMAXD,MMAXD,MMAXD,NATYP), intent(inout) :: ULDAU  !< calculated Coulomb matrix elements (EREFLDAU)
   double complex, dimension(MMAXD,MMAXD,NPOTD), intent(inout) :: DENMATC
   ! .. Local variables
   integer :: IPRINT
   integer :: I1,IT,IPOT,IS,JS,M1,M2,M3,M4,MM,MMAX
   double precision :: FACTOR
   double precision :: DENTOT
   character(len=15) :: STR15
   couble complex :: CSUM,CSUM2
   double precision, dimension(NSPIND) :: DENTOTS
   double precision, dimension(MMAXD,MMAXD,NSPIND) :: DENMAT
   double complex, dimension(MMAXD,MMAXD,NSPIND) :: VLDAU

   !  ..
   data IPRINT /1/
   data FACTOR /1.D0/  ! if this is 1. then: n*(n-1) in Edc and potential
   ! if this is 0. then: n**2    in Edc and potential


   write (1337,'(/,79(1H#),/,16X,A,/,79(1H#))') &
      'LDA+U: Calculating interaction potential VLDAU'
   !----------------------------------------------------------------------------
   do IT = 1,NTLDAU
      I1 = ITLDAU(IT)
      !-------------------------------------------------------------------------
      if ( LOPT(I1).GE.0 ) then
         call RINIT(MMAXD*MMAXD*NSPIND,DENMAT(1,1,1))
         MMAX = 2*LOPT(I1) + 1
         write (1337,99001) I1,LOPT(I1)
         !----------------------------------------------------------------------
         ! Result is in real Ylm basis.
         ! It must be converted to complex Ylm basis:
         !----------------------------------------------------------------------
         if ( IPRINT.GT.1 ) write (1337,99002) 'Occupation matrix in REAL basis:'
         !----------------------------------------------------------------------
         do IS = 1,NSPIN
            IPOT = (I1-1)*NSPIN + IS
            if ( IPRINT.GT.1 ) then
               write (STR15,'(4X,"> ",A,I1)') 'ISPIN = ',IS
               call CMATSTR(STR15,15,DENMATC(1,1,IPOT),MMAXD,MMAX,0,0,0,1d-8,1337)
            end if
            !-------------------------------------------------------------------
            ! Convert DENMATC and DENMAT to complex spherical harmonics.
            !-------------------------------------------------------------------
            call RCLM(1,LOPT(I1),LMAXD,DENMATC(1,1,IPOT))
         end do
         !----------------------------------------------------------------------
         if ( IPRINT.GT.1 )   then
            write (1337,99002) 'Occupation matrix in COMPLEX basis:'
         endif
         DENTOT = 0.D0
         !----------------------------------------------------------------------
         do IS = 1,NSPIN
            IPOT = (I1-1)*NSPIN + IS
            if ( IPRINT.GT.1 ) then
               write (STR15,'(4X,"> ",A,I1)') 'ISPIN = ',IS
               call CMATSTR(STR15,15,DENMATC(1,1,IPOT),MMAXD,MMAX,0,0,0,1d-8,1337)
            end if
            !-------------------------------------------------------------------
            ! DENMAT is real: (imag(denmatc))
            !-------------------------------------------------------------------
            do M2 = 1,MMAX
               do M1 = 1,MMAX
                  ! DENMAT(M1,M2,IS) = DIMAG(DENMATC(M1,M2,IPOT))
                  ! in the new medthod, taking imaginary part included
                  DENMAT(M1,M2,IS) = (DENMATC(M1,M2,IPOT))
               end do
            end do
            !-------------------------------------------------------------------
            ! 2.  Calculate total occupation numbers:
            ! ntot_s = Sum_m n_{m,s,m,s}, ntot = n_1 + n_2
            !-------------------------------------------------------------------
            DENTOTS(IS) = 0.D0
            do MM = 1,MMAX
               DENTOTS(IS) = DENTOTS(IS) + DENMAT(MM,MM,IS)
            end do
            DENTOT = DENTOT + DENTOTS(IS)
            DENTOTS(IS) = DENTOTS(IS)/DFLOAT(3-NSPIN)
         end do
         !----------------------------------------------------------------------
         if ( IPRINT.GT.0 ) then
            write (1337,99002) 'Occupation matrix (real):'
            do IS = 1,NSPIN
               write(1337,99003) IS
               call RWRITE(DENMAT(1,1,IS),MMAXD,MMAX,1337)
               write(1337,99004) 'Trace     =',DENTOTS(IS)
            end do
            write(1337,99005) 'Spins sum =',DENTOT
         end if
         !----------------------------------------------------------------------
         ! In paramagnetic case the spin degeneracy has been accounted
         ! for by the weight DF in tmatrho.
         !----------------------------------------------------------------------
         call CINIT(MMAXD*MMAXD*NSPIND,VLDAU(1,1,1))
         do IS = 1,NSPIN
            !-------------------------------------------------------------------
            ! 3.  Use density matrix and Coulomb matrix ULDAU to calculate the
            ! interaction potential VLDAU
            ! 3a. First part (always diagonal in spin).
            !-------------------------------------------------------------------
            do M2 = 1,MMAX
               do M1 = 1,MMAX
                  CSUM = CZERO
                  do M4 = 1,MMAX
                     do M3 = 1,MMAX
                        CSUM2 = CZERO
                        do JS = 1,NSPIN
                           CSUM2 = CSUM2 + DENMAT(M3,M4,JS)
                        end do
                        CSUM = CSUM + ULDAU(M1,M2,M3,M4,I1)*CSUM2
                     end do
                  end do
                  VLDAU(M1,M2,IS) = VLDAU(M1,M2,IS) + CSUM
               end do
            end do
            !-------------------------------------------------------------------
            ! 3b. Second part (in fully rel. case not diagonal in spin; then this
            ! loop must be changed accordingly).
            !-------------------------------------------------------------------
            do M2 = 1,MMAX
               do M1 = 1,MMAX
                  CSUM = CZERO
                  do M4 = 1,MMAX
                     do M3 = 1,MMAX
                        CSUM = CSUM - ULDAU(M1,M4,M3,M2,I1)*   &
                              DENMAT(M3,M4,IS)/DFLOAT(3-NSPIN)
                     end do
                  end do
                  VLDAU(M1,M2,IS) = VLDAU(M1,M2,IS) + CSUM
               end do
            end do
            !-------------------------------------------------------------------
            ! 3c. Third part (always spin- and m-diagonal).
            !-------------------------------------------------------------------
            do M1 = 1,MMAX
               VLDAU(M1,M1,IS) = VLDAU(M1,M1,IS)- UEFF(I1)*(DENTOT-0.5D0*FACTOR) &
                                 + JEFF(I1)*(DENTOTS(IS)-0.5D0*FACTOR)
            end do
         end do
         !----------------------------------------------------------------------
         ! 4. Calculate total-energy corrections EU and EDC (double-counting).
         ! Then the correction is EU - EDC.
         ! L[LDA+U]=E[LDA]+E[U]-E[DC]
         !> @note EU,EDC initialised outside the routine
         !----------------------------------------------------------------------

         ! Here VLDAU is assumed spin-diagonal (contrary to the spin-orbit case).
         !c              DO M2 = 1,MMAX
         !c                 DO M1 = 1,MMAX
         !c                    EU(I1) = EU(I1) +
         !c    &                    DENMAT(M1,M2,IS) * DREAL(VLDAU(M1,M2,IS))
         !c                 END DO
         !c              END DO

         ! Relativistic case, see the paper H. Ebert et al., Sol. Stat. Comm. 127 (2003) 443
         ! Calculate EDC
         do IS=1,NSPIN
            EDC(I1) = EDC(I1) +JEFF(I1)*DENTOTS(IS)*(DENTOTS(IS)-FACTOR)
         enddo

         EDC(I1) = 0.5D0*(UEFF(I1)*DENTOT*(DENTOT-1.D0)-EDC(I1))

         ! Calculate EU
         do IS = 1,NSPIN
            do JS = 1,NSPIN
               do M4 = 1,MMAX
                  do M3 = 1,MMAX
                     do M2 = 1,MMAX
                        do M1 = 1,MMAX
                           EU(I1) = EU(I1) + DENMAT(M1,M2,IS)* &
                                    ULDAU(M1,M2,M3,M4,I1)*DENMAT(M3,M4,JS)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo

         do IS = 1,NSPIN
            do M4 = 1,MMAX
               do M3 = 1,MMAX
                  do M2 = 1,MMAX
                     do M1 = 1,MMAX
                        EU(I1) = EU(I1) - DENMAT(M1,M2,IS)* &
                                 ULDAU(M1,M4,M3,M2,I1)*DENMAT(M3,M4,IS)
                     enddo
                  enddo
               enddo
            enddo
         enddo

         EU(I1) = 0.5D0*EU(I1)

         !----------------------------------------------------------------------
         IF ( IPRINT.GT.0 ) WRITE (1337,99002)'Interaction potential in COMPLEX basis:'
         !----------------------------------------------------------------------
         do IS = 1,NSPIN
            if ( IPRINT.GT.0 ) then
               write (STR15,'(4X,"> ",A,I1)') 'ISPIN = ',IS
               call CMATSTR(STR15,15,VLDAU(1,1,IS),MMAXD,MMAX,0,0,0,1d-8,1337)
            end if
            !-------------------------------------------------------------------
            ! 5.  Transform VLDAU into real spherical harmonics basis
            !-------------------------------------------------------------------
            call RCLM(2,LOPT(I1),LMAXD,VLDAU(1,1,IS))
            !-------------------------------------------------------------------
            ! Copy transformed VLDAU to real WLDAU
            !
            ! Apply damping to the interaction matrix WLDAU ? Here not.
            !-------------------------------------------------------------------
            do M2 = 1,MMAX
               do M1 = 1,MMAX
                  WLDAU(M1,M2,IS,I1) = DREAL(VLDAU(M1,M2,IS))
               end do
            end do
         end do
         !----------------------------------------------------------------------
         if ( IPRINT.GT.0 ) then
            write (1337,99002) 'Interaction potential in REAL basis:'
            do IS = 1,NSPIN
               write (STR15,'(4X,"> ",A,I1)') 'ISPIN = ',IS
               call CMATSTR(STR15,15,VLDAU(1,1,IS),MMAXD,MMAX,0,0,0,1d-8,1337)
            end do
         end if
         !----------------------------------------------------------------------
         write (1337,99002) 'Interaction potential (real):'
         do IS = 1,NSPIN
            write(1337,99003) IS
            call RWRITE(WLDAU(1,1,IS,I1),MMAXD,MMAX,1337)
         end do
         write(1337,*)
         !----------------------------------------------------------------------
         ! Corrections in total energy:
         ! Write out corrections on energy:
         ! E[LDA+U] = E[LDA] + EU - EDC
         !----------------------------------------------------------------------
         write(1337,99002) 'Corrections to the total energy:'
         write(1337,*)
         write(1337,99004) 'EU  =',EU(I1)
         write(1337,99004) 'Edc =',EDC(I1)
         write(1337,99006) 'E[LDA+U] = E[LDA] + EU - Edc'
      end if
      !-------------------------------------------------------------------------
   end do                    ! I1 = 1,NTLDAU
   !----------------------------------------------------------------------------
   99001 format(/,6X,65(1H=),/,6X,'Atom :',I3,' (l =',I2,')',/,6X,18(1H=))
   99002 format(8X,'* ',A)
   99003 format(/,15X,'> ISPIN =',I1)
   99004 format(10X,A,F10.6)
   99005 format(10X,21(1H-),/,10X,A,F10.6,/,10X,60(1H-),/)
   99006 format(27X,A,/)
end subroutine WMATLDAU

   !*==rwrite.f    processed by SPAG 6.05Rc at 16:58 on 22 Dec 2004
!-------------------------------------------------------------------------------
! SUBROUTINE: RWRITE
!> @brief Auxiliary subroutine to write the entries of the different potentials
!-------------------------------------------------------------------------------
subroutine RWRITE(Z,MMAXD,MMAX,IFILE)

      implicit none
      !.. Input variables
      integer, intent(in) :: IFILE
      integer, intent(in) :: MMAX
      integer, intent(in) :: MMAXD
      double precision, dimension(MMAXD,MMAXD), intent(in) :: Z
      !.. Local variables
      integer :: M1,M2
      !
      write(IFILE,99000)
      do M2 = 1,MMAX
         write (IFILE,99001) (Z(M1,M2),M1=1,MIN(MMAX,7))
      end do
      write(IFILE,99000)
      99000 format(10X,60(1H-))
      99001 format(10X,7F10.6)
end subroutine RWRITE
