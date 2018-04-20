!-------------------------------------------------------------------------------
! SUBROUTINE: PNSQNS
!> @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine PNSQNS(AR,CR,DR,DRDI,EK,ICST,PZ,QZ,FZ,SZ,PNS,QNS,NSRA, &
   VINS,IPAN,IRMIN,IRCUT,CLEB,ICLEB,IEND,LOFLM,LKONV,             & ! Added IRMIN 1.7.2014
   IDOLDAU,LOPT,LMLO,LMHI,WLDAU,WLDAUAV,CUTOFF,MMAXD,LMPOT,IRMIND,&
   LMMAXD,IRM,LMAX)

   use global_variables

   implicit none
   !
   ! .. Input variables
   integer, intent(in) :: IRM       !< Maximum number of radial points
   integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
   integer, intent(in) :: ICST      !< Number of Born approximation
   integer, intent(in) :: IEND      !< Number of nonzero gaunt coefficients
   integer, intent(in) :: IPAN      !< Number of panels in non-MT-region
   integer, intent(in) :: LMLO
   integer, intent(in) :: LMHI
   integer, intent(in) :: LOPT      !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
   integer, intent(in) :: NSRA
   integer, intent(in) :: LKONV
   integer, intent(in) :: IRMIN     !< Max R for spherical treatment
   integer, intent(in) :: MMAXD     !< 2*LMAX+1
   integer, intent(in) :: LMPOT     !< (LPOT+1)**2
   integer, intent(in) :: IRMIND    !< IRM-IRNSD
   integer, intent(in) :: LMMAXD    !< (KREL+KORBIT+1)(LMAX+1)^2
   integer, intent(in) :: IDOLDAU   !< flag to perform LDA+U
   double precision, intent(in) :: WLDAUAV
   double complex, intent(in) :: EK
   integer, dimension(0:IPAND), intent(in)   :: IRCUT   !< R points of panel borders
   integer, dimension(*), intent(in)         :: LOFLM   !< l of lm=(l,m) (GAUNT)
   integer, dimension(NCLEB,4), intent(in)   :: ICLEB   !< Pointer array
   double precision, dimension(IRM), intent(in) :: DRDI  !< Derivative dr/di
   double precision, dimension(IRM), intent(in) :: CUTOFF
   double precision, dimension(NCLEB,2), intent(in)         :: CLEB  !< GAUNT coefficients (GAUNT)
   double precision, dimension(IRMIN:IRM,LMPOT), intent(in) :: VINS  !< Non-spherical part of the potential
   double precision, dimension(MMAXD,MMAXD), intent(in)     :: WLDAU !< potential matrix
   double complex, dimension(IRM,0:LMAX), intent(in)     :: FZ
   double complex, dimension(IRM,0:LMAX), intent(in)     :: QZ
   double complex, dimension(IRM,0:LMAX), intent(in)     :: SZ
   double complex, dimension(IRM,0:LMAX), intent(in)     :: PZ
   double complex, dimension(LMMAXD,LMMAXD), intent(in)  :: DR
   double complex, dimension(LMMAXD,LMMAXD), intent(in)  :: AR
   double complex, dimension(LMMAXD,LMMAXD), intent(in)  :: CR
   double complex, dimension(LMMAXD,LMMAXD,IRMIN:IRM,2), intent(in) :: PNS
   double complex, dimension(LMMAXD,LMMAXD,IRMIN:IRM,2), intent(in) :: QNS
   ! .. Local Scalars
   integer :: I,LM1,LM2,LMMKONV,M1,M2,IR,IRMAX
   ! .. Local Arrays
   double precision, dimension(LMMAXD,LMMAXD,IRMIN:IRM) :: VNSPLL
   double complex, dimension(LMMAXD) :: EFAC
   double complex, dimension(LMMAXD,LMMAXD) :: TMATLL
   double complex, dimension(LMMAXD,LMMAXD,IRMIN:IRM) :: DMAT
   double complex, dimension(LMMAXD,LMMAXD,IRMIN:IRM) :: CMAT
   double complex, dimension(LMMAXD,IRMIN:IRM,2)      :: PZLM
   double complex, dimension(LMMAXD,IRMIN:IRM,2)      :: QZLM
   double complex, dimension(LMMAXD,IRMIN:IRM,2)      :: PZEKDR
   double complex, dimension(LMMAXD,IRMIN:IRM,2)      :: QZEKDR
   ! .. External Subroutines
   external :: IRWNS,REGNS,VLLNS,WFTSCA
   !
   IRMAX = IRCUT(IPAN)                                    ! Added IRMAX 1.7.2014
   call VLLNS(VNSPLL,VINS,CLEB,ICLEB,IEND,IRM,NCLEB,LMPOT,IRMIND,LMMAXD)
   if (LKONV.NE.LMAX) then
      LMMKONV = (LKONV+1)* (LKONV+1)
      do LM1 = 1,LMMAXD
         do LM2 = LMMKONV + 1,LMMAXD
            do I = IRMIND,IRM
               VNSPLL(LM2,LM1,I) = 0.0D0
               VNSPLL(LM1,LM2,I) = 0.0D0
            end do
         end do
      end do
   else
      LMMKONV = LMMAXD
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! LDA+U
   ! Add WLDAU to non-spherical porential VINS in case of LDA+U
   ! Use the average wldau (=wldauav) and the deviation
   ! of wldau from this. Use the deviation in the Born series
   ! for the non-spherical wavefunction, while the average is
   ! used for the spherical wavefunction.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if ( IDOLDAU.EQ.1.AND.LOPT.GE.0 ) then
      do IR = IRMIND,IRM
         !----------------------------------------------------------------------
         ! First add wldau to all elements
         !----------------------------------------------------------------------
         do LM2 = LMLO,LMHI
            M2 = LM2 - LMLO + 1
            do LM1 = LMLO,LMHI
               M1 = LM1 - LMLO + 1
               VNSPLL(LM1,LM2,IR) =  VNSPLL(LM1,LM2,IR)+ WLDAU(M1,M2) * CUTOFF(IR)
            enddo
            !-------------------------------------------------------------------
            ! and then subtract average from diag. elements
            !-------------------------------------------------------------------
            VNSPLL(LM2,LM2,IR) =  VNSPLL(LM2,LM2,IR) - WLDAUAV * CUTOFF(IR)
         enddo
      enddo
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! LDA+U
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !----------------------------------------------------------------------------
   ! Get wfts of same magnitude by scaling with efac
   !----------------------------------------------------------------------------
   call WFTSCA(DRDI,EFAC,PZ,QZ,FZ,SZ,NSRA,PZLM,QZLM,PZEKDR,QZEKDR,&
      EK,LOFLM,IRMIND,IRM,IRMIN,IRMAX,LMAX,LMMAXD)      ! Added IRMIN,IRMAX 1.7.2014
   !----------------------------------------------------------------------------
   ! Determine the irregular non sph. wavefunction
   !----------------------------------------------------------------------------
   call IRWNS(CR,DR,EFAC,QNS,VNSPLL,ICST,IPAN,IRCUT,NSRA,PZLM,QZLM,  &
      PZEKDR,QZEKDR,QNS(1,1,IRMIND,1),CMAT,                          &
      QNS(1,1,IRMIND,2),DMAT,IRMIND,IRM,IRMIN,IRMAX,                 & ! Added IRMIN,IRMAX 1.7.2014
      IPAND,LMMAXD)
   !----------------------------------------------------------------------------
   ! Determine the regular non sph. wavefunction
   !----------------------------------------------------------------------------
   call REGNS(AR,TMATLL,EFAC,PNS,VNSPLL,ICST,IPAN,IRCUT,PZLM,QZLM,&
      PZEKDR,QZEKDR,EK,PNS(1,1,IRMIND,1),CMAT,                    &
      PNS(1,1,IRMIND,2),DMAT,NSRA,IRMIND,IRM,IRMIN,IRMAX,         & ! Added IRMIN,IRMAX 1.7.2014
      IPAND,LMMAXD)
   !
   return

end subroutine PNSQNS
