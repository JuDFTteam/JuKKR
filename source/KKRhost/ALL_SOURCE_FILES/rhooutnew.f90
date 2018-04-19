!-------------------------------------------------------------------------------
! SUBROUTINE: RHOOUTNEW
!> @note -Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine RHOOUTNEW(NSRA,LMMAXD,LMMAXSO,LMAX,GMATLL,EK,LMPOT,&
   DF,NPAN_TOT,NCHEB,CLEB,ICLEB,IEND,                          &
   IRMDNEW,THETASNEW,IFUNM,IMT1,                               &
   LMSP,RLL,RLLLEFT,SLLLEFT,                                   &
   CDEN,CDENLM,CDENNS,RHO2NSC,CORBITAL,                        &
   GFLLE_PART,RPAN_INTERVALL,IPAN_INTERVALL,NTOTD)

   use Constants
   use Profiling
   use global_variables

   implicit none

   integer, intent(in) :: NSRA
   integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
   integer, intent(in) :: IEND      !< Number of nonzero gaunt coefficients
   integer, intent(in) :: IMT1
   integer, intent(in) :: NTOTD
   integer, intent(in) :: NCHEB     !< Number of Chebychev pannels for the new solver
   integer, intent(in) :: LMPOT     !< (LPOT+1)**2
   integer, intent(in) :: LMMAXD    !< (KREL+KORBIT+1)(LMAX+1)^2
   integer, intent(in) :: LMMAXSO   !< 2*LMMAXD)
   integer, intent(in) :: IRMDNEW
   integer, intent(in) :: CORBITAL
   integer, intent(in) :: NPAN_TOT
   double complex, intent(in) :: EK
   double complex, intent(in) :: DF
   integer, dimension(*), intent(in)         :: LMSP !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
   integer, dimension(*), intent(in)         :: IFUNM
   integer, dimension(0:NTOTD), intent(in)   :: IPAN_INTERVALL
   integer, dimension(NCLEB,4), intent(in)   :: ICLEB   !< Pointer array
   double precision, dimension(*), intent(in)         :: CLEB !< GAUNT coefficients (GAUNT)
   double precision, dimension(0:NTOTD), intent(in)   :: RPAN_INTERVALL
   double precision, dimension(NTOTD*(NCHEB+1),NFUND), intent(in) :: THETASNEW
   double complex, dimension(LMMAXSO,LMMAXSO), intent(in) :: GMATLL  !< GMATLL = diagonal elements of the G matrix (system)
   ! Note that SLL is not needed for calculation of density, only needed for calculation of Green function
   double complex, dimension(NSRA*LMMAXSO,LMMAXSO,IRMDNEW), intent(in) :: RLL
   double complex, dimension(NSRA*LMMAXSO,LMMAXSO,IRMDNEW), intent(in) :: RLLLEFT
   double complex, dimension(NSRA*LMMAXSO,LMMAXSO,IRMDNEW), intent(in) :: SLLLEFT

   ! .. Output variables
   double complex, dimension(IRMDNEW,4), intent(out)        :: CDENNS
   double complex, dimension(LMMAXSO,LMMAXSO), intent(out)  :: GFLLE_PART  ! lmlm-dos
   double complex, dimension(IRMDNEW,0:LMAX,4), intent(out) :: CDEN
   double complex, dimension(IRMDNEW,LMMAXD,4), intent(out) :: CDENLM

   ! .. In/Out variables
   double complex, dimension(IRMDNEW,LMPOT,4), intent(inout) :: RHO2NSC

   ! .. Local variables
   integer :: IR,JSPIN,LM1,LM2,LM3,M1,L1,J,IFUN
   integer :: i_stat, i_all
   double precision :: C0LL
   double complex :: CLTDF
   integer, dimension(4) :: LMSHIFT1
   integer, dimension(4) :: LMSHIFT2
   double complex, dimension(LMMAXSO,LMMAXSO,3) :: LOPERATOR

   ! .. Local allocatable arrays
   double complex, dimension(:), allocatable :: CWR      ! lmlm-dos
   double complex, dimension(:,:), allocatable :: QNSI
   double complex, dimension(:,:), allocatable :: PNSI
   double complex, dimension(:,:,:), allocatable :: WR
   double complex, dimension(:,:,:), allocatable :: WR1  ! LDAU
   ! .. External routines
   logical :: TEST,OPT
   external :: TEST,OPT

   allocate(WR(LMMAXSO,LMMAXSO,IRMDNEW),stat=i_stat)
   call memocc(i_stat,product(shape(WR))*kind(WR),'WR','RHOOUTNEW')
   WR=CZERO
   allocate(CWR(IRMDNEW),stat=i_stat)
   call memocc(i_stat,product(shape(CWR))*kind(CWR),'CWR','RHOOUTNEW')
   CWR=CZERO
   allocate(WR1(LMMAXSO,LMMAXSO,IRMDNEW),stat=i_stat)
   call memocc(i_stat,product(shape(WR1))*kind(WR1),'WR1','RHOOUTNEW')
   WR1=CZERO
   allocate(QNSI(LMMAXSO,LMMAXSO),stat=i_stat)
   call memocc(i_stat,product(shape(QNSI))*kind(QNSI),'QNSI','RHOOUTNEW')
   QNSI=CZERO
   allocate(PNSI(LMMAXSO,LMMAXSO),stat=i_stat)
   call memocc(i_stat,product(shape(PNSI))*kind(PNSI),'PNSI','RHOOUTNEW')
   PNSI=CZERO

   ! set LMSHIFT value which is need to construct CDEN
   LMSHIFT1(1)=0
   LMSHIFT1(2)=LMMAXD
   LMSHIFT1(3)=0
   LMSHIFT1(4)=LMMAXD
   LMSHIFT2(1)=0
   LMSHIFT2(2)=LMMAXD
   LMSHIFT2(3)=LMMAXD
   LMSHIFT2(4)=0

   ! for orbital moment
   if (CORBITAL.NE.0) then
      call CALC_ORBITALMOMENT(LMAX,LMMAXSO,LOPERATOR)
   endif

   C0LL=1d0/SQRT(16d0*ATAN(1d0))
   CDEN=CZERO
   CDENLM=CZERO

   do IR = 1,IRMDNEW
      do LM1 = 1,LMMAXSO
         do LM2 = 1,LMMAXSO
            QNSI(LM1,LM2)=SLLLEFT(LM1,LM2,IR)
            !          PNSI(LM1,LM2)=RLL(LM1,LM2,IR)
            PNSI(LM1,LM2)=RLLLEFT(LM1,LM2,IR)
         enddo
      enddo
      !        CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,
      !     +             LMMAXSO,GMATLL,LMMAXSO,EK,QNSI,LMMAXSO)
      call ZGEMM('N','T',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,LMMAXSO,GMATLL,&
         LMMAXSO,EK,QNSI,LMMAXSO)
      do LM1 = 1,LMMAXSO
         do LM2 = 1,LMMAXSO
            PNSI(LM1,LM2)=RLL(LM1,LM2,IR)
         enddo
      enddo
      call ZGEMM('N','T',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,LMMAXSO,QNSI,&
         LMMAXSO,CZERO,WR(1,1,IR),LMMAXSO)
      !
      if (NSRA.EQ.2) then
         do LM1 = 1,LMMAXSO
            do LM2 = 1,LMMAXSO
               !          QNSI(LM1,LM2)=SLLLEFT(LM1+LMMAXSO,LM2,IR)
               QNSI(LM1,LM2)=-SLLLEFT(LM1+LMMAXSO,LM2,IR)
               !          PNSI(LM1,LM2)=RLLLEFT(LM1+LMMAXSO,LM2,IR)
               PNSI(LM1,LM2)=-RLLLEFT(LM1+LMMAXSO,LM2,IR)
            enddo
         enddo
         !        CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,
         !     +             LMMAXSO,GMATLL,LMMAXSO,EK,QNSI,LMMAXSO)
         call ZGEMM('N','T',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,LMMAXSO,GMATLL,&
            LMMAXSO,EK,QNSI,LMMAXSO)
         do LM1 = 1,LMMAXSO
            do LM2 = 1,LMMAXSO
               PNSI(LM1,LM2)=RLL(LM1+LMMAXSO,LM2,IR)
            enddo
         enddo
         call ZGEMM('N','T',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,LMMAXSO,QNSI,&
            LMMAXSO,CONE,WR(1,1,IR),LMMAXSO)
      endif
      !
      ! For orbital moment
      if (CORBITAL.NE.0) then
         call ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,&
            LOPERATOR(1,1,CORBITAL),LMMAXSO,WR(1,1,IR),LMMAXSO,CZERO,PNSI,LMMAXSO)
         do LM1=1,LMMAXSO
            do LM2=1,LMMAXSO
               WR(LM1,LM2,IR)=PNSI(LM1,LM2)
            enddo
         enddo
      endif
      do LM1=1,LMMAXSO
         do LM2=1,LMMAXSO
            WR1(LM1,LM2,IR)=WR(LM1,LM2,IR)
         enddo
      enddo
      do LM1=1,LMMAXSO
         do LM2=1,LM1-1
            WR1(LM1,LM2,IR)=WR1(LM1,LM2,IR)+WR1(LM2,LM1,IR)
         enddo
      enddo
      !
      do JSPIN = 1,4
         do LM1 = 1,LMMAXD
            do LM2 = 1,LM1-1
               WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR)=    &
                  WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR)+ &
                  WR(LM2+LMSHIFT1(JSPIN),LM1+LMSHIFT2(JSPIN),IR)
            enddo
         enddo
      enddo ! JSPIN
   enddo !IR

   ! IF lmdos or LDAU
   if (OPT('lmlm-dos').OR.OPT('LDA+U   ')) then                                  ! lmlm-dos
      ! Integrate only up to muffin-tin radius.                                  ! lmlm-dos
      GFLLE_PART = CZERO                                                         ! lmlm-dos
      do LM2 = 1,LMMAXSO                                                         ! lmlm-dos
         do LM1 = 1,LMMAXSO                                                      ! lmlm-dos
            ! For integration up to MT radius do this:                           ! lmlm-dos
            ! CWR(1:IMT1) = WR(LM1,LM2,1:IMT1)                                   ! lmlm-dos
            ! CWR(IMT1+1:IRMDNEW) = CZERO                                        ! lmlm-dos
            ! CALL INTCHEB_CELL(CWR,GFLLE_PART(LM1,LM2),RPAN_INTERVALL,&         ! lmlm-dos
            !     IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)                         ! lmlm-dos
            ! For full cell integration replace loop content with this:          ! lmlm-dos
            CWR(1:IRMDNEW) = WR1(LM1,LM2,1:IRMDNEW)                              ! lmlm-dos
            ! If LDAU, integrate only up to MT
            do IR=IMT1+1,IRMDNEW
               if (OPT('LDA+U   ')) then
                  CWR(IR)=CZERO  ! LDAU
               else
                  CWR(IR) = CWR(IR)*THETASNEW(IR,1)*C0LL                         ! lmlm-dos
               endif
            enddo
            call INTCHEB_CELL(CWR,GFLLE_PART(LM1,LM2),RPAN_INTERVALL,&
               IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)
         enddo
      enddo
   endif  ! OPT('lmlm-dos').OR.OPT('LDA+U   ')
   !
   !      DO IR = 1,IRMDNEW
   !       DO JSPIN = 1,4
   !        DO LM1 = 1,LMMAXD
   !         DO LM2 = 1,LM1-1
   !          WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR)=
   !    +           WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR)+
   !    +           WR(LM2+LMSHIFT1(JSPIN),LM1+LMSHIFT2(JSPIN),IR)
   !         ENDDO
   !        ENDDO
   !       ENDDO ! JSPIN
   !      ENDDO !IR
   !
   ! First calculate the spherical symmetric contribution
   !
   do L1 = 0,LMAX
      do M1 = -L1,L1
         LM1 = L1*(L1+1)+M1+1
         do IR = 1,IRMDNEW
            do JSPIN=1,4
               CDEN(IR,L1,JSPIN) = CDEN(IR,L1,JSPIN)+&
                  WR(LM1+LMSHIFT1(JSPIN),LM1+LMSHIFT2(JSPIN),IR)
               CDENLM(IR,LM1,JSPIN) =WR(LM1+LMSHIFT1(JSPIN),LM1+LMSHIFT2(JSPIN),IR)
            enddo ! JPSIN
         enddo ! IR
      enddo ! M1
      !
      do JSPIN = 1,4
         do IR = 1,IRMDNEW
            RHO2NSC(IR,1,JSPIN) = RHO2NSC(IR,1,JSPIN)+C0LL*(CDEN(IR,L1,JSPIN)*DF)
         enddo ! IR
         !
         do IR=IMT1+1,IRMDNEW
            CDEN(IR,L1,JSPIN) = CDEN(IR,L1,JSPIN)*THETASNEW(IR,1)*C0LL
            do M1 = -L1,L1
               LM1 = L1*(L1+1)+M1+1
               CDENLM(IR,LM1,JSPIN) = CDENLM(IR,LM1,JSPIN)*THETASNEW(IR,1)*C0LL
            enddo ! M1
         enddo ! IR
      enddo ! JSPIN
   enddo ! L1
   !
   CDENNS=CZERO
   !
   do J = 1,IEND
      LM1 = ICLEB(J,1)
      LM2 = ICLEB(J,2)
      LM3 = ICLEB(J,3)
      CLTDF = DF*CLEB(J)
      do JSPIN = 1,4
         do IR = 1,IRMDNEW
            RHO2NSC(IR,LM3,JSPIN) = RHO2NSC(IR,LM3,JSPIN)+&
               (CLTDF*WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR))
         enddo
         !
         if (LMSP(LM3).GT.0) then
            IFUN = IFUNM(LM3)
            do IR=IMT1+1,IRMDNEW
               CDENNS(IR,JSPIN) = CDENNS(IR,JSPIN)+&
                  CLEB(J)*WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR)*THETASNEW(IR,IFUN)
            enddo
         endif
      enddo ! JSPIN
   enddo ! J

   i_all=-product(shape(WR))*kind(WR)
   deallocate(WR,stat=i_stat)
   call memocc(i_stat,i_all,'WR','RHOOUTNEW')
   i_all=-product(shape(WR1))*kind(WR1)
   deallocate(WR1,stat=i_stat)
   call memocc(i_stat,i_all,'WR1','RHOOUTNEW')
   i_all=-product(shape(CWR))*kind(CWR)
   deallocate(CWR,stat=i_stat)
   call memocc(i_stat,i_all,'CWR','RHOOUTNEW')
   i_all=-product(shape(QNSI))*kind(QNSI)
   deallocate(QNSI,stat=i_stat)
   call memocc(i_stat,i_all,'QNSI','RHOOUTNEW')
   i_all=-product(shape(PNSI))*kind(PNSI)
   deallocate(PNSI,stat=i_stat)
   call memocc(i_stat,i_all,'PNSI','RHOOUTNEW')

end subroutine RHOOUTNEW
