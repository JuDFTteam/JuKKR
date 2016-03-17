       SUBROUTINE RHOOUTNEW(NSRA,LMMAXD,LMMAXSO,LMAX,GMATLL,EK,LMPOTD,
     +                      DF,NPAN_TOT,NCHEB,CLEB,ICLEB,IEND,
     +                      IRMDNEW,NRMAXD,THETASNEW,IFUNM,RNEW,IMT1,
     +                      LMSP,RLL,SLL,RLLLEFT,SLLLEFT,
     +                      CDEN,CDENLM,CDENNS,RHO2NSC,CORBITAL,
     +                      GFLLE_PART,RPAN_INTERVALL,IPAN_INTERVALL)
       IMPLICIT NONE
       include 'inc.p'
       INTEGER NSRA,LMMAXD,LMMAXSO,LMAX,LMPOTD,IEND,CORBITAL
       INTEGER IRMDNEW,NRMAXD,NPAN_TOT,NCHEB,IMT1
       DOUBLE COMPLEX CZERO,CONE
       PARAMETER (CZERO=(0d0,0d0),CONE=(1d0,0d0))
       DOUBLE COMPLEX EK,DF,CLTDF
       DOUBLE PRECISION CLEB(*)
       INTEGER ICLEB(NCLEB,4),IFUNM(*),LMSP(*),IPAN_INTERVALL(0:NTOTD)
       DOUBLE COMPLEX GMATLL(LMMAXSO,LMMAXSO)
       INTEGER IR,JSPIN,LM1,LM2,LM3,M1,L1,J,IFUN
       DOUBLE PRECISION C0LL,RPAN_INTERVALL(0:NTOTD)
       DOUBLE PRECISION RNEW(NRMAXD),
     +                THETASNEW(NTOTD*(NCHEBD+1),NFUND)
       DOUBLE COMPLEX RLL(NSRA*LMMAXSO,LMMAXSO,IRMDNEW),
     +                SLL(NSRA*LMMAXSO,LMMAXSO,IRMDNEW),
     +                RLLLEFT(NSRA*LMMAXSO,LMMAXSO,IRMDNEW),
     +                SLLLEFT(NSRA*LMMAXSO,LMMAXSO,IRMDNEW)
       DOUBLE COMPLEX, ALLOCATABLE :: WR(:,:,:),QNSI(:,:),PNSI(:,:),
     +                CWR(:)                                        ! lmlm-dos
       INTEGER LMSHIFT1(4),LMSHIFT2(4)
       DOUBLE COMPLEX CDEN(IRMDNEW,0:LMAX,4),
     +                CDENLM(IRMDNEW,LMMAXD,4),
     +                RHO2NSC(IRMDNEW,LMPOTD,4),
     +                CDENNS(IRMDNEW,4),
     +                GFLLE_PART(LMMAXSO,LMMAXSO)                         ! lmlm-dos
       DOUBLE COMPLEX LOPERATOR(LMMAXSO,LMMAXSO,3)
       LOGICAL TEST,OPT
       EXTERNAL TEST,OPT
       ALLOCATE(WR(LMMAXSO,LMMAXSO,IRMDNEW))
       ALLOCATE(CWR(IRMDNEW))
       ALLOCATE(QNSI(LMMAXSO,LMMAXSO))
       ALLOCATE(PNSI(LMMAXSO,LMMAXSO))
       WR=CZERO
       CWR=CZERO
       QNSI=CZERO
       PNSI=CZERO
c set LMSHIFT value which is need to construct CDEN
       LMSHIFT1(1)=0
       LMSHIFT1(2)=LMMAXD
       LMSHIFT1(3)=0
       LMSHIFT1(4)=LMMAXD
       LMSHIFT2(1)=0
       LMSHIFT2(2)=LMMAXD
       LMSHIFT2(3)=LMMAXD
       LMSHIFT2(4)=0
       
c for orbital moment
       IF (CORBITAL.NE.0) THEN
        CALL CALC_ORBITALMOMENT(LMAXD,LMMAXSO,LOPERATOR)
       ENDIF      

       C0LL=1d0/SQRT(16d0*ATAN(1d0))
       CDEN=CZERO
       CDENLM=CZERO

       DO IR = 1,IRMDNEW

        DO LM1 = 1,LMMAXSO
         DO LM2 = 1,LMMAXSO
          QNSI(LM1,LM2)=SLLLEFT(LM1,LM2,IR)
c          PNSI(LM1,LM2)=RLL(LM1,LM2,IR)
          PNSI(LM1,LM2)=RLLLEFT(LM1,LM2,IR)
         ENDDO
        ENDDO
c        CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,
c     +             LMMAXSO,GMATLL,LMMAXSO,EK,QNSI,LMMAXSO)
        CALL ZGEMM('N','T',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,
     +             LMMAXSO,GMATLL,LMMAXSO,EK,QNSI,LMMAXSO)
        DO LM1 = 1,LMMAXSO
         DO LM2 = 1,LMMAXSO
          PNSI(LM1,LM2)=RLL(LM1,LM2,IR)
         ENDDO
        ENDDO
        CALL ZGEMM('N','T',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,
     +             LMMAXSO,QNSI,LMMAXSO,CZERO,WR(1,1,IR),LMMAXSO)

        IF (NSRA.EQ.2) THEN
        DO LM1 = 1,LMMAXSO
         DO LM2 = 1,LMMAXSO
c          QNSI(LM1,LM2)=SLLLEFT(LM1+LMMAXSO,LM2,IR)
          QNSI(LM1,LM2)=-SLLLEFT(LM1+LMMAXSO,LM2,IR)
c          PNSI(LM1,LM2)=RLLLEFT(LM1+LMMAXSO,LM2,IR)
          PNSI(LM1,LM2)=-RLLLEFT(LM1+LMMAXSO,LM2,IR)
         ENDDO
        ENDDO
c        CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,
c     +             LMMAXSO,GMATLL,LMMAXSO,EK,QNSI,LMMAXSO)
        CALL ZGEMM('N','T',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,
     +             LMMAXSO,GMATLL,LMMAXSO,EK,QNSI,LMMAXSO)
        DO LM1 = 1,LMMAXSO
         DO LM2 = 1,LMMAXSO
          PNSI(LM1,LM2)=RLL(LM1+LMMAXSO,LM2,IR)
         ENDDO
        ENDDO
        CALL ZGEMM('N','T',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,
     +             LMMAXSO,QNSI,LMMAXSO,CONE,WR(1,1,IR),LMMAXSO)
        ENDIF

c for orbital moment
        IF (CORBITAL.NE.0) THEN
        CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,
     +             LOPERATOR(1,1,CORBITAL),LMMAXSO,WR(1,1,IR),
     +             LMMAXSO,CZERO,PNSI,LMMAXSO)
        DO LM1=1,LMMAXSO
         DO LM2=1,LMMAXSO
          WR(LM1,LM2,IR)=PNSI(LM1,LM2)
         ENDDO
        ENDDO
        ENDIF
 
        DO JSPIN = 1,4
         DO LM1 = 1,LMMAXD
          DO LM2 = 1,LM1-1
           WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR)=
     +           WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR)+
     +           WR(LM2+LMSHIFT1(JSPIN),LM1+LMSHIFT2(JSPIN),IR)
          ENDDO
         ENDDO
        ENDDO ! JSPIN

       ENDDO !IR




       IF (OPT('lmlm-dos')) THEN                                                          ! lmlm-dos
c Integrate only up to muffin-tin radius.                                                 ! lmlm-dos
       GFLLE_PART = CZERO                                                                 ! lmlm-dos
       DO LM2 = 1,LMMAXSO                                                                 ! lmlm-dos
         DO LM1 = 1,LMMAXSO                                                               ! lmlm-dos
c For integration up to MT radius do this:                                                ! lmlm-dos
!              CWR(1:IMT1) = WR(LM1,LM2,1:IMT1)                                             ! lmlm-dos
!              CWR(IMT1+1:IRMDNEW) = CZERO                                                  ! lmlm-dos
!              CALL INTCHEB_CELL(CWR,GFLLE_PART(LM1,LM2),RPAN_INTERVALL,                    ! lmlm-dos
!      +                            IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)                  ! lmlm-dos
c For full cell integration replace loop content with this:                               ! lmlm-dos
            CWR(1:IRMDNEW) = WR(LM1,LM2,1:IRMDNEW)                                        ! lmlm-dos
               DO IR=IMT1+1,IRMDNEW                                                       ! lmlm-dos
                  CWR(IR) = CWR(IR)*THETASNEW(IR,1)*C0LL                                  ! lmlm-dos
               ENDDO                                                                      ! lmlm-dos
             CALL INTCHEB_CELL(CWR,GFLLE_PART(LM1,LM2),RPAN_INTERVALL,                    ! lmlm-dos
     +                            IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)                  ! lmlm-dos
          ENDDO                                                                           ! lmlm-dos
       ENDDO                                                                              ! lmlm-dos
       ENDIF  ! OPT('lmlm-dos')


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


c first calculate the spherical symmetric contribution

       DO L1 = 0,LMAX

        DO M1 = -L1,L1
         LM1 = L1*(L1+1)+M1+1
         DO IR = 1,IRMDNEW
          DO JSPIN=1,4
           CDEN(IR,L1,JSPIN) = CDEN(IR,L1,JSPIN)+
     +          WR(LM1+LMSHIFT1(JSPIN),LM1+LMSHIFT2(JSPIN),IR)
           CDENLM(IR,LM1,JSPIN) =
     +          WR(LM1+LMSHIFT1(JSPIN),LM1+LMSHIFT2(JSPIN),IR)
          ENDDO ! JPSIN
         ENDDO ! IR
        ENDDO ! M1

        DO JSPIN = 1,4
         DO IR = 1,IRMDNEW
          RHO2NSC(IR,1,JSPIN) = RHO2NSC(IR,1,JSPIN)+
     +                               C0LL*(CDEN(IR,L1,JSPIN)*DF)
         ENDDO ! IR
        
         DO IR=IMT1+1,IRMDNEW
          CDEN(IR,L1,JSPIN) = CDEN(IR,L1,JSPIN)*THETASNEW(IR,1)*C0LL
          
          DO M1 = -L1,L1
           LM1 = L1*(L1+1)+M1+1
           CDENLM(IR,LM1,JSPIN) = CDENLM(IR,LM1,JSPIN)
     +                              *THETASNEW(IR,1)*C0LL
          ENDDO ! M1
         ENDDO ! IR

        ENDDO ! JSPIN

       ENDDO ! L1

       CDENNS=CZERO
       
       DO J = 1,IEND  
        LM1 = ICLEB(J,1)
        LM2 = ICLEB(J,2)
        LM3 = ICLEB(J,3)
        CLTDF = DF*CLEB(J)
        
        DO JSPIN = 1,4
         DO IR = 1,IRMDNEW
          RHO2NSC(IR,LM3,JSPIN) = RHO2NSC(IR,LM3,JSPIN) + 
     +          (CLTDF*WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR))
         ENDDO
       
         IF (LMSP(LM3).GT.0) THEN
          IFUN = IFUNM(LM3)
          DO IR=IMT1+1,IRMDNEW
           CDENNS(IR,JSPIN) = CDENNS(IR,JSPIN)+
     +           CLEB(J)*WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR)*
     +           THETASNEW(IR,IFUN)
          if(jspin==1) write(123456789,'(4I5,E15.7,I5,20000E15.7)') 
     &                    j,lm1,lm2,lm3,cleb(j),ifun,thetasnew(ir,ifun),
     &          RLL(LM1+LMMAXSO,LM2,IR),RLLleft(LM2,mod(LM3,lmmaxso),IR)
          ENDDO
         ENDIF
        ENDDO ! JSPIN
       ENDDO ! J   


       DEALLOCATE(WR)
       DEALLOCATE(CWR)
       DEALLOCATE(QNSI) 
       DEALLOCATE(PNSI) 
       END
