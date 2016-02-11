      SUBROUTINE WRITEGREEN(ALAT,NSPIN,MAXMESH,
     +                    IE,IELAST,E,NSHELL,NBX,NBY,NBZ,BBYA,CBYA,
     +                    NAEZ,NATYP,CLS,EQINV,NACLS,RR,RBASIS,EZOA,
     +                    ATOM,RCLS,KAOEZ,ICC,BRAVAIS,RECBV,
     +                    LPOT,YR,WTYR,RIJ,IJEND,IATCONDL,IATCONDR,
     +                    NCONDPAIR,INTERF,INS,GINP,TMATLL,NATOMIMP,
     +                    ATOMIMP,RCLSIMP,IHOST)
!       use omp_lib
      IMPLICIT NONE
     
c-----------------------------------------------------------------
c write down the impurity green function, N. H. Long, Juelich, 05.2013
c-----------------------------------------------------------------
      INCLUDE 'inc.p'
      INCLUDE 'inc.cls'
      INTEGER NSYMAXD
      PARAMETER (NSYMAXD=48)
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ=(LMAXD+1)**2)
      INTEGER LMMAXSO
      PARAMETER (LMMAXSO=NSPD*LMAXSQ)
      INTEGER ALM,ALMSO
      PARAMETER (ALM=LMAXSQ*NAEZD,ALMSO=NAEZD*NSPD*LMAXSQ)
      DOUBLE PRECISION PI
      DOUBLE COMPLEX CONEM,CONE,CZERO,CI,FAC
      PARAMETER (CONEM=(-1d0,0d0),CONE=(1d0,0d0),
     +           CZERO=(0d0,0d0),CI=(0d0,1d0))
      
      DOUBLE COMPLEX E
      DOUBLE PRECISION ALAT,BBYA,CBYA
      INTEGER ICC,IE,IELAST,IJEND,INS,LPOT,MAXMESH,NAEZ,NATYP,
     +        NBX,NBY,NBZ,NSPIN,NSYMAT,ISYMAT
      INTEGER IATCONDL(*),IATCONDR(*),NCONDPAIR,IEGFOUT
      INTEGER ATOM(NACLSD,*),CLS(*),EQINV(*),EZOA(NACLSD,*),KAOEZ(*),
     +        NACLS(*),NSHELL(0:NSHELD)
      INTEGER NATOMIMP,IHOST,ATOMIMP(NATOMIMPD),
     +        ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD),NLAYER       
      DOUBLE PRECISION RBASIS(3,*),RR(3,*),RCLS(3,NACLSD,*),
     +                 BRAVAIS(3,3),RECBV(3,3)
      DOUBLE PRECISION RIJ(IJD,3),WTYR(IJD,*),YR(IJD,*)
      LOGICAL LIRR,INTERF,TEST,OPT

      INTEGER I,IN,L,M,LM,J,NOFKS,IM,JN,LM1,LM2,INFO,I1,IL1,IL2,NDIM
      INTEGER LF(LMAXSQ),NSH1(NSHELD),NSH2(NSHELD),NXYZ(3),
     +        ISYMINDEX(NSYMAXD),IPVT(LMMAXSO),IPVT1(NATOMIMP*LMMAXSO)
      DOUBLE PRECISION VOLBZ,RFCTOR,BZKP(3,KPOIBZ),RSYMAT(64,3,3),
     +                 RATOM(3,NSHELD),RATOMS(3,NSHELD),R1,
     +                 RROT(48,3,NSHELD),VOLCUB(KPOIBZ),RATOMI(3)
      DOUBLE COMPLEX DSYMLL(LMAXSQ,LMAXSQ,NSYMD),
     +               DSYMLL1(LMMAXSO,LMMAXSO,NSYMD),TAUVBZ
      LOGICAL SYMUNIATRY(NSYMAXD)
      CHARACTER*10 ROTNAME(64)
     
      INTEGER INVMOD,ISYM,NS,IC,K,IOFF1,IOFF2,JOFF1,JOFF2,IGF,ILM,JLM 
      DOUBLE COMPLEX GINP(LMAXSQ*NACLSD,LMAXSQ,NCLSD),
     +               TMATLL(LMMAXSO,LMMAXSO,NAEZD),TLL(LMMAXSO,LMMAXSO),
     +               TINVLL(LMMAXSO,LMMAXSO,NAEZD)
      DOUBLE COMPLEX ETAIKR(NSYMAXD,NSHELD),CARG
      DOUBLE COMPLEX, ALLOCATABLE :: GLLKE(:,:),GLLKE0(:,:),
     +               GLLKE1(:,:),GLLKE2(:,:),GS(:,:,:,:),GMATLL(:,:,:)
      
      DOUBLE PRECISION RCLSIMP(3,NATOMIMPD)
      
      DOUBLE COMPLEX,ALLOCATABLE :: GLL(:,:,:,:),GIMP(:,:)
      
!       integer :: thrid, nthrd
      integer :: irec, mu, nscoef, ix, jx
      integer, allocatable :: iatomimp(:)

      WRITE(6,*) 'in writegreen'

      PI=4d0*DATAN(1d0)
      RFCTOR = ALAT/(8.D0*DATAN(1.0D0))

       CALL POINTGRP(RSYMAT,ROTNAME)

       IF (TEST('SETGROUP').OR.TEST('fullBZ  ')) THEN 
        CALL SETGROUP(BRAVAIS,RECBV,RBASIS,RFCTOR,NAEZ,
     &                RSYMAT,ROTNAME,ISYMINDEX,NSYMAT)
       ELSE       
        CALL FINDGROUP(BRAVAIS,RECBV,RBASIS,RFCTOR,NAEZ,
     &                 RSYMAT,ROTNAME,ISYMINDEX,NSYMAT)          
        DO I=1,NSYMAXD
         SYMUNIATRY(I) = .FALSE.
        ENDDO
       END IF
        
       DSYMLL=CZERO
       DSYMLL1=CZERO

       CALL ROTBRILL(DSYMLL,LMAXD,NSYMAT,RSYMAT,
     &               ISYMINDEX,LPOT,YR,WTYR,RIJ,IJEND)
       DO ISYM=1,NSYMAT
        DO LM1=1,LMAXSQ
         DO LM2=1,LMAXSQ
          DSYMLL1(LM1,LM2,ISYM)=DSYMLL(LM1,LM2,ISYM)
          DSYMLL1(LM1+LMAXSQ,LM2+LMAXSQ,ISYM)=DSYMLL(LM1,LM2,ISYM)
         ENDDO
        ENDDO
       ENDDO
c
c Now DSYMLL hold NSYMAT symmetrization matrices
c
       CALL GFSHELLS1(ICC,NSH1,NSH2,NSHELL,NAEZ,NATYP,
     &               RBASIS,BRAVAIS,RATOM,RATOMS,
     &               NSYMAT,ISYMINDEX,RSYMAT,RCLS,CLS,
     &               EQINV,KAOEZ,NACLS,ATOM,RFCTOR,
     &               IATCONDL,IATCONDR,NCONDPAIR,ROTNAME,
     &               NATOMIMP,ATOMIMP,RCLSIMP)

c
c --->  creates difference vectors RROT for BZ integration in KKRMAT1
c
      DO I=1,NSHELL(0)
       DO IN=1,3
        RATOMS(IN,I) =RATOM(IN,I) 
       END DO
      END DO
c
        CALL CRTSTAR(RATOMS,NSHELL(0),RSYMAT,NSYMAT,ISYMINDEX,RROT)
c ------------------------------------------------------------------------
        LM = 0
        DO 100 L = 1,LMAXD + 1
          DO 90 M = 1,L + L - 1
            LM = LM + 1
            LF(LM) = L
   90     CONTINUE
  100   CONTINUE

        IF (NBY.EQ.0) NBY = NBX/BBYA
        IF (NBZ.EQ.0) NBZ = NBX/CBYA

c        DO 230 L=1,MAXMESH
c          IF (L.GT.1) THEN
c            NBX = (NBX)/1.4
c            NBY = (NBY)/1.4
c            NBZ = (NBZ)/1.4
c            IF (NBX.EQ.0) NBX =1
c            IF (NBY.EQ.0) NBY =1
c            IF (NBZ.EQ.0) NBZ =1
c          END IF

c          LIRR=.TRUE.
          NXYZ(1) = NBX
          NXYZ(2) = NBY
          NXYZ(3) = NBZ

          IF (OPT('NO-BREAK')) THEN
          CALL BZIRR3D(NOFKS,NXYZ,KPOIBZ,BZKP,RECBV,BRAVAIS,RFCTOR,
     &                 VOLCUB,VOLBZ,RSYMAT,NSYMAT,ISYMINDEX,LIRR)
c          WRITE(6,*) 'VOLBZ',VOLBZ 
          ENDIF
c 230     END DO                      !  L=1,MAXMESH
    

c       TAUVBZ = 1d0/VOLBZ

c invert tmat
       TINVLL=CZERO

       DO I1=1,NAEZ
        DO LM1=1,LMMAXSO
         TINVLL(LM1,LM1,I1)=CONE
        ENDDO
      
        DO LM1=1,LMMAXSO
         DO LM2=1,LMMAXSO
          TLL(LM1,LM2)=TMATLL(LM1,LM2,I1)
         ENDDO
        ENDDO
        
        CALL ZGETRF(LMMAXSO,LMMAXSO,TLL,LMMAXSO,IPVT,INFO)
        CALL ZGETRS('N',LMMAXSO,LMMAXSO,TLL,LMMAXSO,IPVT,
     +              TINVLL(1,1,I1),LMMAXSO,INFO)
        DO LM1=1,LMMAXSO
         DO LM2=1,LMMAXSO
          TINVLL(LM1,LM2,I1)=TINVLL(LM1,LM2,I1)*RFCTOR
         ENDDO
        ENDDO
       ENDDO ! I1
c calculate green function for the host
       ALLOCATE(GS(LMMAXSO,LMMAXSO,NSYMAT,NSHELL(0)))
       GS=CZERO

       IF (INTERF) INVMOD=1
       
       NLAYER=NAEZ/NPRINCD
       CALL GFMASK(ICHECK,ICC,INVMOD,NSH1,NSH2,NLAYER,NAEZ,IE,
     &       NSHELL,IATCONDL,IATCONDR,NCONDPAIR,IEGFOUT)

       IF (OPT('full inv')) INVMOD=0
       INVMOD=0 ! for FS calculation, we inverse a full matrix anyway

       IF (OPT('NO-BREAK')) THEN
       
        open(8888,file='mu0',form='formatted')
        read(8888,*) mu, nscoef
        allocate(iatomimp(nscoef))
        do i1=1,nscoef
          read(8888,*) iatomimp(i1)
        end do
        close(8888)
        
        if(test('rhoqtest')) then 
           open(9999, access='direct', file='tau0_k', 
     &     form='unformatted', recl=(LMMAXSO*LMMAXSO)*4) ! lm blocks
     
!            write(*,*) 'reading kpts from file'
!            open(8888, file='kpt_in.txt')
!            read(8888,*) nofks
!            write(*,*) nofks
!            do k=1,nofks
!              read(8888,*) (bzkp(ns,k), ns=1,3), volcub(k)
!            end do
!            close(8888)
        end if
!         write(*,*) 'open tau_k unform',LMMAXSO, LMMAXSO*LMMAXSO

! !$omp parallel default(shared) private(ns,i,j,isym,carg,lm1,lm2)
! !$omp& private(etaikr,ic,fac,m,im,jn,gllke,gllke0,gllke1,gllke2)
! !$omp& private(ioff1,ioff2,joff1,joff2,il1,il2,i1,icheck,gs)    
! !$omp& private(irec,ix,jx, nthrd, thrid)
        ALLOCATE(GLLKE(ALMSO,ALMSO))
        ALLOCATE(GLLKE0(LMMAXSO,LMMAXSO))
        ALLOCATE(GLLKE1(ALM,LMAXSQ))
        ALLOCATE(GLLKE2(ALM,ALM))

!         nthrd = omp_get_num_threads()
!         thrid = omp_get_thread_num()

        !print header of statusbar
        write(6,'("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")')
     &      0, 20, 40, 60, 80, 100
        write(6,FMT=190) !beginning of statusbar
c start kpoints loop
!         !$omp do
        DO K=1,NOFKS

c create e^(-ikr)
         DO NS=1,NSHELL(0)
          I = NSH1(NS)
          J = NSH2(NS)
          DO ISYM=1,NSYMAT
           CARG=(BZKP(1,K)*(RROT(ISYM,1,NS)-RBASIS(1,J)+RBASIS(1,I))+
     +           BZKP(2,K)*(RROT(ISYM,2,NS)-RBASIS(2,J)+RBASIS(2,I))+
     +           BZKP(3,K)*(RROT(ISYM,3,NS)-RBASIS(3,J)+RBASIS(3,I)))*
     +           2D0*PI*CI

           ETAIKR(ISYM,NS) = VOLCUB(K)*EXP(-CARG)
          ENDDO ! ISYM
         ENDDO ! NS

         GLLKE=CZERO
         GLLKE0=CZERO
         GLLKE1=CZERO
         GLLKE2=CZERO

         DO I=1,NAEZ
          IC=CLS(KAOEZ(I))
          FAC=CONE
          IF (I.NE.EQINV(I)) FAC = CONEM      
           CALL CINIT(NAEZD*LMAXSQ*LMAXSQ,GLLKE1)
           CALL DLKE1(GLLKE1,ALAT,NACLS(IC),RR,EZOA(1,I),ATOM(1,I),
     +               BZKP(1,K),IE,IC,FAC,GINP(1,1,IC),RCLS(1,1,IC))   
           DO M=1,LMAXSQ
            IM=(I-1)*LMAXSQ+M
            DO JN=1,LMAXSQ*NAEZ
             GLLKE2(JN,IM) = GLLKE2(JN,IM)+ GLLKE1(JN,M)
            ENDDO
           ENDDO
          ENDDO ! atom loop
c  for spin-orbit coupling G_LL'(k) double size
          IF (NSPD == 2) THEN
           DO I=1,NAEZ
            DO J=1,NAEZ
             IOFF1=LMMAXSO*(I-1)
             JOFF1=LMAXSQ*(I-1)

             IOFF2=LMMAXSO*(J-1)
             JOFF2=LMAXSQ*(J-1)
             DO LM1=1,LMAXSQ
              DO LM2=1,LMAXSQ
               GLLKE(IOFF1+LM1,IOFF2+LM2)=GLLKE2(JOFF1+LM1,JOFF2+LM2)
               GLLKE(IOFF1+LM1+LMAXSQ,IOFF2+LM2+LMAXSQ)=
     +                     GLLKE2(JOFF1+LM1,JOFF2+LM2)
              
              END DO
             END DO
            ENDDO
           ENDDO
          END IF !NSPD
         
          DO I1=1,NAEZ
           DO LM1=1,LMMAXSO
            DO LM2=1,LMMAXSO
             IL1=LMMAXSO*(I1-1)+LM1
             IL2=LMMAXSO*(I1-1)+LM2
             GLLKE(IL1,IL2)=GLLKE(IL1,IL2)-TINVLL(LM1,LM2,I1)
            ENDDO
           ENDDO
          ENDDO
        
          CALL INVERSION_SO(GLLKE,INVMOD,ICHECK,LMMAXSO)
          DO NS=1,NSHELL(0)
           I=NSH1(NS)
           J=NSH2(NS)
           ILM=LMMAXSO*(I-1)+1
           JLM=LMMAXSO*(J-1)
           DO LM1=1,LMMAXSO
            CALL ZCOPY(LMMAXSO,GLLKE(ILM,JLM+LM1),1,GLLKE0(1,LM1),1)
           ENDDO
           DO ISYM=1,NSYMAT
            DO LM1=1,LMMAXSO
             DO LM2=1,LMMAXSO
              GS(LM1,LM2,ISYM,NS)=GS(LM1,LM2,ISYM,NS)+
     +                           ETAIKR(ISYM,NS)*GLLKE0(LM1,LM2)
             ENDDO
            ENDDO
           ENDDO
           if(test('rhoqtest')) then 
!              write(*,*) 'writing gllke',shape(gllke0),ns,nshell(0),K
!              if(k==1) write(*,*) '',i,j, mu, iatomimp, 
!      &        (I==mu) , (J==mu) , any(I==iatomimp) ,any(J==iatomimp), 
!      &              ((I==mu) .or. (J==mu) .or. 
!      &            any(I==iatomimp) .or. any(J==iatomimp))
             irec = (nscoef*2)*(K-1) + (nscoef*2)*NOFKS*(IE-1)
             if( ((I==mu) .and. any(J==iatomimp)) .or. ((J==mu) .and. 
     &            any(I==iatomimp)) ) then
               ix=0
               lm1 = 1
               do while (ix==0)
                 if(iatomimp(lm1)==i) ix = lm1
                 lm1 = lm1 + 1
               end do
               jx=0
               lm1 = 1
               do while(jx==0)
                 if(iatomimp(lm1)==j) jx = lm1
                 lm1 = lm1 + 1
               end do
               if(j==mu) then
                 irec = irec + nscoef + ix
               else
                 irec = irec + jx
               endif
!                if (k==1 .and. ns==1) write(*,*) ix,jx,nscoef,irec,k,
!                write(*,'(10I9)') i,j,mu,ix,jx,irec,k, ie,
!      &                                          nofks,nscoef
               write(9999,rec=irec) GLLKE0(1:LMMAXSO,1:LMMAXSO)
             end if ! i==mu ...
           end if ! test('rhoqtest')
          ENDDO ! NS

      !update statusbar
         if(mod(K,NOFKS/50)==0) write(6,FMT=200)
!          if(thrid==0.and.mod(K,NOFKS/50/nthrd)==0) write(6,FMT=200)
         ENDDO ! K loop
!          !$omp end do
         
190      FORMAT('                 |'$)   ! status bar
200      FORMAT('|'$)                    ! status bar
         write(6,*)                      ! status bar
!          if(thrid==0) write(6,*)                      ! status bar
         DEALLOCATE(GLLKE)
         DEALLOCATE(GLLKE0)
         DEALLOCATE(GLLKE1)
         DEALLOCATE(GLLKE2)
         WRITE(6,*) 'After k-points integration'
!          if(thrid==0) WRITE(6,*) 'After k-points integration'
!          !$omp end parallel
         
         if(test('rhoqtest')) then
           ! close tau_k file
           close(9999)
           ! save kpoints
           open(8888, file='kpts.txt', form='formatted')
           write(8888,'(I9)') NOFKS
           CALL GETVOLBZ(RECBV,BRAVAIS,VOLBZ)
           write(8888,'(E16.7)') VOLBZ
           do k=1,nofks
             write(8888,'(4E16.7)') (BZKP(i1,K), i1=1,3), VOLCUB(K)
           end do
           close(8888)
           ! save shell info
           open(8888, file='shellinfo.txt')
           write(8888,'(2I9)') NSHELL(0), mu, nscoef
           write(8888,'(1000I9)') (iatomimp(i1), i1=1,nscoef)
           write(8888,'(20000I9)') NSH1(1:NSHELL(0)), NSH2(1:NSHELL(0))
           close(8888)
         end if

       ELSE ! OPT('NO-BREAK')
        IF (OPT('BREAK-1 ')) THEN
         ! write to files
          IF(IE==1)THEN
            open(unit=180,file='greentrans.dat',form='formatted',
     +           action='write')

            write(180,'(100I8)') NSHELD, NATYPD, NAEZD, NEMBD, NCLSD,
     +                           NACLSD, NRD, LMMAXSO, LMAXSQ, NAEZ,
     +                           NSPD
            write(180,'(100I8)') NSHELL(0), NSYMAT, NXYZ(1:3)
            write(180,'(5000I8)') NSH1(1:NSHELD), NSH2(1:NSHELD),
     +              CLS(1:NATYPD),
     +              KAOEZ(1:(NAEZD+NEMBD)), EQINV(1:NAEZD),
     +              NACLS(1:NCLSD), ATOM(1:NACLSD,1:NAEZD),
     +              EZOA(1:NACLSD,1:NAEZD)
            write(180,'(50000ES25.16)') ALAT, RROT(1:48,1:3,1:NSHELD),
     +              RBASIS(1:3,1:(NAEZD+NEMBD)), RR(1:3,1:NRD+1),
     +              RCLS(1:3,1:NACLSD,1:NCLSD), RECBV(1:3,1:3)

          END IF!IE==1

          write(180,'(10000ES25.16)')
     +                 TINVLL(1:LMMAXSO,1:LMMAXSO,1:NAEZD),
     +                 GINP(1:LMAXSQ*NACLSD,1:LMAXSQ,1:NCLSD)

          IF(IE==3) close(180)

        ELSEIF(OPT('BREAK-2 '))THEN
          !read in GS
          write(*,*) 'read in GS, IE=',IE
          if(IE==1)
     +      open(unit=182,file='greentrans_back.dat',form='formatted',
     +           action='read')
          read(182,'(2ES25.16)') GS
          if(IE==3) close(182)
        ELSE!OPT('BREAK-1 ')
          !rise error
          stop 'either BREAK-1 or BREAK-2 or NO-BREAK must be given '
        END IF!OPT('BREAK-1 ')

       END IF

       IF (OPT('NO-BREAK').OR.OPT('BREAK-2 ')) THEN           
        ALLOCATE(GMATLL(LMMAXSO,LMMAXSO,NSHELL(0)))
        GMATLL=CZERO
        CALL GETVOLBZ(RECBV,BRAVAIS,VOLBZ)
        WRITE(6,*) 'VOLBZ',VOLBZ
        TAUVBZ=1d0/VOLBZ
        CALL TAUTOG1_SO(GS,TINVLL,DSYMLL1,NSHELL(0),RFCTOR,
     +                   GMATLL,IGF,TAUVBZ,NSYMAT,NSH1,NSH2,RATOM)
        DEALLOCATE(GS)
c expand GLL for all pair
        ALLOCATE(GLL(LMMAXSO,LMMAXSO,NATOMIMP,NATOMIMP))
        ALLOCATE(GLLKE0(LMMAXSO,LMMAXSO))
        GLL=CZERO
        GLLKE0=CZERO
        DO 10 I=1,NATOMIMP
         DO 10 J=1,NATOMIMP
          DO I1=1,NSHELL(0)
           DO ISYM=1,NSYMAT
            ISYMAT=ISYMINDEX(ISYM)
            DO K=1,3
             RATOMI(K)=RSYMAT(ISYMAT,K,1)*RATOM(1,I1)+
     +                 RSYMAT(ISYMAT,K,2)*RATOM(2,I1)+
     +                 RSYMAT(ISYMAT,K,3)*RATOM(3,I1)
            ENDDO
            IF (ATOMIMP(I).EQ.NSH1(I1).AND.ATOMIMP(J).EQ.NSH2(I1).OR.
     +          ATOMIMP(I).EQ.NSH2(I1).AND.ATOMIMP(J).EQ.NSH1(I1)) THEN
             R1=(RCLSIMP(1,J)-RCLSIMP(1,I)-RATOMI(1))**2+
     +          (RCLSIMP(2,J)-RCLSIMP(2,I)-RATOMI(2))**2+
     +          (RCLSIMP(3,J)-RCLSIMP(3,I)-RATOMI(3))**2
             IF (R1.LT.1d-10) THEN
              CALL ZGEMM('T','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,
     +                   DSYMLL1(1,1,ISYM),LMMAXSO,GMATLL(1,1,I1),
     +                   LMMAXSO,CZERO,GLLKE0,LMMAXSO)   
              CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,
     +                   GLLKE0,LMMAXSO,DSYMLL1(1,1,ISYM),
     +                   LMMAXSO,CZERO,GLL(1,1,I,J),LMMAXSO)
              GO TO 10
             ENDIF 
            ENDIF
           ENDDO ! ISYM
          ENDDO ! I1
  10     CONTINUE
         DEALLOCATE(GMATLL)
         DEALLOCATE(GLLKE0)
         ALLOCATE(GIMP(NATOMIMP*LMMAXSO,NATOMIMP*LMMAXSO))
         GIMP=CZERO
         WRITE(58,'(2e17.9)') E
         DO J=1,NATOMIMP
          DO LM2=1,LMMAXSO
           JLM=(J-1)*LMMAXSO+LM2
           DO I=1,NATOMIMP
            DO LM1=1,LMMAXSO
             ILM=(I-1)*LMMAXSO+LM1
             GIMP(ILM,JLM)=GLL(LM1,LM2,I,J)
             WRITE(58,'((2I5),(2e17.9))') JLM,ILM,GIMP(ILM,JLM)
            ENDDO
           ENDDO
          ENDDO
         ENDDO 
        DEALLOCATE(GLL)
        DEALLOCATE(GIMP)
        WRITE(6,*) 'After host Green function'
       ENDIF ! OPT('NO-BREAK').OR.OPT('BREAK-2 ')
       
       END
