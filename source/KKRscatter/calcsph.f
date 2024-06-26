       SUBROUTINE CALCSPH(NSRA,IRMDNEW,LMAX,NSPIN,Z,C,E,LMPOTD,LMMAXSO,
     +                   RNEW,VINS,NCHEB,NPAN_TOT,RPAN_INTERVALL,
     +                   JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR,TMAT,
     +                   USE_SRATRICK)
       IMPLICIT NONE
c construct wavefunctions for spherical potentials
       INTEGER NSRA,IRMDNEW,NSPIN,LMAX,LMPOTD,NCHEB,NPAN_TOT,LMMAXSO
       INTEGER USE_SRATRICK
       DOUBLE PRECISION Z,C
       DOUBLE COMPLEX E,GMATPREFACTOR
       DOUBLE PRECISION RNEW(IRMDNEW),RPAN_INTERVALL(0:NPAN_TOT)
       DOUBLE PRECISION VINS(IRMDNEW,LMPOTD,NSPIN)
       DOUBLE COMPLEX HLK(4*(LMAX+1),IRMDNEW),
     +                JLK(4*(LMAX+1),IRMDNEW),
     +                HLK2(4*(LMAX+1),IRMDNEW),
     +                JLK2(4*(LMAX+1),IRMDNEW)
       INTEGER JLK_INDEX(2*LMMAXSO)

c local
       INTEGER LMSIZE,LMSIZE2,NVEC
       INTEGER IVEC,LVAL,IR,ISPIN,LSPIN,LSRA,I,L1,M1,LM1
       INTEGER, ALLOCATABLE :: JLK_INDEXTEMP(:)
       DOUBLE COMPLEX, ALLOCATABLE :: VLL0(:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: VLL(:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: RLLTEMP(:,:,:),SLLTEMP(:,:,:),
     +                                HLKTEMP(:,:),JLKTEMP(:,:),
     +                                HLK2TEMP(:,:),JLK2TEMP(:,:),
     +                                HLKNEW(:,:),JLKNEW(:,:)
       DOUBLE COMPLEX TMATTEMP(1,1)
       DOUBLE COMPLEX TMAT(2*(LMAX+1))
       LMSIZE=1
       IF (NSRA.EQ.2) THEN
        LMSIZE2=2
        NVEC=2
       ELSE
        LMSIZE2=1
        NVEC=1
       ENDIF
       ALLOCATE (RLLTEMP(LMSIZE2,LMSIZE,IRMDNEW))
       ALLOCATE (SLLTEMP(LMSIZE2,LMSIZE,IRMDNEW))
       ALLOCATE (HLKTEMP(NVEC,IRMDNEW))
       ALLOCATE (JLKTEMP(NVEC,IRMDNEW))
       ALLOCATE (HLK2TEMP(NVEC,IRMDNEW))
       ALLOCATE (JLK2TEMP(NVEC,IRMDNEW))
       ALLOCATE (JLK_INDEXTEMP(LMSIZE2))
       ALLOCATE (HLKNEW(NVEC*NSPIN*(LMAX+1),IRMDNEW))
       ALLOCATE (JLKNEW(NVEC*NSPIN*(LMAX+1),IRMDNEW))

       DO IVEC=1,NVEC
        JLK_INDEXTEMP(IVEC)=IVEC
       ENDDO
       ALLOCATE(VLL0(LMSIZE,LMSIZE,IRMDNEW))   
       IF (NSRA.EQ.2) THEN
         ALLOCATE(VLL(2*LMSIZE,2*LMSIZE,IRMDNEW))   
       ELSE
         ALLOCATE(VLL(LMSIZE,LMSIZE,IRMDNEW))   
       ENDIF
c spin loop
       DO ISPIN=1,NSPIN
 
        LSPIN=(LMAX+1)*(ISPIN-1)
        LSRA=(LMAX+1)*NVEC
c each value of l, the Lippmann-Schwinger equation is solved using
c the free-potential wavefunctions and potentials corresponding to l-value
        DO LVAL=0,LMAX

         DO IR=1,IRMDNEW
          VLL0(LMSIZE,LMSIZE,IR)=VINS(IR,1,ISPIN)-2d0*Z/RNEW(IR)
         ENDDO
          
         IF (NSRA.EQ.2) THEN
          CALL VLLMATSRA(VLL0,VLL,RNEW,LMSIZE,IRMDNEW,
     +                  E,C,LMAX,LVAL,'Ref=0')
         ELSE
          VLL(:,:,:)=VLL0(:,:,:)
         ENDIF

         JLKTEMP(1,:)=JLK(LVAL+1,:)
         HLKTEMP(1,:)=HLK(LVAL+1,:)
         JLK2TEMP(1,:)=JLK2(LVAL+1,:)
         HLK2TEMP(1,:)=HLK2(LVAL+1,:)
         IF (NSRA.EQ.2) THEN
          JLKTEMP(2,:)=JLK(LMAX+LVAL+2,:)
          HLKTEMP(2,:)=HLK(LMAX+LVAL+2,:)
          JLK2TEMP(2,:)=JLK2(LMAX+LVAL+2,:)
          HLK2TEMP(2,:)=HLK2(LMAX+LVAL+2,:)
         ENDIF
         CALL RLLSLL(RPAN_INTERVALL,RNEW,VLL,RLLTEMP,SLLTEMP,TMATTEMP,
     +               NCHEB,NPAN_TOT,LMSIZE,LMSIZE2,NVEC,IRMDNEW,NVEC,
     +               JLK_INDEXTEMP,HLKTEMP,JLKTEMP,HLK2TEMP,JLK2TEMP,
     +               GMATPREFACTOR,'1','1','0',USE_SRATRICK)
         
         DO IR=1,IRMDNEW
          HLKNEW(LSPIN+LVAL+1,IR)=SLLTEMP(1,1,IR)/RNEW(IR)
          JLKNEW(LSPIN+LVAL+1,IR)=RLLTEMP(1,1,IR)/RNEW(IR)
         ENDDO
         IF (NSRA.EQ.2) THEN
          DO IR=1,IRMDNEW
           HLKNEW(LSPIN+LSRA+LVAL+1,IR)=SLLTEMP(2,1,IR)/RNEW(IR)
           JLKNEW(LSPIN+LSRA+LVAL+1,IR)=RLLTEMP(2,1,IR)/RNEW(IR)
          ENDDO
         ENDIF
         TMAT(LSPIN+LVAL+1)=TMATTEMP(1,1)          
        ENDDO ! LMAX
       ENDDO ! NSPIN
       
       LM1=1
       DO IVEC=1,NVEC
        DO I=1,2
         DO L1=0,LMAX
          DO M1=-L1,L1
           JLK_INDEX(LM1)=L1+(IVEC-1)*NSPIN*(LMAX+1)+(I-1)*(LMAX+1)+1
           LM1=LM1+1
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       DO IR=1,IRMDNEW
        DO L1=1,NVEC*(LMAX+1)*NSPIN
         HLK(L1,IR)=HLKNEW(L1,IR)
         JLK(L1,IR)=JLKNEW(L1,IR)
        ENDDO
       ENDDO
       IF (NSRA.EQ.2) THEN
        DO IR=1,IRMDNEW
         DO L1=1,(LMAX+1)*NSPIN
          HLK2(L1,IR)=-HLKNEW(L1+LMAX+1,IR)
          JLK2(L1,IR)=-JLKNEW(L1+LMAX+1,IR)
         ENDDO
         DO L1=NSPIN*(LMAX+1)+1,NVEC*(LMAX+1)*NSPIN
          HLK2(L1,IR)=HLKNEW(L1-(LMAX+1)*NSPIN,IR)
          JLK2(L1,IR)=JLKNEW(L1-(LMAX+1)*NSPIN,IR)
         ENDDO
        ENDDO
       ELSE
        DO IR=1,IRMDNEW
         DO L1=1,NVEC*(LMAX+1)*NSPIN
          HLK2(L1,IR)=-HLKNEW(L1,IR)
          JLK2(L1,IR)=-JLKNEW(L1,IR)
         ENDDO
        ENDDO
       ENDIF

       DEALLOCATE (RLLTEMP)
       DEALLOCATE (SLLTEMP)
       DEALLOCATE (HLKTEMP)
       DEALLOCATE (JLKTEMP)
       DEALLOCATE (HLK2TEMP)
       DEALLOCATE (JLK2TEMP)
       DEALLOCATE (JLK_INDEXTEMP)
       DEALLOCATE (HLKNEW)
       DEALLOCATE (JLKNEW)
       DEALLOCATE (VLL0)
       DEALLOCATE (VLL)
       END 
           
