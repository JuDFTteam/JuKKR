SUBROUTINE calcmvec(nfilcbwf,splitss,iepath,nepath,  &
    irel,iprint,nt,nl,mezz,mezj,taut,tsst, iqat,nkmq,nkm,iecurr,netab,igrid,we,  &
    mvevdl0,mvevil,bmvevdl0,bmvevil, r2drdi,jrws,imt,amemvec,  &
    ikmllim1,ikmllim2,imkmtab,ntmax,nlmax,nmuemax,  &
    nqmax,nkmmax,nmmax,nmvecmax,nrmax)
!   ********************************************************************
!   *                                                                  *
!   *                                                                  *
!   ********************************************************************
use mod_types, only: t_inc
IMPLICIT COMPLEX*16(A-H,O-Z)

! PARAMETER definitions

COMPLEX*16 C0,C1,CI
PARAMETER (C0=(0.0D0,0.0D0),C1=(1.0D0,0.0D0),CI=(0.0D0,1.0D0))
REAL*8 PI
PARAMETER (PI=3.141592653589793238462643D0)
COMPLEX*16 CPRE
PARAMETER (CPRE=-C1/PI)

! Dummy arguments

COMPLEX*16 WE
INTEGER IECURR,IEPATH,IGRID,IPRINT,IREL,NEPATH,NETAB,NFILCBWF,NKM, &
        NKMMAX,NL,NLMAX,NMMAX,NMUEMAX,NMVECMAX,NQMAX,NRMAX,NT, &
        NTMAX
LOGICAL SPLITSS
REAL*8 AMEMVEC(NKMMAX,NKMMAX,3,NMVECMAX),R2DRDI(NRMAX,NMMAX)
COMPLEX*16 BMVEVDL0(NLMAX,NTMAX,3,NMVECMAX), &
           BMVEVIL(NLMAX,NTMAX,3,NMVECMAX), &
           MEZJ(NKMMAX,NKMMAX,NTMAX,NMVECMAX), &
           MEZZ(NKMMAX,NKMMAX,NTMAX,NMVECMAX), &
           MVEVDL0(NLMAX,NTMAX,3,NMVECMAX), &
           MVEVIL(NLMAX,NTMAX,3,NMVECMAX), &
           TAUT(NKMMAX,NKMMAX,NTMAX),TSST(NKMMAX,NKMMAX,NTMAX)
INTEGER IKMLLIM1(NKMMAX),IKMLLIM2(NKMMAX),IMKMTAB(NKMMAX), &
        IMT(NTMAX),IQAT(NQMAX,NTMAX),JRWS(NMMAX),NKMQ(NQMAX)

! Local variables
!
!F77--------------------------------------------------------------------
!ccc      COMPLEX*16 JF(NRMAX,2,NKMMAX),JG(NRMAX,2,NKMMAX),
!ccc     &           ZF(NRMAX,2,NKMMAX),ZG(NRMAX,2,NKMMAX),
!F77--------------------------------------------------------------------
!F90--------------------------------------------------------------------
      COMPLEX*16 JF(:,:,:),JG(:,:,:),ZF(:,:,:),ZG(:,:,:)
      ALLOCATABLE JF,JG,ZF,ZG
!F90--------------------------------------------------------------------
COMPLEX*16 BMVEVD(NTMAX,3,NMVECMAX), &
           BMVEVDL(NLMAX,3,NMVECMAX), &
           BMVEVDM(NLMAX,NMUEMAX,3,NMVECMAX), &
           BMVEVI(NTMAX,3,NMVECMAX),CWGT, &
           MEIRR(NKMMAX,NKMMAX,3,NMVECMAX), &
           MEREG(NKMMAX,NKMMAX,3,NMVECMAX), &
           MVEVD(NTMAX,3,NMVECMAX),MVEVDL(NLMAX,3,NMVECMAX), &
           MVEVDM(NLMAX,NMUEMAX,3,NMVECMAX), &
           MVEVI(NTMAX,3,NMVECMAX),W1(NKMMAX,NKMMAX), &
           W2(NKMMAX,NKMMAX),W3(NKMMAX,NKMMAX), &
           ZFJF(2,2),ZFZF(2,2),ZGJG(2,2),ZGZG(2,2)
LOGICAL CHECK
INTEGER I,I0,IKM,IKM1,IKM2,IKMCB(2,NKMMAX),IKMT1,IKMT2, &
        IL,IM,IMKM1,IMKM2,IMV,IMVEC,IPOL,IT,ITI,J,J1,J2,JTOP,K, &
        K1,K2,KAPCB,L,LI,LMAX,M,MUE,MUETOP,N,NMVEC,NOSWF,NPOL, &
        NSOL,NSOLCB(NKMMAX)
REAL*8 MJ,SUM
CHARACTER*3 STR3

! index 3:  ipol= 1,2,3  ==  (+),(-),(z)


check = .true.
check = .false.
npol = 3
nmvec = MIN(4,nmvecmax)

IF ( iecurr == 1 .AND. iepath == 1 ) THEN
  
  DO imv = 1,nmvec
    DO ipol = 1,npol
      DO it = 1,nt
        DO il = 1,nl
          mvevdl0(il,it,ipol,imv) = c0
          mvevil(il,it,ipol,imv) = c0
          bmvevdl0(il,it,ipol,imv) = c0
          bmvevil(il,it,ipol,imv) = c0
        END DO
      END DO
    END DO
  END DO
  
END IF

noswf = nt*nkm
nmvec = 3

!F90--------------------------------------------------------------------
allocate (jf(nrmax,2,nkmmax),jg(nrmax,2,nkmmax),stat=it)
IF ( it /= 0 ) STOP '      < CALCMVEC > : allocate JF/JG '
allocate (zf(nrmax,2,nkmmax),zg(nrmax,2,nkmmax),stat=it)
IF ( it /= 0 ) STOP '      < CALCMVEC > : allocate ZF/ZG '
!F90--------------------------------------------------------------------
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
DO it = 1,nt
  
  m = nkmmax
  n = nkmq(iqat(1,it))
  im = imt(it)
  im = 1
  jtop = jrws(im)
  
  lmax = nl - 1
  
  CALL cinit(nkmmax*nkmmax*3*nmvecmax,mereg)
  CALL cinit(nkmmax*nkmmax*3*nmvecmax,meirr)
  
!=======================================================================
!                      calculate matrix elements
!=======================================================================
  
  DO ikm = 1,n
    READ (nfilcbwf,REC=ikm+(it-1)*nkm) iti,li,mj,nsol,str3,  &
        (kapcb,ikmcb(k,ikm),(zg(i,k,ikm),zf(i,k,ikm),i=1,jtop) ,k=1,nsol)
    IF ( it /= iti ) STOP ' IT(INI) <> IT  in <MENABIRR>'
    IF ( str3 /= 'REG' ) STOP 'WFT(INI) <> REG in <MENABIRR>'
    
    READ (nfilcbwf,REC=ikm+(it-1)*nkm+noswf) iti,li,mj,nsol, str3,  &
        (kapcb,ikmcb(k,ikm),(jg(i,k,ikm),jf(i,k,ikm),i=1,jtop) ,k=1,nsol)
    nsolcb(ikm) = nsol
    IF ( it /= iti ) STOP ' IT(INI) <> IT  in <MENABIRR>'
    IF ( str3 /= 'IRR' ) STOP 'WFT(INI) <> IRR in <MENABIRR>'
    IF ( iprint > 3 .AND. (t_inc%i_write>0))  &
        WRITE (1337,*) iti,li,mj,nsol,str3,kapcb
  END DO
  
  
  DO ikmt2 = 1,n
    
    DO ikmt1 = ikmllim1(ikmt2),ikmllim2(ikmt2)
      
      DO imvec = 1,nmvec
        DO k2 = 1,nsolcb(ikmt2)
          j2 = ikmcb(k2,ikmt2)
          DO k1 = 1,nsolcb(ikmt1)
            j1 = ikmcb(k1,ikmt1)
            DO ipol = 1,npol
              IF ( ABS(amemvec(j1,j2,ipol,imvec)) > 1E-8 ) GO TO 10
            END DO
          END DO
        END DO
      END DO
! -------------------------------------- all angular matrix elements = 0
      GO TO 20
! ---------------------------------- non-0 angular matrix elements found
! ------------------------------------- calculate radial matrix elements
      
      10            CONTINUE
      CALL cintabr(zg(1,1,ikmt1),zg(1,1,ikmt2),zgzg,  &
          zf(1,1,ikmt1),zf(1,1,ikmt2),zfzf,  &
          r2drdi(1,im),nsolcb(ikmt1),nsolcb(ikmt2), jtop,nrmax)
      
      CALL cintabr(zg(1,1,ikmt1),jg(1,1,ikmt2),zgjg,  &
          zf(1,1,ikmt1),jf(1,1,ikmt2),zfjf,  &
          r2drdi(1,im),nsolcb(ikmt1),nsolcb(ikmt2), jtop,nrmax)
      
! -------------------------------------- calculate total matrix elements
      
      DO k2 = 1,nsolcb(ikmt2)
        ikm2 = ikmcb(k2,ikmt2)
        imkm2 = imkmtab(ikm2)
        
        DO k1 = 1,nsolcb(ikmt1)
          ikm1 = ikmcb(k1,ikmt1)
          imkm1 = imkmtab(ikm1)
          
          DO imv = 1,nmvec
            DO ipol = 1,npol
              mereg(ikmt1,ikmt2,ipol,imv) = mereg(ikmt1,ikmt2,ipol,imv)  &
                  + amemvec(ikm1,ikm2,ipol,imv)*zgzg(k1,k2)
            END DO
          END DO
          
        END DO
        
        DO imv = 1,nmvec
          DO ipol = 1,npol
            mereg(ikmt2,ikmt2,ipol,imv) = mereg(ikmt2,ikmt2,ipol,imv)  &
                - amemvec(imkm2,imkm2,ipol,imv)*zfzf(k2,k2)
          END DO
        END DO
        
      END DO
      
      IF ( ikmt1 == ikmt2 ) THEN
        DO k2 = 1,nsolcb(ikmt2)
          ikm2 = ikmcb(k2,ikmt2)
          imkm2 = imkmtab(ikm2)
          DO k1 = 1,nsolcb(ikmt1)
            ikm1 = ikmcb(k1,ikmt1)
            imkm1 = imkmtab(ikm1)
            
            DO imv = 1,nmvec
              DO ipol = 1,npol
                meirr(ikmt1,ikmt2,ipol,imv) = meirr(ikmt1,ikmt2,ipol,imv)  &
                    + amemvec(ikm1,ikm2,ipol,imv) *zgjg(k1,k2)
              END DO
            END DO
            
          END DO
          
          DO imv = 1,nmvec
            DO ipol = 1,npol
              meirr(ikmt2,ikmt2,ipol,imv) = meirr(ikmt2,ikmt2,ipol,imv)  &
                  - amemvec(imkm2,imkm2,ipol,imv) *zfjf(k2,k2)
            END DO
          END DO
          
        END DO
      END IF
      
    20         END DO
    
  END DO
  
  IF ( check ) THEN
    DO i = 1,nkm
      DO j = ikmllim1(i),ikmllim2(i)
        sum = ABS(mezz(i,j,it,1)) + ABS(mezj(i,j,it,1))
        sum = sum + ABS(mereg(i,j,3,1)) + ABS(meirr(i,j,3,1))
        IF ( sum > 1D-8 .AND. (t_inc%i_write>0)) THEN
          WRITE (1337,*) ' spin '
          WRITE (1337,'(2i3,2e17.8,2x,2e17.8)') i,j,  &
              mezz(i,j,it,2),mezj(i,j,it,2)
          WRITE (1337,'(6x,2e17.8,2x,2e17.8)') mereg(i,j,3,1)  &
              ,meirr(i,j,3,1), (mezz(i,j,it,2)-mereg(i,j,3,1)),  &
              (mezj(i,j,it,2)-meirr(i,j,3,1))
          WRITE (1337,*) ' orb '
          WRITE (1337,'(2i3,2e17.8,2x,2e17.8)') i,j,  &
              mezz(i,j,it,3),mezj(i,j,it,3)
          WRITE (1337,'(6x,2e17.8,2x,2e17.8)') mereg(i,j,3,2)  &
              ,meirr(i,j,3,2), (mezz(i,j,it,3)-mereg(i,j,3,2)),  &
              (mezj(i,j,it,3)-meirr(i,j,3,2))
        END IF
      END DO
    END DO
  END IF
  
!=======================================================================
  DO imv = 1,nmvec
    DO ipol = 1,npol
      
      mvevd(it,ipol,imv) = 0.0D0
      mvevi(it,ipol,imv) = 0.0D0
      
      IF ( .NOT.splitss ) THEN
        bmvevd(it,ipol,imv) = 0.0D0
        bmvevi(it,ipol,imv) = 0.0D0
      END IF
      
      CALL zgemm('N','N',n,n,n,cpre,mereg(1,1,ipol,imv),m,  &
          taut(1,1,it),m,c0,w1,m)
      cwgt = -1D0
      DO j = 1,n
        CALL zaxpy(n,-cpre,meirr(1,j,ipol,imv),1,w1(1,j),1)
        CALL zcopy(n,taut(1,j,it),1,w2(1,j),1)
        CALL zaxpy(n,cwgt,tsst(1,j,it),1,w2(1,j),1)
      END DO
      CALL zgemm('N','N',n,n,n,cpre,mereg(1,1,ipol,imv),m,w2,m, c0,w3,m)
      
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      DO l = 0,lmax
        il = l + 1
        
        IF ( irel <= 1 ) THEN
          i0 = l*l
          DO mue = 1,2*l + 1
            mvevdm(il,mue,ipol,imv) = w1(i0+mue,i0+mue)
            bmvevdm(il,mue,ipol,imv) = w3(i0+mue,i0+mue)
          END DO
        ELSE
          i0 = 2*l*l + l + l
          muetop = 2*l + 2
          DO mue = 1,muetop
            mvevdm(il,mue,ipol,imv) = w1(i0+mue,i0+mue)
            bmvevdm(il,mue,ipol,imv) = w3(i0+mue,i0+mue)
          END DO
          i0 = 2*(l-1)*l + l + l
          DO mue = 2,muetop - 1
            mvevdm(il,mue,ipol,imv) = mvevdm(il,mue,ipol,imv)  &
                + w1(i0+mue-1,i0+mue-1)
            bmvevdm(il,mue,ipol,imv) = bmvevdm(il,mue,ipol,imv)  &
                + w3(i0+mue-1,i0+mue-1)
          END DO
        END IF
        
        mvevdl(il,ipol,imv) = 0.0D0
        
        IF ( .NOT.splitss ) bmvevdl(il,ipol,imv) = 0.0D0
        
        IF ( irel > 1 ) THEN
          muetop = 2*l + 2
        ELSE
          muetop = 2*l + 1
        END IF
        
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
        DO mue = 1,muetop
          mvevdl(il,ipol,imv) = mvevdl(il,ipol,imv) + mvevdm(il,mue,ipol,imv)
          
          IF ( .NOT.splitss ) bmvevdl(il,ipol,imv) = bmvevdl(il,ipol,imv)  &
              + bmvevdm(il,mue,ipol,imv)
        END DO
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
        
        mvevil(il,it,ipol,imv) = mvevil(il,it,ipol,imv)  &
            + we*mvevdl(il,ipol,imv)
        mvevdl0(il,it,ipol,imv) = mvevdl(il,ipol,imv)
        
        IF ( .NOT.splitss ) THEN
          bmvevil(il,it,ipol,imv) = bmvevil(il,it,ipol,imv)  &
              + we*bmvevdl(il,ipol,imv)
          bmvevdl0(il,it,ipol,imv) = bmvevdl(il,ipol,imv)
        END IF
        
        IF ( igrid /= 4 .OR. iecurr <= netab ) THEN
          
          mvevd(it,ipol,imv) = mvevd(it,ipol,imv) + mvevdl(il,ipol,imv)
          mvevi(it,ipol,imv) = mvevi(it,ipol,imv) + mvevil(il,it,ipol,imv)
          
          IF ( .NOT.splitss ) THEN
            bmvevd(it,ipol,imv) = bmvevd(it,ipol,imv) + bmvevdl(il,ipol,imv)
            bmvevi(it,ipol,imv) = bmvevi(it,ipol,imv)  &
                + bmvevil(il,it,ipol,imv)
          END IF
        END IF
        
      END DO
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    END DO
  END DO
  
END DO
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

!F90--------------------------------------------------------------------
deallocate (jf,jg,zf,zg,stat=it)
IF ( it /= 0 ) STOP '      < CALCMVEC > : deallocate JF/JG/ZF/ZG'
!F90--------------------------------------------------------------------
IF ( splitss .AND. ((iepath == 1) .AND. (iecurr == netab)) ) THEN
  DO imv = 1,nmvec
    DO ipol = 1,npol
      DO it = 1,nt
        DO il = 1,nl
          bmvevdl0(il,it,ipol,imv) = mvevdl0(il,it,ipol,imv)
        END DO
      END DO
    END DO
  END DO
END IF

!=======================================================================
IF ( (igrid >= 6) .OR. (iecurr /= netab) .OR. (iepath /= nepath) ) RETURN

!     this part of the original Munich subroutine has been moved to
!     < mvecglobal >
!     main2 --> tbkkr2 --> mvecglobal -- see makefile2

END SUBROUTINE calcmvec
