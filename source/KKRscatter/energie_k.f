c 28.9.99 ***************************************************************
      SUBROUTINE KKR_ENERGY_K(VOLCUB,E,LF,
     +                  TINVLL,TMATLL,DELTALL,RROT,DSYMLL,
     +                  NSHELL,RFCTOR,IE,IELAST,ALAT,NSYMAT,
     +                  NAEZ,CLS,EQINV,NACLS,RR,EZOA,ATOM,KAOEZ,
     +                  NSH1,NSH2,GINP,ICC,RBASIS,RCLS,FACTINV,
     &                  TAUVBZ,RATOM,IHANDLE1,IGF,
     &                  IEGFOUT,REFPOT)
      implicit none
c ************************************************************************
c It calculates the band structure along symmetry directions.
c (optionally calculates also complex band structure)
c use OPT('COMPLEX')
c ------------------------------------------------------------------------
c     .. parameters ..
      include 'inc.p'
      include 'inc.cls'
      INTEGER LMAX
      PARAMETER (LMAX=LMAXD)
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAX+1)**2)
      INTEGER ALM
      PARAMETER (ALM = NAEZD*LMAXSQ)
      INTEGER NKPOID
      PARAMETER(NKPOID=1700000)
      INTEGER NAUX
      PARAMETER(NAUX = 2*ALM**2+5*ALM)
      INTEGER NSPBLOCK,NDIM
      PARAMETER (NSPBLOCK=1,NDIM = NSPBLOCK*LMAXSQ)   ! can be put to nprincd, change inc.p also
      DOUBLE PRECISION ZERO
      DOUBLE COMPLEX CI,CZERO,CONE
      PARAMETER (ZERO  = 0.0D0)
      PARAMETER (CI=(0.D0,1.D0),CZERO=(0.D0,0.D0),CONE=(1.D0,0.D0))
      
c     ..
c     .. scalar arguments ..
      DOUBLE COMPLEX E,TAUVBZ
      DOUBLE PRECISION ALAT,RFCTOR
      INTEGER ICC,IE,IELAST,IUMAX,IVMAX,NAEZ,NOFKS,NSHELL,NZ,NSYMAT
      INTEGER NLBASIS,NRBASIS,IGF,IEGFOUT
c     ..
c     .. array arguments ..
      DOUBLE COMPLEX
     +     DELTALL(LMAXSQ,LMAXSQ,*),
     +     DSYMLL(LMAXSQ,LMAXSQ,*),
     +     GINP(LMAXSQ*NACLSD,LMAXSQ,*),
     +     TINVLL(LMAXSQ,LMAXSQ,*),
     +     TMATLL(LMAXSQ,LMAXSQ,*)
      DOUBLE PRECISION
     +     RROT(48,3,*),
     +     VOLCUB(*),
     +     RBASIS(3,*),             ! position of atoms in the unit cell
                                    ! in units of bravais vectors
     +     RR(3,0:NRD),
     +     RCLS(3,NACLSD,*),RATOM(3,NSHELD),
     +     DAUX(2*ALM),WUP,WDO
      INTEGER
     +     ATOM(NACLSD,*),
     +     CLS(*),
     +     EQINV(*),
     +     EZOA(NACLSD,*),
     +     KAOEZ(*),
     +     NACLS(*),
     +     LF(*),
     +     NSH1(*),NSH2(*),
     &     ICOUPLE(NAEZD,NAEZD)
C     ..
C     .. LOCAL SCALARS ..
      DOUBLE COMPLEX CARG,GTEMP,g1,g2,FACI
      INTEGER DI,DJ,DIM,
     +        I,I1,I2,I3,I4,ID,IETA,II,IICC,
     +        ILM,ILM1,IMAX,INCNP,INFO,IO,IOPT,IROUND,ISYM,IU,
     +        J,J1,JJ,JLM,JLM2,IL1,IL2,IP1,IP2,
     +        K,II1,II2,LDI1,LDI2,
     +        LL,LM,LM1,LM2,LMJ,
     +        NI,NL,NLAYER,NROUND,NS,NSTAU,NUMGSP,NUMTAU,
     +        IDECI,INVMOD,IER,IL,IHOST,IHANDLE1,INV,
     +        IRF,JRF,RF,INS,NK
      INTEGER L1,l2,ip1t,ip2t,m1,m2,ldi1t,ldi2t,ISTEP1,ISTEP2,ILT1,ILT2
     +        NKPOINTS,NUM,NKPOINTS,KSUM,POS,KSUM_STORE
      INTEGER status_read,status_open
      DOUBLE PRECISION QDECI,TPI,VOL,BOUND,FACTINV,test1
C     ..
C     .. LOCAL ARRAYS ..
      DOUBLE COMPLEX
     +     AUX(NAUX),
     +     G(LMAXSQ,LMAXSQ),GTRANS(LMAXSQ,LMAXSQ),
     +     TAU(NDIM,NDIM),
     +     TESTLL(LMAXSQ,LMAXSQ),
     +     LVDUMMY(1,ALM),RVDUMMY(1,ALM),
     +     W(ALM),TLL(LMAXSQ,LMAXSQ)
      double complex GSP(NDIM,NDIM,NAUXSPD)
      DOUBLE PRECISION BZKPK(6),WEIGHT
      DOUBLE PRECISION BZKP(6,NKPOID),DK(6)
      INTEGER IPVT(ALM),IPVT1(LMAXSQ),IPIV(NDIM),
     +     INGSP(NLSPD,NLSPD),
     +     INGSP0(NLSPD,NLSPD),
     +     INTAU(NLSPD,NLSPD),
     +     ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD),ITERMAX,ICHCK,
     +     IND(ALM),
     +     REFPOT(*)
      LOGICAL LGSP,LIO,TEST,OPT,TLAYPOS,LINTERFACE
      CHARACTER*80 UIO
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL CINIT,DLKEB,DLKE0,OPT,RCSTOP,
     +         TEST,TLAYPOS,
     +         ZGETRF,ZGETRS,ZXMYPZ
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC DATAN,EXP
      DATA BOUND/1.0D-10/
C
c     .. arrays in common
C
      DOUBLE COMPLEX GLLKE(ALM,ALM)
c      
      DATA LGSP,LIO /2*.TRUE./
c     ..
c ------------------------------------------------------------------------
      if(test('flow     '))
     +     write(6,*) '>>> kkrmateigen: loop over k-points'
c
      IF (OPT('sparse   ')) THEN
       WRITE(6,*) ' Do not use option  -> sparse <-'
       STOP 'Error stop from sub:   kkrmateigen'
      END IF 
       TPI = 8.D0*DATAN(1.D0)    ! = 2*PI 

c
c     reading from the inputcard the points defining the symmetry
c     directions for a band structure calculation.

c      CALL BANDINPUT(IE,ksum,BZKP,NKPOID) 


c     open(unit=178, file="kpoints_upper",status="old",
c    +              action="read",iostat=status_open)

c     if ( status_open /= 0 ) then
c       stop "problem with file kpoints_upper"   
c     end if

c     ksum=0      
c     einlese_schleife: do
c       ksum=ksum+1
c       read (178,"(3(E17.9,X))",iostat=status_read)
c    +            BZKP(1,ksum),BZKP(2,ksum),BZKP(3,ksum)
c       do j=4,6
c         BZKP(j,ksum)=0.d0
c       end do
c       if ( status_read /= 0 ) exit   
c     end do einlese_schleife

c     If(status_read < 0) then
c       ksum=ksum-1
c     END IF

c     close(178)


      open(unit=178, file="kpoints_lower",status="old",
     +              action="read",iostat=status_open)

      if ( status_open /= 0 ) then
        stop "problem with file kpoints_lower"   
      end if

      ksum=0      
      einlese_schleife: do
        ksum=ksum+1
        read (178,"(3(E17.9,X))",iostat=status_read)
     +            BZKP(1,ksum),BZKP(2,ksum),BZKP(3,ksum)
        do j=4,6
          BZKP(j,ksum)=0.d0
        end do
        if ( status_read /= 0 ) exit   
      end do einlese_schleife

      If(status_read < 0) then
        ksum=ksum-1
      END IF

      close(178)



c     open(unit=178, file="kpoints_lower",status="old",
c    +              action="read",iostat=status_open)

c     if ( status_open /= 0 ) then
c       stop "problem with file kpoints_lower"   
c     end if

c     do nk=ksum+1,2*ksum
c        write (6,*) "ksum:",ksum
c       read (178,"(3(E17.9,X))",iostat=status_read)
c    +            BZKP(1,nk),BZKP(2,nk),BZKP(3,nk)
c       do j=4,6
c         BZKP(j,nk)=0.d0
c       end do
c       if ( status_read /= 0 ) exit   
c     end do 

c     ksum=2*ksum
      write (6,*) "Number of kpoints on the FS:",ksum
c     close(178)

      CALL IoInput('INS       ' ,UIO,1,7,IER)   
      READ (UNIT=UIO,FMT=*) INS

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     If we want just real eigenvalues we have to transform the
c     KKR matrix in an hermitian complex matrix.
c     See the PHD thesis of P. Zahn
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     this is the P = {delta_LL'*exp(i*delta_LL')} transformation.
      
c         write (6,*) 't-matrices BEFORE the transformation'
c         do lm1=1,lmaxsq
c            write (6,*) lm1,tinvll(lm1,lm1,1)
c         enddo
c         
         ! write (6,*) 'phase shift'
         ! do lm1=1,lmaxsq
         !    write (6,*) lm1,deltall(lm1,lm1,1)
         ! enddo
      
      
      DO 40 I1 = 1,NAEZ   
         
         INV = KAOEZ(I1)
         RF  = REFPOT(INV)
         
         IF (INS.EQ.0) THEN
            DO 470 LM1 = 1,LMAXSQ
               TMATLL(LM1,LM1,I1) = TMATLL(LM1,LM1,I1)
     &              /DELTALL(LM1,LM1,RF)/DELTALL(LM1,LM1,RF)
 470        END DO

            CALL CINIT(LMAXSQ*LMAXSQ,TINVLL(1,1,I1))

            DO 350 LM1 = 1,LMAXSQ
              TINVLL(LM1,LM1,I1) = RFCTOR/TMATLL(LM1,LM1,I1)   
 350        CONTINUE

         ELSE
            DO 471 LM1 = 1,LMAXSQ
               DO 472 LM2 = 1,LMAXSQ
                  TMATLL(LM1,LM2,I1) = TMATLL(LM1,LM2,I1)
     +                 /DELTALL(LM1,LM1,RF)/DELTALL(LM2,LM2,RF)
 472           END DO
 471        END DO

            CALL CINIT(LMAXSQ*LMAXSQ,TINVLL(1,1,I1))
            
            DO 360 LM1 = 1,LMAXSQ
               TINVLL(LM1,LM1,I1) = RFCTOR
 360        CONTINUE
c     
            CALL ZCOPY(LMAXSQ*LMAXSQ,TMATLL(1,1,I1),1,TLL,1)
c     
c     ---> invert t-matrix and count negative eigenvalues (real part)
c     
            CALL ZGETRF(LMAXSQ,LMAXSQ,TLL,LMAXSQ,IPVT1,INFO)
            CALL ZGETRS('N',LMAXSQ,LMAXSQ,TLL,LMAXSQ,
     +           IPVT1,TINVLL(1,1,I1),LMAXSQ,INFO)
         END IF
                  
 40   CONTINUE 
         
      
      ! write (6,*) 't-matrices AFTER the transformation'
      ! do lm1=1,lmaxsq
      ! write (6,*) lm1,tinvll(lm1,lm1,1)
      ! enddo    
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Starts the Dyson part with the Fourier transform.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     
      write (39,*) KSUM,IE ,'  K-POINTS AND ENERGY '
      
      DO 300 K = 1,KSUM         ! K-POINT-LOOP
         ! write(6,*)  ' k-point ',K        
         DO J=1,6
            BZKPK(J) = BZKP(J,K)
         ENDDO


         if (ie.eq.1) then
            if (k.eq.1) then 
            write(6,*) 
     &'       Real       and         imaginary parts of k-point '
            end if
            write (6,9090)k,(BZKPK(J),j=1,3), (BZKPK(J),j=4,6)
         endif
 9090    FORMAT(I5,3F10.5,3X,3F10.5)
c     
c ---> fourier transformation
c

c         do i1=1,20
c           do i2=1,16
c            write(6,*) i1,i2,GINP(i1,i2,1)
c           end do
c        end do
c        write (6,*) ' 000000GINP GINP 0000000 '


        CALL DLKE0(GLLKE,GSP,INGSP,NSPBLOCK,ALAT,NAEZ,
     +       CLS,EQINV,NACLS,RR,EZOA,ATOM,BZKPK,IE,KAOEZ,RCLS,GINP)
c        write(6,*) 'After dlke0 '
c        do i1=1,16
c           do i2=1,16
c            write(6,*) i1,i2,GLLKE(i1,i2)
c           end do
c        end do
c        write (6,*) ' 0000000000000 '
c
c     here the Green function is transformed via the I transformation
c     (see PHD thesis of P.Zahn), in order to have a final hermitian
c     complex KKR matrix.


        
        FACI = CI
        DO 90 I = 1,NAEZ
           IRF = REFPOT(KAOEZ(I))
           DO 80 J = 1,NAEZ  
              JRF = REFPOT(KAOEZ(J))
              DO 70 LM1 = 1,LMAXSQ
                 ILM1=LMAXSQ*(I-1) + LM1
                 DO 60 LM2 = 1,LMAXSQ
                    JLM2=LMAXSQ*(J-1) + LM2
c     
                    GLLKE(ILM1,JLM2) =  
     +                   GLLKE(ILM1,JLM2)
     +                   *(FACI**(LF(LM1)-LF(LM2)))
     +                   *DELTALL(LM1,LM1,IRF)*DELTALL(LM2,LM2,JRF)
c     
 60              END DO
 70           END DO
 80        END DO
 90     END DO
            
c     It constructs the matrix M=[-(t)^-1 + G^r] and it is stored
c     in the same matrix GLLKE where G^r was stored.

        DO I1=1,NAEZ
           DO LM1=1,LMAXSQ
              DO LM2=1,LMAXSQ
                 IL1=LMAXSQ*(I1-1)+LM1
                 IL2=LMAXSQ*(I1-1)+LM2
                 GLLKE(IL1,IL2)=GLLKE(IL1,IL2)-TINVLL(LM1,LM2,I1)
              ENDDO
           ENDDO
           
        ENDDO
        
c     eigenvalues determination
c        do i1=1,20
c           do i2=1,20
c            write(6,*) i1,i2,GLLKE(i1,i2)
c           end do
c        end do
c        write (6,*) 'before zgeev '
        CALL ZGEEV('N','N',ALM,GLLKE,ALM,W,LVDUMMY,1,RVDUMMY,1,
     +              AUX,NAUX,DAUX,INFO)
c        write (6,*) 'after zgeev '
c        DO I=1,ALM
c           write (6,2000) W(I)
c        ENDDO

        CALL ZSORT(W,IND,ALM,POS)
        NUM = 0
        
        DO  I = 1,ALM
           II = IND(I)
           IF (DREAL(W(II)).lt.zero) NUM = NUM + 1
        END DO
        info = 0
        if (info.gt.1) then
         write(39,2000) num,(bzkpk(j),j=1,6)
        DO I=1,ALM
           write (39,2001) W(IND(I))
        ENDDO
        end if
        write (39,2002) ie,num,(bzkpk(j),j=1,6)
 2000   FORMAT(I5,6D15.6)
 2001   format(2D15.6,'   eigenval')
 2002   FORMAT(2I5,6D15.6)
c     test purpose!!!!!!!!!!!!

c        IF (IE.eq.1) THEN
c           DO I=1,ALM
c              write (6,*) IE,K,W(IND(I))
c           ENDDO
c        ENDIF
        

 300  END DO                    ! K = 1,KSUM

      RETURN

 9000 FORMAT(20I4)
 9010 FORMAT(3f12.5)
 9011 FORMAT('K = ',i4,3f12.5)
 9030 FORMAT(A1,4x,':',I6,2d12.3,2f10.4)
 9040 FORMAT(' NUMTAU =',I6)
 9050 FORMAT(' NUM_LR =',I6)

      END















