C **********************************************************************
      SUBROUTINE GENERALPOT(IFILE,NATPS,NATYP,NSPIN,Z,ALAT,RMT,RMTNEW,
     +                 RWS,ITITLE,R,DRDI,VM2Z,IRWS,A,B,TXC,KXC,INS,IRNS,
     +                 LPOT,VINS,QBOUND,IRC,KSHAPE,EFERMI,VBC,ECORE,
     +                 LCORE,NCORE,LMPOTD,IRMD,IRMIND)
c **************************************************
c * The subroutine writes out the potential cards
c * in a standard r-mesh that can be read in and 
c * interpolated to a different r-mesh from subroutine
c * start No shape function information is needed
c * and all nessecery data are stored in the potential 
c * card.
c *                                      ver. 18.5.2000
c ***************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER LMPOTD,IRMD,IRMIND
      DOUBLE PRECISION ALAT,QBOUND
      INTEGER IFILE,INS,KSHAPE,KXC,LPOT,NATPS,NATYP,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),B(*),DRDI(IRMD,*),ECORE(20,*),EFERMI,
     +                 R(IRMD,*),RMT(*),RMTNEW(*),RWS(*),VBC(2),
     +                 VINS(IRMIND:IRMD,LMPOTD,*),VM2Z(IRMD,*),Z(*)
      INTEGER IRC(*),IRNS(*),IRWS(*),ITITLE(20,*),LCORE(20,*),NCORE(*)
      CHARACTER*24 TXC(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,B1,RMAX,RMT1,RMTNW1,RV,SUM,Z1,PARSUM,
     &                 PARSUMDERIV,R0,RINTER,DR,MAXA,VCON
      INTEGER I,ICORE,IH,IP,IR,IRMIN,IRNS1,IS,ISAVE,J,LM,LMNR,
     &        LMPOT,NCORE1,NR,NZ1,NR_U,IRMIN_U,IRNS_U,
     &        IMT1,LM1,IRNSTOT
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DRADI(IRMD),ECORE1(20),RA(IRMD),VM2ZA(IRMD),
     &                 RR_U(IRMD),DRDI_U(IRMD)
      DOUBLE PRECISION VM2ZB(IRMD),VM2Z_U(IRMD),
     &            VINS_U(IRMIND:IRMD,LMPOTD),VINSA(IRMIND:IRMD,LMPOTD),
     &            VINSB(IRMIND:IRMD,LMPOTD)
      INTEGER LCORE1(20)
       CHARACTER*4 ELEMNAME(0:113)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
c        1      2      3      4      5      6      7      8      9    
      DATA ELEMNAME/'VAC',
     &  'H   ','He  ','Li  ','Be  ','B   ','C   ','N   ','O   ','F   ',
     &  'Ne  ',
     &  'Na  ','Mg  ','Al  ','Si  ','P   ','S   ','Cl  ','Ar  ','K   ',
     &  'Ca  ',
     &  'Sc  ','Ti  ','V   ','Cr  ','Mn  ','Fe  ','Co  ','Ni  ','Cu  ',
     &  'Zn  ',
     &  'Ga  ','Ge  ','As  ','Se  ','Br  ','Kr  ','Rb  ','Sr  ','Y   ',
     &  'Zr  ',
     &  'Nb  ','Mo  ','Tc  ','Ru  ','Rh  ','Pd  ','Ag  ','Cd  ','In  ',
     &  'Sn  ',
     &  'Sb  ','Te  ','I   ','Xe  ','Cs  ','Ba  ','La  ','Ce  ','Pr  ',
     &  'Nd  ',
     &  'Pm  ','Sm  ','Eu  ','Gd  ','Tb  ','Dy  ','Ho  ','Er  ','Tm  ',
     &  'Yb  ',
     &  'Lu  ','Hf  ','Ta  ','W   ','Re  ','Os  ','Ir  ','Pt  ','Au  ',
     &  'Hg  ',
     &  'Tl  ','Pb  ','Bi  ','Po  ','At  ','Rn  ','Fr  ','Ra  ','Ac  ',
     &  'Th  ',
     &  'Pa  ','U   ','Np  ','Pu  ','Am  ','Cm  ','Bk  ','Cf  ','Es  ',
     &  'Fm  ',
     &  'Md  ','No  ','Lr  ','Rf  ','Db  ','Sg  ','Bh  ','Hs  ','Mt  ',
     &  'Uun ','Uuu ','Uub ','NoE '/

C     ..
      ISAVE = 1
      LMPOT = (LPOT+1)* (LPOT+1)
   
      DO 60 IH = 1,NATYP
         DO 50 IS = 1,NSPIN
          DO LM =1,LMPOTD
          DO IR =IRMIND,IRMD
          VINSA(IR,LM) = 0.d0
          VINSB(IR,LM) = 0.d0
          ENDDO
          ENDDO
          IP = NSPIN* (IH-1) + IS
          RMT1 = RMT(IH)
          RMTNW1 = RMTNEW(IH)
          Z1 = Z(IH)
          RMAX = RWS(IH)
          IF (KSHAPE.EQ.0) THEN
            NR = IRWS(IH)
            IRNS1 = 0
          ELSE
            NR = IRC(IH)
            IRNS1 = IRNS(IH)
          END IF

          IRMIN = NR - IRNS1
          A1 = A(IH)
          B1 = B(IH)
          NCORE1 = NCORE(IP)
c
          DO 10 J = 1,NR
            RA(J) = R(J,IH)
            DRADI(J) = DRDI(J,IH)
            VM2ZA(J) = VM2Z(J,IP) 
   10     CONTINUE
          DO LM1=1,LMPOT
            DO J=IRMIND,IRMD
            VINSA(J,LM1) = VINS(J,LM1,IP)
            END DO
          END DO 
c
          IF (NCORE1.GE.1) THEN
c
            DO 20 J = 1,NCORE1
              LCORE1(J) = LCORE(J,IP)
              ECORE1(J) = ECORE(J,IP)
   20       CONTINUE
          END IF
c
c  Generate uniform mesh RUNI
c     
          NR_U = NR
          IRNS_U = IRNS1
          IRMIN_U = NR_U
          if (ins.gt.0) IRMIN_U = NR_U - IRNS_U
  
          IF (INS.EQ.0) THEN
             DO I=1,NR_U
                RR_U(I) = RA(I)
                DRDI_U(I) = DRADI(I)
             END DO 
             IMT1 = 0
          ELSE
             IMT1 = ANINT(LOG(RMTNW1/B1+1.0D0)/A1) + 1
             DO I=1,IMT1
                RR_U(I) = RA(I)
                DRDI_U(I) = DRADI(I)
             END DO
             RINTER =  RMAX - RMTNW1
             DR = RINTER/FLOAT(NR-IMT1) 
             Do I=1,NR-IMT1
                DRDI_U(IMT1+I) = DR
                RR_U(IMT1+I) = RR_U(IMT1) + DR*FLOAT(I)
             END DO
             CALL DOUBLERAUS1(NR,IRMIN,LMPOT,RA,DRADI,VM2ZA,VINSA,
     &            IRMD,IRMIND,LMPOTD)
c     
c     After this sub the arrays are rearanged and nr is not 
c     the same anymore in the case of FP. If ins.eq.0 there is
c     no nead for doubleraus1. IRMIN should remain the same
c  
          END IF
c ----------------------------------------------------------------
c Now the new mesh is generated
c     
c test     
        !  write(6,*) nr_u,imt1,irns_u
        !  do i=1,nr_u
        !    write(6,*) i,ra(i),rr_u(i)  
        !  end do
c test
c
          maxa = 1.d35
          CALL SPLINE(IRMD,RA,VM2ZA,NR,maxa,maxa,VM2ZB)
          IF (INS.GT.0) THEN
             DO LM1=1,LMPOT 
                IRNSTOT = NR - IRMIN ! nr has changed irmin is the same
                !write(6,*) ' Testing ',nr,irmin,irnstot,irmind 
                CALL SPLINE(IRMD-IRMIND,RA(IRMIND),
     &                      VINSA(IRMIND,LM1),IRNSTOT,maxa,maxa,
     &                      VINSB(IRMIND,LM1))
             ENDDO              ! LM1
          END IF
c     
c OK with spline 
c
          DO IR = 1,NR_U
             R0 = RR_U(IR)
             CALL SPLINT(RA,VM2ZA,VM2ZB,NR,R0,PARSUM,PARSUMDERIV)
             VM2Z_U(IR) = PARSUM
          END DO
          IF (INS.GT.0) THEN
             !IRNSTOT = NR_U - IRMIN_U 
             DO LM1=1,LMPOT
                DO IR = IRMIN_U,NR_U
                   R0 = RR_U(IR)
                   CALL SPLINT(RA(IRMIND),VINSA(IRMIND,LM1),
     &                         VINSB(IRMIND,LM1),IRNSTOT,R0,
     &                         PARSUM,PARSUMDERIV)
                   VINS_U(IR,LM1) = PARSUM 
                END DO
             end do
          END IF 
          !write(6,*) ' All interpolation ok now write'             
c     --------------------------------------------------------------
          WRITE (IFILE,FMT=8000)
          NZ1 = Z1
          IF (NSPIN.EQ.1) THEN
          WRITE (IFILE,FMT=8010) ELEMNAME(NZ1),Z1
          ELSEIF (IS.EQ.1) THEN 
          WRITE (IFILE,FMT=8012) ELEMNAME(NZ1),Z1
          ELSEIF (IS.EQ.2) THEN 
          WRITE (IFILE,FMT=8011) ELEMNAME(NZ1),Z1
          END IF
          WRITE (IFILE,FMT=8020)
c          write (ifile,*) ALAT,RMAX,RMTNW1,RMT1
          WRITE (IFILE,FMT=8030) ALAT,RMAX,RMTNW1,RMT1
          WRITE (IFILE,FMT=8040) NR_U,IMT1,IRNS1
          WRITE (IFILE,FMT=8050) A1,B1
          WRITE (IFILE,FMT=8060) EFERMI,VBC(IS)!,VCON
          WRITE (IFILE,FMT=8070) NCORE1,LMPOT
          IF (NCORE1.GE.1) WRITE (IFILE,FMT=9040) (LCORE1(ICORE),
     +         ECORE1(ICORE),ICORE=1,NCORE1)
          
          IF (INS.EQ.0 .OR. (IH.LT.NATPS.AND.INS.LE.2)) THEN
c     
c---  >       store only the spherically averaged potential 
c     (in mt or as - case) 
c     this is done always for the host
c     
             WRITE (IFILE,FMT=9051) (VM2Z_U(IR),IR=1,NR_U)
          ELSE
c     
c---  >     store the full potential , but the non spherical contribution
c     only from irns1 up to irws1 ;
c     remember that the lm = 1 contribution is multiplied
c     by a factor 1/sqrt(4 pi)
c     
             WRITE (IFILE,FMT=9060) NR_U,IRNS1,LMPOT,ISAVE
             WRITE (IFILE,FMT=9070) (VM2Z_U(IR),IR=1,NR_U)
             IF (LPOT.GT.0) THEN
                LMNR = 1
                DO 40 LM = 2,LMPOT
                   SUM = 0.0D0
                   DO 30 IR = IRMIN,NR_U
                      RV = VINS_U(IR,LM)*RR_U(IR)
                      SUM = SUM + RV*RV*DRADI(IR)
 30                CONTINUE
                   
                   IF (SQRT(SUM).GT.QBOUND) THEN
                      LMNR = LMNR + 1
                    WRITE (IFILE,FMT=9060) LM
                    WRITE (IFILE,FMT=9070) (VINS_U(IR,LM),IR=IRMIN,NR_U)
                   END IF
                   
 40             CONTINUE
c     
c---  >         write a one to mark the end
c     
                IF (LMNR.LT.LMPOT) WRITE (IFILE,FMT=9060) ISAVE
             END IF
             
          END IF
          
 50    CONTINUE
 60   CONTINUE
      
      
 8000 format (' GENERAL POTENTIAL MESH             exc:')
 8010 Format ('#  ',A4,'POTENTIAL             Z = ',F8.3)
 8011 format ('#  ',A4,'POTENTIAL SPIN UP     Z=  ',F8.3)
 8012 format ('#  ',A4,'POTENTIAL SPIN DOWN   Z=  ',F8.3)
 8020 format ('#')
 8030 format (4f12.8, '   # alat, rmax, rmaxlog, rmt')
 8040 format (1p,3I6,31X,'  # IRWS, IRMT, IRNS ')
 8050 format (2D15.8,19X,'  # A , B ')
 8060 format (3f12.8,13X,'  # Ef, vbc, vcon ')
 8070 format (1p,2I5,39X,'  # NCORE, LMPOT' )
 9000 FORMAT (7a4,6x,'  exc:',a24,3x,a10)
 9010 FORMAT (3f12.8)
 9020 FORMAT (f10.5,/,f10.5,2f15.10)
 9030 FORMAT (i3,/,2d15.8,/,2i2)
 9040 FORMAT (i5,1p,d20.11)
 9050 FORMAT (1p,2d15.6,1p,d15.8)
 9051 FORMAT (1p,4d20.12)
 9060 FORMAT (10i5)
 9070 FORMAT (1p,4d20.13)
      END
C **********************************************************************
C
c***********************************************************************
      SUBROUTINE SPLINE(NMAX,X,Y,N,YP1,YPN,Y2) 
      IMPLICIT NONE
      INTEGER N,NMAX 
      DOUBLE PRECISION YP1,YPN,X(NMAX),Y(NMAX),Y2(NMAX) 
c Given arrays x(1:n) and  y(1:n) containing a tabulated function, 
c i.e., y i = f(xi), with x1<x2<...<xN , and given values yp1 and ypn 
c for the rst derivative of the interpolating function at points 
c 1 and n, respectively, this routine returns an array y2(1:n) of 
c length n which contains the second derivatives of the interpolating 
c function at the tabulated points xi. 
c If yp1 and/or ypn are equal to 1.e30 or larger, the routine is
c signaled to set the corresponding boundary condition for a natural
c spline, with zero second derivative on that boundary. 
c Parameter: NMAX is the largest anticipated value of n. 
c Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
      INTEGER I,K 
      DOUBLE PRECISION P,QN,SIG,UN,U(NMAX) 

      IF (N.GT.NMAX) STOP 'SPLINE: N > NMAX.'
      IF (YP1.GT.0.99D30) THEN
c The lower boundary condition is set either to be "natural" 
         Y2(1) = 0.D0
         U(1) = 0.D0
      ELSE
c or else to have a specified first derivative. 
         Y2(1) = -0.5D0
         U(1)=(3.D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1) 
      ENDIF 
C
      DO I = 2,N-1  
c This is the decomposition loop of the tridiagonal algorithm. y2 and u 
c are used for temporary storage of the decomposed factors. 
         SIG = (X(I)-X(I-1)) / (X(I+1)-X(I-1))
         P = SIG * Y2(I-1) + 2.D0 
         Y2(I) = (SIG-1.D0)/P
         U(I)=(6.D0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) 
     &        /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1)) / P 
      ENDDO
C
      IF (YPN.GT..99D30) THEN   
c The upper boundary condition is set either to be "natural"
         QN = 0.D0
         UN = 0.D0
      ELSE
c or else to have a specified 1rst derivative. 
         QN = 0.5D0
         UN = (3.D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N) = (UN-QN*U(N-1)) / (QN*Y2(N-1)+1.D0) 
      DO K = N-1,1,-1 
c This is the backsubstitution loop of the tridiagonal algorithm. 
         Y2(K)=Y2(K)*Y2(K+1)+U(K) 
      ENDDO
C
      RETURN 
      END
C **********************************************************************
C
C **********************************************************************
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y,YDERIV)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X,Y,YDERIV,XA(*),YA(*),Y2A(*)
c Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
c function (with the xai's in order), and given the array y2a(1:n), which
c is the output from spline above, and given a value of x, this routine
c returns a cubic-spline interpolated value y and the derivative yderiv.
c Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
      INTEGER K,KHI,KLO
      DOUBLE PRECISION A,B,H
c We will find the right place in the table by means of bisection.
c This is optimal if sequential calls to this routine are at random
c values of x. If sequential calls are in order, and closely
c spaced, one would do better to store previous values of
c klo and khi and test if they remain appropriate on the
c next call.
      KLO=1
      KHI=N
 1    IF (KHI-KLO.GT.1) THEN
         K=(KHI+KLO)/2
         IF(XA(K).GT.X)THEN
            KHI=K
         ELSE
            KLO=K
         ENDIF
         GOTO 1
      ENDIF
c klo and khi now bracket the input value of x.
      H=XA(KHI)-XA(KLO)
c The xa's must be distinct.
      IF (H.EQ.0.D0) PAUSE 'Bad XA input in SPLINT'
c Cubic spline polynomial is now evaluated.
      A = (XA(KHI)-X)/H
      B = (X-XA(KLO))/H
      Y = A*YA(KLO) + B*YA(KHI) +
     &     ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI)) * (H**2)/6.D0
      YDERIV = (YA(KHI)-YA(KLO))/H - 
     &     ((3.D0*A*A-1.D0)*Y2A(KLO) - (3.D0*B*B-1.D0)*Y2A(KHI))*H/6.D0
C
      RETURN
      END
C **********************************************************************
C
C **********************************************************************
c************************************************************************
      SUBROUTINE DOUBLERAUS1(IRMAX,IRMIN,LMPOT,RR,DRDI,VPOT,VINS,
     &                       IRMD,IRMIND,LMPOTD)
c Gets rid of the double-points in the radial mesh, i.e. the points
c where RR(I) = RR(I-1). Returns the "new" mesh in the same array,
c and rearranges accordingly the WAVEF defined at the same mesh.
c IRMAX is also altered to the new value.
      IMPLICIT NONE
      INTEGER IRMD,LMPOTD,IRMIND
      INTEGER NCOUNTMAX
      PARAMETER(NCOUNTMAX=500)
c Input and output:
      INTEGER IRMAX
      DOUBLE PRECISION RR(IRMD),DRDI(IRMD),VPOT(IRMD),
     &                 VINS(IRMIND:IRMD,LMPOTD)
c Inside:
      INTEGER IR,ICOUNT,NC,NCOUNT,LMPOT,IRMIN
      INTEGER LM1
      INTEGER IDOUBLE(NCOUNTMAX)

c Find double-points:
      NCOUNT = 0
      DO IR = 2,IRMAX
         IF ((RR(IR)-RR(IR-1)).LT.1.D-20) THEN
            NCOUNT = NCOUNT + 1
            IDOUBLE(NCOUNT) = IR
         ENDIF
      ENDDO
      IF (NCOUNT+1.GT.NCOUNTMAX) 
     &                STOP 'DOUBLERAUS2: Too many double-points.'
      IDOUBLE(NCOUNT+1) = IRMAX+1   ! To be used below.
      
c Rearrange the arrays.
      DO ICOUNT = 1,NCOUNT
         DO IR = IDOUBLE(ICOUNT)-ICOUNT+1,IDOUBLE(ICOUNT+1)-ICOUNT
            IF((IR+ICOUNT).LE.IRMAX) THEN
               RR(IR) = RR(IR+ICOUNT)
               DRDI(IR) = DRDI(IR+ICOUNT)
               VPOT(IR) = VPOT(IR+ICOUNT)
            END IF
         ENDDO
      ENDDO
      IRMAX = IRMAX - NCOUNT
      NCOUNT = 0
      DO NC = 1,NCOUNTMAX
       IDOUBLE(NC) = 0 
      ENDDO

      DO IR = IRMIN,IRMAX
         IF ((RR(IR)-RR(IR-1)).LT.1.D-20) THEN
            NCOUNT = NCOUNT + 1
            IDOUBLE(NCOUNT) = IR
         ENDIF
      ENDDO
c Rearrange the arrays.
      DO ICOUNT = 1,NCOUNT
         DO IR = IDOUBLE(ICOUNT)-ICOUNT+1,IDOUBLE(ICOUNT+1)-ICOUNT
            DO LM1=1,LMPOT
            VINS(IR,LM1) = VINS(IR+ICOUNT,LM1)
            END DO
         ENDDO
      ENDDO      
      RETURN
      END
C **********************************************************************
