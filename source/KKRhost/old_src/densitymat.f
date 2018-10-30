C*==densitymat.f    processed by SPAG 6.05Rc at 15:54 on 21 Dec 2004
      SUBROUTINE DENSITYMAT(DF,PZ,QZ,PNS,QNS,AR,CR,
     &                      DR,GMATLL,IPAN,IRCUT,DRDI,EK,
     &                      IRMIN,LOPT,MMAX,LMSTART,LMEND,PHI,DENMATC
     &     ,den,ie) ! test fivos 19.9.08
C **********************************************************************
C *                                                                    *
C * Calculation of density matrix needed in evaluating the Coulomb     *
C * interaction potential in LDA+U                                     *
C * non-relativistic case -- otherwise matrices DENMAT and VLDAU must  *
C *                          have double dimension                     *
C *                                                                    *
C * The density matrix n is calculated using the Green function and    *
C * the reference functions Phi as                                     *
C *                                                                    *
C *  n_{m,s,m',s'} = -1/pi Int dE                                      *
C *    [ Sum_{L''L'''} (Phi_L,R_L'') G_{L''L'''}(E) (R_L''',Phi_L') +  *
C *                                                                    *
C *      + Sum_L'' (Phi_L,R_L'')(H_L'',Phi_L') ]                       *
C *                                                                    *
C * This is expressed in complex spherical harmonics basis. It is then *
C * transformed into the real spherical harmonics basis - subroutine   *
C * WMATLDAU                                                           *
C *                                                                    *
C * Here, only the l-th subblock of G is used. (Is this correct? qldau)*
C *                                                                    *
C * The integration is implicit: this routine is called within an      *
C * energy loop and the result is summed up.                           *
C *                                                                    *
C *                  ph. mavropoulos, h.ebert munich/juelich 2002-2004 *
C **********************************************************************
      IMPLICIT NONE
      INCLUDE 'inc.p'
C
C PARAMETER definitions
C
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (KREL+1) * (LMAXD+1)**2)
      INTEGER MMAXD
      PARAMETER ( MMAXD = 2*LMAXD + 1 )
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD) 
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C
      DOUBLE COMPLEX CZERO,CONE
      PARAMETER (CZERO=(0.0D0,0.0D0),CONE=(1.D0,0.D0))
C
C Dummy arguments
C
      DOUBLE COMPLEX DF,EK
      INTEGER IRMIN,LOPT
      INTEGER LMSTART,LMEND,MMAX
      DOUBLE COMPLEX AR(LMMAXD,LMMAXD),CR(LMMAXD,LMMAXD),
     &           DENMATC(MMAXD,MMAXD),DR(LMMAXD,LMMAXD),
     &           GMATLL(LMMAXD,LMMAXD),PHI(IRMD),
     &           PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),PZ(IRMD,0:LMAXD),
     &           QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),QZ(IRMD,0:LMAXD)
      DOUBLE PRECISION DRDI(IRMD)
      INTEGER IPAN,IRCUT(0:IPAND)
C
C Local variables
C
      DOUBLE COMPLEX DENMATC2(MMAXD,MMAXD),GTEMP(MMAXD,MMAXD),
     &               PHIQ(MMAXD,MMAXD),PHIR(MMAXD,MMAXD)
      INTEGER LM1,LM2,M1,M2
      INTEGER LMAXD1,IE ! test fivos 19.9.08
      PARAMETER (LMAXD1= LMAXD+1) ! test fivos 19.9.08
      DOUBLE COMPLEX DEN(0:LMAXD1,IEMXD*(1+KREL)) ! test fivos 19.9.08
      EXTERNAL CINIT,OVERLAP,RINIT
C
      CALL CINIT(MMAXD*MMAXD,DENMATC2(1,1))
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C 1.  Within implicit energy loop:
C Calculate density matrix.
C
C 1a. Calculate overlap integral (inner product) with wavefunctions of
C current energy (Phi,R) (result: PHIR) and (Phi,Q) (result:PHIQ)
C Result is in real Ylm basis.
C
C
      CALL CINIT(MMAXD*MMAXD,PHIR)
      CALL OVERLAP(PHIR,PHI,PZ,QZ,PNS,AR,DR,.FALSE.,IPAN,IRCUT,DRDI,
     &             IRMIN,LOPT,IPAND,LMAXD,LMMAXD,MMAXD,
     &             LMPOTD,IRMIND,IRMD)
C
      CALL CINIT(MMAXD*MMAXD,PHIQ)
      CALL OVERLAP(PHIQ,PHI,PZ,QZ,QNS,CR,DR,.TRUE.,IPAN,IRCUT,DRDI,
     &             IRMIN,LOPT,IPAND,LMAXD,LMMAXD,MMAXD,
     &             LMPOTD,IRMIND,IRMD)
C
C
C 1b. Use overlap integrals and Green function matrix to find density 
C matrix (implicit integration)
C Result is in real Ylm basis.
C
C Copy l-th subblock of G into Gtemp (Is this correct? qldau)
      DO LM2 = LMSTART,LMEND
         M2 = LM2 - LMSTART + 1
         DO LM1 = LMSTART,LMEND
            M1 = LM1 - LMSTART + 1
            GTEMP(M1,M2) = GMATLL(LM1,LM2)
         END DO
      END DO
C
C First step: PHIQ = G*PHIR+EK*PHIQ.
C (If phi=pz, the trace of this should be similar to the dos).
C
      CALL ZGEMM('N','N',MMAX,MMAX,MMAX,CONE,GTEMP,MMAXD,PHIR,MMAXD,EK,
     &           PHIQ,MMAXD)
C Second step: DENMATC2 = PHIR*PHIQ
      CALL ZGEMM('T','N',MMAX,MMAX,MMAX,CONE,PHIR,MMAXD,PHIQ,MMAXD,
     &           CZERO,DENMATC2,MMAXD)
C
C Third step: Integration step: DENMATC = DF*DENMATC2 + DENMATC
C
      CALL CINIT(MMAXD*MMAXD,DENMATC2(1,1)) ! test fivos 19.9.08
      do m1 = 1,mmax                         ! test fivos 19.9.08
      denmatc2(m1,m1) = den(lopt,ie)/dfloat(mmax) ! test fivos 19.9.08
      enddo                                 ! test fivos 19.9.08
      DO M2 = 1,MMAX
         DO M1 = 1,MMAX
            DENMATC(M1,M2) = DENMATC(M1,M2) + DF*DENMATC2(M1,M2)
         END DO
      END DO
c test fivos
c         write(*,9001) denmatc(1,1),denmatc(2,2),denmatc(3,3),
c    &        denmatc(4,4),denmatc(5,5),denmatc(6,6),denmatc(7,7)
c        write(*,9001) -denmatc2(1,1)/3.14159,
c    &                 -denmatc2(2,2)/3.14159,
c    &                 -denmatc2(3,3)/3.14159,
c    &                 -denmatc2(4,4)/3.14159,
c    &                 -denmatc2(5,5)/3.14159,
c    &                 -denmatc2(6,6)/3.14159,
c    &                 -denmatc2(7,7)/3.14159
c9001 format(14e12.4,' test fivos')
C
C Energy loop ends
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      END
