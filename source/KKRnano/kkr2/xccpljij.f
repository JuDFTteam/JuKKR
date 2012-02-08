      SUBROUTINE XCCPLJIJ(
     >                    INFO,I1,IER,WGTE,
     >                    RXIJ,NXIJ,IXCP,RXCCLS,
     >                    GMATXIJ,DTIXIJ,
     >                    LMPIC,LCOMM,
     >                    MYRANK,EMPIC,EMYRANK,
     <                    JXCIJINT,ERESJIJ,
C                         new input parameters after inc.p removal
     &                    naez, lmax, nxijd, nspind,
     &                    lmpid, smpid, empid)
C   ********************************************************************
C   *                                                                  *
C   *  calculates the site-off diagonal  XC-coupling parameters  J_ij  *
C   *  according to  Lichtenstein et al. JMMM 67, 65 (1987)            *
C   *                                                                  *
C   *  adopted for TB-KKR code from Munich SPR-KKR package Sep 2004    *
C   *  adopted for KKRnano, Jun 2009                                   *
C   ********************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'mpif.h'

      INTEGER naez
      INTEGER lmax
      INTEGER nxijd
      INTEGER nspind
      INTEGER lmpid
      INTEGER smpid
      INTEGER empid

C     .. Parameters
      DOUBLE COMPLEX   CONE,CZERO
      PARAMETER        ( CONE  = (1D0,0D0) )
      PARAMETER        ( CZERO = (0D0,0D0) )

C     INTEGER          LMMAXD
C     PARAMETER        (LMMAXD= (LMAXD+1)**2)
C     ..
C     .. Scalar arguments
      DOUBLE COMPLEX   WGTE
      INTEGER          I1,IER,NXIJ,
     +                 LMPIC
      CHARACTER        INFO
      LOGICAL          ERESJIJ
C     ..
C     .. Array arguments
C     DOUBLE COMPLEX   GMATXIJ(LMMAXD,LMMAXD,NXIJD,NSPIND),
C    &                 DTIXIJ(LMMAXD,LMMAXD)
C     DOUBLE PRECISION RXIJ(NXIJD),RXCCLS(3,NXIJD)
C     INTEGER          IXCP(NXIJD),LCOMM(LMPID*SMPID*EMPID)

      DOUBLE COMPLEX   GMATXIJ((LMAX+1)**2,(LMAX+1)**2,NXIJD,NSPIND)
      DOUBLE COMPLEX   DTIXIJ((LMAX+1)**2,(LMAX+1)**2)
      DOUBLE PRECISION RXIJ(NXIJD)
      DOUBLE PRECISION RXCCLS(3,NXIJD)
      INTEGER          IXCP(NXIJD)
      INTEGER          LCOMM(LMPID*SMPID*EMPID)

C     ..
C     .. Local scalars
      INTEGER          XIJ,ISPIN,LM1,LM2,
     &                 D1,D10,D100,D1000,
     &                 SEND,RECV
      DOUBLE COMPLEX   CSUM,JSCAL
      DOUBLE PRECISION JOUT
      CHARACTER*12     FNAME
      LOGICAL          LSAME
C     ..
C     .. Local arrays
      INTEGER          OFF(3)

C     DOUBLE COMPLEX   JXCIJINT(NXIJD),
C    &                 JXCE(NXIJD,EMPID),
C    &                 JRECV(NXIJD),
C    &                 DTNXIJ(LMMAXD,LMMAXD,NAEZD),
C    &                 GMIJ(LMMAXD,LMMAXD),
C    &                 GMJI(LMMAXD,LMMAXD),
C    &                 W1(LMMAXD,LMMAXD),W2(LMMAXD,LMMAXD),
C    &                 W3(LMMAXD,LMMAXD)

      DOUBLE COMPLEX   JXCIJINT(NXIJD)
      DOUBLE COMPLEX   JXCE(NXIJD,EMPID)
      DOUBLE COMPLEX   JRECV(NXIJD)
      DOUBLE COMPLEX   GMIJ((LMAX+1)**2,(LMAX+1)**2)
      DOUBLE COMPLEX   GMJI((LMAX+1)**2,(LMAX+1)**2)
      DOUBLE COMPLEX   W1((LMAX+1)**2,(LMAX+1)**2)
      DOUBLE COMPLEX   W2((LMAX+1)**2,(LMAX+1)**2)
      DOUBLE COMPLEX   W3((LMAX+1)**2,(LMAX+1)**2)

C     large local array
C     DOUBLE COMPLEX   DTNXIJ((LMAX+1)**2,(LMAX+1)**2,NAEZ)
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: DTNXIJ

C     .. MPI ..
      INTEGER, dimension(MPI_STATUS_SIZE) :: STATUS
      INTEGER          IERR,MAPBLOCK
C     .. N-MPI
      INTEGER          MYRANK
C     .. E-MPI
      INTEGER          EMYRANK(EMPID,NAEZ*LMPID*SMPID),
     &                 EMPI,EMPIR,EMPIC
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC        SQRT
C     ..
C     .. External Subroutines ..
      EXTERNAL         CINIT,CMATMUL,ZCOPY,LSAME
C     ..

      INTEGER memory_stat
      LOGICAL memory_fail

      INTEGER          LMMAXD

      LMMAXD= (LMAX+1)**2
      memory_fail = .false.

C=======================================================================
C energy parallelization requires splitting of this routine into
C part 'R' = running over energy points and into part 'F' the 
C output part
C=======================================================================
C
      IF (LSAME(INFO,'R')) THEN
C
        JSCAL = CONE/4D0
C
!        IF (LMPIC.EQ.1) THEN

C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C ==>                   IE.EQ.1 -- initialisation step --            <==
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

C     Allocate array

      allocate(DTNXIJ(LMMAXD,LMMAXD,NAEZ), stat = memory_stat)
      if (memory_stat /= 0) memory_fail = .true.

      if (memory_fail .eqv. .true.) then
        write(*,*) "XCCPLJIJ: FATAL Error, failure to allocate memory."
        write(*,*) "       Probably out of memory."
      stop
      end if

        IF ( IER.EQ.1 ) THEN
          DO XIJ = 2, NXIJ
          JXCIJINT(XIJ) = CZERO
          ENDDO

          IF (ERESJIJ) THEN
            D1 = mod(I1,10)
            D10 = int( (mod(I1,100) + 0.5)/10 )
            D100 = int( (mod(I1,1000) + 0.5)/100 )
            D1000 = int( (mod(I1,10000) + 0.5)/1000 )

            OFF(1) = iachar('1')-1
            OFF(2) = iachar('1')-1
            OFF(3) = iachar('1')-1

            IF ( D10.GE.10 ) OFF(1) = iachar('7')
            IF ( D100.GE.100 ) OFF(2) = iachar('7')
            IF ( D1000.GE.1000 ) OFF(3) = iachar('7')

            FNAME='Eij.'
     +       //achar(D1000+OFF(3))
     +       //achar(D100+OFF(2))
     +       //achar(D10+OFF(1))
     +       //achar(D1+iachar('1')-1)
     +       //'.dat'

            OPEN(75,FILE=FNAME,FORM='formatted')
          ENDIF
        ENDIF

C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C ==>                   INITIALISATION END                           <==
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
C XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
C
C ==>  get DTJXIJ = T(UP) - T(DOWN) for all atoms
C

       CALL MPI_ALLGATHER(DTIXIJ,LMMAXD*LMMAXD,MPI_DOUBLE_COMPLEX,
     +                    DTNXIJ,LMMAXD*LMMAXD,MPI_DOUBLE_COMPLEX,
     +                    LCOMM(LMPIC),IERR)


        DO XIJ = 2, NXIJ ! loop XIJ = 1, NXIJ(I1)
C
C ==>  get the off-diagonal Green function matrix Gij(UP) and Gji(DOWN)
C
         DO ISPIN = 1,2
           IF ( ISPIN.EQ.1 ) THEN
             CALL ZCOPY(LMMAXD*LMMAXD,GMATXIJ(1,1,XIJ,ISPIN),
     +                  1,GMIJ,1)
           ELSE
             DO LM2 = 1,LMMAXD
               DO LM1 = 1,LMMAXD
C
C -> use Gji = Gij^T
C
                  GMJI(LM1,LM2) = GMATXIJ(LM2,LM1,XIJ,ISPIN)

               ENDDO
             ENDDO
           ENDIF
         ENDDO

C ----------------------------------------------------------------------
C
C ==> calculate the exchange coupling constant J_ij via Eq. (19)
C     modified for G instead of tau:
C          J_ij ~ Trace [ (t_i(D)-t_i(U)) * Gij(U)
C                       * (t_j(D)-t_j(U)) * Gji(D)]
C
C -------------------------------------------------- loop over occupants
C --> Delta_j * Gjt,it
C
         CALL CMATMUL(LMMAXD,LMMAXD,DTNXIJ(1,1,IXCP(XIJ)),GMJI,W2)
C
C --> Delta_i * Git,jt
C
         CALL CMATMUL(LMMAXD,LMMAXD,DTIXIJ,GMIJ,W3)
C
C --> Delta_i * Git,jt * Delta_j * Gjt,it
C
         CALL CMATMUL(LMMAXD,LMMAXD,W3,W2,W1)
C
         CSUM = CZERO
         DO LM1 = 1,LMMAXD
           CSUM = CSUM + W1(LM1,LM1)
         ENDDO
C
         JOUT = -DIMAG(WGTE*CSUM*JSCAL)
C
         IF (ERESJIJ) THEN
           WRITE(75,73002) 
     &     IER,XIJ,RXIJ(XIJ),JOUT,
C     &     IER,XIJ,RXIJ(XIJ),-DIMAG(WGTE*CSUM),
     &     RXCCLS(1,XIJ),RXCCLS(2,XIJ),RXCCLS(3,XIJ),IXCP(XIJ)
         ENDIF
C
         JXCIJINT(XIJ) = JXCIJINT(XIJ) - WGTE*CSUM
C
C                  -------> perform substraction instead of addition
C                           because WGTE ~ -1/pi
C ======================================================================
        ENDDO             ! loop XIJ = 1, NXIJ
C XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
!        ENDIF ! LMPIC = 1
C

      deallocate(DTNXIJ)

      if (ERESJIJ) close(75)

      ENDIF ! called with INFO='R'
C=======================================================================
C finished subroutine if called with INFO = 'R'
C=======================================================================
C
C
C
C
C
C
C
C
C
C=======================================================================
C energy parallelization requires splitting of this routine into
C part 'R' = running over energy points and into part 'F' the 
C output part
C=======================================================================
C
      IF (LSAME(INFO,'F')) THEN
C
        IF (EMPID.GT.1) THEN
C
          DO EMPI = 1, EMPID
C
            SEND  =
     +      EMYRANK(MAPBLOCK(EMPI,1,EMPID,1,0,EMPID-1)+1,EMPIC)
C
            DO EMPIR = 1, EMPID
C
              RECV =
     +        EMYRANK(MAPBLOCK(EMPIR,1,EMPID,1,0,EMPID-1)+1,EMPIC)
C
              IF (RECV.EQ.SEND.AND.RECV.EQ.MYRANK) THEN
C
                DO XIJ=1,NXIJ
                  JXCE(XIJ,EMPI) = JXCIJINT(XIJ)
                ENDDO
C
              ELSE
C
                IF (MYRANK.EQ.SEND) THEN
C
                  CALL MPI_SEND(JXCIJINT,NXIJD,
     +                          MPI_DOUBLE_COMPLEX,
     +                          RECV,93,MPI_COMM_WORLD,IERR)
C
                ENDIF
C
                IF (MYRANK.EQ.RECV) THEN
C
                  CALL MPI_RECV(JRECV,NXIJD,
     +                          MPI_DOUBLE_COMPLEX,
     +                          SEND,93,MPI_COMM_WORLD,STATUS,IERR)
C
                  DO XIJ=1,NXIJ
                    JXCE(XIJ,EMPI) = JRECV(XIJ)
                  ENDDO
C
                ENDIF
C
              ENDIF
C
            ENDDO
C
          ENDDO
C
          DO XIJ = 1, NXIJ
            JXCIJINT(XIJ) = CZERO
          ENDDO
C
          DO EMPI = 1, EMPID
            DO XIJ = 1, NXIJ
              JXCIJINT(XIJ) = JXCIJINT(XIJ) + JXCE(XIJ,EMPI)
            ENDDO
          ENDDO
C
        ENDIF
C
C
        IF (LMPIC.EQ.1) THEN
C
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        JSCAL = CONE/4D0
C .. ...........................................................
C write Jij's to file Jij.I1.dat
C ..
        D1 = mod(I1,10)
        D10 = int( (mod(I1,100) + 0.5)/10 )
        D100 = int( (mod(I1,1000) + 0.5)/100 )
        D1000 = int( (mod(I1,10000) + 0.5)/1000 )

        OFF(1) = iachar('1')-1
        OFF(2) = iachar('1')-1
        OFF(3) = iachar('1')-1

        IF ( D10.GE.10 ) OFF(1) = iachar('7')
        IF ( D100.GE.100 ) OFF(2) = iachar('7')
        IF ( D1000.GE.1000 ) OFF(3) = iachar('7')

        FNAME='Jij.'
     +   //achar(D1000+OFF(3))
     +   //achar(D100+OFF(2))
     +   //achar(D10+OFF(1))
     +   //achar(D1+iachar('1')-1)
     +   //'.dat'

        OPEN(73,FILE=FNAME,FORM='formatted')

        WRITE(73,73000) I1

        DO XIJ = 2, NXIJ

          JXCIJINT(XIJ) = JSCAL*JXCIJINT(XIJ)
          WRITE(73,73001) 
     &    XIJ,RXIJ(XIJ),DIMAG(JXCIJINT(XIJ)),
     &    RXCCLS(1,XIJ),RXCCLS(2,XIJ),RXCCLS(3,XIJ),IXCP(XIJ)

        ENDDO

        CLOSE(73)
C
        ENDIF ! LMPIC = 1
C
      ENDIF ! called with INFO='R'
C=======================================================================
C finished subroutine if called with INFO = 'F'
C=======================================================================
C
C
73000 FORMAT("# off-diagonal exchange-coupling constants Jij ",/,
     &       "# for atom i = ",I1,/,
     &       "# j    R_ij( ALAT )   J_ij( Ry )      RXCCLS      ",
     &       "             IXCP")
73001 FORMAT(I3,3X,F9.5,6X,D9.3,6X,3(1X,F7.4),I5)
73002 FORMAT(I3,1X,I3,3X,F9.5,6X,D9.3,6X,3(1X,F7.4),I5)

C
      END





      SUBROUTINE CMATMUL(N,M,A,B,C)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix operation           C = A * B       *
C   *                                                                  *
C   *   A,B,C   complex  SQUARE  N x N - matrices                      *
C   *   N       dimension of A, B and C                                *
C   *   M       array size of A, B, C with M >= N                      *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT DOUBLE COMPLEX(A-H,O-Z)
C
C PARAMETER definitions
C
      DOUBLE COMPLEX C0
      PARAMETER (C0=(0.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER M,N
      DOUBLE COMPLEX A(M,M),B(M,M),C(M,M)
C
C Local variables
C
      DOUBLE COMPLEX BLJ
      INTEGER I,J,L
C
      DO J = 1,N
         DO I = 1,N
            C(I,J) = C0
         END DO
      END DO
C
      DO J = 1,N
         DO L = 1,N
            BLJ = B(L,J)
            IF ( BLJ.NE.C0 ) THEN
               DO I = 1,N
                  C(I,J) = C(I,J) + A(I,L)*BLJ
               END DO
            END IF
         END DO
      END DO
C
      END
