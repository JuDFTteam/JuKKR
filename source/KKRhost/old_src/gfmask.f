      SUBROUTINE GFMASK(LINTERFACE,ICHECK,ICC,INVMOD,NSH1,NSH2,NAEZ,
     &                  NSHELL,NAEZD,NPRINCD)
C **********************************************************************
C *                                                                    *
C * This subroutine prepares the ICHECK matrix that is used for        *
C * calculating the proper off-diagonal GF matrix elements ( e.g.      *
C * impurity) in case of no full inversion algorithm                   *
C *                                                                    *
C * ICHECK(I,J) points to the block (I,J) of the GF matrix having the  *
C * size NPRINCD                                                       *
C *                                                                    *
C *                                            29.02.2000              *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar arguments
      INTEGER NAEZD,NPRINCD
      INTEGER ICC,INVMOD,NLAYER,NAEZ,NSHELL
      LOGICAL LINTERFACE
C     ..
C     .. Array arguments
      INTEGER ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD)    
      INTEGER NSH1(*),NSH2(*)
C     .. Local variables
      INTEGER ICOUPLE(NAEZD,NAEZD)
      INTEGER I,J,K,II,ISTEP1,ILT1,ISTEP2,ILT2,IL2,IL1,LFCHK
      CHARACTER*80 FMTCHK
      CHARACTER*35 INVALG(0:3)  ! GODFRIN
C     ..
C     .. External functions
      LOGICAL OPT,TEST
      EXTERNAL OPT,TEST
C     ..
C     .. Data statements
      DATA INVALG /'FULL MATRIX                        ',
     &             'BANDED MATRIX (slab)               ',
     &             'BANDED + CORNERS MATRIX (supercell)',
     &             'godfrin module                     ' /   ! GODFRIN

      WRITE (1337,99000)
C
C --> set default inversion to SUPERCELL mode = banded matrix + corners
C
      INVMOD = 2
C
C --> LINTERFACE = use band diagonal mode
C
      IF (LINTERFACE) INVMOD = 1
C
C --> full inversion is performed ONLY BY EXPLICIT request
C
      IF ( OPT('full inv') )  INVMOD = 0
!
! ----------------------------------------------------------------------
      if (opt('godfrin ')) invmod = 3  ! GODFRIN
! ----------------------------------------------------------------------
C        
C 21.10.2014 GODFRIN Flaviano 
      IF ( ( INVMOD.NE.0 ).AND.(INVMOD.NE.3).AND.
     &     ( MOD(NAEZ,NPRINCD).NE.0 ) ) THEN
         WRITE(6,99001) NAEZ,NPRINCD
         STOP
      END IF
C
      WRITE (1337,99002) INVALG(INVMOD)
C
      NLAYER=NAEZ/NPRINCD  
C ----------------------------------------------------------- INVMOD = 1
C                                                   band-diagonal matrix
      IF (INVMOD.EQ.1) THEN
         DO I=1,NLAYER
            DO J=1,NLAYER
               IF (I.EQ.J) THEN
                  ICHECK(I,J)=1
               ELSE
                  ICHECK(I,J)=0
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C ----------------------------------------------------------- INVMOD = 2
C              band-diagonal matrix with corners (slab periodic along z)
      IF (INVMOD.EQ.2) THEN
         DO I=1,NLAYER
            DO J=1,NLAYER
               IF ( ( I.EQ.J ).OR. ( J.EQ.NLAYER ) 
     &                        .OR. ( I.EQ.NLAYER ) ) THEN
                  ICHECK(I,J)=1
               ELSE
                  ICHECK(I,J)=0
               ENDIF
            ENDDO
         ENDDO
      ENDIF                     
C ================================================= INVMOD = 1, ICC <> 0
C                   band-diagonal matrix, off-diagonal G elements needed
C
C --> prepare the matrix ICOUPLE which has 1 in all nn' blocks 
C     (atomic sites) that are needed
C
C ======================================================================
      IF ( ( ICC.NE.0 ) .AND. ( INVMOD.EQ.1 ) ) THEN
         DO I=1,NAEZ
            DO J=1,NAEZ
               ICOUPLE(I,J) = 0
C
               DO II=1,NSHELL
                  IF ( ( ( NSH1(II).EQ.I ) .AND. ( NSH2(II).EQ.J ) ) 
     &            .OR. ( ( NSH1(II).EQ.J ) .AND. ( NSH2(II).EQ.I ) )
     &                 ) ICOUPLE(I,J) = 1 
               ENDDO
            ENDDO
         ENDDO
CcccC ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
CcccC                                               conductivity calculation
Cccc         IF (OPT('CONDUCT ')) THEN 
Cccc            DO I=1,NLAYER
Cccc               DO J=1,NLAYER
Cccc                  ICHECK(I,J)=0
Cccc               ENDDO
Cccc            ENDDO
Cccc            DO I=1,NAEZ
Cccc               DO J=1,NAEZ
Cccc                  ICOUPLE(I,J) = 0
Cccc               END DO
Cccc            END DO
Cccc            DO II=1,NCONDPAIR
Cccc               I = IATCONDL(II)
Cccc               J = IATCONDR(II)
Cccc               ICOUPLE(I,J) = 1       
Cccc               ICOUPLE(J,I) = 1
Cccc            ENDDO
Cccc         END IF                 ! Conductivity calculation
CcccC ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C ----------------------------------------------------------------------
C Now given the matrix ICOUPLE prepare the matrix ICHECK which has 1 in
C all principal-layer blocks that we need -- this will be used in the 
C matrix inversion
C ----------------------------------------------------------------------
         ISTEP1=0
         ILT1=1
C ----------------------------------------------------------------------
         DO IL1=1,NAEZ
            ISTEP1=ISTEP1+1
C
            IF ( ISTEP1.GT.NPRINCD ) THEN
               ILT1=ILT1+1
               ISTEP1=1
            ENDIF
C
            ILT2=1
            ISTEP2=0
C ......................................................................
            DO IL2=1,NAEZ 
               ISTEP2=ISTEP2+1
C
               IF ( ISTEP2.GT.NPRINCD ) THEN
                  ILT2=ILT2+1
                  ISTEP2=1
               ENDIF
C
               IF ( ICOUPLE(IL1,IL2).EQ.1 ) ICHECK(ILT1,ILT2)=1 
            END DO
C ......................................................................
         END DO
C ----------------------------------------------------------------------
C in the case of calculation of single blocks it has to put the correct 
C value to ICHECK in order to calculate all the elements also necessary 
C to calculate that single block          ?????
C ----------------------------------------------------------------------
         DO J=1,NLAYER
C
C --> loop over the element ICHECK(I,J) with fixed J and I < J
C
            IF (J.NE.1) THEN
               DO I=1,J-1
                  IF (ICHECK(I,J).EQ.1) THEN
                     DO K=I+1,J
                        ICHECK(K,J)=1
                     ENDDO
                     DO K=J,NLAYER
                        ICHECK(K,K)=1
                     END DO
                  ENDIF
               ENDDO
            ENDIF
C
            IF ( .NOT.OPT('CONDUCT ') ) THEN 
C
C --> loop over the element ICHECK(I,J) with fixed J and I > J
C
                IF (J.NE.NLAYER) THEN
                    DO I=NLAYER,J+1,-1
                        IF (ICHECK(I,J).EQ.1) THEN
                            DO K=I-1,J,-1
                                ICHECK(K,J)=1
                            ENDDO
                        ENDIF
                    ENDDO
                ENDIF
            END IF
        ENDDO
C ----------------------------------------------------------------------
      ENDIF      
C ======================================================================
C
      IF ( TEST('ICHECK  ') ) THEN 
C
         FMTCHK=' '
         LFCHK = 1
         DO I = 1,MIN(35,NLAYER)
            FMTCHK=FMTCHK(1:LFCHK)//'--'
            LFCHK = LFCHK+2
         END DO
C
         WRITE (1337,'(8X,A,/,8X,A)') 'ICHECK matrix :',FMTCHK(1:LFCHK)
         DO I=1,NLAYER
            WRITE (1337,'(9X,35I2)') (ICHECK(I,J),J=1,MIN(35,NLAYER))
         ENDDO
         WRITE (1337,'(8X,A,/)') FMTCHK(1:LFCHK)
      END IF
C
99000 FORMAT (5X,'< GFMASK > : set KKR matrix inversion algorithm',/)
99001 FORMAT(6X,"ERROR: Number of sites (NAEZ) =",I3,
     &     " not an integer multiplier",/,6X,
     &     "of principal layers (NPRINCD) =",I3,/,6X,
     &     "Use ONLY  full inversion in this case")
99002 FORMAT (8X,'INVERSION algorithm used : ',A35,/)
      END
