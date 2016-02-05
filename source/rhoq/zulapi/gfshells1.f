      SUBROUTINE GFSHELLS1(ICC,NSH1,NSH2,NSHELL,NAEZ,NATYP,
     &                    RBASIS,BRAVAIS,RATOM,RATOMS,
     &                    NSYMAT,ISYMINDEX,RSYMAT,RCLS,CLS,
     &                    EQINV,KAOEZ,NACLS,ATOM,RFCTOR,
     &                    IATCONDL,IATCONDR,NCONDPAIR,ROTNAME,
     &                    NATOMIMP,ATOMIMP,RCLSIMP)
      implicit none
c     ****************************************************
c     * This subroutine constructs mainly the index arrays
c     * NSH1,NSH1,NSHELL etc to be used to write out the 
c     * impurity Green's function
c     ****************************************************
      include 'inc.p'
      include 'inc.cls'
      INTEGER NSYMAXD
      PARAMETER (NSYMAXD=48)
      INTEGER ICC,NSH1(*),NSH2(*),IATCONDL(*),IATCONDR(*)
      INTEGER NSYMAT,ISYMINDEX(NSYMAXD),NAEZ,NCONDPAIR
      INTEGER NATOMIMP,ATOMIMP(NATOMIMPD)
      DOUBLE PRECISION RCLSIMP(3,NATOMIMPD)
      DOUBLE PRECISION RATOM(3,NSHELD),RATOMS(3,NSHELD),
     &     RBASIS(3,*),BRAVAIS(3,3),RSYMAT(64,3,3)
      DOUBLE PRECISION RCLSNEW(3,NATOMIMPD),
     &                 VEC1(3),VEC2(3,NAEZD),RCLS(3,NACLSD,*)
      INTEGER CLS(*),EQINV(*),KAOEZ(*),NACLS(*),
     +        NSHELL(0:NSHELD),
     &        NSH1S(NSHELD),NSH2S(NSHELD),NSHELLS(NSHELD)
      INTEGER ATOM(NACLSD,*)
      DOUBLE PRECISION RFCTOR,RSORT(NSHELD)
      DOUBLE PRECISION R,diff
c     
      INTEGER N,I,J,NATYP,POS,II,IC,N1,N2,N3
      INTEGER NS,IN,NMAX,K,IAT
      INTEGER ISORT(NSHELD)
      LOGICAL OPT,EXIST(NATYPD),TEST
      CHARACTER*10 ROTNAME(64)
c
c
      EXTERNAL DSORT,SHELLGEN2K,OPT,TEST
c --------------------------------------------------------

c--->   construction of ratom, nsh1 and nsh2 for a self-consistent
c       calculation
c     
       NSHELL(0) = NATYP
      
        DO 170 I=1,NSHELL(0)
           RATOM(1,I) = 0.0D0
           RATOM(2,I) = 0.0D0
           RATOM(3,I) = 0.0D0
           NSHELL(I) = 0

           DO 190 J=1,NAEZ
              IF (KAOEZ(J).EQ.I) THEN
                 NSHELL(I) = NSHELL(I) + 1
                 IF (NSHELL(I).EQ.1) THEN
                    NSH1(I) = J
                    NSH2(I) = J                    
                 END IF
              END IF               
 190       END DO

           IF (NSHELL(I).EQ.0) THEN
              WRITE(6,*) 'THERE ARE SOME INCONSISTENCIES ',
     +             'IN THE KAOEZ ARRAY.'
              WRITE(6,*) 'NOT ALL ATOMS DEFINED BY NATYP FOUND.'
              STOP
           END IF
 170    CONTINUE
        CALL SHELLGEN2K(RCLSIMP(1,1),ATOMIMP(1),NATOMIMP,RATOM(1,1),
     +                  NSHELL,NSH1,NSH2,RSYMAT,NSYMAT,ISYMINDEX,
     +                  RFCTOR,RBASIS,NAEZ)      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        N = 0
        DO NS=1,NSHELL(0)
          N = N + NSHELL(NS)     
        END DO
c
        IF (  NSHELL(0).GT.NSHELD ) THEN
          WRITE (6,*) 'Please change the parameter NSHELD in ',
     +         'inc.p to ',NSHELL(0)
          CALL RCSTOP('gfshells')
        END IF
 9000 FORMAT(3F12.5,F15.8)
 9010 FORMAT(2F12.6)
 9020 FORMAT(' NMESH : ',I4,'  NOFKS : ',I7,'  VOLBZ :',f14.8)
 9030 FORMAT(I3,I7,I5,F14.3,2F11.3,I8,F8.3)
       END









