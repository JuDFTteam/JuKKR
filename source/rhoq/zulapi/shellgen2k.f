c 23.2.2000 ***************************************************************
      SUBROUTINE SHELLGEN2K(RCLS,ATOM,NATOM,RATOM,NSHELL,
     +                    NSH1,NSH2,ND,NROT,ISYMINDEX,RFCTOR,
     &                    rbasis,naez)
c ************************************************************************
      implicit none
C     .. Parameters ..
      include 'inc.p'
      include 'inc.cls'
c
c     .. scalar arguments
      INTEGER 
     +     NATOM,                   ! number of atoms in cluster
     +     NROT                     ! number of symmetry operations in ND()
      DOUBLE PRECISION 
     +     RFCTOR                   ! 2*PI/ALATA
c     .. array arguments
      INTEGER ATOM(*),              ! number of position in EZ     
     +     NSH1(*),NSH2(*),
     +     NSHELL(0:NSHELD)         ! NSHELL(0) number of shells
                                    ! NSHELL(1..NSHELL(0)) number of atoms
                                    ! or pairs in shell (not really used)
      INTEGER ISYMINDEX(*)

      DOUBLE PRECISION 
     +     ND(64,3,*),
     +     RATOM(3,*),              ! difference vector of pair
     +     RCLS(3,*)                ! position of atom in cluster 
                                    ! (rel. to atom 1)
      double precision rbasis(3,*)
      integer naez
      INTEGER AI,AJ,FOUND,I,IROT,J,NS,NS0,NSS,ID,ISYM,II
      DOUBLE PRECISION R1,SMALL
      DOUBLE PRECISION RI(3),RJ(3)

      INTEGER NSHELL0
      PARAMETER(NSHELL0 =5500)
      INTEGER NSH1I(NSHELL0),NSH2I(NSHELL0),NSHELLI(NSHELL0)
      DOUBLE PRECISION RATOMI(3,NSHELL0)

      LOGICAL TEST
      EXTERNAL TEST

      DATA SMALL /  1.0D-10/
c ------------------------------------------------------------------------
      IF(TEST('flow    ')) write(6,*) '>>> SHELLGEN'
      WRITE(6,*) 'Number of cites to calculate Green Function ',natom
      WRITE(6,*) 'Number of rotation for this system          ',nrot
c      write(6,*) 'NATOM,NROT : ',natom,nrot

c      write(6,FMT='(3f12.4  )') ((rcls(j,i),j=1,3),i=1,natom)
c      write(6,FMT='(10I4)') (atom(i),i=1,natom)
c      write(6,FMT='((3(3I6,/),/)  )')
c     +     (((nd(n,i,j),i=1,3),j=1,3),n=1,nrot)

c
c ---> search shells for test (number.lt.NSHELD)
c

      NSS=0      
    

      DO 11 I = 1,NATOM
        AI = ATOM(I)
        
        IF (AI.LT.0) GOTO 11
        DO 21 J = 1,NATOM

c     it avoids to consider again the diagonal elements
c     already prepared (the shells) in the gll2k subroutine

          IF (J.EQ.I) GOTO 21

          AJ = ATOM(J)
         
          IF (AJ.LT.0) GOTO 21
          FOUND = 0                 ! found a symmetric equivalent 
                                    ! pair of atoms

          IF (NSHELL(0)+NSS.GE.1) THEN
            
            DO 31 ID = 1,NROT
              ISYM = ISYMINDEX(ID)
              DO 61 II = 1,3
                RI(II) = ND(ISYM,II,1)*RCLS(1,I) +
     +               ND(ISYM,II,2)*RCLS(2,I) + 
     +               ND(ISYM,II,3)*RCLS(3,I)
                RJ(II) = ND(ISYM,II,1)*RCLS(1,J) +
     +               ND(ISYM,II,2)*RCLS(2,J) + 
     +               ND(ISYM,II,3)*RCLS(3,J)
 61           CONTINUE

c     if ICC.gt.0 then I'm calculating just the cluster 
c     around the impurity. nshell(0) is equal to zero so 
c     the following loop is used just in case if I'm calculating
c     all the shell in the upper or lower diagonal.

              IF (NSHELL(0).GT.0) THEN
                DO 51 NS = 1,NSHELL(0) ! loop over shells found before
                                       ! in other clusters
                  IF (FOUND.EQ.0) THEN
                    IF ( (AI.EQ.NSH1(NS) .AND. 
     +                    AJ.EQ.NSH2(NS)) .OR.
     +                   (AI.EQ.NSH2(NS) .AND. 
     +                    AJ.EQ.NSH1(NS))  ) THEN
                      R1 = (RI(1)-RJ(1)+RATOM(1,NS))**2  +
     +                     (RI(2)-RJ(2)+RATOM(2,NS))**2  +
     +                     (RI(3)-RJ(3)+RATOM(3,NS))**2
                      IF (R1.LT.SMALL) THEN
                        FOUND = 1
                        NSHELL(NS) = NSHELL(NS) + 1
                        
        write (6,*) 'entered the loop NSHELL(0).GT.0 and
     +                       found a previous shell'
                        
                      END IF        ! (R1.LT.SMALL)
                    END IF
                  END IF            ! (FOUND.EQ.0)
 51             END DO              ! NS = 1,NSHELL(0)
              END IF

cccccccc

c     when the rotation and the representative pair that
c     identify a pair of atoms is found then FOUND=1 and 
c     the search for a different pair of atoms starts.

              IF (NSS.GT.0 .AND. FOUND.EQ.0) THEN
                DO 52 NS = 1,NSS    ! loop over shells found before
                                    ! in this cluster
                  IF (FOUND.EQ.0) THEN
                    IF ( (AI.EQ.NSH1I(NS) .AND. 
     +                    AJ.EQ.NSH2I(NS)) .OR.
     +                   (AI.EQ.NSH2I(NS) .AND. 
     +                    AJ.EQ.NSH1I(NS))  ) THEN
                      R1 = (RI(1)-RJ(1)+RATOMI(1,NS))**2  +
     +                     (RI(2)-RJ(2)+RATOMI(2,NS))**2  +
     +                     (RI(3)-RJ(3)+RATOMI(3,NS))**2
                      IF (R1.LT.SMALL) THEN
                        FOUND = 1
                        NSHELLI(NS) = NSHELLI(NS) + 1
                      END IF        ! (R1.LT.SMALL)
                    END IF
                  END IF            ! (FOUND.EQ.0)
 52             END DO              ! NS = 1,NSS
              END IF

 31         END DO                  ! ID = 1,NROT
          END IF                    ! (NSHELL(0)+NSS.GE.1)
            
          IF (FOUND.EQ.0) THEN
            NSS = NSS + 1
            NSH1I(NSS) = AI
            NSH2I(NSS) = AJ
            NSHELLI(NSS) = 1
            DO 41 II=1,3
              RATOMI(II,NSS) = RCLS(II,J)-RCLS(II,I)
 41         END DO

            

          END IF                    ! (FOUND.EQ.0)
         

          IF (NSS.GT.NSHELL0) THEN
            write(6,*) 
     +           'Please increase the parameter NSHELL0 in ',
     +           'subroutine SHELLGEN in ''str.f'''
            STOP 'Dimension error.'
          END IF
 21     END DO                      ! J = 1,NATOM
 11   END DO                        ! I = 1,NATOM
c
c ---> test number of shells
c
      IF (NSS+NSHELL(0).GT.NSHELD) THEN
        write(6,*) 
     +       'Please change the parameter NSHELD in ''inc.p'' to ',
     +       NSS,NSHELL(0),NSS+NSHELL(0)
        STOP 'Dimension error.'
      END IF

      IF (NSS.GT.0) THEN
        DO 10 NS = 1,NSS
          NS0 = NSHELL(0) + NS
          NSH1(NS0) = NSH1I(NS) 
          NSH2(NS0) = NSH2I(NS)


          
          NSHELL(NS0) = NSHELLI(NS)
          DO 40 II=1,3
            RATOM(II,NS0) = RATOMI(II,NS)
c test test test 10.12.2001 nikos
c            RATOM(II,NS0) = RATOMI(II,NS) + RBASIS(II,NSH2(NS0)) -
c     &                                      RBASIS(II,NSH1(NS0)) 
c test test test 10.12.2001 nikos
 40       END DO
 10     END DO                      ! NS = 1,NSS
        NSHELL(0) = NSHELL(0) + NSS
      END IF                        ! (NSS.GT.0)

      IF(TEST('flow    ')) write(6,*) '<<< SHELLGEN'

      RETURN
      END                           ! SUBROUTINE SHELLGEN




