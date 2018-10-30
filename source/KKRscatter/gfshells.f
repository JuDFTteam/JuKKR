      SUBROUTINE GFSHELLS(ICC,NATOMIMP,NSH1,NSH2,NSHELL,NAEZ,NATYP,
     &                    RBASIS,BRAVAIS,RATOM,RATOMS,RCLSIMP,
     &                    NSYMAT,ISYMINDEX,RSYMAT,RCLS,CLS,
     &                    EQINV,KAOEZ,NACLS,ATOM,ATOMIMP,RFCTOR,
     &                    IATCONDL,IATCONDR,NCONDPAIR,ROTNAME,HOSTIMP)
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
      DOUBLE PRECISION RATOM(3,NSHELD),RATOMS(3,NSHELD),
     &     RBASIS(3,*),BRAVAIS(3,3),RSYMAT(64,3,3)
      DOUBLE PRECISION RCLSIMP(3,NATOMIMPD),RCLSNEW(3,NATOMIMPD),
     &                 VEC1(3),VEC2(3,NAEZD),RCLS(3,NACLSD,*)
      INTEGER CLS(*),EQINV(*),KAOEZ(*),NACLS(*),
     +        NSHELL(0:NSHELD),ATOMIMP(NATOMIMPD),
     &        NSH1S(NSHELD),NSH2S(NSHELD),NSHELLS(NSHELD),
     &        IMPHOST(natypd),hostimp(0:natypd)
      INTEGER ATOM(NACLSD,*)
      DOUBLE PRECISION RFCTOR,RSORT(NSHELD)
      DOUBLE PRECISION R,diff
c     
      INTEGER N,I,J,NATOMIMP,NATYP,POS,II,IC,N1,N2,N3
      INTEGER NS,IN,NMAX,K,IAT,NHOST,nrep,ndim,natper
      INTEGER ISORT(NSHELD),NSHOUT(NATOMIMPD)
      LOGICAL OPT,exist(NATYPD),test
      CHARACTER*10 ROTNAME(64)
c
c
      EXTERNAL DSORT,SHELLGEN2K,OPT,test
c --------------------------------------------------------
      if (icc.gt.0) WRITE(6,*) 'ICC     = ',ICC
      IF (OPT('CONDUCT ')) THEN
          WRITE(6,*) '::::::::::::::::::::::::::::::::::::::::'
          WRITE(6,*) 'SUBROUTINE GFSHELLS IN CONDUCTANCE MODE '
          WRITE(6,*) '::::::::::::::::::::::::::::::::::::::::'
          ICC = 1
      END IF
c     
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
        
c            WRITE(6,*) 'NSHELL    = ',NSHELL(0)
c            WRITE(6,*) 'NSH1(I)   : ',(NSH1(I)  ,I=1,NSHELL(0))
c            WRITE(6,*) 'NSH2(I)   : ',(NSH2(I)  ,I=1,NSHELL(0))
c            WRITE(6,*) 'NSHELL(I) : ',(NSHELL(I),I=1,NSHELL(0))
        IF (ICC.NE.0) THEN      ! (ICC.NE.0)
           
           IF (ICC.GT.0) THEN
              

              IF (OPT('CONDUCT ')) THEN
c In case of conductance calculation set up the rclsimp
                N = 0 
                DO I=1,NCONDPAIR
                   IAT = IATCONDL(I)
                   N = N + 1
                   DO J=1,3
                   RCLSIMP(J,N) = RBASIS(J,IAT)
                   END DO
c
                   IAT = IATCONDR(I)
                   N = N + 1
                   DO J=1,3
                   RCLSIMP(J,N) = RBASIS(J,IAT)
                   END DO
                EnD DO
                NATOMIMP = N
              ELSE
c     it reads the cluster coordinates from an external file
              REWIND 25
              READ (25,FMT=*) NATOMIMP
              write(6,*) 'Reading impurity cites .....',natomimp
              DO I=1,NATOMIMP
                 READ (25,FMT=*) (RCLSIMP(J,I),J=1,3),ATOMIMP(I)
                 write(6,*) i,(RCLSIMP(J,I),J=1,3),ATOMIMP(I)
                 ATOMIMP(I) = ATOMIMP(I) +ICC -1 
              ENDDO
      
              END IF

              IF (NATOMIMP.GT.NATOMIMPD ) THEN
                 WRITE (6,FMT=*) ' ERROR IN GLL99: DIMENSIONS TOO SMALL'
                 WRITE (6,FMT=*) 
     &                    'NATOMIMP=',  NATOMIMP,'NATOMIMPD=',NATOMIMPD
                 CALL RCSTOP('gfshells')
              END IF
              
c     
c--->  shells around atom icc are prepared for storing the 
c      cluster-gf in subroutine kkrmat in GMATLL(LMAXSQ,LMAXSQ,*) 
c     
               if (TEST('impurity')) then
      write(6,*) 'This checks if the impurity atoms are ok'
      write(6,*) 'CAUTION the TB-clusters are used for this test'
      write(6,*) 'This will give wrong results if tb-cluster too small'
              do i=1,natomimp
                 
                 do j=1,3
                    rclsnew(j,i) = 
     &                   rclsimp(j,i) + rbasis(j,icc)
                 enddo
                 nmax = 7
                 do n1=-nmax,nmax
                    do n2=-nmax,nmax
                       do n3=-nmax,nmax
                          do j=1,3
                             vec1(j) = n1*bravais(j,1)+n2*bravais(j,2)+
     +                            n3*bravais(j,3)
                          enddo
                          do k=1,naez
                             do j=1,3
                                vec2(j,k) = vec1(j) + rbasis(j,k)
                             enddo                              
                             diff=sqrt((rclsnew(1,i)-vec2(1,k))**2 +
     +                            (rclsnew(2,i)-vec2(2,k))**2 +
     +                            (rclsnew(3,i)-vec2(3,k))**2 )
                             
                             if (diff.le.(1.d-5)) then
                                if (k.ne.atomimp(i)) THEN
                              WRITE(6,*) 'ATOM TYPE OF ATOM ',I,' WAS ',
     &                        atomimp(i),' AND IS SET TO ',k
                             END IF
                                atomimp(i) = k
                             endif
                             
                          enddo !atoms in the cell
                       enddo
                    enddo
                 enddo
                  
              enddo             ! natomimp
              end if            ! test(impurity)
              do i=1,natyp
                 exist(i) = .false.
              end do

              write(6,8050)              
              do i=1,natomimp
                 write (6,111) i,atomimp(i)
 111             format('The impurity atom  ',I5,
     +      '  corresponds to atom  ',I5,'  of the host ')
              if (atomimp(i).gt.natyp) STOP ' Error atomoimp, GFSHELLS'
                 exist(atomimp(i)) = .true.
              end do  
              
              nhost = 0
              do i=1,natyp
                 imphost(i) = 0
                 if (exist(i)) THEN
                    nhost = nhost + 1
                    imphost(i) = nhost
                    hostimp(nhost) = i
                    write(6,*) 'Host type ',i,
     &                   '  will be impurity-host type',nhost 
                 end if
              end do
              hostimp(0) = nhost 
              if (nhost.ne.natyp) Write(6,8040)
c              
              CALL SHELLGEN2K(RCLSIMP(1,1),ATOMIMP(1),
     +             NATOMIMP,RATOM(1,1),NSHELL,NSH1,NSH2,
     +             RSYMAT,NSYMAT,ISYMINDEX,RFCTOR,rbasis,
     &             naez)            ! changed
           
C     
           ELSE                 ! (ICC.GT.0)
c
c --->        all shells are prepared
c
              DO 250 I=1,NAEZ
                 IC=CLS(EQINV(I))
                 write(6,*) 'I,IC:',I,IC
                 CALL SHELLGEN2K(RCLS(1,1,IC),ATOM(1,I),
     +                NACLS(IC),RATOM(1,1),NSHELL,NSH1,NSH2,
     +                RSYMAT,NSYMAT,ISYMINDEX,RFCTOR,rbasis,
     &                naez)   
                 write(6,*) 'NSHELL(0) :',NSHELL(0)
 250          END DO            ! I=1,NAEZ
c     
c --->        interchange NSH1 and NSH2 if NSH1.gt.NSH2
c
              DO 103 I=1,NSHELL(0)
                 IF (NSH1(I).GT.NSH2(I)) THEN
                    IN = NSH2(I)
                    NSH2(I) = NSH1(I)
                    NSH1(I) = IN
                    DO 102 IN=1,3
                       RATOM(IN,I) = -RATOM(IN,I)
 102                END DO
                 END IF
 103          END DO
c
c --->        sort the shells due to distance
c
              DO 104 I=1,NSHELL(0)
                 DO 101 IN=1,3
                    RATOMS(IN,I) = RATOM(IN,I)
 101             END DO
                 NSH1S(I)=NSH1(I)
                 NSH2S(I)=NSH2(I)
                 NSHELLS(I)=NSHELL(I)
                 R=RATOM(1,I)**2+RATOM(2,I)**2+RATOM(3,I)**2
                 RSORT(I) = SQRT(R)+1.D-7*NSH1(I)
 104          END DO
              CALL DSORT(RSORT,ISORT,NSHELL(0),POS)
              DO 105 I=1,NSHELL(0)
                 POS=ISORT(I)
                 DO 106 IN=1,3
                    RATOM(IN,I) = RATOMS(IN,POS)
 106             END DO
                 NSH1(I)=NSH1S(POS)
                 NSH2(I)=NSH2S(POS)
                 NSHELL(I)=NSHELLS(POS)
 105          END DO
              
           END IF               ! (ICC.GT.0)
           
c           DO 210 NS=1,NSHELL(0)
c              DO 220 I=1,3
c                 JSHELL(I,NS) = RATOM(I,NS)*2.0D0
c 220          END DO
c              JSHELL(4,NS) = NSH1(NS)
c              JSHELL(5,NS) = NSH2(NS)
c 210       END DO
           
        END IF                  ! (ICC.EQ.0)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (icc.gt.0) then
        write(6,*) 'Different shells for GF calculation NSHELL :',
     &                                              NSHELL(0)
        write(6,*) 'SHELL CITE1 CITE2   '   
        write(6,FMT=9030) 
     +       (NS,NSH1(NS),NSH2(NS),(RATOM(II,NS),II=1,3),
     +       NSHELL(NS),
     +       SQRT(RATOM(1,NS)**2+RATOM(2,NS)**2+RATOM(3,NS)**2),
     +       NS=1,NSHELL(0))
        end if
        N = 0
        DO NS=1,NSHELL(0)
          N = N + NSHELL(NS)     
        END DO
        write(6,*) 'number of block elements to be calculated :',N
c
c now write out the tucoef file for the impurity calculation
c
              NREP = 1
              NDIM = 1
              NATPER = natomimp
              WRITE(6,*) NATPER,NATOMIMP
              DO i=1,natper
                 NSHOUT(i) = 1
                WRITE(6,*) I,NSHOUT(I)
              end do
              if (icc.ne.0) THEN
              open(58,file='impurity.coefs',form='formatted') 
              write(58,8001) NREP,natomimp,lmaxd,NATPER,
     &             (NSHOUT(I),I=1,NATPER)
              write(58,8010)
              do i=1,natomimp
              if (natyp.eq.nhost) THEN 
              write(58,8002) 
     &                (RCLSIMP(J,I),J=1,3),atomimp(i),i,i,
     & SQRT(RCLSIMP(1,I)**2+ RCLSIMP(2,I)**2+RCLSIMP(3,I)**2),
     &                atomimp(i)
              else
              write(58,8002) 
     &                (RCLSIMP(J,I),J=1,3),imphost(atomimp(i)),i,i,
     & SQRT(RCLSIMP(1,I)**2+ RCLSIMP(2,I)**2+RCLSIMP(3,I)**2),
     &                atomimp(i)
              end if
              end do
              write(58,8001) natomimp*(lmaxd+1)**2,N
              write(58,8001) NDIM
              write(58,8005) nhost
              write(58,8006) (hostimp(i),i=1,nhost) 
              write(58,8007) NSYMAT
              write(58,8003) (rotname(isymindex(i)),i=1,nsymat)
              write(58,8004)
              write(58,8001) NSHELL(0)
              write(58,8030) 
     +       (NS,NSH1(NS),NSH2(NS),(RATOM(II,NS),II=1,3),
     +       NSHELL(NS),
     +       SQRT(RATOM(1,NS)**2+RATOM(2,NS)**2+RATOM(3,NS)**2),
     +       NS=1,NSHELL(0))  
              close(58)
             end if
 8001  format(11i5)
 8002  FORMAT (3f12.8,i4,i4,i5,F10.6,I5)
 8003  format (5A10)
 8004  format ('Shells in the reduced format')
 8005  format ('Host order, no of host sites: ',I5)
 8006  format (12I4)
 8007  format (I5,'    Symmetries for the Bulk')
 8010  format ('     Position of Impurity            Host Imp Shell',
     &         '   Dist     Host id in Bulk')
 8030  FORMAT(3I5,3F12.8,I8,F10.5)
 8040  format(' >>> Some host atoms are missing in the impurity <<<'/
     &        '     INDEXING WILL BE CHANGED:  PLEASE CHECK ALSO  ') 
 8050  format(' ******************************************'/
     &        ' *  Prepearing indexing for impurity GF   *'/
     &        ' *  Unsymmetrized GF is written out in    *'/
     &        ' *  this version                          *'/
     &        ' *  Producing:   green,                   *'/
     &        ' *               impurity.coefs,          *'/
     &        ' *               intercell_ref            *'/
     &        ' *                            v. 20.09.01 *'/
     &        ' ******************************************')
c
c
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









