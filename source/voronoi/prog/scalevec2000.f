        SUBROUTINE SCALEVEC2000(RBASIS,ABASIS,BBASIS,CBASIS,NLBASIS,
     &     NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT, 
     &     TLEFT,TRIGHT,LINTERFACE,NAEZ,NEMB,BRAVAIS,KAOEZ)
      implicit none
      include 'inc.geometry'
      INTEGER  NRBASIS,NRIGHT,NLBASIS,NLEFT,NAEZ,NEMB
      REAL*8        ABASIS,BBASIS,CBASIS
      REAL*8        RBASIS(3,*),ZPERLEFT(3),ZPERIGHT(3),
     &                 TLEFT(3,*),TRIGHT(3,*),BRAVAIS(3,3)
       
      REAL*8        TEMP(3),RBASIS1(3,NAEZD+NEMBD)
      INTEGER I,J,I1,IER
      INTEGER KAOEZ(NAEZD+NEMBD)
      REAL*8        TX,TY,TZ    
      CHARACTER*200 UIO    
      LOGICAL LCARTESIAN,LINTERFACE   
c
c---->  normalization of basis vectors
c
      DO 1 I=1,NAEZ+NEMB
        RBASIS1(1,I)=RBASIS(1,I)/ABASIS
        RBASIS1(2,I)=RBASIS(2,I)/BBASIS
        RBASIS1(3,I)=RBASIS(3,I)/CBASIS
 1    CONTINUE
ccccccccccccccccccccccccccccccccccccc added 11.10.99
      IF (LINTERFACE) THEN
         DO I=1,NLBASIS
            TLEFT(1,I)=TLEFT(1,I)/ABASIS
            TLEFT(2,I)=TLEFT(2,I)/BBASIS
            TLEFT(3,I)=TLEFT(3,I)/CBASIS
         ENDDO
         ZPERLEFT(1)=ZPERLEFT(1)/ABASIS
         ZPERLEFT(2)=ZPERLEFT(2)/BBASIS
         ZPERLEFT(3)=ZPERLEFT(3)/CBASIS
      
         DO I=1,NRBASIS
            TRIGHT(1,I)=TRIGHT(1,I)/ABASIS
            TRIGHT(2,I)=TRIGHT(2,I)/BBASIS
            TRIGHT(3,I)=TRIGHT(3,I)/CBASIS
         ENDDO      
         ZPERIGHT(1)=ZPERIGHT(1)/ABASIS
         ZPERIGHT(2)=ZPERIGHT(2)/BBASIS
         ZPERIGHT(3)=ZPERIGHT(3)/CBASIS
      ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c ---> normalization of atomic positions in the unit cell
c
      write(6,*) 'position of atoms in unit cell :'
      LCARTESIAN= .false.   ! defalt is false then bravais lattice
                           ! coordinates are used...
      CALL IoInput('CARTESIAN ',UIO,1,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) LCARTESIAN
      WRITE(*,*) 'SCALEVEC2000: CARTESIAN=',LCARTESIAN
c    
c if lcartesian is true cartesian coordinates are used
c else the basis are in units of the 
              
              IF (.NOT.LINTERFACE) THEN
                 IF (.NOT.LCARTESIAN) THEN  ! Rescale lattice
                    DO 2 I=1,NAEZ+NEMB
                       DO 3 J=1,3
                          RBASIS(J,I)=( RBASIS1(1,I)*BRAVAIS(J,1)
     +                         +RBASIS1(2,I)*BRAVAIS(J,2)
     +                         +RBASIS1(3,I)*BRAVAIS(J,3) ) ! /ALATC
 3                     CONTINUE
                       WRITE(6,2025) (RBASIS(J,I), J=1,3),KAOEZ(I) ! 1.11.99
 2                  CONTINUE
                 ELSE           !.NOT.LCARTESIAN
c     changed by v.Bellini 21/10/99
                    DO J=1,NAEZ+NEMB
                       DO I=1,3
                          RBASIS(I,J) = RBASIS1(I,J)
                       END DO
                       WRITE(6,2025) (RBASIS(I,J), I=1,3),KAOEZ(J) ! 1.11.99
                    ENDDO
c     end of the change
                 END IF
              ELSE              ! .NOT.LINTERFACE
                 IF (.NOT.LCARTESIAN) THEN
                    DO I=1,NAEZ+NEMB
                       DO J=1,2
                          RBASIS(J,I)=( RBASIS1(1,I)*BRAVAIS(J,1)
     +                         +RBASIS1(2,I)*BRAVAIS(J,2) ) 
                       END DO

                       RBASIS(3,I) = RBASIS1(3,I) ! added 11.10.99      

                       WRITE(6,2025) (RBASIS(J,I), J=1,3),KAOEZ(I) ! 1.11.99
                    END DO
c     -------------------------------------------------------
c     Do the same for the boundary vectors
c     
                    DO I=1,NLBASIS
                       DO I1=1,2
                          TEMP(I1) = TLEFT(I1,I)
                       END DO
                       DO J=1,2
                          TLEFT(J,I) =( TEMP(1)*BRAVAIS(J,1)
     +                         +TEMP(2)*BRAVAIS(J,2) ) 
                       END DO
                    END DO
                    DO I1=1,2
                       TEMP(I1) = ZPERLEFT(I1)
                    END DO
                    DO J=1,2
                       ZPERLEFT(J) =( TEMP(1)*BRAVAIS(J,1)
     +                      +TEMP(2)*BRAVAIS(J,2) ) 
                    END DO
c     Now right
                    DO I=1,NRBASIS
                       DO I1=1,2
                          TEMP(I1) = TRIGHT(I1,I)
                       END DO
                       DO J=1,2
                          TRIGHT(J,I) =( TEMP(1)*BRAVAIS(J,1)
     +                         +TEMP(2)*BRAVAIS(J,2) ) 
                       END DO
                    END DO
c     
                    DO I1=1,2
                       TEMP(I1) = ZPERIGHT(I1)
                    END DO
                    DO J=1,2
                       ZPERIGHT(J) =( TEMP(1)*BRAVAIS(J,1)
     +                      +TEMP(2)*BRAVAIS(J,2) ) 
                    END DO
c     -------------------------------------------------------
c     
                 ELSE
                    
                    DO I=1,3
                       DO J=1,NAEZ+NEMB
                          RBASIS(I,J) = RBASIS1(I,J)
                       ENDDO
                    ENDDO
                    
                 END IF         !  (.NOT.LCARTESIAN) 
                 WRITE(6,9470)
                 DO I=NLEFT,1,-1
                    DO I1=NLBASIS,1,-1
                       tx = TLEFT(1,i1) + (I-1)*ZPERLEFT(1)
                       ty = TLEFT(2,i1) + (I-1)*ZPERLEFT(2)
                       tz = TLEFT(3,i1) + (I-1)*ZPERLEFT(3)
                       WRITE(6,9420) (I-1)*NLBASIS+i1,tx,ty,tz,
     &                                                 kaoez(NAEZ+i1)
                    END DO 
                 END DO
                 WRITE(6,9475)
                 DO I=1,NAEZ
                    WRITE(6,9420) I, (RBASIS(I1,I),I1=1,3),kaoez(i)
                 END DO
                 WRITE(6,9480)
                 DO I=1,NRIGHT
                    DO I1=1,NRBASIS
                       tx = TRIGHT(1,i1) + (I-1)*ZPERIGHT(1)
                       ty = TRIGHT(2,i1) + (I-1)*ZPERIGHT(2)
                       tz = TRIGHT(3,i1) + (I-1)*ZPERIGHT(3) 
                       WRITE(6,9420) (I-1)*NRBASIS+i1,tx,ty,tz,
     &                                         kaoez(NAEZ+NLBASIS+i1)       
                    END DO 
                 END DO   
                 
              END IF            !   .NOT.LINTERFACE 
 9420            format(I5,3F12.6,I5)
                 
 9470            format('--------------- Left  Host -------------- ')
 9475            format('---------------   S L A B  -------------- ')
 9480            format('--------------- Right Host -------------- ')    
 2025  FORMAT((3F15.8,I6))
c
c NOW !!! RBASIS are the basis vectors in units of au/alat in (xyz) 
c         reference
c
          RETURN
          END









