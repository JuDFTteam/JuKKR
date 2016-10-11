      SUBROUTINE WRITE_TBKKR_FILES(LMAX,NEMB,NCLS,NATYP,NAEZ,IELAST,INS,
     +                             ALAT,BRAVAIS,RECBV,RBASIS,CLS,NACLS,
     +                             RCLS,EZOA,ATOM,RR)
      use mod_version_info, only: serialnr
      IMPLICIT NONE
      INCLUDE 'inc.p'
!interface
      INTEGER LMAX,NEMB,NCLS,NATYP,NAEZ,IELAST,INS
      DOUBLE PRECISION ALAT,BRAVAIS(3,3),RECBV(3,3),
     +                 RBASIS(3,NAEZD+NEMBD),RCLS(3,NACLSD,NCLSD),
     +                 RR(3,0:NRD)
      INTEGER CLS(NAEZD),NACLS(NCLSD),EZOA(NACLSD,NAEZD),
     +        ATOM(NACLSD,NAEZD)
!local
      INTEGER I1,I2,J
C     .. External Functions ..
      LOGICAL OPT
      EXTERNAL OPT


          OPEN(934,FILE='TBkkr_params.txt',FORM='formatted')          
          WRITE(934,'(A,A)') '#FILEVERSION= 2'//'   # serial: ',serialnr
          WRITE(934,'(I8,4X,A)') LMAXD, 'lmaxd',   LMAX, 'lmax',
     +                           KORBIT, 'korbit',
     +                           NSPIND, 'nspind', NRD, 'nrd',        
     +                           NEMBD, 'nembd',   NEMB, 'nemb',      
     +                           NCLSD, 'nclsd',   NCLS, 'ncls',      
     +                           NATYPD, 'natypd', NATYP, 'natyp',    
     +                           NAEZD, 'naezd',   NAEZ, 'naez',      
     +                           NACLSD, 'naclsd', IELAST, 'ielast',  
     +                           INS, 'ins'                           
          CLOSE(934)                                                  

          OPEN(935,FILE='TBkkr_container.txt',FORM='formatted')
          WRITE(935,'(A,A)') '#FILEVERSION= 2'//'   # serial: ',serialnr
          !write out lattice information
          write(935,'(A)') 'alat:'
          write(935,'(ES25.16)') ALAT
          write(935,'(A)') 'bravais:'
          write(935,'(3ES25.16)') ((BRAVAIS(I1,I2),I1=1,3),I2=1,3)
          write(935,'(A)') 'recbv:'
          write(935,'(3ES25.16)') ((RECBV(I1,I2),I1=1,3),I2=1,3)
          write(935,'(A)') 'RBASIS:'
          write(935,'(3ES25.16)') ((RBASIS(J,I1), J=1,3),
     +                              I1=1,NAEZD+NEMBD)

          !write out cluster information
          write(935,'(A)') 'CLS:'
          write(935,'(1I8)') (CLS(I1),I1=1,NATYPD)
          write(935,'(A)') 'NACLS:'
          write(935,'(1I8)') (NACLS(I1), I1=1,NCLSD)
          write(935,'(A)') 'RCLS:'
          do I2=1,NCLSD
           do I1=1,NACLSD
             write(935,'(3ES25.16)') RCLS(:,I1,I2)
           end do
          end do
          write(935,'(A)') 'EZOA:'
          write(935,'(1I8)') ((EZOA(I1,I2),I1=1,NACLSD),I2=1,NAEZD)
          write(935,'(A)') 'ATOM:'
          write(935,'(1I8)') ((ATOM(I1,I2),I1=1,NACLSD),I2=1,NAEZD)
          write(935,'(A)') 'RR:'
          do I1=0,NRD
            write(935,'(3ES25.16)') RR(:,I1)
          end do

          CLOSE(935)



      END SUBROUTINE WRITE_TBKKR_FILES
