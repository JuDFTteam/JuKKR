      subroutine writekkrflex(natomimp,nspin,ielast,lmpot,lmmaxd, alat,
     +                  natyp,kshape,VBC,ATOMIMP,HOSTIMP,NOQ,ZAT,KAOEZ,
     +                  CONC,CMOM,CMINST,VINTERS)
     
      use mod_types, only: t_tgmat
      use mod_wunfiles, only: t_params, read_angles
      use mod_version_info
      use mod_md5sums
     
      implicit none
      include 'inc.p'
      INTEGER LMPOTD
      PARAMETER (LMPOTD=(LPOTD+1)**2)
      INTEGER LMGF0D
      PARAMETER (LMGF0D= (LMAXD+1)**2)
!interface
      integer          :: natomimp
      integer          :: nspin
      integer          :: ielast
      integer          :: lmpot
      integer          :: lmmaxd
      double precision :: alat
      integer       :: natyp
      integer       :: kshape
      DOUBLE PRECISION VBC(2)
      INTEGER ATOMIMP(NATOMIMPD)
      INTEGER HOSTIMP(0:NATYPD)
      integer       :: NOQ(NAEZD)
      double precision :: ZAT(NATYPD)
      integer          :: KAOEZ(NATYPD,NAEZD+NEMBD)
      DOUBLE PRECISION CONC(NATYPD)
      DOUBLE PRECISION CMOM(LMPOTD,NATYPD),CMINST(LMPOTD,NATYPD)
      DOUBLE PRECISION VINTERS(LMPOTD,NAEZD)
!local
      integer :: ispin,ie,i1,iatom,irec,i,lm
      DOUBLE COMPLEX TMAT0(LMMAXD,LMMAXD)
      DOUBLE PRECISION THETA(NATYPD), PHI(NATYPD)
C     .. External Functions ..
      LOGICAL OPT
      EXTERNAL OPT

      write(1337,*) 'KKRFLEX WRITEOUT'
      write(1337,*) OPT('KKRFLEX ')

      IF ( OPT('KKRFLEX ') ) THEN
        OPEN (6699,FILE='kkrflex_tmat',STATUS='unknown')
        call version_print_header(6699,
     &           '; '//md5sum_potential//'; '//md5sum_shapefun)
        write(6699,*) '#',NATOMIMP,NSPIN,IELAST,LMMAXD,KORBIT
        if (t_tgmat%tmat_to_file) then
           OPEN (69,ACCESS='direct',RECL=WLENGTH*4*LMMAXD*LMMAXD,
     &             FILE='tmat',  FORM='unformatted')
        end if

        !read in non-collinear angles
        call read_angles(t_params,NATYP,THETA,PHI)

        DO IATOM = 1,NATOMIMP
          I1=ATOMIMP(IATOM)
          IF (KORBIT.EQ.0) THEN
          DO ISPIN=1,NSPIN
            DO IE=1,IELAST
              IF (I1<=NATYP) THEN 
                  IREC = IE+IELAST*(ISPIN-1)+IELAST*NSPIN*(I1-1) 
                  if (t_tgmat%tmat_to_file) then
                     READ (69,REC=IREC) TMAT0
                  else
                       stop 'WRONG tmat_to_file for KKRFLEX writeout!!'
                       !not correctly read in with this option, to fix this communication is needed
!                      tmat0(:,:) = t_tgmat%tmat(:,:,irec)
                  end if
              ELSE
                 TMAT0=(0.0D0,0.0D0)
              END IF 
              WRITE(6699,'(4I,50000E)') IATOM,ISPIN,IE,0,TMAT0
          END DO !IE=1,IELAST
        END DO !ISPIN=1,NSPIN
        ELSEIF (KORBIT.EQ.1) THEN
         ISPIN=1
         DO IE=1,IELAST
          IF (I1<=NATYP) THEN
           IREC=IE+IELAST*(I1-1)
            READ(69,REC=IREC) TMAT0
            !perform here a transformation from the local to the global
            !spin-frame of reference
            CALL ROTATEMATRIX(TMAT0,THETA(I1),PHI(I1),LMGF0D,0)
          ELSE
           TMAT0=(0d0,0d0)
          ENDIF
          WRITE(6699,'(4I,50000E)') IATOM,ISPIN,IE,0,TMAT0         
         ENDDO
        ENDIF
        END DO
        CLOSE(69)
        CLOSE(6699)

        OPEN (91,FILE='kkrflex_intercell_ref',STATUS='unknown')
        call version_print_header(91,
     &           '; '//md5sum_potential//'; '//md5sum_shapefun)
        WRITE(91,*) '# Intercell potential of each atom'
        WRITE(91,*) '# '
        WRITE(91,*) '# NATOMIMP',NATOMIMP
        WRITE(91,*) '# lmpot',lmpot
        WRITE(91,*) '# KSHAPE',KSHAPE
        WRITE(91,*) '# NATOMIMP, lmpot, ALAT VBC(1), VBC(2)'
        WRITE(91,'(2I,10F)') NATOMIMP,lmpot,ALAT,VBC(1),VBC(2)
        DO IATOM = 1,NATOMIMP      ! Bauer 2011-10-11
          I=ATOMIMP(IATOM)         !
         write(1337,*) 'ac2',I,HOSTIMP(I),lmpot,
     +                (VINTERS(LM,I),LM=1,lmpot)
          WRITE(91,'(5000G)') (VINTERS(LM,I),LM=1,lmpot)
        END DO
        CLOSE(91)

        OPEN (91,FILE='kkrflex_intercell_cmoms',STATUS='unknown')
        call version_print_header(91,
     &           '; '//md5sum_potential//'; '//md5sum_shapefun)
        WRITE(91,*) '# Charge moments of each atom in the unit cell'
        WRITE(91,*) '# Values given are CMOM + CMINST'
        WRITE(91,*) '# First colums is the core charge other'
        WRITE(91,*) '# other colums the charge moments'
        WRITE(91,*) '# NATOMIMP',NATOMIMP
        WRITE(91,*) '# lmpot',lmpot
        WRITE(91,*) '# KSHAPE',KSHAPE

        DO IATOM = 1,NATOMIMP
          I=ATOMIMP(IATOM)
          I1=KAOEZ(1,I)
          WRITE(1337,*) 'NOQ',I,NOQ(I)
          IF (NOQ(I)/=1 .and. NOQ(I)/=0) 
     +      stop '[vmadelblk] VIRATOMS: NOQ/=1'
          IF (NOQ(I)==0) then
            WRITE(91,'(5000G)') 0.0D0,(0.0D0,LM=1,lmpot)
          ELSE
            IF ( KSHAPE.NE.0 ) THEN
              WRITE(91,'(5000G)') ZAT(I1), 
     +                          ((CMOM(LM,I1) 
     +             + CMINST(LM,I1))*CONC(I1),LM=1,lmpot)
            ELSE
              WRITE(91,'(5000G)') ZAT(I1), 
     +            (CMOM(LM,I1)*CONC(I1),LM=1,lmpot)
            END IF 
          END IF
        END DO
        CLOSE(91)


      END IF



      end subroutine writekkrflex
