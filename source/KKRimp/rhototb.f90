MODULE MOD_RHOTOTB

  CONTAINS
    !-------------------------------------------------------------------------
    !> Summary: Calculation of total charge density
    !> Category: physical-observables, KKRimp
    !> Get total charge density (as sum of valence and core charge)
    !> and output atomic charge and spin-moments.
    !>
    !>   add core and valence density expanded in spherical harmonics
    !>       ( convention see subroutine rholm )
    !>   in the paramagnetic case (nspin=1) the core valence charge times
    !>       r**2 is add to the valence charge density times r**2
    !>       then only rho2ns(irmd,lmxtsq,natypd,1) is used .
    !>   in the spin-polarized case (nspin=2) the spin-splitted core
    !>       charge density times r**2 is converted into core charge
    !>        density times r**2 and core spin density times r**2 .
    !>        then these parts are added to corresponding parts of
    !>        the valence densities times r**2 , that are rho2ns(...,1)
    !>        which contains the charge density  and rho2ns(...,2) which
    !>        contains in that case the spin density .
    !>             (see notes by b.drittler)
    !>
    !>    attention : the core density is spherically averaged and multi-
    !>                plied by 4 pi. therefore the core density is only
    !>                added to l=0 part .
    !>
    !>                              b.drittler   nov. 1989
    !>
    !>    total orbital moment within the WS sphere is also calculated
    !>    in the relativistic case; orbital density is normalised in the
    !>    same way as the charge density. v.popescu march 2002
    !-------------------------------------------------------------------------
    SUBROUTINE RHOTOTB(ITSCF,NATOM,NSPIN,LMAXATOM,density, &
                         ZATOM,CELL,SHAPEFUN,INS)

      use nrtype
!       use configparams, only: ins
      use type_density
      use type_cell
      use type_shapefun
      use mod_simp3
      use mod_simpk
      use mod_config, only: config_testflag
      use global_variables, only: ipand
      implicit none
!C     .. Parameters ..
!       include 'inc.p'
!       INTEGER LMPOTD
!       PARAMETER (LMPOTD= (LPOTD+1)**2)
      integer                      :: itscf
      integer                      :: natom
      integer                      :: nspin
      integer                      :: lmaxatom(natom)
      type(density_type)           :: density(natom)
      real(kind=dp)                :: zatom(natom)
      type(cell_type)              :: cell(natom)
      type(shapefun_type)          :: shapefun(natom)
      integer                      :: ins

!local
      integer                      :: ispindn,ispinup
      real(kind=dp)                :: rho(cell(1)%nrmaxd)
      real(kind=dp)                :: chrgnt
      real(kind=dp)                :: catom(nspin,natom)

      real(kind=dp)                :: DIFF,FACTOR,RFPI,SUM,TOTSMOM
      integer                      :: IR,IATOM,IFUN,IPAN1,IRC1,IRS1,ISPIN, &
                                      LM,LMPOT
!       INTRINSIC ATAN,SQRT
!       write(*,*) 'nspin',nspin,'natom',natom
      RFPI = SQRT(16.0D0*ATAN(1.0D0))
      TOTSMOM=0.0D0
      write(*,*) ''
      DO IATOM = 1, NATOM 
        LMPOT = (2*LMAXATOM(IATOM)+1)**2

          IF (NSPIN.EQ.2) THEN
            ISPINDN = 1 !2*IATOM - 1
            ISPINUP = 2 !*IATOM
            FACTOR = 1.0D0
          ELSE
            ISPINDN = 1 !2*IATOM - 1
            ISPINUP = 1 !*IATOM
!             IPOTD = IATOM
!             IPOTU = IATOM
            FACTOR = 0.5D0
          END IF

          IF (INS.NE.0) THEN
            IPAN1 = CELL(IATOM)%NPAN !(IATOM)
            IRS1  = CELL(IATOM)%NRCUT(1)
            IRC1  = CELL(IATOM)%NRCUT(IPAN1)
          ELSE
            IRS1 = CELL(IATOM)%NRMAX!IRWS(IATOM)
            IRC1 = IRS1
          END IF
!

!-----------------------------------------------------------------------
!--->     convert core density
!-----------------------------------------------------------------------
          DO IR = 2,IRS1
            SUM  = (DENSITY(IATOM)%RHOC(IR,ISPINUP) + DENSITY(IATOM)%RHOC(IR,ISPINDN))*FACTOR/RFPI
            DIFF = (DENSITY(IATOM)%RHOC(IR,ISPINUP) - DENSITY(IATOM)%RHOC(IR,ISPINDN))/RFPI
!--->       add this to the lm=1 component of rho2ns
            DENSITY(IATOM)%RHO2NS(IR,1,1)     = DENSITY(IATOM)%RHO2NS(IR,1,1) + SUM
            DENSITY(IATOM)%RHO2NS(IR,1,NSPIN) = DENSITY(IATOM)%RHO2NS(IR,1,NSPIN) + DIFF
          END DO
!-----------------------------------------------------------------------


          if (config_testflag('spinsplit')) then
            if (itscf==1) then
              write(*,*) 'shifting the density in the first iteration'
              density(iatom)%RHO2NS(1:cell(iatom)%nrcut(1),1,2)=20.0D0 / ( 4.0D0/3.0D0*cell(iatom)%rmesh(cell(iatom)%nrcut(1))**3 )
            end if
          end if


!-----------------------------------------------------------------------
!--->   calculate  charge and moment of the atom
!-----------------------------------------------------------------------
          DO ISPIN = 1,NSPIN
            IF (INS==0) THEN
!--->         integrate over wigner seitz sphere - no shape correction
              CALL SIMP3(DENSITY(IATOM)%RHO2NS(1,1,ISPIN), &
                         SUM,1,IRS1,CELL(IATOM)%DRMESHDI(:))
!--->         the result has to be multiplied by sqrt(4 pi)
!             (4 pi for integration over angle and 1/sqrt(4 pi) for
!             the spherical harmonic y(l=0))
              SUM = SUM*RFPI
            ELSE  !INS==0
!--->         convolute charge density with shape function to get the
!             charge in the exact cell - if kshape .gt. 0
              DO IR = 1,IRS1
                RHO(IR) = DENSITY(IATOM)%RHO2NS(IR,1,ISPIN)*RFPI
              END DO
!
              DO IR = IRS1 + 1,IRC1
                RHO(IR) = 0.0D0
              END DO
!
              DO IFUN = 1,SHAPEFUN(IATOM)%NLMSHAPE
                LM = SHAPEFUN(IATOM)%index2lm(IFUN)
!                 write(*,*) 'lm',lm
                IF (LM.LE.LMPOT) THEN
                  DO IR = IRS1 + 1,IRC1
                    RHO(IR) = RHO(IR) + DENSITY(IATOM)%RHO2NS(IR,LM,ISPIN)* &
                             SHAPEFUN(IATOM)%THETAS(IR-IRS1,IFUN)
                  END DO 
                END IF
              END DO

!--->         integrate over circumscribed sphere
              ipand = ipan1
              CALL SIMPK(RHO,SUM,IPAN1,CELL(IATOM)%NRCUT(:),CELL(IATOM)%DRMESHDI(:)) !,ipan1)

            END IF !INS==0

            CATOM(ISPIN,IATOM) = SUM

            IF (ISPIN.NE.1) THEN
               IF (INS.NE.0) THEN
                  WRITE (1337,FMT=9010) SUM
               ELSE
                  WRITE (1337,FMT=9050) SUM
               END IF
            ELSE                      ! (ISPIN.NE.1)
              IF (INS.NE.0) THEN
                WRITE (1337,FMT=9000) IATOM,SUM
              ELSE
                WRITE (1337,FMT=9040) IATOM,SUM
              END IF
            END IF                    ! (ISPIN.NE.1)
          END DO                      ! ISPIN = 1,NSPIN

!-----------------------------------------------------------------------
!           WRITE(6,'(2X,77(1H-))')
!       WRITE(6,*)


      CHRGNT = CATOM(1,IATOM) - ZATOM(IATOM)
      WRITE(1337,'(79(1H+))')
      WRITE (1337,FMT=9020) ITSCF,CHRGNT


      IF (NSPIN.EQ.2) THEN
        TOTSMOM = TOTSMOM+CATOM(NSPIN,IATOM) !*CONC(I1)
        WRITE (1337,FMT=9031) CATOM(NSPIN,IATOM) !TOTSMOM
        WRITE (*,FMT=9034) IATOM,CATOM(1,IATOM),CATOM(NSPIN,IATOM)!TOTSMOM
      END IF
      WRITE (1337,*)
      END DO !iatom
      WRITE (1337,FMT=9030) TOTSMOM
      WRITE (*,FMT=9030) TOTSMOM




      RETURN

 9000 FORMAT ('  Atom ',I4,' charge in wigner seitz cell =',f12.6)
 9010 FORMAT (7X,'spin moment in wigner seitz cell =',f12.6)
 9040 FORMAT ('  Atom ',I4,' charge in wigner seitz sphere =',f12.6)
 9050 FORMAT (7X,'spin moment in wigner seitz sphere =',f12.6)
 9051 FORMAT (7X,'orb. moment in wigner seitz sphere =',f12.6)
 9052 FORMAT (7X,'total magnetic moment in WS sphere =',f12.6)
 9020 FORMAT ('      ITERATION',I4,&
           ' charge neutrality in unit cell = ',f12.6)
 9030 FORMAT ('                   ',&
           ' TOTAL mag. moment in unit cell = ',f12.6)
 9031 FORMAT ('                   ',&
           '           spin magnetic moment = ',f12.6)
 9032 FORMAT ('                   ',&
           '        orbital magnetic moment = ',f12.6)
 9034 FORMAT ('Atom ',i5,':  ',&
           '        charge in WS-cell = ',f12.6,'    spin moment =',f12.6)
 9071 FORMAT ('      Site ',i3,' total charge =',f10.6)
 9072 FORMAT ('         ',' total spin moment =',f10.6)
 9073 FORMAT ('         ',' total orb. moment =',f10.6)
 9074 FORMAT ('      total magnetic moment =',f10.6)

    END SUBROUTINE RHOTOTB

END MODULE MOD_RHOTOTB
