module mod_tbxccpljij

use mod_DataTypes

IMPLICIT NONE

private
public tbxccpljij

contains


SUBROUTINE tbxccpljij(iftmat,ielast,ez,wez,nspin,ncpa,naez, natyp,  &
    noq,itoq,iqat,nshell,natomimp,atomimp,ratom,  &
    nofgij,nqcalc,iqcalc,ijtabcalc, ijtabsym,ijtabsh,ish,jsh,dsymll,iprint,  &
    natypd,nsheld,lmmaxd,npol)
!   ********************************************************************
!   *                                                                  *
!   *  calculates the site-off diagonal  XC-coupling parameters  J_ij  *
!   *  according to  Lichtenstein et al. JMMM 67, 65 (1987)            *
!   *                                                                  *
!   *  adopted for TB-KKR code from Munich SPR-KKR package Sep 2004    *
!   *                                                                  *
!   *  for mpi-parallel version: moved energy loop from main1b into    *
!   *   here.                                B. Zimmermann, Dez 2015   *
!   *                                                                  *
!   ********************************************************************

#ifdef CPP_MPI
use mpi
#ENDIF
use mod_types, only: t_tgmat, t_mpi_c_grid, t_cpa
use mod_mympi, only: myrank, master
use mod_version_info
use mod_md5sums

IMPLICIT NONE
!.
!. Parameters
DOUBLE COMPLEX CONE,CZERO
PARAMETER ( CONE  = (1D0,0D0) )
PARAMETER ( CZERO = (0D0,0D0) )
INTEGER NIJMAX
PARAMETER ( NIJMAX = 200 )
!.
!. Scalar arguments
INTEGER IELAST,NSPIN,IFTMAT,IPRINT,LMMAXD,NAEZ,NATOMIMP, &
        NATYP,NATYPD,NCPA,NOFGIJ,NQCALC,NSHELD,NPOL
DOUBLE COMPLEX EZ(*), WEZ(*)
!.
!. Array arguments
INTEGER ATOMIMP(*),IJTABCALC(*),IJTABSH(*),IJTABSYM(*), &
        IQAT(*),IQCALC(*),ISH(NSHELD,*),ITOQ(NATYPD,*), &
        JSH(NSHELD,*),NOQ(*),NSHELL(0:NSHELD)
DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,*)
DOUBLE PRECISION RATOM(3,*)
!.
!. Local scalars
INTEGER I1,IA,IFGMAT,IFMCPA,IQ,IREC,ISPIN,ISYM,IT,J1,JA, &
        JQ,JT,L1,LM1,LM2,LSTR,NS,NSEFF,NSHCALC,NSMAX,NTCALC, &
        IE, ie_end, ie_num
#ifdef CPP_MPI
INTEGER ie_start, IERR
#endif
DOUBLE COMPLEX CSUM, WGTEMP
#ifdef CPP_MPI
DOUBLE COMPLEX XINTEGDTMP
#endif
CHARACTER*8 FMT1
CHARACTER*22 FMT2
CHARACTER*8 JFBAS
CHARACTER*13 JFNAM
CHARACTER*26 JFNAM2
CHARACTER*80 STRBAR,STRTMP
!.
!. Local arrays
INTEGER NIJCALC(:),KIJSH(:,:),JIJDONE(:,:,:)
DOUBLE COMPLEX JXCIJINT(:,:,:)
#ifndef CPP_MPI
DOUBLE COMPLEX, ALLOCATABLE :: XINTEGD(:,:,:)
#else
DOUBLE COMPLEX, ALLOCATABLE :: csum_store(:,:,:,:), &
                               csum_store2(:,:,:,:)
#endif
ALLOCATABLE NIJCALC,JIJDONE,KIJSH,JXCIJINT
DOUBLE COMPLEX DELTSST(LMMAXD,LMMAXD,NATYP), &
               DMATTS(LMMAXD,LMMAXD,NATYP,NSPIN), &
               DTILTS(LMMAXD,LMMAXD,NATYP,NSPIN), &
               GMIJ(LMMAXD,LMMAXD), &
               GMJI(LMMAXD,LMMAXD),GS(LMMAXD,LMMAXD,NSPIN), &
               TSST(LMMAXD,LMMAXD,NATYP,2),W1(LMMAXD,LMMAXD), &
               W2(LMMAXD,LMMAXD),W3(LMMAXD,LMMAXD)
DOUBLE PRECISION RSH(NSHELD),PI
INTEGER JTAUX(NATYP)
!..
!.. Intrinsic Functions ..
INTRINSIC MAX,SQRT
!..
!.. External Subroutines ..
EXTERNAL CINIT,CMATMUL,INITABJIJ,ZCOPY
LOGICAL TEST
EXTERNAL TEST
!..
!.. Save statement
SAVE IFGMAT,IFMCPA,JIJDONE,JXCIJINT,KIJSH,NIJCALC
#ifndef CPP_MPI
SAVE XINTEGD
#endif
SAVE NSHCALC,NSMAX
!..
!.. Data statement
DATA JFBAS/'Jij.atom'/
DATA PI /3.14159265358979312D0/
!     ..
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! ==>                   -- initialisation step --                    <==
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
!          open(22,STATUS='unknown',FILE='integrand1.dat',
!     &                         FORM='formatted')
!           open(44,STATUS='unknown',FILE='integrand2.dat',
!     &                         FORM='formatted')
!      write(*,*) 'test brahim 2'

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
IF(iprint>0)  WRITE(1337,99000)
#ifdef CPP_MPI
! this part of the code does not take much time, thus only the energy
! loop is parallel and not the second level over atoms, only the master
! in every energy loop does the calculation of the Jijs
IF(t_mpi_c_grid%myrank_ie==0) THEN
  
  ie_end = t_mpi_c_grid%ntot_pt2(t_mpi_c_grid%myrank_at)
#ELSE
  ie_end = ielast
#ENDIF
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

!         open(22,STATUS='unknown',FILE='integrand.dat',
!     &                         FORM='formatted')
ifgmat = iftmat + 1
ifmcpa = ifgmat + 1
nsmax = MAX(naez,natyp)
nshcalc = nshell(0) - nsmax
!ccc         IF ( NSHLOC.LT.NSHELD) STOP ' NSHLOC'
!ccc         IF ( NTLOC.LT.NATYP) STOP ' NTLOC'

allocate (kijsh(nijmax,nshell(0)),stat=lm1)
IF ( lm1 /= 0 ) THEN
  WRITE(6,99001) 'KIJSH'
  STOP
endif
allocate (nijcalc(nshell(0)),jijdone(natyp,natyp,nshell(0)), stat=lm1)
IF ( lm1 /= 0 ) THEN
  WRITE(6,99001) 'JIJDONE/NIJCALC'
  STOP
endif

!        The array CSUM_STORE could become quite large -> store to MPI-IO-file in future????
!        Only allocate it for MPI-usage. If MPI is not used, two smaller array XINTEGD arrays is needed.
#ifdef CPP_MPI
allocate (csum_store(natyp,natyp,nshell(0),ielast),stat=lm1)
IF ( lm1 /= 0 ) THEN
  WRITE(6,99001) 'csum_store'
  STOP
endif
#ELSE
!        ALLOCATE (XINTEGD1(NATYP,NATYP,IELAST),STAT=LM1)
allocate (xintegd(natyp,natyp,nshell(0)),stat=lm1)
IF ( lm1 /= 0 ) THEN
  WRITE(6,99001) 'XINTEGD'
  STOP
endif
#ENDIF
allocate (jxcijint(natyp,natyp,nshell(0)),stat=lm1)
IF ( lm1 /= 0 ) THEN
  WRITE(6,99001) 'JXCIJINT'
  STOP
endif
!..................
#ifdef CPP_MPI
DO ie=1,ielast
  DO ns = 1,nshell(0)
    DO jt = 1,natyp
      DO it = 1,natyp
        csum_store(it,jt,ns,ie)= czero
      END DO
    END DO
  END DO
END DO!IE
#ELSE
DO ns = 1,nshell(0)
  DO jt = 1,natyp
    DO it = 1,natyp
      xintegd(it,jt,ns)= czero
!                  XINTEGD1(IT,JT,IE)= CZERO
    END DO
  END DO
END DO
#ENDIF
DO ns = 1,nshell(0)
  nijcalc(ns) = 0
  DO jt = 1,natyp
    DO it = 1,natyp
      jijdone(it,jt,ns) = 0
      jxcijint(it,jt,ns) = czero
    END DO
  END DO
END DO

CALL initabjij(iprint,naez,natyp,natomimp,nofgij,nqcalc,  &
    nsmax,nshell,iqcalc,atomimp,ish,jsh,ijtabcalc,  &
    ijtabsh,ijtabsym,nijcalc(1),kijsh(1,1), nijmax,nshell(0),nsheld)
! --------------------------------------------------------------------
DO ns = nsmax + 1,nshell(0)
  nseff = ns - nsmax
  DO i1 = 1,nijcalc(ns)
    ia = ish(ns,kijsh(i1,ns))
    ja = jsh(ns,kijsh(i1,ns))
    lm1 = (ia-1)*natomimp + ja
    iq = atomimp(ia)
    jq = atomimp(ja)
    
    DO l1 = 1,noq(jq)
      jt = itoq(l1,jq)
      DO j1 = 1,noq(iq)
        it = itoq(j1,iq)
        jijdone(it,jt,nseff) = 1
      END DO
    END DO
  END DO
END DO
! ----------------------------------------------------------------------

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
IF ( iprint > 0 ) THEN
  DO it = 1,natyp
    DO jt = 1,natyp
      jtaux(jt) = 0
      DO ns = nsmax + 1,nshell(0)
        nseff = ns - nsmax
        IF ( jijdone(it,jt,nseff) /= 0 ) jtaux(jt) = jtaux(jt) + 1
      END DO
    END DO
!ccc               WRITE (6,99012) IT,(JTAUX(JT),JT=1,MIN(25,NATYP))
  END DO
  WRITE (1337,99013)
endif
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! ==>                   INITIALISATION END                           <==
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII


! ==>  read in the single site matrices (TSST) in the LOCAL frame

! ==>  set up the Delta t_i matrices DELTSST = TSST(UP) - TSST(DOWN)

! ==>  read in projection matrices DMAT and DTIL used to get
!                     ij            ij    _
!                    G     =  D  * G    * D
!                     ab       a    CPA    b
!     with a/b the atom of type a/b sitting on site i/j
!     for an atom having occupancy 1, DMAT/DTIL = unit matrix


! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES
#ifdef CPP_MPI
ie_start = t_mpi_c_grid%ioff_pt2(t_mpi_c_grid%myrank_at)
ie_end   = t_mpi_c_grid%ntot_pt2(t_mpi_c_grid%myrank_at)

DO ie_num=1,ie_end
  ie = ie_start+ie_num
#ELSE
  DO ie = 1,ielast
    ie_num = ie
    ie_end = ielast
#ENDIF
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS SPIN
  DO ispin = 1,nspin
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NATYP
    DO it = 1,natyp
!            write(*,*) 'test brahim 3'
      IF (t_tgmat%tmat_to_file) THEN
        irec = ie + ielast*(ispin-1) + ielast*2*(it-1)
        READ (iftmat,REC=irec) w1
      ELSE
        ie_end = t_mpi_c_grid%ntot2
        irec = ie_num + ie_end*(ispin-1) + ie_end*2*(it-1)
        w1(:,:) = t_tgmat%tmat(:,:,irec)
      endif
      DO j1 = 1,lmmaxd
        CALL zcopy(lmmaxd,w1(1,j1),1,tsst(1,j1,it,ispin),1)
      END DO
! ----------------------------------------------------------------------
      IF ( ispin == 2 ) THEN
        DO j1 = 1,lmmaxd
          DO i1 = 1,lmmaxd
            deltsst(i1,j1,it) = tsst(i1,j1,it,2) - tsst(i1,j1,it,1)
          END DO
        END DO
      endif
! ----------------------------------------------------------------------
      IF ( ncpa == 0 ) THEN
        CALL cinit(lmmaxd*lmmaxd,dmatts(1,1,it,ispin))
        DO i1 = 1,lmmaxd
          dmatts(i1,i1,it,ispin) = cone
        END DO
        CALL zcopy(lmmaxd*lmmaxd,dmatts(1,1,it,ispin),1,  &
            dtilts(1,1,it,ispin),1)
      ELSE!NCPA.EQ.0
        IF(t_cpa%dmatproj_to_file)THEN
          WRITE(*,*) 'test read proj'
          irec = ie + ielast*(ispin-1) + ielast*2*(it-1)
          READ (ifmcpa,REC=irec) w1,w2
          DO j1 = 1,lmmaxd
            DO i1 = 1,lmmaxd
              dmatts(i1,j1,it,ispin) = w1(i1,j1)
              dtilts(i1,j1,it,ispin) = w2(i1,j1)
            END DO
          END DO
        ELSE!t_cpa%dmatproj_to_file
          irec = ie_num + ie_end*(ispin-1)
          dmatts(:,:,it,ispin) = t_cpa%dmatts(:,:,it,irec)
          dtilts(:,:,it,ispin) = t_cpa%dtilts(:,:,it,irec)
          
        endif!t_cpa%dmatproj_to_file
      endif!NCPA.EQ.0
! ----------------------------------------------------------------------
    END DO
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  END DO
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  IF ( iprint > 1 ) THEN
    WRITE(1337,'(8X,60(1H-),/,10X,  &
        "Single-site and projection matrices read in", " for IT=1,",I3,/)') natyp
    DO ispin = 1,nspin
      WRITE(1337,'(8X,60(1H+),/,30X, " ISPIN = ",I1,/,8X,60(1H+),/)') ispin
      DO it = 1,natyp
        WRITE(1337,'(12X," IE = ",I2," IT =",I3)') ie,it
        CALL cmatstr(' T MAT ',7,tsst(1,1,it,ispin),  &
            lmmaxd,lmmaxd,0,0,0,1D-8,6)
        CALL cmatstr(' D MAT ',7,dmatts(1,1,it,ispin),  &
            lmmaxd,lmmaxd,0,0,0,1D-8,6)
        CALL cmatstr(' D~ MAT',7,dtilts(1,1,it,ispin),  &
            lmmaxd,lmmaxd,0,0,0,1D-8,6)
        IF ( it /= natyp) WRITE(1337,'(8X,60(1H-),/)')
      END DO
      WRITE(1337,'(8X,60(1H+),/)')
    END DO
    WRITE(1337,'(8X,60(1H-),/,10X,  &
        "Delta_t = t(it,DN) - t(it,UP) matrices for IT=1,", I3,/)') natyp
    DO it = 1,natyp
      WRITE(1337,'(12X," IE = ",I2," IT =",I3)') ie,it
      CALL cmatstr(' DEL T ',7,deltsst(1,1,it), lmmaxd,lmmaxd,0,0,0,1D-8,6)
      IF ( it /= natyp) WRITE(1337,'(8X,60(1H-),/)')
    END DO
    WRITE(1337,'(8X,60(1H-),/)')
  endif
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  
! ***************************************************** loop over shells
  DO ns = nsmax + 1,nshell(0)
    
! ======================================================================
    
! ==>  get the off-diagonal Green function matrix Gij(UP) and Gji(DOWN)
!      using the appropiate rotation (similar to ROTGLL routine)
!      step one: read in GS for the representative shell
! ======================================================================
    DO ispin = 1,nspin
      IF (t_tgmat%gmat_to_file) THEN
        irec = ie + ielast*(ispin-1) + ielast*nspin*(ns-1)
        READ (ifgmat,REC=irec) w1
      ELSE
        ie_end = t_mpi_c_grid%ntot2
        irec = ie_num + ie_end*(ispin-1) + ie_end*nspin*(ns-1)
        w1(:,:) = t_tgmat%gmat(:,:,irec)
      endif
      CALL zcopy(lmmaxd*lmmaxd,w1,1,gs(1,1,ispin),1)
    END DO
! ======================================================================
!      step two: scan all I,J pairs needed out of this shell
!                transform with the appropriate symmetry to get Gij
    
! ====================================================== loop over pairs
    nseff = ns - nsmax
    DO l1 = 1,nijcalc(ns)
      
      ia = ish(ns,kijsh(l1,ns))
      ja = jsh(ns,kijsh(l1,ns))
      
      lm1 = (ia-1)*natomimp + ja
      isym = ijtabsym(lm1)
! ----------------------------------------------------------------------
      DO ispin = 1,nspin
        CALL zgemm('C','N',lmmaxd,lmmaxd,lmmaxd,  &
            cone,dsymll(1,1,isym),lmmaxd, gs(1,1,ispin),lmmaxd,czero,w2,lmmaxd)
        
        CALL zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,  &
            cone,w2,lmmaxd,dsymll(1,1,isym),lmmaxd, czero,w1,lmmaxd)
        
        IF ( ispin == 1 ) THEN
          CALL zcopy(lmmaxd*lmmaxd,w1,1,gmij,1)
        ELSE
          DO lm2 = 1,lmmaxd
            DO lm1 = 1,lmmaxd
              
! -> use Gji = Gij^T
              
              gmji(lm1,lm2) = w1(lm2,lm1)
            END DO
          END DO
        endif
      END DO
! ----------------------------------------------------------------------
      
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      IF ( iprint > 2 ) THEN
        WRITE(1337,'(8X,60(1H-),/,10X,  &
            " G_ij(DN) and G_ji(UP) matrices I =",I3," J =",I3,  &
            " IE =",I3,/,8X,60(1H-))') ia,ja,ie
        CALL cmatstr(' Gij DN',7,gmij,lmmaxd,lmmaxd,0,0,0,1D-8,6)
        CALL cmatstr(' Gji UP',7,gmji,lmmaxd,lmmaxd,0,0,0,1D-8,6)
      endif
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      
      
! ----------------------------------------------------------------------
      
! ==> calculate the exchange coupling constant J_ij via Eq. (19)
!     modified for G instead of tau:
!          J_ij ~ Trace [ (t_i(D)-t_i(U)) * Gij(U)
!                       * (t_j(D)-t_j(U)) * Gji(D)]
!       in case of alloy system: perform projection on atom types
      
! ----------------------------------------------------------------------
      iq = atomimp(ia)
      jq = atomimp(ja)
! -------------------------------------------------- loop over occupants
      DO j1 = 1,noq(jq)
        jt = itoq(j1,jq)
        DO i1 = 1,noq(iq)
          it = itoq(i1,iq)
          
! --> Gjq,iq is projected on jt,it ==> Gjt,it
          
          CALL cmatmul(lmmaxd,lmmaxd,gmji,dtilts(1,1,it,2),w2)
          CALL cmatmul(lmmaxd,lmmaxd,dmatts(1,1,jt,2),w2,w1)
          
! --> Delta_j * Gjt,it
          
          CALL cmatmul(lmmaxd,lmmaxd,deltsst(1,1,jt),w1,w2)
          
! --> Giq,jq is projected on it,jt ==> Git,jt * Delta_j * Gjt,it
          
          CALL cmatmul(lmmaxd,lmmaxd,gmij,dtilts(1,1,jt,1),w3)
          CALL cmatmul(lmmaxd,lmmaxd,dmatts(1,1,it,1),w3,w1)
          
! --> Delta_i * Git,jt
          
          CALL cmatmul(lmmaxd,lmmaxd,deltsst(1,1,it),w1,w3)
          
! --> Delta_i * Git,jt * Delta_j * Gjt,it
          
          CALL cmatmul(lmmaxd,lmmaxd,w3,w2,w1)
          
          csum = czero
          DO lm1 = 1,lmmaxd
            csum = csum + w1(lm1,lm1)
          END DO
          
#ifdef CPP_MPI
!BZ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BZ! to ensure the correct order of energy points, the csum=values !!
!BZ! are stored for the MPI version, and the evaluation of JXCIJINT!!
!BZ! and XINTEGD is performed later                                !!
!BZ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          csum_store(it,jt,nseff,ie)=csum
#ELSE
          
          jxcijint(it,jt,nseff) = jxcijint(it,jt,nseff)  &
              - wez(ie)*csum/DBLE(nspin)
          
          xintegd(it,jt,nseff)= csum/(pi*4.d0)
          
!                  -------> perform substraction instead of addition
!                           because WGTE ~ -1/pi (WGTE = WEZ(IE)/NSPIN)
! Write out energy-resorved integrand and integral
! Phivos Mavropoulos 24.10.2012
          IF(npol==0 .OR. test('Jijenerg'))THEN
            fmt2 = '(A,I5.5,A,I5.5,A,I5.5)'
            WRITE (jfnam2,fmt2) 'Jij_enrg.',it,'.',jt,'.',ns
            IF (ie == 1) THEN
              OPEN (499,FILE=jfnam2,STATUS='UNKNOWN')
              CALL version_print_header(499,  &
                  '; '//md5sum_potential//'; '//md5sum_shapefun)
              WRITE(499,FMT='(A)') '# Energy Re,Im ; j(E) Re,Im; J(E) Re,Im '
              WRITE(499,FMT='(3(A,I5))') '# IT=',it,' JT=',jt,' SHELL=',ns
              WRITE(499,FMT='(A,I6)') '#ENERGIES: ',ielast
            ELSE
              OPEN (499,FILE=jfnam2,STATUS='OLD', POSITION='APPEND')
            endif
            WRITE (499,FMT='(6E12.4)') ez(ie),xintegd(it,jt,nseff),  &
                jxcijint(it,jt,nseff)/4.d0
            CLOSE(499)
          endif!(npol==0 .or. test('Jijenerg'))
#ENDIF
        
      END DO !I1
    END DO   !J1, loop over occupants
! ----------------------------------------------------------------------
  END DO   !L1 = 1,NIJCALC(NS)
  
! ======================================================================
END DO  !loop over shells
! **********************************************************************
!       write(22,*)DIMAG(XINTEGD(1,1,1)),DIMAG(JXCIJINT(1,1,1))
!       write(44,*)DIMAG(XINTEGD(2,2,1)),DIMAG(XINTEGD(2,3,1))

END DO!IE
! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC MPI Communication
#ifdef CPP_MPI
!allocate
allocate(csum_store2(natyp,natyp,nshell(0),ielast),stat=lm1)
IF ( lm1 /= 0 ) THEN
  WRITE(6,99001) 'csum_store2'
  STOP
endif

!initialize with zero
DO ie=1,ielast
  DO ns = 1,nshell(0)
    DO jt = 1,natyp
      DO it = 1,natyp
        csum_store2(it,jt,ns,ie)= czero
      END DO
    END DO
  END DO
END DO!IE

!perform communication and collect resutls in csum_store2
lm1 = natyp*natyp*nshell(0)*ielast
CALL mpi_reduce(csum_store,csum_store2,lm1,mpi_double_complex,  &
    mpi_sum,master,t_mpi_c_grid%mympi_comm_at,ierr)
IF(ierr/=mpi_success)THEN
  WRITE(*,*) 'Problem in MPI_REDUCE(csum_store)'
  STOP 'TBXCCPLJIJ'
endif ! IERR/=MPI_SUCCESS

!copy array back to original one
IF(myrank==master)THEN
  DO ie=1,ielast
    DO ns = 1,nshell(0)
      DO jt = 1,natyp
        DO it = 1,natyp
          csum_store(it,jt,ns,ie)= csum_store2(it,jt,ns,ie)
        END DO
      END DO
    END DO
  END DO!IE
  
  deallocate(csum_store2)
endif!myrank==master



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!~~~~~~~ NEW WRITEOUT: integrand and energy-resolved integral ~~~~~~~~~~
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
IF(myrank==master)THEN
  DO ns = nsmax + 1,nshell(0)
    nseff = ns - nsmax
    DO l1 = 1,nijcalc(ns)
      ia = ish(ns,kijsh(l1,ns))
      ja = jsh(ns,kijsh(l1,ns))
      iq = atomimp(ia)
      jq = atomimp(ja)
      DO j1 = 1,noq(jq)
        jt = itoq(j1,jq)
        DO i1 = 1,noq(iq)
          it = itoq(i1,iq)
!                  -------> perform substraction instead of addition
!                           because WGTE ~ -1/pi (WGTE = WEZ(IE)/NSPIN)
! Write out energy-resorved integrand and integral
! Phivos Mavropoulos 24.10.2012
          IF(npol==0 .OR. test('Jijenerg'))THEN
            fmt2 = '(A,I5.5,A,I5.5,A,I5.5)'
            WRITE (jfnam2,fmt2) 'Jij_enrg.',it,'.',jt,'.',ns
            OPEN (499,FILE=jfnam2,STATUS='UNKNOWN')
            CALL version_print_header(499,  &
                '; '//md5sum_potential//'; '//md5sum_shapefun)
            WRITE(499,FMT='(A)') '# Energy Re,Im ; j(E) Re,Im; J(E) Re,Im '
            WRITE(499,FMT='(3(A,I5))') '# IT=',it,' JT=',jt,' SHELL=',ns
            WRITE(499,FMT='(A,I6)') '#ENERGIES: ',ielast
          endif!(NPOL==0 .OR. TEST('Jijenerg'))then
          
          DO ie=1,ielast
            jxcijint(it,jt,nseff) = jxcijint(it,jt,nseff)  &
                - wez(ie)*csum_store(it,jt,nseff,ie)/DBLE(nspin)
            xintegdtmp= csum_store(it,jt,nseff,ie)/(pi*4.d0)
            IF(npol==0 .OR. test('Jijenerg'))THEN
              WRITE (499,FMT='(6E12.4)')  &
                  ez(ie),xintegdtmp,jxcijint(it,jt,nseff)/4.d0
            endif!(NPOL==0 .OR. TEST('Jijenerg'))then
            
          END DO!IE
          
          IF(npol==0 .OR. test('Jijenerg')) CLOSE(499)
          
        END DO !I1
      END DO   !J1, loop over occupants
    END DO   !L1 = 1,NIJCALC(NS)
  END DO  !NS, loop over shells
endif!myrank==master
#ENDIF

IF(myrank==master)THEN
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  wgtemp = cone/4D0
!     -------> here factor 1/pi omitted since it is included in WGTE
  DO ns = 1,nshcalc
    DO jt = 1,natyp
      DO it = 1,natyp
        jxcijint(it,jt,ns) = wgtemp*jxcijint(it,jt,ns)
      END DO
    END DO
    nseff = ns + nsmax
    rsh(ns) = 0D0
    DO i1 = 1,3
      rsh(ns) = rsh(ns) + ratom(i1,nseff)*ratom(i1,nseff)
    END DO
    rsh(ns) = SQRT(rsh(ns))
  END DO
  
  lm2 = 0
  ntcalc = 0
  DO it = 1,natyp
    l1 = 0
    DO ns = 1,nshcalc
      lm1 = 0
      DO jt = 1,natyp
        IF ( jijdone(it,jt,ns) /= 0 ) lm1 = lm1 + 1
      END DO
      lm2 = MAX(lm1,lm2)
      l1 = MAX(l1,lm1)
    END DO
    IF ( l1 > 0 ) THEN
      ntcalc = ntcalc + 1
      jtaux(ntcalc) = it
    endif
  END DO
  WRITE (strbar,'(19(1H-))')
  lstr = 19
  DO i1 = 1,lm2
    WRITE(strtmp,'(A,15(1H-))') strbar(1:lstr)
    lstr = lstr+15
    strbar(1:lstr)=strtmp(1:lstr)
  END DO
  
  WRITE(1337,99002) strbar(1:lstr),strbar(1:lstr)
  DO i1 = 1,ntcalc
    it = jtaux(i1)
    WRITE (1337,99003, advance='no') it,iqat(it)
    l1 = 0
    DO ns = 1,nshcalc
      lm1 = 0
      DO jt = 1,natyp
        IF ( jijdone(it,jt,ns) /= 0 ) lm1 = lm1 + 1
      END DO
      IF ( lm1 /= 0 ) THEN
        lm2 = 0
        IF ( l1 == 0 ) THEN
          WRITE(1337,99005, advance='no') rsh(ns)
          l1 = 1
        ELSE
          WRITE (1337,99006, advance='no') rsh(ns)
        endif
        DO jt = 1,natyp
          IF ( jijdone(it,jt,ns) /= 0 ) THEN
            lm2 = lm2 + 1
            IF ( lm2 == 1 ) THEN
              WRITE (1337,99007, advance='no') iqat(jt), aimag(jxcijint(it,jt,ns))*1D3,jt
              WRITE(1337,*) ns+nsmax,' shell'
            ELSE
              WRITE (1337,99008, advance='no') aimag(jxcijint(it,jt,ns))*1D3,jt
              WRITE(1337,*) ns+nsmax,' shell'
            endif
            IF ( lm2 == lm1 ) WRITE (1337,*)
          endif
        END DO
      endif
    END DO
    WRITE (1337,99004) strbar(1:lstr)
  END DO !I1 = 1,NTCALC
  WRITE (1337,*)
! ----------------------------------------------------------------------
! --> prepare output files
  
  DO i1 = 1,ntcalc
    it = jtaux(i1)
    DO ns = 1,nshcalc
      l1 = 1
      DO jt = 1,natyp
        IF ( jijdone(it,jt,ns) /= 0 ) THEN
          jijdone(it,l1,ns) = jt
          jxcijint(it,l1,ns) = jxcijint(it,jt,ns)
          l1 = l1 + 1
        endif
      END DO
      DO jt = l1,natyp
        jijdone(it,jt,ns) = 0
      END DO
    END DO
  END DO
  
  DO i1 = 1,ntcalc
    it = jtaux(i1)
!         FMT1 = '(A,I1)'
!         IF ( IT.GE.10 ) THEN
!            FMT1 = '(A,I2)'
!            IF ( IT.GE.100 ) FMT1 = '(A,I3)'
!         endif
    fmt1 = '(A,I5.5)'
    WRITE (jfnam,fmt1) jfbas,it
    
!         write(JFNAME,FMT1)JFINTEG, IT
    OPEN (49,FILE=jfnam)
    CALL version_print_header(49,  &
        '; '//md5sum_potential//'; '//md5sum_shapefun)
    WRITE (49,99009) it,iqat(it)
!         open(22,FILE=JFNAME)
!          write(22,99014)IT,IQAT(IT)
    DO l1 = 1,natyp
      lm1 = 0
      DO ns = 1,nshcalc
        IF ( jijdone(it,l1,ns) /= 0 ) THEN
          lm1 = lm1 + 1
          WRITE (49,99010) rsh(ns), aimag(jxcijint(it,l1,ns)),  &
              jijdone(it,l1,ns),ns+nsmax ! fivos added NS+NSMAX
!              do IE=1,IELAST
!                XINTEGD(IT,JT,IE)= XINTEGD(IT,JT,IE)*(1.D0/(PI*4.D0))
          
!                write(22,99015)XINTEGD(IT,JT,IE),JIJDONE(IT,L1,NS)
!             end do
        endif
      END DO
!            end do
      IF ( (lm1 /= 0) .AND. (l1 /= natyp) ) WRITE (49,*) '&'
!           IF ( (LM1.NE.0) .AND. (L1.NE.NATYP) ) WRITE (22,*) '&'
    END DO !l1=1,natyp
    CLOSE (49)
!         Close(22)
  END DO !i1=1,ntcalc
  WRITE (1337,99011)
! ----------------------------------------------------------------------
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
endif!myrank==master
#ifdef CPP_MPI
endif ! t_mpi_c_grid%myrank_ie==0
#ENDIF

deallocate (kijsh,nijcalc,jijdone,jxcijint,stat=lm1)
99000 FORMAT(79(1H=),/,10X, "TBXCCPLJIJ : Off-diagonal exchange coupling",  &
    " constants J_ij",/,79(1H=),/)
99001 FORMAT(6X,"ERROR: Cannot allocate array(s) :",a,/)
99002 FORMAT(4X,a,4(1H-),/,5X," IT/IQ ",3X,"R_IQ,JQ",2X,"JQ ",  &
    " (J_IT,JT  JT)",/,15X," [ ALAT ] ",4X,"[ mRy ]",/, 4X,a,4(1H-))
99003 FORMAT(5X,i3,1X,i3)
99004 FORMAT(4X,a,4(1H-))
99005 FORMAT(f10.6)
99006 FORMAT(12X,f10.6)
99007 FORMAT(i4,f12.8,i3)
99008 FORMAT(f12.8,i3)
99009 FORMAT("# off-diagonal exchange coupling constants ",/,  &
    "# for atom IT = ",i3," on site IQ = ",i3,/,  &
    "# R_IQ,JQ      J_IT,JT       JT",/, "# ( ALAT )       ( Ry )",/,"#      ")
99010 FORMAT(f12.8,2X,e15.8,2X,i3,i5)
99011 FORMAT(6X,"Output written into the files Jij.atomX",/,  &
    6X,"  X = atom index in the input file")
99012       FORMAT (10X,i3,3X,25(i3))
99013       FORMAT (/,8X,60('-'),/)

99014       FORMAT("# Integrand  ",/,  &
    "# for atom IT = ",i3," on site IQ = ",i3,/, "#      J_IT,JT       JT",/,  &
    "#        ( Ry )",/,"#      ")
99015 FORMAT(f12.8)
END SUBROUTINE tbxccpljij

END module mod_tbxccpljij
