      subroutine coreshell(ih,i,z,ncore,lcore,ecore)
C     The number of shells in NCORE
C     LCORE = 0 means s-electrons, 1 -- p, 2 -- d, 3 -- f
      double precision z(100),ecore(20,200)
      integer lcore(20,200),icore,ncore(200),ih,i,LQNTAB(15),NQNTAB(15)
      DATA NQNTAB/1,2,2,3,3,3,4,4,4,5,5,4,5,6,6/
      DATA LQNTAB/0,0,1,0,1,2,0,1,2,0,1,3,2,0,1/
C
      NCORE(I) = 0
      IF ( INT(Z(IH)).GT.2 )  NCORE(I) = 1
      IF ( INT(Z(IH)).GT.10 ) NCORE(I) = 3
      IF ( INT(Z(IH)).GT.18 ) NCORE(I) = 5
      IF ( INT(Z(IH)).GT.30 ) NCORE(I) = 6
      IF ( INT(Z(IH)).GT.36 ) NCORE(I) = 8
      IF ( INT(Z(IH)).GT.48 ) NCORE(I) = 9
      IF ( INT(Z(IH)).GT.54 ) NCORE(I) = 11
      IF ( INT(Z(IH)).GT.70 ) NCORE(I) = 12
      IF ( INT(Z(IH)).GT.80 ) NCORE(I) = 13
      IF ( INT(Z(IH)).GT.86 ) NCORE(I) = 15
      DO ICORE=1,NCORE(I)
          L = LQNTAB(ICORE)
          N = NQNTAB(ICORE)
          LCORE(ICORE,I)=L
          ECORE(ICORE,i)=-Z(IH)**2/(2.D0*N*N)
      END DO
      end
C
C
C
      FUNCTION TABCHSYM(Z)
C   ********************************************************************
C   *                                                                  *
C   *     return chemical symbol for atomic number Z                   *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER Z
      CHARACTER*2 TABCHSYM
C
C Local variables
C
      CHARACTER*2 TABLE(0:110)
C
C*** End of declarations rewritten by SPAG
C
      DATA TABLE/'Vc','H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     &     'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti',
     &     'V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se',
     &     'Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd',
     &     'Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce',
     &     'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     &     'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb',
     &     'Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U ','Np','Pu',
     &     'Am','Cm','Bk','Cf','Es','Fm','Md','No','Lw','Kt','??','??',
     &     '??','??','??','??'/
C
      IF ( Z.GE.0 .AND. Z.LE.104 ) THEN
         TABCHSYM = TABLE(Z)
      ELSE
         WRITE (*,*) ' atomic number = ',Z
         STOP ' in <TABCHSYM>'
      END IF
      END
      
