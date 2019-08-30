      subroutine elementdatabase(elem_z,elem_name)
      implicit none
!#@# KKRtags: VORONOI visualization core-electrons
!#@# KKRmerge: can be simplified strongly
! -----------------------------------------------------------
!  This is an element batabase
!  It takes as input the atomic number ELEM_Z
!  and returns the element name
! ------------------------------------------------------------
! ...Input 
      REAL*8      ELEM_Z    
! ...Output      
      CHARACTER*3 ELEM_NAME
! ...Local variables
      CHARACTER*3 DATA1(0:113)
      INTEGER IN_Z
      DATA DATA1/'H1 ',
     &     'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',  
     &     'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ', 
     &     'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &     'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &     'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &     'Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &     'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &     'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &     'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &     'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &     'Md ','No ','Lr ','Rf ','Db ','Sg ','Bh ','Hs ','Mt ','Uun',
     &     'Uuu','Uub','NoE'/
! --------------------------------------------------------------------
! Get the integer part of the atomic number
      IN_Z = ELEM_Z
      IF (ABS(FLOAT(IN_Z)-ELEM_Z).GT.1.D-6) IN_Z = 113
      IF ((IN_Z.LT.0).OR.(IN_Z.GT.113)) THEN             
         WRITE (6,*) ' Ooops, This atomic number does not exist!'
         WRITE (6,*) ' ERROR in elementbatabase'
         STOP
      END IF
      ELEM_NAME = DATA1(IN_Z)
 1000 format(2A8,I3,1X,3F8.4,1X,A24)
      RETURN
      END

