
      SUBROUTINE AMN2010(JATOM,A,B,ALAT,GAUNTCOEFF,LMAX,NATOM,RATOM)
      use type_gauntcoeff
       implicit none

!       INTEGER LMMAXD,L3D,LM3D
!       PARAMETER (LMMAXD= (LPOTD+1)**2,L3D=2*LPOTD,LM3D= (L3D+1)**2)

      INTEGER JATOM
      REAL*8 A(NTPERD,LMMAXD,*),B(NTPERD,*)
      REAL*8 ALAT
      type(gauntcoeff_type),intent(in)      ::   gauntcoeff
      INTEGER LMMAX,LMAX,L3MAX
      INTEGER NATOM
      INTEGER RATOM(3,NATOM)
      REAL*8 EPI,FPI,PI
      REAL*8 R,R1,R2,R3,S

      REAL*8 CLEB(NCLEB,2),DFAC(0:2*LMAX,0:2*LMAX),Y(4*LMAX)

      INTEGER LX, LY,I,L
      INTEGER IATOM,JATOM
      INTEGER LM1,L1
!       INTEGER LOFLM((4*LMAX+1)**2)
      REAL*8 Y((4*LMAX+1)**2)


      FPI = 4.0D0*PI
      EPI = 8.0D0*PI
      L3MAX = 2*LMAX
      LMMAX = (LMAX+1)**2



      A=0.0D0
      B=0.0D0


      DFAC(0,0) = 1.D0
      DO 40 LX = 1,LMAX
        DFAC(LX,0) = DFAC(LX-1,0)*REAL(2*LX-1)/REAL(2*LX+1)
        DFAC(0,LX) = DFAC(LX,0)
        DO 30 LY = 1,LX
          DFAC(LX,LY) = DFAC(LX,LY-1)*REAL(2* (LX+LY)-1)/REAL(2*LY+1)
          DFAC(LY,LX) = DFAC(LX,LY)
   30   CONTINUE
   40 CONTINUE


       DO IATOM=1,NATOM

              IF (NATOM.NE.IATOM) THEN
                R1 = RATOM(1,JATOM) - RATOM(1,IATOM)
                R2 = RATOM(2,JATOM) - RATOM(2,IATOM)
                R3 = RATOM(3,JATOM) - RATOM(3,IATOM)
c
                CALL YMY(R1,R2,R3,R,Y,2*LMAX)
c
                R = R*ALAT

                DO 150 LM1 = 1,LMMAX
                  L1 = gauntcoeff.LOFLM(LM1)
                  B(IATOM,LM1) = B(IATOM,LM1) -
     +                           EPI/REAL(2*L1+1)*Y(LM1)/ (R** (L1+1))
  150           CONTINUE



               DO 170 I = 1,gauntcoeff.IEND
                  LM1 = gauntcoeff.ICLEB(I,1)
                  LM4 = gauntcoeff.ICLEB(I,2)
                  LM3 = gauntcoeff.ICLEB(I,3)
                  L1  = gauntcoeff.LOFLM(LM1)
                  L2  = gauntcoeff.LOFLM(LM4)
                    A(IATOM,LM1,LM4) = A(IATOM,LM1,LM4) +
     +                                 EPI*FPI*DFAC(L1,L2)*Y(LM3)*
     +                                 gauntcoeff.CLEB(I,1)/
     +                                 (R** (L1+L2+1))
  170           CONTINUE

       END DO !IATOM=1,NATOM
