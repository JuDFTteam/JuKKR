      SUBROUTINE CSIMP3(CF,CFOUT,IRMIN,IRM,DRDI)
c ************************************************************************
c     this subroutine does an integration up to rcut of an
c     complex function cf with an extended 3-point-simpson, r is real
c
c                                  rmax
c                      cfout(ll) = { cf(ll,r') dr'
c                                  rmin
c
c     attention : input cf is destroyed !
c     modified by Long, Juelich 1/2013
c
c-----------------------------------------------------------------------
      INTEGER IRMIN,IRM,I
      DOUBLE COMPLEX CF(IRM),CFINT(IRM),CFOUT
      DOUBLE PRECISION DRDI(*)
      DOUBLE PRECISION A1,A2

      A1 = 1.0D0/3.0D0
      A2 = 4.0D0/3.0D0
c
c
        DO I = IRMIN,IRM
          CF(I) = CF(I)*DRDI(I)
        ENDDO
c
          CFINT(IRMIN) = 0d0
          CFINT(IRMIN+1)=(CF(IRMIN+3)-5D0*CF(IRMIN+2)+
     +                    19D0*CF(IRMIN+1)+9D0*CF(IRMIN))/24D0

c---> calculate with an extended 3-point-simpson
c
        DO I=IRMIN+2,IRM
         CFINT(I) = CFINT(I-2)+A1*CF(I-2)+
     +                   A2*CF(I-1)+A1*CF(I)
        ENDDO
         CFOUT=CFINT(IRM)
c
      END
