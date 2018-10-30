!-------------------------------------------------------------------------------
!> Summary: Outward integration of multiple functions with ext. 3-point Simpson
!> Author: B. Drittler
!> Date: 1989
!> 
!> This subroutine does an outwards integration of \(l_{lmax}\)
!> functions f with an extended 3-point-Simpson:
!>
!> \begin{equation}
!> f_{int}\left(ll,i\right) = \int_{ir}^{ir_{max}}f\left(ll,i'\right)di'           
!> \end{equation}
!>
!> The starting value for this integration at \(ir_{min}+1\) is determined by
!> a 4 point lagrangian integration, coefficients given by
!> m. abramowitz and i.a. stegun, handbook of mathematical functions,
!> nbs applied mathematics series 55 (1968)
!-------------------------------------------------------------------------------
!> @warning In case of radial integration :
!> the weights drdi have to be multiplied before calling this
!> subroutine.
!> @endwarning
!>
!> B. Drittler Mar. 1989
!>
!> Modified for functions with kinks - at each kink the integration
!> is restarted
!>
!> @warning It is supposed that \(ir_{min} + 3\) is less than imt! 
!> @endwarning
!>
!> B. Drittler july 1989
!-------------------------------------------------------------------------------
      MODULE mod_csout
      CONTAINS
!-------------------------------------------------------------------------------
!> Summary: Outward integration of multiple functions with ext. 3-point Simpson
!> Author: B. Drittler
!> Date: 1989
!> Category: KKRimp, radial-grid, numerical-tools
!> Deprecated: False 
!> 
!> This subroutine does an outwards integration of $l_{lmax}$
!> functions f with an extended 3-point-Simpson:
!>
!> \begin{equation}
!> f_{int}\left(ll,i\right) = \int_{ir}^{ir_{max}}f\left(ll,i'\right)di'           
!> \end{equation}
!>
!> The starting value for this integration at irmin+1 is determined by
!> a 4 point lagrangian integration, coefficients given by
!> m. abramowitz and i.a. stegun, handbook of mathematical functions,
!> nbs applied mathematics series 55 (1968)
!-------------------------------------------------------------------------------
      SUBROUTINE CSOUT(F,FINT,LMMSQD,IRMIND,IRMD,IPAN,IRCUT)
C     .. Parameters ..
      DOUBLE PRECISION A1,A2
      PARAMETER (A1=1.D0/3.D0,A2=4.D0/3.D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER IPAN,IRMD,IRMIND,LMMSQD
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX F(LMMSQD,IRMIND:IRMD),FINT(LMMSQD,IRMIND:IRMD)
      INTEGER IRCUT(0:IPAN)
C     ..
C     .. Local Scalars ..
      INTEGER I,IEN,IP,IST,LL
C     ..
c
c---> loop over kinks
c

      DO 50 IP = 1,IPAN
        IEN = IRCUT(IP)
        IST = IRCUT(IP-1) + 1
c
        IF (IP.EQ.1) THEN
          IST = IRMIND
          DO 10 LL = 1,LMMSQD
            FINT(LL,IST) = 0.0D0
c---> integrate fint(ist+1) with a 4 point lagrangian
            FINT(LL,IST+1) = (F(LL,IST+3)-5.0D0*F(LL,IST+2)+
     +                       19.0D0*F(LL,IST+1)+9.0D0*F(LL,IST))/24.0D0
   10     CONTINUE

        ELSE
          DO 20 LL = 1,LMMSQD
            FINT(LL,IST) = FINT(LL,IST-1)
c---> integrate fint(ist+1) with a 4 point lagrangian
            FINT(LL,IST+1) = FINT(LL,IST-1) +
     +                       (F(LL,IST+3)-5.0D0*F(LL,IST+2)+
     +                       19.0D0*F(LL,IST+1)+9.0D0*F(LL,IST))/24.0D0
   20     CONTINUE
        END IF

c
c---> calculate fint with an extended 3-point-simpson
c
        DO 40 I = IST + 2,IEN
          DO 30 LL = 1,LMMSQD
            FINT(LL,I) = ((FINT(LL,I-2)+F(LL,I-2)*A1)+F(LL,I-1)*A2) +
     +                   F(LL,I)*A1
   30     CONTINUE
   40   CONTINUE
   50 CONTINUE
c
      END SUBROUTINE
      END MODULE mod_CSOUT
