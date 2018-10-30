!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module mod_dlke

  implicit none

    private
    public :: dlke0, dlke0dk

  contains

      !-------------------------------------------------------------------------------
      !> Summary: Driver for Fourier transform of structure constants
      !> Author: B. Zimmermann
      !> Date: 01.08.2010
      !> Category: PKKprime, k-points, reference-system, structural-greensfunction
      !> Deprecated: False ! This needs to be set to True for deprecated subroutines
      !>
      !> @note Simplified copy of dlke0 from host code @endnote
      !>
      !> modified by Bernd Zimmermann on 01.08.2010:
      !>   - deleted sparse matrix option
      !>   - deleted complex k-point option -> length from array reduced from 6 to 3
      !>   - symmetrization commented
      !>   - transformed to fortran90 standard
      !> modified by Bernd Zimmermann on 13.10.2014:
      !>   - deleted EQINV and KAOEZ-arrays
      !-------------------------------------------------------------------------------
      SUBROUTINE DLKE0(GLLKE, inc, alat, cls, nacls, rr, ezoa, atom, bzkp, rcls, ginp)

      use type_inc

      implicit none
!     .. Parameters ..
      DOUBLE COMPLEX CZERO,CI
      PARAMETER (CZERO= (0.0D0,0.0D0),CI= (0.0D0,1.0D0))
      DOUBLE COMPLEX CONE
      PARAMETER (CONE= (1.0D0,0.0D0))

!     .. Scalar Arguments ..
      type(inc_type),   intent(in) :: inc
      double precision, intent(in) :: alat

!     .. Array Arguments ..
      INTEGER,          INTENT(IN)  :: ATOM(inc%naclsd,inc%naezd), CLS(inc%natypd), EZOA(inc%naclsd,inc%naezd), NACLS(inc%nclsd)
      DOUBLE PRECISION, INTENT(IN)  :: BZKP(3), RR(3,0:inc%nrd), RCLS(3,inc%naclsd,inc%nclsd)
      DOUBLE COMPLEX,   INTENT(IN)  :: GINP(inc%lmmax*inc%naclsd,inc%lmmax,inc%nclsd)
      DOUBLE COMPLEX,   INTENT(OUT) :: GLLKE(inc%alm,inc%alm)

!     .. Local Scalars ..
      INTEGER I,IC,IM,M,JN
      DOUBLE COMPLEX FAC

!     .. Local Arrays ..
      DOUBLE COMPLEX GLLKE1(inc%alm,inc%lmmax)
      DOUBLE PRECISION KP(3)

!     .. External Subroutines ..
      EXTERNAL ZAXPY

      GLLKE=CZERO
      DO I=1,inc%naez

        KP = BZKP
        IC  = CLS(I)

        CALL DLKE1(GLLKE1,inc,alat,nacls(IC),rr,ezoa(:,I),atom(:,I),kp,GINP(:,:,IC),RCLS(:,:,IC))
             ! Changed on 22.03.2000 ! forgoten I 3.11.2000
             ! Correction for fourier trans.

         DO M=1,inc%lmmax
            IM=(I-1)*inc%lmmax+M
            DO JN=1,inc%lmmax*inc%naez
               GLLKE(JN,IM) = GLLKE(JN,IM)+ GLLKE1(JN,M)
            END DO
         END DO

      END DO! I=1,NAEZ
      
      RETURN

      END SUBROUTINE DLKE0

      !-------------------------------------------------------------------------------
      !> Summary: Calculates Fourier transform of structure constants
      !> Author: B. Zimmermann
      !> Category: PKKprime, k-points, reference-system, structural-greensfunction
      !> Deprecated: False ! This needs to be set to True for deprecated subroutines
      !>
      !> @note Simplified copy of dlke1 from host code @endnote
      !>
      !> Fourier transformation of the cluster Greens function GINP
      !-------------------------------------------------------------------------------
      SUBROUTINE DLKE1(GLLKE,inc,alat,nacls,rr,ezoa,atom,bzkp,ginp,rcls)

      use type_inc
      use mod_mathtools, only: tpi
      implicit none

!     .. Parameters ..
      DOUBLE COMPLEX CZERO,CI
      PARAMETER (CZERO= (0.0D0,0.0D0),CI= (0.0D0,1.0D0))
      DOUBLE COMPLEX CONE
      PARAMETER (CONE=(1.0D0,0.0D0))

!     .. Scalar Arguments ..
      type(inc_type),   intent(in) :: inc
      INTEGER,          INTENT(IN) :: nacls
      DOUBLE PRECISION, INTENT(IN) :: alat

!     .. Array Arguments ..
      INTEGER,          INTENT(IN) :: atom(inc%naclsd), ezoa(inc%naclsd)
      DOUBLE COMPLEX,   INTENT(IN) :: ginp(inc%lmmax*inc%naclsd,inc%lmmax)
      DOUBLE PRECISION, INTENT(IN) :: bzkp(3),rr(3,0:inc%nrd),rcls(3,inc%naclsd)

      DOUBLE COMPLEX, INTENT(OUT):: GLLKE(inc%alm,inc%lmmax)

!     .. Local Scalars ..
      DOUBLE PRECISION CONVPU
      INTEGER AM,IM,LM2,M
      DOUBLE COMPLEX  EIKR,TT

!     .. Local Arrays ..
      DOUBLE COMPLEX ARG(3)
!     ..
!     .. External Subroutines ..
      EXTERNAL ZAXPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS,DATAN,EXP
!     ..
!     .. COMMON BLOCK
!     ..
!     .. Data statements ..
! ------------------------------------------------------------------------

      GLLKE=CZERO

      CONVPU = ALAT/tpi   !=alat / (2*pi)

      DO M = 1,NACLS

! --->  for option 'WIRE': avoid artifical couplings in the
!       structural Greens Function in in-plane-direction (perp. to c-axis)
        IF (ATOM(M).LT.0) CYCLE
!     
!     Here we do   --                  nn'
!                  \                   ii'          ii'
!                  /  exp(+ik(x  -x ))G   (E)  =   G   (k,E)
!                  --          n'  n   LL'          LL'
!                  n'
!  Be carefull a minus sign must be included here. RR is not
!  symmetric around each atom. The minus comes from the fact that
!  the repulsive potential GF is calculated for 0n and not n0!          
!  and that is why we nead a minus sign extra!
!  
        ARG(1) = -CI*TPI*RR(1,EZOA(M))
        ARG(2) = -CI*TPI*RR(2,EZOA(M))
        ARG(3) = -CI*TPI*RR(3,EZOA(M))

        TT = BZKP(1)*ARG(1)+BZKP(2)*ARG(2)+BZKP(3)*ARG(3)
        EIKR = EXP(TT) * CONVPU    ! convert to p.u.

        IM = 1 + (M-1)      *inc%lmmax
        AM = 1 + (ATOM(M)-1)*inc%lmmax
        DO LM2 = 1,inc%lmmax
          CALL ZAXPY(inc%lmmax,EIKR,GINP(IM,LM2),1,GLLKE(AM,LM2),1)
        END DO

      END DO! M = 1,NACLS

      RETURN
      END SUBROUTINE DLKE1





      !-------------------------------------------------------------------------------
      !> Summary: Driver for k-derivative of structure constants
      !> Author: B. Zimmermann
      !> Date: 25.11.2014
      !> Category: PKKprime, k-points, reference-system, structural-greensfunction
      !> Deprecated: False ! This needs to be set to True for deprecated subroutines
      !>
      !> Calculate the ananlytical derivative of the Fourier transformed structure constants with respect to k
      !> @warning Might be buggy. Needs testing! @endwarning
      !-------------------------------------------------------------------------------
      SUBROUTINE DLKE0DK(GLLKE, inc, alat, cls, nacls, rr, ezoa, atom, bzkp, rcls, ginp,ixyz)

      use type_inc

      implicit none
!     .. Parameters ..
      DOUBLE COMPLEX CZERO,CI
      PARAMETER (CZERO= (0.0D0,0.0D0),CI= (0.0D0,1.0D0))
      DOUBLE COMPLEX CONE
      PARAMETER (CONE= (1.0D0,0.0D0))

!     .. Scalar Arguments ..
      type(inc_type),   intent(in) :: inc
      integer,          intent(in) :: ixyz
      double precision, intent(in) :: alat

!     .. Array Arguments ..
      INTEGER,          INTENT(IN)  :: ATOM(inc%naclsd,inc%naezd), CLS(inc%natypd), EZOA(inc%naclsd,inc%naezd), NACLS(inc%nclsd)
      DOUBLE PRECISION, INTENT(IN)  :: BZKP(3), RR(3,0:inc%nrd), RCLS(3,inc%naclsd,inc%nclsd)
      DOUBLE COMPLEX,   INTENT(IN)  :: GINP(inc%lmmax*inc%naclsd,inc%lmmax,inc%nclsd)
      DOUBLE COMPLEX,   INTENT(OUT) :: GLLKE(inc%alm,inc%alm)

!     .. Local Scalars ..
      INTEGER I,IC,IM,M,JN
      DOUBLE COMPLEX FAC

!     .. Local Arrays ..
      DOUBLE COMPLEX GLLKE1(inc%alm,inc%lmmax)
      DOUBLE PRECISION KP(3)

!     .. External Subroutines ..
      EXTERNAL ZAXPY

      GLLKE=CZERO
      DO I=1,inc%naez

        KP = BZKP
        IC  = CLS(I)

        CALL DLKE1DK(GLLKE1,inc,alat,nacls(IC),rr,ezoa(:,I),atom(:,I),kp,GINP(:,:,IC),RCLS(:,:,IC),ixyz)
             ! Changed on 22.03.2000 ! forgoten I 3.11.2000
             ! Correction for fourier trans.

       DO M=1,inc%lmmax
         IM=(I-1)*inc%lmmax+M
         DO JN=1,inc%lmmax*inc%naez
           GLLKE(JN,IM) = GLLKE(JN,IM)+ GLLKE1(JN,M)
         END DO!JN
       END DO!M

      END DO! I=1,NAEZ
      
      RETURN

      END SUBROUTINE DLKE0DK




      !-------------------------------------------------------------------------------
      !> Summary: Calculates the k-derivative of structure constants
      !> Author: B. Zimmermann
      !> Date: 25.11.2014
      !> Category: PKKprime, k-points, reference-system, structural-greensfunction
      !> Deprecated: False ! This needs to be set to True for deprecated subroutines
      !>
      !> Calculate the ananlytical derivative of the Fourier transformed structure constants with respect to k
      !> @warning Might be buggy. Needs testing! @endwarning
      !>
      !>     Here we do   --                  nn'
      !>                  \                   ii'          ii'
      !>                  /  exp(+ik(x  -x ))G   (E)  =   G   (k,E)
      !>                  --          n'  n   LL'          LL'
      !>                  n'
      !>  Be carefull a minus sign must be included here. RR is not
      !>  symmetric around each atom. The minus comes from the fact that
      !>  the repulsive potential GF is calculated for 0n and not n0!          
      !>  and that is why we nead a minus sign extra!
      !-------------------------------------------------------------------------------
      SUBROUTINE DLKE1DK(GLLKE,inc,alat,nacls,rr,ezoa,atom,bzkp,ginp,rcls,ixyz)

      use type_inc
      use mod_mathtools, only: tpi
      implicit none
!     .. Parameters ..
      DOUBLE COMPLEX CZERO,CI
      PARAMETER (CZERO= (0.0D0,0.0D0),CI= (0.0D0,1.0D0))
      DOUBLE COMPLEX CONE
      PARAMETER (CONE=(1.0D0,0.0D0))

!     .. Scalar Arguments ..
      type(inc_type),   intent(in) :: inc
      INTEGER,          INTENT(IN) :: nacls, ixyz
      DOUBLE PRECISION, INTENT(IN) :: alat

!     .. Array Arguments ..
      INTEGER,          INTENT(IN) :: atom(inc%naclsd), ezoa(inc%naclsd)
      DOUBLE COMPLEX,   INTENT(IN) :: ginp(inc%lmmax*inc%naclsd,inc%lmmax)
      DOUBLE PRECISION, INTENT(IN) :: bzkp(3),rr(3,0:inc%nrd),rcls(3,inc%naclsd)

      DOUBLE COMPLEX, INTENT(OUT):: GLLKE(inc%alm,inc%lmmax)

!     .. Local Scalars ..
      DOUBLE PRECISION CONVPU
      INTEGER AM,IM,LM2,M
      DOUBLE COMPLEX  TT

!     .. Local Arrays ..
      DOUBLE COMPLEX ARG(3), EIKR
!     ..
!     .. External Subroutines ..
      EXTERNAL ZAXPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS,DATAN,EXP
!     ..
!     .. COMMON BLOCK
!     ..
!     .. Data statements ..
! ------------------------------------------------------------------------

      GLLKE=CZERO

      CONVPU = ALAT/tpi   !=alat / (2*pi)

      DO M = 1,NACLS

! --->  for option 'WIRE': avoid artifical couplings in the
!       structural Greens Function in in-plane-direction (perp. to c-axis)
        IF (ATOM(M).LT.0) CYCLE
!     
!     Here we do   --                  nn'
!                  \                   ii'          ii'
!                  /  exp(+ik(x  -x ))G   (E)  =   G   (k,E)
!                  --          n'  n   LL'          LL'
!                  n'
!  Be carefull a minus sign must be included here. RR is not
!  symmetric around each atom. The minus comes from the fact that
!  the repulsive potential GF is calculated for 0n and not n0!          
!  and that is why we nead a minus sign extra!
!  
        ARG(1) = -CI*TPI*RR(1,EZOA(M))
        ARG(2) = -CI*TPI*RR(2,EZOA(M))
        ARG(3) = -CI*TPI*RR(3,EZOA(M))

        TT = BZKP(1)*ARG(1)+BZKP(2)*ARG(2)+BZKP(3)*ARG(3)
!       the following block is modified to give the derivative with respect to k
!                                                         b.zimmermann, Nov.2014
        EIKR = ARG(ixyz) * EXP(TT) !* CONVPU    ! convert to p.u.

        IM = 1 + (M-1)      *inc%lmmax
        AM = 1 + (ATOM(M)-1)*inc%lmmax
        DO LM2 = 1,inc%lmmax
          CALL ZAXPY(inc%lmmax,EIKR,GINP(IM,LM2),1,GLLKE(AM,LM2),1)
        END DO

      END DO! M = 1,NACLS

      RETURN
      END SUBROUTINE DLKE1DK





end module mod_dlke
