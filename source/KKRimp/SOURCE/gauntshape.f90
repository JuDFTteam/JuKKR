   !-------------------------------------------------------------------------------
   !> Summary: this is module is used to generat the Gauntcoefficients convoluted with shapefunction 
   !> Author: Who wrote this subroutine
   !> Category: TAGS for the code they must be written as TAG1, TAG2, ..., TAGN
   !> Deprecated: False ! This needs to be set to True for deprecated subroutines
   !> A More detailed explanation with the math, concepts, etc necessary to understand the routine
   !-------------------------------------------------------------------------------
   !> @note Notes on the code
   !> @endnote
   !> @todo things that must be checked
   !> @endtodo
   !> @warning Important precautions
   !> @endwarning
   !> @bug If nasty things are found
   !> @endbug
   !-------------------------------------------------------------------------------
MODULE MOD_GAUNTSHAPE
CONTAINS
   !-------------------------------------------------------------------------------
   !> Summary: This subroutine generats the Gauntcoefficients convoluted with shapefunction 
   !> Author: 
   !> Category: Shapefun, Potential, kkrimp 
   !> Deprecated: False ! This needs to be set to True for deprecated subroutines
   !> A More detailed explanation with the math, concepts, etc necessary to understand the routine
   !-------------------------------------------------------------------------------
   !> @note Notes on the code
   !> @endnote
   !> @todo things that must be checked
   !> @endtodo
   !> @warning Important precautions
   !> @endwarning
   !> @bug If nasty things are found
   !> @endbug
   !-------------------------------------------------------------------------------

     SUBROUTINE GEN_GAUNTSHAPE(lmaxd,NATOM,LMAXATOM,gauntshape,shapefun,ins)
use type_gauntshape
use type_shapefun
use mod_gauntharmonics, only: gauntcoeff
implicit none
!interface
integer                                    :: lmaxd
integer                                    :: natom
integer                                    :: lmaxatom(natom)
type(gauntshape_type),allocatable          :: gauntshape(:)
type(shapefun_type)                        :: shapefun(natom)
!local
! integer :: lmaxbounds(2)
integer,parameter                          :: NGSHD=4*13079 !bauer march2012
integer                                    :: ins,lval
!     GET ALL LMAX VALUES
! call getlmaxbounds(lmaxatom,natom,lmaxbounds)

 allocate(gauntshape( lmaxd ) )

 do lval = 2, lmaxd

  allocate( gauntshape(lval)%gsh(NGSHD), gauntshape(lval)%ilm(NGSHD,3), gauntshape(lval)%imaxsh(0:(2*lmaxd+1)**2) )

gauntshape(lval)%gsh=0.0D0
gauntshape(lval)%ilm=0
gauntshape(lval)%imaxsh=0

if (ins/=0) then 

  call shape(2*lval,natom,gauntshape(lval)%gsh,gauntshape(lval)%ilm,gauntshape(lval)%imaxsh, &
              shapefun, gauntcoeff(lval)%wg,gauntcoeff(lval)%yrg,4*lval,(2*lval+1)**2,NGSHD )
  gauntshape%ngshd=NGSHD

end if
 end do !ival

!    allocate typegauntshape

!    call shape() for all lmax values

     END SUBROUTINE GEN_GAUNTSHAPE

   !-------------------------------------------------------------------------------
   !> Summary: Prepares shape corrections using gaussian quadrature as given by
   !>          m. abramowitz and i.a. stegun, handbook of mathematical functions
   !>          nbs applied mathematics series 55 (1968), pages 887 and 916
   !> Author: Who wrote this subroutine
   !> Category: Shapefun, Potential 
   !> Deprecated: False ! This needs to be set to True for deprecated subroutines
   !> A More detailed explanation with the math, concepts, etc necessary to understand the routine
   !-------------------------------------------------------------------------------
   !> @note the parameter LASSLD has to be chosen such that l1+l2+l3 .le. 2*LASSLD 
   !> @endnote
   !> @todo things that must be checked
   !> @endtodo
   !> @warning Important precautions
   !> @endwarning
   !> @bug If nasty things are found
   !> @endbug
   !-------------------------------------------------------------------------------


      SUBROUTINE SHAPE(LPOT,NATOM,GSH,ILM,IMAXSH,SHAPEFUN,W,YR, &
                      LASSLD,LMPOTD,NGSHD)
!C **********************************************************************
!C *  Prepares shape corrections using gaussian quadrature as given by  *
!C *  m. abramowitz and i.a. stegun, handbook of mathematical functions *
!C *  nbs applied mathematics series 55 (1968), pages 887 and 916       *
!C *                                                                    *
!C *  the parameter LASSLD has to be chosen such that                   *
!C *                        l1+l2+l3 .le. 2*LASSLD                      *
!C *                                                                    *
!C **********************************************************************
      USE TYPE_SHAPEFUN
      IMPLICIT NONE
!C     ..
!C     .. Scalar Arguments ..
      INTEGER LASSLD,LMPOTD,NGSHD
      INTEGER LPOT,NATOM
!C     ..
!C     .. Array Arguments ..
      DOUBLE PRECISION GSH(NGSHD),W(LASSLD),YR(LASSLD,0:LASSLD,0:LASSLD)
      INTEGER ILM(NGSHD,3),IMAXSH(0:LMPOTD)
      TYPE(SHAPEFUN_TYPE)       :: shapefun(natom)
! LMSP(:,:)
!C     ..
!C     .. Local Scalars ..
      DOUBLE PRECISION FACTOR,GAUNT,S
      INTEGER I,IAT,ISUM,J,L1,L2,L3,LM1,LM2,LM3,M1,M1A,M1S,M2,M2A, &
             M2S,M3,M3A,M3S
!       LOGICAL TRIANGLE
!C     ..
!C     .. Intrinsic Functions ..
!       INTRINSIC ABS,DBLE,SIGN
!C     ..
!C     .. External Subroutines ..
!       EXTERNAL RCSTOP,TRIANGLE
!C     ..
!C
!C -> set up of the gaunt coefficients with an index field
!C    so that  c(lm,lm',lm'') is mapped to c(i)
      I = 1
!C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      DO L1 = 0,LPOT
!C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
         DO M1 = -L1,L1
!C
            LM1 = L1*L1 + L1 + M1 + 1
            IMAXSH(LM1-1) = I - 1
!C llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
            DO L3 = 0,LPOT*2
!C mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
               DO M3 = -L3,L3
!C
                  LM3 = L3*L3 + L3 + M3 + 1
                  ISUM = 0
!C     
                  DO IAT = 1,NATOM
!                      ICELL = NTCELL(IAT)
                IF (LM3<=UBOUND(SHAPEFUN(IAT)%LMUSED,1)) ISUM = ISUM + SHAPEFUN(IAT)%LMUSED(LM3)
!                 IF (lm3<=ubound(lmsp,1)>lm3) ISUM = ISUM + LMSP(IAT,LM3)
!                      ISUM = ISUM + LMSP(IAT,LM3)
                  END DO
!C     
!C ======================================================================
                  IF ( ISUM.GT.0 ) THEN
                     DO L2 = 0,LPOT
!C ----------------------------------------------------------------------
                        IF ( TRIANGLE(L1,L2,L3) ) THEN
                           DO M2 = -L2,L2
!C
                              LM2 = L2*L2 + L2 + M2 + 1
!C
!C -> use the m-conditions for the gaunt coefficients not to be 0
!C
                              M1S = SIGN(1,M1)
                              M2S = SIGN(1,M2)
                              M3S = SIGN(1,M3)
!C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                              IF ( M1S*M2S*M3S.GE.0 ) THEN
                                 M1A = ABS(M1)
                                 M2A = ABS(M2)
                                 M3A = ABS(M3)
                                 FACTOR = 0.0D0
!C
                                 IF (M1A+M2A.EQ.M3A) FACTOR = FACTOR + &
                                         DBLE(3*M3S+SIGN(1,-M3))/8.0D0
!C
                                 IF (M1A-M2A.EQ.M3A) FACTOR = FACTOR + &
                                         DBLE(M1S)/4.0D0
!C
                                 IF (M2A-M1A.EQ.M3A) FACTOR = FACTOR + &
                                         DBLE(M2S)/4.0D0
!C ......................................................................
                                 IF (FACTOR.NE.0.0D0) THEN
!C
                                    IF ( M1S*M2S.NE.1 .OR. M2S*M3S.NE.1  &
                                     .OR.M1S*M3S.NE.1 ) FACTOR = -FACTOR
!C
                                    S = 0.0D0
                                    DO J = 1,LASSLD
                                       S = S + W(J) * YR(J,L1,M1A) &
                                           * YR(J,L2,M2A) * YR(J,L3,M3A) 
                                    END DO
!C
                                    GAUNT = S*FACTOR
                                    IF ( ABS(GAUNT).GT.1D-10 ) THEN
                                       GSH(I) = GAUNT
                                       ILM(I,1) = LM1
                                       ILM(I,2) = LM2
                                       ILM(I,3) = LM3
                                       I = I + 1
                                    END IF
                                 END IF
!C ......................................................................
                              END IF
!C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                           END DO
                        END IF
!C ----------------------------------------------------------------------
                     END DO
                  END IF
!C ======================================================================
               END DO
!C mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
            END DO
!C llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
         END DO
!C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      END DO
!C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
!C
      IMAXSH(LM1) = I - 1
      WRITE (1337,FMT=9000) IMAXSH(LM1),NGSHD
      IF ( IMAXSH(LM1).GT.NGSHD ) STOP 'SHAPE   '
!C
 9000 FORMAT(' >>> SHAPE : IMAXSH(',I6,'),NGSHD :',2I6)
!C
      END SUBROUTINE SHAPE
!C

   !-------------------------------------------------------------------------------
   !> Summary: This function when used insures to have only non zero Gauntcoffecients
   !>           The Gauant coefficients are different from zero only if the L1, L2 and L3
   !>           satisfie the equation: abs(L3-L2)<= L1 >= (L3+L2)
   !>          
   !> Author: Who wrote this function
   !> Category: Potential
   !> Deprecated: False ! This needs to be set to True for deprecated functions
   !> A More detailed explanation with the math, concepts, etc necessary to understand the function
   !-------------------------------------------------------------------------------
   !> @note Notes on the code
   !> @endnote
   !> @todo things that must be checked
   !> @endtodo
   !> @warning Important precautions
   !> @endwarning
   !> @bug If nasty things are found
   !> @endbug
   !-------------------------------------------------------------------------------



      FUNCTION TRIANGLE(L1,L2,L3)
      IMPLICIT NONE
      INTEGER L1,L2,L3
      LOGICAL TRIANGLE
      INTRINSIC MOD
!C     ..
      TRIANGLE = (L1.GE.ABS(L3-L2)) .AND. (L1.LE.(L3+L2)) &
           .AND. (MOD((L1+L2+L3),2).EQ.0)
      END FUNCTION TRIANGLE

   !-------------------------------------------------------------------------------
   !> Summary: To get the lmax bounds for every atom
   !> Author: Who wrote this subroutine
   !> Category: TAGs
   !> Deprecated: False ! This needs to be set to True for deprecated subroutines
   !> A More detailed explanation with the math, concepts, etc necessary to understand the routine
   !-------------------------------------------------------------------------------
   !> @note Notes on the code
   !> @endnote
   !> @todo things that must be checked
   !> @endtodo
   !> @warning Important precautions
   !> @endwarning
   !> @bug If nasty things are found
   !> @endbug
   !-------------------------------------------------------------------------------


subroutine getlmaxbounds(lmaxatom,natom,lmaxbounds)
    implicit none
    integer, intent(in)       :: lmaxatom(natom)
    integer, intent(in)       :: natom
    integer, intent(inout)    :: lmaxbounds(2)
    integer                   :: iatom
  lmaxbounds=0
  do iatom=1,natom
    if (lmaxbounds(1)==0 .and. lmaxbounds(2)==0) lmaxbounds=lmaxatom(iatom)
    if (lmaxatom(iatom)<lmaxbounds(1)) lmaxbounds(1)=lmaxatom(iatom)
    if (lmaxatom(iatom)>lmaxbounds(2)) lmaxbounds(2)=lmaxatom(iatom)
  end do !iatom
end subroutine getlmaxbounds


END MODULE MOD_GAUNTSHAPE
