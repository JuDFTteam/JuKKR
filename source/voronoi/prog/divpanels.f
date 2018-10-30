      SUBROUTINE DIVPANELS(
     > NFUN,NMIN,NSMALL,     
     X THETAS,XRN,DRN,NPAN,MESHN,NM)
c Divide the already existing panels of shape functions into smaller subpanels.
c Keep radial points, weights and calculated shape functions unchanged.
c Original arrays THETAS,XRN,DRN,NPAN,MESHN,NM contain new values on output.
c#@# KKRtags: VORONOI geometry shape-functions single-site radial-grid
      IMPLICIT NONE
      include 'inc.geometry'
      INTEGER IBMAXD
      PARAMETER (IBMAXD=(LMAXD1+1)*(LMAXD1+1))
c Input shapes:
      REAL*8 THETAS(IRID,IBMAXD),XRN(IRID),DRN(IRID)
      INTEGER NPAN,MESHN,NM(NPAND),NFUN,NMIN
c New (output) shapes:
      REAL*8 THETASNEW(IRID,IBMAXD),XRNNEW(IRID),DRNNEW(IRID)
      INTEGER NPANNEW,MESHNNEW,NMNEW(NPAND),NSMALL
c Criterion for panel width
      REAL*8 WIDTHMAX

c Local:
      INTEGER IP,NDIVIDE,IDIVIDE,NADDED,IR,IRNEW,IRUN
      INTEGER IFUN,NSMALL1,NMOD
      REAL*8 WIDTH
      

c Initialize
      THETASNEW(1:IRID,1:IBMAXD) = 0.D0
      XRNNEW(1:IRID) = 0.D0
      DRNNEW(1:IRID) = 0.D0
      NMNEW(1:NPAND) = 0
      NPANNEW = 0
      MESHNNEW = 0

c Take care that lower limit of points is not violated
      NSMALL1 = MAX(NMIN,NSMALL)
c Apply criterion to each panel
      NPANNEW = 0
      IR = 0
      IRNEW = 0
      MESHNNEW = MESHN
      DO IP = 1,NPAN
c Find number of (sub)panels in new division


         NMOD = MOD(NM(IP)-1,NSMALL1-1)
         NDIVIDE = (NM(IP)-1)/(NSMALL1-1)

c        NDIVIDE = (NM(IP)-1) / NSMALL1  ! Integer. Divide in panels of NSMALL1 points each.
         NADDED = NDIVIDE - 1            ! New points have to be added at the end of each new sub-panel
         DO IDIVIDE = 1,NDIVIDE - 1      ! Loop over new panels except last
            NPANNEW = NPANNEW + 1
            NMNEW(NPANNEW) = NSMALL1
c Copy shapes into sub-panel. The radial points, weights & shape values are the same as before.
            DO IRUN = 1,NMNEW(NPANNEW)
               IRNEW = IRNEW + 1
               IR = IR + 1
               XRNNEW(IRNEW) = XRN(IR)
               DRNNEW(IRNEW) = DRN(IR)
               DO IFUN = 1,NFUN
                  THETASNEW(IRNEW,IFUN) = THETAS(IR,IFUN)
               ENDDO
            ENDDO
            IR = IR - 1  ! In the old mesh, the end-point of the new panel was not taken double.
         ENDDO
         NPANNEW = NPANNEW + 1     ! Last sub-panel of original panel IP can have NMNEW > NSMALL1.
         NMNEW(NPANNEW) = NSMALL1 + NMOD  ! Total number of pts minus already-used pts
         IF (NDIVIDE.EQ.0) NMNEW(NPANNEW) = NM(IP)
         DO IRUN = 1,NMNEW(NPANNEW)
            IRNEW = IRNEW + 1
            IR = IR + 1
            XRNNEW(IRNEW) = XRN(IR)
            DRNNEW(IRNEW) = DRN(IR)
            DO IFUN = 1,NFUN
               THETASNEW(IRNEW,IFUN) = THETAS(IR,IFUN)
            ENDDO
         ENDDO
         ! Here the  IR = IR - 1 is not applied, because the end of the original panel is a double-point.
      ENDDO ! IP = 1,NPAN
      MESHNNEW = IRNEW


c Copy to original arrays.
      THETAS(1:IRID,1:IBMAXD) = THETASNEW(1:IRID,1:IBMAXD)
      XRN(1:IRID) = XRNNEW(1:IRID)
      DRN(1:IRID) = DRNNEW(1:IRID)
      NM(1:NPAND) = NMNEW(1:NPAND)
      NPAN = NPANNEW
      MESHN = MESHNNEW


      RETURN
      END
