      SUBROUTINE APPBLCKCIRC(VECS,GLLHBLCK,
C                            new after inc.p removal - arg. INFO removed
     &                       naez,lmmaxd,nthrds,
     &                       natbld, xdim, ydim, zdim)
C
      IMPLICIT NONE

      INCLUDE 'fftw3.f'

      INTEGER naez
      INTEGER lmmaxd
C     number of OpenMP threads
      INTEGER nthrds
C     number of atoms per preconditioning block
      INTEGER natbld
C     number of preconditioning blocks in each direction
      INTEGER xdim
      INTEGER ydim
      INTEGER zdim


      DOUBLE COMPLEX CONE
      PARAMETER      (CONE  = ( 1.0D0,0.0D0))
      DOUBLE COMPLEX CZERO
      PARAMETER      (CZERO = ( 0.0D0,0.0D0))

      DOUBLE COMPLEX VECS(NAEZ*LMMAXD,LMMAXD)
      DOUBLE COMPLEX GLLHBLCK(NATBLD*LMMAXD,
     &                        XDIM*YDIM*ZDIM*NATBLD*LMMAXD)

C     local arrays
      DOUBLE COMPLEX TBLCK(NATBLD*LMMAXD,NATBLD*LMMAXD)
      DOUBLE COMPLEX TXK(NATBLD*LMMAXD)
      DOUBLE COMPLEX TYK(NATBLD*LMMAXD)

C     local arrays
      DOUBLE COMPLEX X(XDIM,YDIM,ZDIM)
      DOUBLE COMPLEX XK(NATBLD*LMMAXD,XDIM,YDIM,ZDIM)
      DOUBLE COMPLEX Y(XDIM,YDIM,ZDIM)
      DOUBLE COMPLEX YK(NATBLD*LMMAXD,XDIM,YDIM,ZDIM)

C ..
C local scalars ..
      DOUBLE PRECISION FAC
      INTEGER        LM1,LMATBL,IX,IY,IZ,I,J
      INTEGER*8      FFTWPLAN

      LOGICAL        LSAME
      INTEGER        MYTHRD
!$    INTEGER        OMP_GET_THREAD_NUM

C..
C
      EXTERNAL       LSAME
C
C
C=======================================================================
C outer loop over all columns LM1=1, LMMAXD       begin
C=======================================================================
C

      FAC = 1/FLOAT(XDIM*YDIM*ZDIM)
C
      DO LM1=1,LMMAXD
C
C=======================================================================
C
C
C
C-----------------------------------------------------------------------
C perform Fourier-transform backward              begin
C-----------------------------------------------------------------------
C
        CALL DFFTW_PLAN_DFT_3D(FFTWPLAN,XDIM,YDIM,ZDIM,X,X,
     +                              FFTW_BACKWARD,FFTW_ESTIMATE)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!$      CALL OMP_SET_NUM_THREADS(NTHRDS)
!$OMP  PARALLEL PRIVATE (LMATBL,IX,IY,IZ,X,MYTHRD)
        MYTHRD = 0
!$      MYTHRD = OMP_GET_THREAD_NUM()
C
        DO LMATBL = 1,LMMAXD*NATBLD
        IF (MOD(LMATBL,NTHRDS).EQ.MYTHRD) THEN
C
          DO IX =1,XDIM
            DO IY =1,YDIM
              DO IZ =1,ZDIM

              X(IX,IY,IZ) = VECS(((IX-1)+(IY-1)*XDIM+
     +                      (IZ-1)*XDIM*YDIM)*LMMAXD*NATBLD+LMATBL,LM1)

              ENDDO
            ENDDO
          ENDDO

          CALL DFFTW_EXECUTE_DFT(FFTWPLAN,X,X)

          DO IX =1,XDIM
            DO IY =1,YDIM
              DO IZ =1,ZDIM
                XK(LMATBL,IX,IY,IZ) = X(IX,IY,IZ)*FAC
              ENDDO
            ENDDO
          ENDDO
C
        ENDIF
        ENDDO

!$OMP END PARALLEL
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
        CALL DFFTW_DESTROY_PLAN(FFTWPLAN)
C
C-----------------------------------------------------------------------
C perform Fourier-transform backward              end
C-----------------------------------------------------------------------
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!$      CALL OMP_SET_NUM_THREADS(NTHRDS)
!$OMP  PARALLEL PRIVATE (IX,IY,IZ,TBLCK,TXK,TYK,MYTHRD)
        MYTHRD = 0
!$      MYTHRD = OMP_GET_THREAD_NUM()

        DO IZ =1,ZDIM
        IF (MOD(IZ,MIN(ZDIM,NTHRDS)).EQ.MYTHRD) THEN
          DO IY =1,YDIM
            DO IX =1,XDIM
C
              DO I=1,LMMAXD*NATBLD
                DO J=1,LMMAXD*NATBLD
                  TBLCK(I,J) = GLLHBLCK(I,LMMAXD*NATBLD*
     +                    ((IX-1)+(IY-1)*XDIM+(IZ-1)*XDIM*YDIM)+J)
                ENDDO
                TXK(I) = XK(I,IX,IY,IZ)
              ENDDO
C
              CALL ZGEMV('N',LMMAXD*NATBLD,LMMAXD*NATBLD,
     +                   CONE,TBLCK,
     +                   LMMAXD*NATBLD,TXK,
     +                   1,CZERO,TYK,1)

              DO I=1,LMMAXD*NATBLD
                YK(I,IX,IY,IZ) = TYK(I)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ENDDO
C
!$OMP END PARALLEL
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C
C-----------------------------------------------------------------------
C perform Fourier-transform forward               begin
C-----------------------------------------------------------------------
C
        CALL DFFTW_PLAN_DFT_3D(FFTWPLAN,XDIM,YDIM,ZDIM,Y,Y,
     +                              FFTW_FORWARD,FFTW_ESTIMATE)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!$      CALL OMP_SET_NUM_THREADS(NTHRDS)
!$OMP  PARALLEL PRIVATE (LMATBL,IX,IY,IZ,Y,MYTHRD)
        MYTHRD = 0
!$      MYTHRD = OMP_GET_THREAD_NUM()
C
        DO LMATBL = 1,LMMAXD*NATBLD
        IF (MOD(LMATBL,NTHRDS).EQ.MYTHRD) THEN
C
          DO IX =1,XDIM
            DO IY =1,YDIM
              DO IZ =1,ZDIM
                Y(IX,IY,IZ) = YK(LMATBL,IX,IY,IZ)
              ENDDO
            ENDDO
          ENDDO

          CALL DFFTW_EXECUTE_DFT(FFTWPLAN,Y,Y)

          DO IX =1,XDIM
            DO IY =1,YDIM
              DO IZ =1,ZDIM

        VECS(((IX-1)+(IY-1)*XDIM+
     +  (IZ-1)*XDIM*YDIM)*LMMAXD*NATBLD+LMATBL,LM1) = Y(IX,IY,IZ)

              ENDDO
            ENDDO
          ENDDO
C
        ENDIF
        ENDDO
C
!$OMP END PARALLEL
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
        CALL DFFTW_DESTROY_PLAN(FFTWPLAN)
C
C
C-----------------------------------------------------------------------
C perform Fourier-transform forward               end
C-----------------------------------------------------------------------
C
C
C=======================================================================
C
      ENDDO
C
C=======================================================================
C outer loop over all columns LM1=1, LMMAXD       end
C=======================================================================
C
      END
