      SUBROUTINE APPBLCKCIRC(VECS,GLLHBLCK,
C                            new after inc.p removal - arg. INFO removed
     &                       naez,lmax,nthrds,
     &                       natbld, xdim, ydim, zdim)
C
      IMPLICIT NONE
C
C      INCLUDE 'inc.p'
      INCLUDE 'fftw3.f'

      INTEGER naez
      INTEGER lmax
C     number of OpenMP threads
      INTEGER nthrds
C     number of atoms per preconditioning block
      INTEGER natbld
C     number of preconditioning blocks in each direction
      INTEGER xdim
      INTEGER ydim
      INTEGER zdim


      INTEGER        LMMAXD

      DOUBLE COMPLEX CONE
      PARAMETER      (CONE  = ( 1.0D0,0.0D0))
      DOUBLE COMPLEX CZERO
      PARAMETER      (CZERO = ( 0.0D0,0.0D0))
C

C
C     DOUBLE COMPLEX VECS(NAEZD*LMMAXD,LMMAXD),
C     +               GLLHBLCK(LMMAXD*NATBLD,
C     +                        LMMAXD*NATBLD*XDIM*YDIM*ZDIM),
C     +               BLCK(LMMAXD*NATBLD,LMMAXD*NATBLD),
C     +               TBLCK(LMMAXD*NATBLD,LMMAXD*NATBLD),
C     +               TXK(LMMAXD*NATBLD),TYK(LMMAXD*NATBLD)

      DOUBLE COMPLEX VECS(NAEZ*(LMAX+1)**2,(LMAX+1)**2),
     +               GLLHBLCK(NATBLD*(LMAX+1)**2,
     +                        XDIM*YDIM*ZDIM*NATBLD*(LMAX+1)**2),
     +               BLCK(NATBLD*(LMAX+1)**2,NATBLD*(LMAX+1)**2),
     +               TBLCK(NATBLD*(LMAX+1)**2,NATBLD*(LMAX+1)**2),
     +               TXK(NATBLD*(LMAX+1)**2),TYK(NATBLD*(LMAX+1)**2)

C
C      DOUBLE COMPLEX X(XDIM,YDIM,ZDIM),
C     +               XK(NATBLD*LMMAXD,XDIM,YDIM,ZDIM),
C     +               Y(XDIM,YDIM,ZDIM),
C     +               YK(NATBLD*LMMAXD,XDIM,YDIM,ZDIM)

      DOUBLE COMPLEX X(XDIM,YDIM,ZDIM),
     +               XK(NATBLD*(LMAX+1)**2,XDIM,YDIM,ZDIM),
     +               Y(XDIM,YDIM,ZDIM),
     +               YK(NATBLD*(LMAX+1)**2,XDIM,YDIM,ZDIM)

C ..
C local scalars ..
      DOUBLE COMPLEX TRACE
      DOUBLE PRECISION FAC
      INTEGER        LM1,LMATBL,IX,IY,IZ,I,J
      INTEGER*8      FFTWPLAN

      LOGICAL        LSAME
      INTEGER        MYTHRD,OMP_GET_THREAD_NUM

C..
C
      EXTERNAL       LSAME
C
C
C=======================================================================
C outer loop over all columns LM1=1, LMMAXD       begin
C=======================================================================
C
      LMMAXD= (LMAX+1)**2

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
        CALL OMP_SET_NUM_THREADS(NTHRDS)
!$OMP  PARALLEL PRIVATE (LMATBL,IX,IY,IZ,X,MYTHRD)
        MYTHRD = OMP_GET_THREAD_NUM()
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
        CALL OMP_SET_NUM_THREADS(NTHRDS)
!$OMP  PARALLEL PRIVATE (IX,IY,IZ,TBLCK,TXK,TYK,MYTHRD)
        MYTHRD = OMP_GET_THREAD_NUM()

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
C              CALL ZGEMV_SRC(LMMAXD*NATBLD,LMMAXD*NATBLD,
C     +                       CONE,TBLCK,LMMAXD*NATBLD,TXK,TYK)
C
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
        CALL OMP_SET_NUM_THREADS(NTHRDS)
!$OMP  PARALLEL PRIVATE (LMATBL,IX,IY,IZ,Y,MYTHRD)
        MYTHRD = OMP_GET_THREAD_NUM()
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
C
C
C
C
      SUBROUTINE ZGEMV_SRC(M,N,ALPHA,A,LDA,X,Y)
      DOUBLE COMPLEX ALPHA
      INTEGER        LDA,M,N
      DOUBLE COMPLEX A(LDA,*),X(*),Y(*)
      DOUBLE COMPLEX ZERO
      PARAMETER     (ZERO= (0.0D+0,0.0D+0))
      DOUBLE COMPLEX TEMP
      INTEGER        I,J,JX,JY
C     ..
      DO 10 I = 1,M
          Y(I) = ZERO
   10 CONTINUE
      JX = 1
      DO 60 J = 1,N
          IF (X(JX).NE.ZERO) THEN
              TEMP = ALPHA*X(JX)
              DO 50 I = 1,M
                  Y(I) = Y(I) + TEMP*A(I,J)
   50         CONTINUE
          END IF
          JX = JX + INCX
   60 CONTINUE
      RETURN
      END
