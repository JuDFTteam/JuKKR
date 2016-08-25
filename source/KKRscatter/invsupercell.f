c ************************************************************************
      SUBROUTINE INVSUPERCELL(M2,M1,M3,GIN,ICHECK)
c ************************************************************************
c
c ---> ALGORITM FOR SUPERCELL GEOMETRY
c
c ------------------------------------------------------------------------
c
c ---> factorization D ^-1 = (prod L) * M * (prod U)
c
c      see notes R. Zeller
c
c ------------------------------------------------------------------------
      implicit none

c     .. Parameters ..
      include 'inc.p'
c

      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAXD+1)**2)
      INTEGER NDIM
      PARAMETER (NDIM= NPRINCD*LMAXSQ)
      INTEGER ALMD
      PARAMETER (ALMD= NAEZD*LMAXSQ)
c     .. array arguments
      DOUBLE COMPLEX
     +     M1(NDIM,NDIM,NLAYERD),
     +     M2(NDIM,NDIM,NLAYERD),
     +     M3(NDIM,NDIM,NLAYERD),
     +     GIN(ALMD,ALMD)
C
c     .. local scalars
      INTEGER N,LM,INFO,
     +        IROW,ICOL,NL
c
c     .. local arrays
      INTEGER IPVT(NDIM),ICHECK(NLAYERD,NLAYERD)
      DOUBLE COMPLEX A(NDIM,NDIM),
     +               B(NDIM,NDIM,NLAYERD),
     +               C(NDIM,NDIM,NLAYERD),
     +               D(NDIM,NDIM,NLAYERD),
     +               E(NDIM,NDIM),
     +               F(NDIM,NDIM),
     +               G(NDIM,NDIM),
     +               CUNIT(NDIM,NDIM)
c ------------------------------------------------------------------------
c     .. Data Statements ..
      DOUBLE COMPLEX CZERO,CONE
      PARAMETER (CZERO=(0.0D0,0.0D0),CONE=(1.0D0,0.0D0))
c
c     .. External Subroutines ..
      EXTERNAL CINIT,ZCOPY,ZGEMM,ZGETRF,ZGETRS,BTOM
      INTRINSIC ABS,DIMAG


c
c ---> START OF THE FACTORIZATION L * M * U
c

c ---> N =1
c
c     initialize all the matricex

      call cinit(ndim*ndim,a)
      call cinit(ndim*ndim*nlayerd,b)
      call cinit(ndim*ndim*nlayerd,c)
      call cinit(ndim*ndim*nlayerd,d)
      call cinit(ndim*ndim,e)
      call cinit(ndim*ndim,f)
      call cinit(ndim*ndim,g)
      call cinit(ndim*ndim,cunit)

c ------------------------------------------------------------------------
c
c ---> cunit = complex unity matrix of order NDIM
c
      DO 10 N = 1,NDIM
         CUNIT(N,N) = CONE
 10   CONTINUE


      CALL ZCOPY(NDIM*NDIM,M2(1,1,1),1,E(1,1),1)
      CALL ZCOPY(NDIM*NDIM,CUNIT,1,D(1,1,1),1)
      CALL ZGETRF(NDIM,NDIM,E(1,1),NDIM,IPVT,INFO)
      CALL ZGETRS('N',NDIM,NDIM,E(1,1),NDIM,IPVT,D(1,1,1),NDIM,INFO)


      NL=NLAYERD


      IF (NL.EQ.1) GOTO 980

      CALL ZCOPY(NDIM*NDIM,M2(1,1,NL),1,A(1,1),1)
      CALL ZCOPY(NDIM*NDIM,M1(1,1,NL),1,B(1,1,1),1)
      CALL ZCOPY(NDIM*NDIM,M3(1,1,NL),1,C(1,1,1),1)


c ------------------------------------------------------------------------
c
c ---> 2 <= N < NL-1
c
      IF (NL.EQ.2) GOTO 970

      IF (NL.EQ.3) GOTO 960

      DO 20 N = 2,NL-2

c
c ---> E = D(N-1) * C(N-1)
c

         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,N-1),NDIM,
     +        C(1,1,N-1),NDIM,CZERO,E(1,1),NDIM)

c
c ---> F = D(N-1) * M1(N-1)
c

         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,N-1),NDIM,
     +        M1(1,1,N-1),NDIM,CZERO,F(1,1),NDIM)

c
c ---> A = A - B(N-1)*D(N-1)*C(N-1)
c

         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,B(1,1,N-1),NDIM,
     +        E(1,1),NDIM,CONE,A(1,1),NDIM)

c
c ---> B(N) = - B(N-1)*D(N-1)*M1(N-1)
c
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,B(1,1,N-1),NDIM,
     +        F(1,1),NDIM,CZERO,B(1,1,N),NDIM)

c
c ---> C(N) = - M3(N-1)*D(N-1)*C(N-1)
c

         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,M3(1,1,N-1),NDIM,
     +        E(1,1),NDIM,CZERO,C(1,1,N),NDIM)

c
c ---> D(N) = [ M2(N) - M3(N-1)*D(N-1)*M1(N-1) ]^-1
c

         CALL ZCOPY(NDIM*NDIM,M2(1,1,N),1,E(1,1),1)
         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,M3(1,1,N-1),NDIM,
     +        F(1,1),NDIM,CONE,E(1,1),NDIM)
         CALL ZCOPY(NDIM*NDIM,CUNIT(1,1),1,D(1,1,N),1)
         CALL ZGETRF(NDIM,NDIM,E(1,1),NDIM,IPVT,INFO)
         CALL ZGETRS('N',NDIM,NDIM,E(1,1),NDIM,IPVT,D(1,1,N),NDIM,INFO)


 20   CONTINUE

c ------------------------------------------------------------------------
c
c ---> N = NL - 1
c
c

 960  N = NL - 1

c
c ---> E = D(N-1) * C(N-1)
c

      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,N-1),NDIM,
     +     C(1,1,N-1),NDIM,CZERO,E(1,1),NDIM)

c
c ---> F = D(N-1) * M1(N-1)
c

      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,N-1),NDIM,
     +     M1(1,1,N-1),NDIM,CZERO,F(1,1),NDIM)

c
c ---> A = A - B(N-1)*D(N-1)*C(N-1)
c

      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,B(1,1,N-1),NDIM,
     +     E(1,1),NDIM,CONE,A(1,1),NDIM)

c
c ---> B(N) = - B(N-1)*D(N-1)*M1(N-1) + M3(N)
c

      CALL ZCOPY(NDIM*NDIM,M3(1,1,N),1,B(1,1,N),1)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,B(1,1,N-1),NDIM,
     +     F(1,1),NDIM,CONE,B(1,1,N),NDIM)

c
c ---> C(N) = - M3(N-1)*D(N-1)*C(N-1) + M1(N)
c

      CALL ZCOPY(NDIM*NDIM,M1(1,1,N),1,C(1,1,N),1)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,M3(1,1,N-1),NDIM,
     +     E(1,1),NDIM,CONE,C(1,1,N),NDIM)

c
c ---> D(N) = [ M2(N) - M3(N-1)*D(N-1)*M1(N-1) ]^-1
c

      CALL ZCOPY(NDIM*NDIM,M2(1,1,N),1,E(1,1),1)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,M3(1,1,N-1),NDIM,
     +     F(1,1),NDIM,CONE,E(1,1),NDIM)
      CALL ZCOPY(NDIM*NDIM,CUNIT(1,1),1,D(1,1,N),1)
      CALL ZGETRF(NDIM,NDIM,E(1,1),NDIM,IPVT,INFO)
      CALL ZGETRS('N',NDIM,NDIM,E(1,1),NDIM,IPVT,D(1,1,N),NDIM,INFO)


c ------------------------------------------------------------------------
c
c ---> N = NL
c

 970  N = NL

c
c ---> D(NL) = (A - B(NL-1)*D(NL-1)*C(NL-1))^-1
c
c
c ---> E = D(NL-1) * C(NL-1)
c

      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,NL-1),NDIM,
     +     C(1,1,NL-1),NDIM,CZERO,E(1,1),NDIM)

c
c ---> A = A - B(NL-1) * E
c

      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,B(1,1,NL-1),NDIM,
     +     E(1,1),NDIM,CONE,A(1,1),NDIM)

c
c ---> D(NL) = (A)^-1
c

      CALL ZCOPY(NDIM*NDIM,CUNIT(1,1),1,D(1,1,NL),1)
      CALL ZGETRF(NDIM,NDIM,A(1,1),NDIM,IPVT,INFO)
      CALL ZGETRS('N',NDIM,NDIM,A(1,1),NDIM,IPVT,D(1,1,NL),NDIM,INFO)


 980  CONTINUE                  ! jump label for NL=1

c ------------------------------------------------------------------------
c
c --->  END OF FACTORIZATION
c
c ------------------------------------------------------------------------

c
c ---> HERE IT STARTS LOOP OVER THE DIAGONAL ELEMENT.
c ---> THE PROGRAM CHECKS ID ICHECK(N,N) = 1 AND, IF SO,
c ---> IT CALCULATES THE DIAGONAL BLOCK (N,N)
c ---> THEN IT MAKES TWO LOOPS, ONE OVER A ROW INDEX `IROW`
c ---> AND THE OTHER OVER A COLUMN INDEX `ICOL`, AND CHECKS
c ---> WHICH ARE THE ELEMENTS THAT HAS TO BE CALCULATED.
c ---> (THE ONES FOR WHICH ICHECK = 1)

c
c ---> IT STARTS THE LOOP OVER N
c

      DO 100 N=NL,1,(-1)        ! START OF THE LOOP OVER THE DIAGONAL

         IF (N.EQ.NL) THEN

c     write (6,*) 'it calculates the element ','(',nl,',',nl,')'

c
c ---> GTOT(NL,NL) = D(NL)
c

            CALL ZCOPY(NDIM*NDIM,D(1,1,NL),1,E(1,1),1)

            CALL BTOM(NL,NL,E,NDIM,GIN,ALMD,.FALSE.)


         ELSE IF (N.EQ.NL-1) THEN

            IF (ICHECK(NL-1,NL-1).EQ.1) THEN
c
c ---> GTOT(NL-1,NL-1) = D(NL-1) + D(NL-1)*C(NL-1)*D(NL)*B(NL-1)*D(NL-1)
c

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,B(1,1,NL-1),NDIM,
     +              D(1,1,NL-1),NDIM,CZERO,E(1,1),NDIM)

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,NL),NDIM,
     +              E(1,1),NDIM,CZERO,F(1,1),NDIM)

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,C(1,1,NL-1),NDIM,
     +              F(1,1),NDIM,CZERO,E(1,1),NDIM)

               DO LM = 1,NDIM
                  E(LM,LM) = CONE + E(LM,LM)
               ENDDO

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,NL-1),NDIM,
     +              E(1,1),NDIM,CZERO,F(1,1),NDIM)

               CALL BTOM(NL-1,NL-1,F,NDIM,GIN,ALMD,.FALSE.)

c     write (6,*) 'it calculates the element ','(',nl-1,',',nl-1,')'

            ENDIF

         ELSE

            IF (ICHECK(N,N).EQ.1) THEN

c
c ---> GTOT(N,N) = D(N) + D(N)*( M(N,N+1)*GTOT(N+1,N+1) +
c                  + C(N)*Z(NL,N+1) )*M(N+1,N)*D(N) -
c                  - Z(N,NL)*B(N)*D(N)

               CALL BOFM(NL,N+1,F,NDIM,GIN,ALMD)

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,C(1,1,N),NDIM,
     +              F(1,1),NDIM,CZERO,E(1,1),NDIM)

               CALL BOFM(N+1,N+1,F,NDIM,GIN,ALMD)

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,M1(1,1,N),NDIM,
     +              F(1,1),NDIM,CONE,E(1,1),NDIM)

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,M3(1,1,N),NDIM,
     +              D(1,1,N),NDIM,CZERO,F(1,1),NDIM)

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,E(1,1),NDIM,
     +              F(1,1),NDIM,CZERO,G(1,1),NDIM)

               DO LM = 1,NDIM
                  G(LM,LM) = CONE + G(LM,LM)
               ENDDO

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,N),NDIM,
     +              G(1,1),NDIM,CZERO,E(1,1),NDIM)

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,B(1,1,N),NDIM,
     +              D(1,1,N),NDIM,CZERO,F(1,1),NDIM)

               CALL BOFM(N,NL,G,NDIM,GIN,ALMD)

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,G(1,1),NDIM,
     +              F(1,1),NDIM,CONE,E(1,1),NDIM)

               CALL BTOM(N,N,E,NDIM,GIN,ALMD,.FALSE.)

c     write (6,*) 'it calculates the element ','(',n,',',n,')'

            ENDIF

         ENDIF


         DO IROW = (N-1),1,(-1) ! LOOP OVER THE ROW FOR THE COLUMN N

            IF (ICHECK(IROW,N).EQ.1) THEN

               IF ((N.EQ.NL).AND.(IROW.EQ.(NL-1))) THEN

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,NL-1),NDIM,
     +              C(1,1,NL-1),NDIM,CZERO,F(1,1),NDIM)

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,F(1,1),NDIM,
     +              D(1,1,NL),NDIM,CZERO,G(1,1),NDIM)

               CALL BTOM(NL-1,NL,G,NDIM,GIN,ALMD,.FALSE.)

               ELSE

c     M(I,I+1) * Z(I+1,J)

               CALL BOFM(IROW+1,N,E,NDIM,GIN,ALMD)

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,M1(1,1,IROW),NDIM,
     +              E(1,1),NDIM,CZERO,G(1,1),NDIM)

c     M(I,I+1) * Z(I+1,J) + C(I) * Z(N,J)

               CALL BOFM(NL,N,E,NDIM,GIN,ALMD)

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,C(1,1,IROW),NDIM,
     +              E(1,1),NDIM,CONE,G(1,1),NDIM)

c     -D(I) * ( M(I,I+1)*Z(I+1,J)+C(I)*Z(N,J) )

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,D(1,1,IROW),NDIM,
     +              G(1,1),NDIM,CZERO,E(1,1),NDIM)

               CALL BTOM(IROW,N,E,NDIM,GIN,ALMD,.FALSE.)

            ENDIF

c     write (6,*) 'it calculates the element ','(',irow,',',n,')'

         ENDIF

      ENDDO                     ! LOOP OVER THE ROW FOR THE COLUMN N


      DO ICOL = (N-1),1,(-1)    ! LOOP OVER THE COLUMN FOR THE ROW N

         IF (ICHECK(N,ICOL).EQ.1) THEN

            IF ((N.EQ.NL).AND.(ICOL.EQ.(NL-1))) THEN

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,B(1,1,NL-1),NDIM,
     +              D(1,1,NL-1),NDIM,CZERO,E(1,1),NDIM)

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,D(1,1,NL),NDIM,
     +              E(1,1),NDIM,CZERO,G(1,1),NDIM)

               CALL BTOM(NL,NL-1,G,NDIM,GIN,ALMD,.FALSE.)

            ELSE

c     Z(I,J+1) * M(J+1,J)

               CALL BOFM(N,ICOL+1,E,NDIM,GIN,ALMD)

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,E(1,1),NDIM,
     +              M3(1,1,ICOL),NDIM,CZERO,G(1,1),NDIM)

c     Z(I,J+1) * M(J+1,J) + Z(I,N) * B(J)

               CALL BOFM(N,NL,E,NDIM,GIN,ALMD)

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,E(1,1),NDIM,
     +              B(1,1,ICOL),NDIM,CONE,G(1,1),NDIM)

c     -( Z(I,J+1) * M(J+1,J)+Z(I,N) * B(J) ) * D(J)

               CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,G(1,1),NDIM,
     +              D(1,1,ICOL),NDIM,CZERO,E(1,1),NDIM)

               CALL BTOM(N,ICOL,E,NDIM,GIN,ALMD,.FALSE.)

            ENDIF

c     write (6,*) 'it calculates the element ','(',n,',',icol,')'

         ENDIF

      ENDDO                     ! LOOP OVER THE COLUMN FOR THE ROW N


 100  CONTINUE                  ! END OF THE LOOP OVER THE DIAGONAL


C
      RETURN
C
      END
