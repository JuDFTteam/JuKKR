C ************************************************************************
      SUBROUTINE SPATPR(A,B,C,V)
C ************************************************************************
C SPATPR COMPUTES THE SPATIAL PRODUCT OF THREE VECTORS A,B AND C
C RETURNING IT INTO V: V=AXB.C.
C ------------------------------------------------------------------------
      REAL*8        A(*), B(*), C(*), V
      V=0.0D0
      V=V+C(1)*(A(2)*B(3)-A(3)*B(2))
      V=V+C(2)*(A(3)*B(1)-A(1)*B(3))
      V=V+C(3)*(A(1)*B(2)-A(2)*B(1))
      RETURN
      END
