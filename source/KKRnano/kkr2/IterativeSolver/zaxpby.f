C
C**********************************************************************
**********************************************************************
C
C     Copyright (C) 1992  Roland W. Freund and Noel M. Nachtigal
C     All rights reserved.
C
C     This code is part of a copyrighted package.  For details, see the
C     file `cpyrit.doc' in the top-level directory.
C
C     *****************************************************************
C     ANY USE OF  THIS CODE CONSTITUTES ACCEPTANCE OF  THE TERMS OF THE
C                             COPYRIGHT NOTICE
C     *****************************************************************
C
C**********************************************************************
C
      subroutine zaxpby (n,zz,za,zx,zb,zy)
c
c     purpose:
c     this subroutine computes zz = za * zx + zb * zy.  several special
c     cases are handled separately:
c        za =  0.0, zb =  0.0 => zz = 0.0
c        za =  0.0, zb =  1.0 => zz = zy  (this is copy)
c        za =  0.0, zb = -1.0 => zz = -zy
c        za =  0.0, zb =   zb => zz = zb * zy  (this is scal)
c        za =  1.0, zb =  0.0 => zz = zx  (this is copy)
c        za =  1.0, zb =  1.0 => zz = zx + zy
c        za =  1.0, zb = -1.0 => zz = zx - zy
c        za =  1.0, zb =   zb => zz = zx + zb * zy (this is axpy)
c        za = -1.0, zb =  0.0 => zz = -zx
c        za = -1.0, zb =  1.0 => zz = -zx + zy
c        za = -1.0, zb = -1.0 => zz = -zx - zy
c        za = -1.0, zb =   zb => zz = -zx + zb * zy
c        za =   za, zb =  0.0 => zz = za * zx  (this is scal)
c        za =   za, zb =  1.0 => zz = za * zx + zy  (this is axpy)
c        za =   za, zb = -1.0 => zz = za * zx - zy
c        za =   za, zb =   zb => zz = za * zx + zb * zy
c     zz may be the same as zx or zy.
c
c     parameters:
c     n  = the dimension of the vectors (input).
c     zz = the vector result (output).
c     za = scalar multiplier for zx (input).
c     zx = one of the vectors (input).
c     zb = scalar multiplier for zy (input).
c     zy = the other vector (input).
c
c     noel m. nachtigal
c     march 23, 1993
c
c**********************************************************************
c
      intrinsic dimag, dreal
c
      integer n
      double complex za, zb, zx(n), zy(n), zz(n)
c
c     local variables.
c
      integer i
      double precision dai, dar, dbi, dbr
c
      if (n.le.0) return
c
      dai = dimag(za)
      dar = dreal(za)
      dbi = dimag(zb)
      dbr = dreal(zb)
      if ((dar == 0.d0).and.(dai == 0.d0)) then
         if ((dbr == 0.d0).and.(dbi == 0.d0)) then
c           za = 0.0, zb = 0.0 => zz = 0.0.
            do 10 i = 1, n
               zz(i) = (0.d0,0.d0)
 10         continue
         else if ((dbr == 1.d0).and.(dbi == 0.d0)) then
c           za = 0.0, zb = 1.0 => zz = zy (this is copy).
            do 20 i = 1, n
               zz(i) = zy(i)
 20         continue
         else if ((dbr == -1.d0).and.(dbi == 0.d0)) then
c           za = 0.0, zb = -1.0 => zz = -zy.
            do 30 i = 1, n
               zz(i) = -zy(i)
 30         continue
         else
c           za = 0.0, zb = zb => zz = zb * zy (this is scal).
            do 40 i = 1, n
               zz(i) = zb * zy(i)
 40         continue
         end if
      else if ((dar == 1.d0).and.(dai == 0.d0)) then
         if ((dbr == 0.d0).and.(dbi == 0.d0)) then
c           za = 1.0, zb = 0.0 => zz = zx (this is copy).
            do 50 i = 1, n
               zz(i) = zx(i)
 50         continue
         else if ((dbr == 1.d0).and.(dbi == 0.d0)) then
c           za = 1.0, zb = 1.0 => zz = zx + zy.
            do 60 i = 1, n
               zz(i) = zx(i) + zy(i)
 60         continue
         else if ((dbr == -1.d0).and.(dbi == 0.d0)) then
c           za = 1.0, zb = -1.0 => zz = zx - zy.
            do 70 i = 1, n
               zz(i) = zx(i) - zy(i)
 70         continue
         else
c           za = 1.0, zb = zb => zz = zx + zb * zy (this is axpy).
            do 80 i = 1, n
               zz(i) = zx(i) + zb * zy(i)
 80         continue
         end if
      else if ((dar == -1.d0).and.(dai == 0.d0)) then
         if ((dbr == 0.d0).and.(dbi == 0.d0)) then
c           za = -1.0, zb = 0.0 => zz = -zx
            do 90 i = 1, n
               zz(i) = -zx(i)
 90         continue
         else if ((dbr == 1.d0).and.(dbi == 0.d0)) then
c           za = -1.0, zb = 1.0 => zz = -zx + zy
            do 100 i = 1, n
               zz(i) = -zx(i) + zy(i)
 100        continue
         else if ((dbr == -1.d0).and.(dbi == 0.d0)) then
c           za = -1.0, zb = -1.0 => zz = -zx - zy.
            do 110 i = 1, n
               zz(i) = -zx(i) - zy(i)
 110        continue
         else
c           za = -1.0, zb = zb => zz = -zx + zb * zy
            do 120 i = 1, n
               zz(i) = -zx(i) + zb * zy(i)
 120        continue
         end if
      else
         if ((dbr == 0.d0).and.(dbi == 0.d0)) then
c           za = za, zb = 0.0 => zz = za * zx (this is scal).
            do 130 i = 1, n
               zz(i) = za * zx(i)
 130        continue
         else if ((dbr == 1.d0).and.(dbi == 0.d0)) then
c           za = za, zb = 1.0 => zz = za * zx + zy (this is axpy)
            do 140 i = 1, n
               zz(i) = za * zx(i) + zy(i)
 140        continue
         else if ((dbr == -1.d0).and.(dbi == 0.d0)) then
c           za = za, zb = -1.0 => zz = za * zx - zy.
            do 150 i = 1, n
               zz(i) = za * zx(i) - zy(i)
 150        continue
         else
c           za = za, zb = zb => zz = za * zx + zb * zy.
            do 160 i = 1, n
               zz(i) = za * zx(i) + zb * zy(i)
 160        continue
         end if
      end if
c
      return

      end
