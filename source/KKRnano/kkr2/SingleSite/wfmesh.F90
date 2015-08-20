      subroutine wfmesh(e, ek, cvlight, nsra, z, r, s, rs, irm, irmd, lmaxd)
      double complex, intent(in) :: e
      double complex, intent(out) :: ek
      double precision, intent(in) :: cvlight, z
      integer, intent(in) :: irm, irmd, lmaxd, nsra
      double precision, intent(in) :: r(irmd)
      double precision, intent(out) :: rs(irmd,0:lmaxd), s(0:lmaxd)

      double precision :: s1
      integer :: ir,l

      if (nsra == 2) then
        ek = sqrt(e+e*e/ (cvlight*cvlight))
      else
        ! assume(nsra == 1)
        ek = sqrt(e)
      endif
        
      do l = 0, lmaxd

        if (nsra == 2) then
          s1 = sqrt(dble(l*l+l+1) - 4.d0*z*z/(cvlight*cvlight))
          if (z == 0.d0) s1 = dble(l)
        else
          s1 = dble(l)
        end if
        s(l) = s1
        rs(1,l) = 0.d0
        do ir = 2, irm
          rs(ir,l) = r(ir)**s1
        enddo ! ir
        do ir = irm+1, irmd
           rs(ir,l) = 0.d0
        enddo ! ir

      enddo ! l
      
      endsubroutine ! wfmesh
