! ************************************************************************
    Subroutine inversion(gllke, invmod, icheck)
      Use mod_datatypes, Only: dp
! ************************************************************************
! This subroutine calculates the inversion of a matrix
! in 4 different ways depending on the form of the matrix

!     INVMOD = 0  ----> total inversion scheme
!     INVMOD = 1  ----> band matrix inversion scheme
!     INVMOD = 2  ----> corner band matrix inversion scheme
!     INVMOD = 3  ----> sparse matrix inversion scheme

! ------------------------------------------------------------------------
      Implicit None

!     .. parameters ..
      Include 'inc.p'
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************


! changed 3.11.99

      Integer :: lmmaxd
      Parameter (lmmaxd=(krel+korbit+1)*(lmaxd+1)**2)
      Integer :: almd, ndim
      Parameter (almd=naezd*lmmaxd, ndim=nprincd*lmmaxd)
      Complex (Kind=dp) :: ci, czero, cone
      Parameter (ci=(0.E0_dp,1.E0_dp), czero=(0.E0_dp,0.E0_dp), &
        cone=(1.E0_dp,0.E0_dp))

      Complex (Kind=dp) :: gllke(almd, almd), gdi(ndim, ndim, nlayerd), &
        gup(ndim, ndim, nlayerd), gdow(ndim, ndim, nlayerd)
      Complex (Kind=dp), Allocatable :: gtemp(:, :)
      Integer :: i, i1, ip1, ii1, il1, ldi1, ip2, ii2, il2, ldi2, j, invmod
      Integer :: lm1, lm2, info, ipvt(almd), nlayer
      Integer :: icheck(naezd/nprincd, naezd/nprincd)
! total matrix inversion
      External :: zgetrf, zgetrs, zcopy, invslab

      Allocate (gtemp(almd,almd))

      nlayer = naezd/nprincd

      If (invmod==0) Then

!     write (6,*) '-------full inversion calculation--------'

        Do i = 1, almd
          Do j = 1, almd
            gtemp(i, j) = czero
            If (i==j) Then
              gtemp(i, j) = cone
            End If
          End Do
        End Do




        Call zgetrf(almd, almd, gllke, almd, ipvt, info)
        Call zgetrs('N', almd, almd, gllke, almd, ipvt, gtemp, almd, info)
! slab or supercell
        Call zcopy(almd*almd, gtemp, 1, gllke, 1)
! inversion


      Else If ((invmod>=1) .And. (invmod<=2)) Then
!     this part now is correct also for    ! changes 20/10/99
!     supercell geometry : 20/10/99
!---> upper linear part
        Do i1 = 1, nlayerd
          Do ip1 = 1, nprincd
            Do ip2 = 1, nprincd
              ii1 = (i1-1)*nprincd + ip1
              ii2 = (i1-1)*nprincd + ip2
              Do lm1 = 1, lmmaxd
                Do lm2 = 1, lmmaxd
                  ldi1 = lmmaxd*(ip1-1) + lm1
                  il1 = lmmaxd*(ii1-1) + lm1
                  ldi2 = lmmaxd*(ip2-1) + lm2
                  il2 = lmmaxd*(ii2-1) + lm2
                  gdi(ldi1, ldi2, i1) = gllke(il1, il2)
                End Do
              End Do
            End Do
          End Do
        End Do


!---> lower linear part

        Do i1 = 1, nlayerd
          Do ip1 = 1, nprincd
            Do ip2 = 1, nprincd
              Do lm1 = 1, lmmaxd
                Do lm2 = 1, lmmaxd
                  ldi1 = lmmaxd*(ip1-1) + lm1
                  ldi2 = lmmaxd*(ip2-1) + lm2
                  If (i1<=(nlayerd-1)) Then
                    ii1 = (i1-1)*nprincd + ip1
                    ii2 = i1*nprincd + ip2
                    il1 = lmmaxd*(ii1-1) + lm1
                    il2 = lmmaxd*(ii2-1) + lm2
                    gup(ldi1, ldi2, i1) = gllke(il1, il2)
                  Else
                    ii1 = ip1
                    ii2 = (nlayerd-1)*nprincd + ip2
                    il1 = lmmaxd*(ii1-1) + lm1
                    il2 = lmmaxd*(ii2-1) + lm2
                    gdow(ldi1, ldi2, i1) = gllke(il1, il2)
                  End If
                End Do
              End Do
            End Do
          End Do
        End Do
!     end of the corrected part  20/10/99


        Do i1 = 1, nlayerd
          Do ip1 = 1, nprincd
            Do ip2 = 1, nprincd
              Do lm1 = 1, lmmaxd
                Do lm2 = 1, lmmaxd
                  ldi1 = lmmaxd*(ip1-1) + lm1
                  ldi2 = lmmaxd*(ip2-1) + lm2
                  If (i1<=(nlayerd-1)) Then
                    ii1 = i1*nprincd + ip1
                    ii2 = (i1-1)*nprincd + ip2
                    il1 = lmmaxd*(ii1-1) + lm1
                    il2 = lmmaxd*(ii2-1) + lm2
                    gdow(ldi1, ldi2, i1) = gllke(il1, il2)
                  Else
                    ii1 = (nlayerd-1)*nprincd + ip1
                    ii2 = ip2
                    il1 = lmmaxd*(ii1-1) + lm1
                    il2 = lmmaxd*(ii2-1) + lm2
                    gup(ldi1, ldi2, i1) = gllke(il1, il2)
                  End If
                End Do
              End Do
            End Do
          End Do
        End Do

!          write (6,*) '-------slab calculation--------'

        If (invmod==1) Then
! supercell geometry inversion
          Call invslab(gdi, gup, gdow, gllke, icheck)


!          write (6,*) '-------supercell calculation--------'
        Else If (invmod==2) Then

          Call invsupercell(gdi, gup, gdow, gllke, icheck)

! sparse matrix inversion

        End If
!     NOT YET IMPLEMENTED!!!!!!!!!

      Else




      End If

! ************************************************************************
      Deallocate (gtemp)
! ************************************************************************
! This subroutine calculates the inversion of a matrix
      Return
! in 4 different ways depending on the form of the matrix
    End Subroutine
