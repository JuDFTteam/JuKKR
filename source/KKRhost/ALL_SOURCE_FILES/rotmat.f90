    Subroutine rotmat(iopt, li, nrot, symopm, vecg)
      Use mod_datatypes, Only: dp
!- Converts rotation/rotoinversion matrix <-> (nrot,vecg,li)
! ----------------------------------------------------------------------
!i Inputs:
!i   iopt  := -1 to convert (nrot,vecg,li) in symopm
!i          =  1 to convert symopm in (nrot,vecg,li)
!i Inputs/Outputs:
!io  li    :if T: inversion or rotoinversion
!io  nrot  :rotation angle = 2*pi/nrot
!io  symopm:symmetry operation matrix
!io  vecg  :rotation axis
! ----------------------------------------------------------------------
      Implicit None
! Passed parameters:                                                    
      Integer :: iopt, nrot
      Real (Kind=dp) :: vecg(3), symopm(3, 3)
      Logical :: li
! Local parameters:                                                     
      Integer :: i, idamax, in, j
      Real (Kind=dp) :: costbn, detop, ddet33, dnrm2, omcos, sintbn, sinpb3, &
        tiny, twopi, vfac
      Character (Len=144) :: messg
      Parameter (twopi=6.28318530717958648E0_dp)
      Parameter (tiny=1.0E-3_dp)
! External calls:                                                       
      External :: daxpy, dcopy, ddet33, rinit, dnrm2, dscal, errmsg, idamax, &
        nrmliz
! Intrinsic functions:                                                  
      Intrinsic :: abs, acos, cos, max, sign, sin, sqrt, abs, nint

      If (iopt==-1) Then
        Call rinit(symopm, 9)
        in = abs(nrot)
        If (in==1) Then
          Call dcopy(3, 1.E0_dp, 0, symopm, 4)
        Else If (in==2 .Or. in==3 .Or. in==4 .Or. in==6) Then
          sintbn = sin(twopi/nrot)
          costbn = cos(twopi/nrot)
          omcos = 1.E0_dp - costbn
          If (dnrm2(3,vecg,1)<tiny) Call errmsg( &
            ' ROTMAT: zero rotation vector.$', 4)
          Call nrmliz(1, vecg, vecg)
          symopm(1, 1) = omcos*vecg(1)*vecg(1) + costbn
          symopm(1, 2) = omcos*vecg(1)*vecg(2) - sintbn*vecg(3)
          symopm(1, 3) = omcos*vecg(1)*vecg(3) + sintbn*vecg(2)
          symopm(2, 1) = omcos*vecg(2)*vecg(1) + sintbn*vecg(3)
          symopm(2, 2) = omcos*vecg(2)*vecg(2) + costbn
          symopm(2, 3) = omcos*vecg(2)*vecg(3) - sintbn*vecg(1)
          symopm(3, 1) = omcos*vecg(3)*vecg(1) - sintbn*vecg(2)
          symopm(3, 2) = omcos*vecg(3)*vecg(2) + sintbn*vecg(1)
          symopm(3, 3) = omcos*vecg(3)*vecg(3) + costbn
        Else
          Call errmsg(' ROTMAT: bad nrot.$', 3)
        End If
        If (li) Call dscal(9, -1.E0_dp, symopm(1,1), 1)

      Else If (iopt==1) Then
! ----- First calculate determinant.
        detop = ddet33(symopm)
        If (abs(abs(detop)-1.0E0_dp)>tiny) Call errmsg( &
          ' ROTMAT: determinant is not +/- 1$', 4)
        detop = sign(1.E0_dp, detop)
        li = detop < 0.E0_dp
! ----- multiply operation symopm with detop
        Call dscal(9, detop, symopm(1,1), 1)
! ----- For the rotation angle we have due to the normalization of v:
! ----- sum_i symopm(i,i) = sum_i (1-cos) v_i*v_i+3*cos = 1 + 2 * cos,
        costbn = -0.5E0_dp
        Call daxpy(3, 0.5E0_dp, symopm(1,1), 4, costbn, 0)
        If (abs(costbn-1.E0_dp)<tiny) Then
          nrot = 1
          Call rinit(vecg, 3)
        Else
          nrot = nint(twopi/acos(max(-1.E0_dp,costbn)))
! ------- for nrot > 2 the matrix is non-symmetric and the rotation
! ------- axis can be calculated from the antisymmetric part.
! ------- for nrot = 2 this not possible. However, the squared vector
! ------- components are given by:  mat(i,i) = 2 v_i * v_i - 1.
! ------- This is used for the largest component. The others are taken
! ------- from: mat(i,j) = 2 v_i * v_j for i ne j. This way we also
! ------- get the right phases between the components.
          If (nrot==2) Then
            Do i = 1, 3
              vecg(i) = 0.5E0_dp*(symopm(i,i)+1.0E0_dp)
            End Do
            j = idamax(3, vecg, 1)
            If (vecg(j)<0.0E0_dp) Then
              Write (messg, 100) j, symopm(j, j)
              Call errmsg(messg, 4)
            End If
            vecg(j) = sqrt(vecg(j))
            vfac = 0.5E0_dp/vecg(j)
            Do i = 1, 3
              If (i/=j) vecg(i) = vfac*symopm(i, j)
            End Do
          Else
            vecg(1) = symopm(3, 2) - symopm(2, 3)
            vecg(2) = symopm(1, 3) - symopm(3, 1)
            vecg(3) = symopm(2, 1) - symopm(1, 2)
          End If
! ------- next renormalize at least one component to 1 in order to
! ------- allow for abbreviations as 'D', 'X', 'Y' or 'Z'
          sinpb3 = sqrt(.75E0_dp)
          If (abs((sinpb3-abs(vecg(1)))*(sinpb3-abs(vecg(2)))*(sinpb3- &
            abs(vecg(3))))>tiny) Then
            Do j = 3, 1, -1
              vfac = abs(vecg(j))
              If (vfac>tiny) Call dscal(3, 1.E0_dp/vfac, vecg, 1)
            End Do
          End If
        End If
        Call dscal(9, detop, symopm(1,1), 1)
      End If

100   Format (' ROTMAT: Bad component ', I1, ' of operation ', &
        '. Diagonal element =', F9.5, '$')
    End Subroutine
