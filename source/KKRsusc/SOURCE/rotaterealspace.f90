module rotaterealspace 

integer,parameter :: verbose = 0

contains

subroutine rotaterealspace_vector(vec1,theta,phi,mode1)
use mod_mathtools
implicit none
double precision :: vec1(3)
double precision :: theta,phi
character(len=*) :: mode1
!local
double precision :: rotmat(3,3),vec2(3)


if (mode1=='loc->glob') then
  rotmat=rotaterealspace_createrotmatrix(theta,phi)
elseif (mode1=='glob->loc') then
  rotmat=rotaterealspace_createrotmatrix(theta,phi)
  call transpose_dm(rotmat)
else
  stop '[rotaterealspace] mode1 error'
end if

vec2=matvec_dmdm(rotmat,vec1)

vec1=vec2

end subroutine

subroutine rotaterealspace_matrix(mat1,theta,phi,mode1)
use mod_mathtools
implicit none
double precision :: mat1(3,3)
double precision :: theta,phi
character(len=*) :: mode1

!local
double precision :: rotmat(3,3),mat2(3,3)

if (mode1=='loc->glob') then
  rotmat=rotaterealspace_createrotmatrix(theta,phi)
elseif (mode1=='glob->loc') then
  rotmat=rotaterealspace_createrotmatrix(theta,phi)
  call transpose_dm(rotmat)
else
  stop '[rotaterealspace] mode1 error'
end if

mat2=matmat1T_dmdm(mat1,rotmat)
mat1=matmat_dmdm(rotmat,mat2)

end subroutine

subroutine rotaterealspace_matrix2(mat1,theta,phi,theta2,phi2,mode1)
use mod_mathtools
implicit none
double precision :: mat1(3,3)
double precision :: theta,phi,theta2,phi2
character(len=*) :: mode1

!local
double precision :: rotmat(3,3),mat2(3,3)
double precision :: rotmat2(3,3)

if (mode1=='loc->glob') then
  rotmat = rotaterealspace_createrotmatrix(theta,phi)
  rotmat2= rotaterealspace_createrotmatrix(theta2,phi2)

elseif (mode1=='glob->loc') then
  rotmat = rotaterealspace_createrotmatrix(theta,phi)
  rotmat2= rotaterealspace_createrotmatrix(theta2,phi2)

  call transpose_dm(rotmat)
  call transpose_dm(rotmat2)
else
  stop '[rotaterealspace] mode1 error'
end if

mat2=matmat1T_dmdm(mat1,rotmat2)
mat1=matmat_dmdm(rotmat,mat2)


end subroutine

function rotaterealspace_createrotmatrix(theta,phi)
implicit none
double precision :: theta,phi
!local
double precision :: rotaterealspace_createrotmatrix(3,3)

 rotaterealspace_createrotmatrix(1,1) =  cos(phi)* cos(theta)
 rotaterealspace_createrotmatrix(1,2) = -sin(phi)
 rotaterealspace_createrotmatrix(1,3) =  cos(phi)* sin(theta)
 rotaterealspace_createrotmatrix(2,1) =  sin(phi)* cos(theta)
 rotaterealspace_createrotmatrix(2,2) =  cos(phi)
 rotaterealspace_createrotmatrix(2,3) =  sin(phi)* sin(theta)
 rotaterealspace_createrotmatrix(3,1) =             -sin(theta)
 rotaterealspace_createrotmatrix(3,2) =  0.0D0
 rotaterealspace_createrotmatrix(3,3) =              cos(theta)

if (verbose=='1') then
  print *,'-----------------------------------------------------------------'
  print *, rotaterealspace_createrotmatrix(1,1),rotaterealspace_createrotmatrix(1,2),rotaterealspace_createrotmatrix(1,3)
  print *, rotaterealspace_createrotmatrix(2,1),rotaterealspace_createrotmatrix(2,2),rotaterealspace_createrotmatrix(2,3)
  print *, rotaterealspace_createrotmatrix(3,1),rotaterealspace_createrotmatrix(3,2),rotaterealspace_createrotmatrix(3,3)
  print *,'-----------------------------------------------------------------'
end if

end function rotaterealspace_createrotmatrix



function rotaterealspace_createrotmatrix_nvec(alpha,nvec)
implicit none
double precision :: alpha
double precision :: nvec(3)
!local
double precision :: rotaterealspace_createrotmatrix_nvec(3,3)
double precision :: cosa,sina

 cosa = cos(alpha)
 sina = sin(alpha)

 rotaterealspace_createrotmatrix_nvec(1,1) = nvec(1)*nvec(1)*(1.0D0-cosa) +         cosa
 rotaterealspace_createrotmatrix_nvec(1,2) = nvec(1)*nvec(2)*(1.0D0-cosa) - nvec(3)*sina
 rotaterealspace_createrotmatrix_nvec(1,3) = nvec(1)*nvec(3)*(1.0D0-cosa) + nvec(2)*sina
 rotaterealspace_createrotmatrix_nvec(2,1) = nvec(2)*nvec(1)*(1.0D0-cosa) + nvec(3)*sina
 rotaterealspace_createrotmatrix_nvec(2,2) = nvec(2)*nvec(2)*(1.0D0-cosa) +         cosa
 rotaterealspace_createrotmatrix_nvec(2,3) = nvec(2)*nvec(3)*(1.0D0-cosa) - nvec(1)*sina
 rotaterealspace_createrotmatrix_nvec(3,1) = nvec(3)*nvec(1)*(1.0D0-cosa) - nvec(2)*sina
 rotaterealspace_createrotmatrix_nvec(3,2) = nvec(3)*nvec(2)*(1.0D0-cosa) + nvec(1)*sina
 rotaterealspace_createrotmatrix_nvec(3,3) = nvec(3)*nvec(3)*(1.0D0-cosa) +         cosa

if (verbose=='1') then
  print *,'-----------------------------------------------------------------'
  print *, rotaterealspace_createrotmatrix_nvec(1,1),rotaterealspace_createrotmatrix_nvec(1,2),rotaterealspace_createrotmatrix_nvec(1,3)
  print *, rotaterealspace_createrotmatrix_nvec(2,1),rotaterealspace_createrotmatrix_nvec(2,2),rotaterealspace_createrotmatrix_nvec(2,3)
  print *, rotaterealspace_createrotmatrix_nvec(3,1),rotaterealspace_createrotmatrix_nvec(3,2),rotaterealspace_createrotmatrix_nvec(3,3)
  print *,'-----------------------------------------------------------------'
end if

end function rotaterealspace_createrotmatrix_nvec

function rotaterealspace_createdrotmatrixdalpha_nvec(alpha,nvec)
implicit none
double precision :: alpha
double precision :: nvec(3)
!local
double precision :: rotaterealspace_createdrotmatrixdalpha_nvec(3,3)
double precision :: cosa,sina

 cosa = cos(alpha)
 sina = sin(alpha)

 rotaterealspace_createdrotmatrixdalpha_nvec(1,1) = nvec(1)*nvec(1)*(1.0D0+sina) -         sina
 rotaterealspace_createdrotmatrixdalpha_nvec(1,2) = nvec(1)*nvec(2)*(1.0D0+sina) - nvec(3)*cosa
 rotaterealspace_createdrotmatrixdalpha_nvec(1,3) = nvec(1)*nvec(3)*(1.0D0+sina) + nvec(2)*cosa
 rotaterealspace_createdrotmatrixdalpha_nvec(2,1) = nvec(2)*nvec(1)*(1.0D0+sina) + nvec(3)*cosa
 rotaterealspace_createdrotmatrixdalpha_nvec(2,2) = nvec(2)*nvec(2)*(1.0D0+sina) -         sina
 rotaterealspace_createdrotmatrixdalpha_nvec(2,3) = nvec(2)*nvec(3)*(1.0D0+sina) - nvec(1)*cosa
 rotaterealspace_createdrotmatrixdalpha_nvec(3,1) = nvec(3)*nvec(1)*(1.0D0+sina) - nvec(2)*cosa
 rotaterealspace_createdrotmatrixdalpha_nvec(3,2) = nvec(3)*nvec(2)*(1.0D0+sina) + nvec(1)*cosa
 rotaterealspace_createdrotmatrixdalpha_nvec(3,3) = nvec(3)*nvec(3)*(1.0D0+sina) -         sina

if (verbose=='1') then
  print *,'-----------------------------------------------------------------'
  print *, rotaterealspace_createdrotmatrixdalpha_nvec(1,1),rotaterealspace_createdrotmatrixdalpha_nvec(1,2),rotaterealspace_createdrotmatrixdalpha_nvec(1,3)
  print *, rotaterealspace_createdrotmatrixdalpha_nvec(2,1),rotaterealspace_createdrotmatrixdalpha_nvec(2,2),rotaterealspace_createdrotmatrixdalpha_nvec(2,3)
  print *, rotaterealspace_createdrotmatrixdalpha_nvec(3,1),rotaterealspace_createdrotmatrixdalpha_nvec(3,2),rotaterealspace_createdrotmatrixdalpha_nvec(3,3)
  print *,'-----------------------------------------------------------------'
end if

end function rotaterealspace_createdrotmatrixdalpha_nvec


function rotaterealspace_created2rotmatrixdalpha2_nvec(alpha,nvec)
implicit none
double precision :: alpha
double precision :: nvec(3)
double precision :: rotaterealspace_created2rotmatrixdalpha2_nvec(3,3)

rotaterealspace_created2rotmatrixdalpha2_nvec = -1.0D0 * rotaterealspace_createrotmatrix_nvec(alpha,nvec)

end function rotaterealspace_created2rotmatrixdalpha2_nvec


end module rotaterealspace