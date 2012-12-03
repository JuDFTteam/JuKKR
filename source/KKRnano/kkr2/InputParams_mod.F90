!> Module that defines a datastructure that contains input parameters.
!> TODO: add OPTC TESTC ??

module InputParams_mod

  type InputParams
    double precision  :: RCUTJIJ
    double precision  :: ALAT
    double precision  :: FCM
    double precision  :: MIXING
    double precision  :: RMAX
    double precision  :: GMAX
    integer  :: SCFSTEPS
    integer  :: IMIX
    integer  :: KPRE
    integer  :: KTE
    integer  :: KXC
    integer  :: KFORCE
    logical  :: JIJ
    logical  :: LDAU
    double precision  :: QMRBOUND
    integer  :: ICST
    integer  :: NSRA
    ! TODO: get rid of the following parameters - not really input params
    integer  :: NCLS
    integer  :: NREF
    integer  :: IELAST
    integer  :: NSYMAT
    integer  :: MAXMESH
    integer  :: NR

  end type InputParams

  CONTAINS
!-------------------------------------------------------------------------------
!> Returns RCUTJIJ.
!> @param[in] self InputParams object
double precision  function getRCUTJIJ(self)
  implicit none
  type(InputParams), intent(in) :: self

  getRcutjij = self%RCUTJIJ
end function

!-------------------------------------------------------------------------------
!> Returns ALAT.
!> @param[in] self InputParams object
double precision  function getALAT(self)
  implicit none
  type(InputParams), intent(in) :: self

  getAlat = self%ALAT
end function

!-------------------------------------------------------------------------------
!> Returns FCM.
!> @param[in] self InputParams object
double precision  function getFCM(self)
  implicit none
  type(InputParams), intent(in) :: self

  getFcm = self%FCM
end function

!-------------------------------------------------------------------------------
!> Returns MIXING.
!> @param[in] self InputParams object
double precision  function getMIXING(self)
  implicit none
  type(InputParams), intent(in) :: self

  getMixing = self%MIXING
end function

!-------------------------------------------------------------------------------
!> Returns RMAX.
!> @param[in] self InputParams object
double precision  function getRMAX(self)
  implicit none
  type(InputParams), intent(in) :: self

  getRmax = self%RMAX
end function

!-------------------------------------------------------------------------------
!> Returns GMAX.
!> @param[in] self InputParams object
double precision  function getGMAX(self)
  implicit none
  type(InputParams), intent(in) :: self

  getGmax = self%GMAX
end function

!-------------------------------------------------------------------------------
!> Returns SCFSTEPS.
!> @param[in] self InputParams object
integer  function getSCFSTEPS(self)
  implicit none
  type(InputParams), intent(in) :: self

  getScfsteps = self%SCFSTEPS
end function

!-------------------------------------------------------------------------------
!> Returns IMIX.
!> @param[in] self InputParams object
integer  function getIMIX(self)
  implicit none
  type(InputParams), intent(in) :: self

  getImix = self%IMIX
end function

!-------------------------------------------------------------------------------
!> Returns KPRE.
!> @param[in] self InputParams object
integer  function getKPRE(self)
  implicit none
  type(InputParams), intent(in) :: self

  getKpre = self%KPRE
end function

!-------------------------------------------------------------------------------
!> Returns KTE.
!> @param[in] self InputParams object
integer  function getKTE(self)
  implicit none
  type(InputParams), intent(in) :: self

  getKte = self%KTE
end function

!-------------------------------------------------------------------------------
!> Returns KXC.
!> @param[in] self InputParams object
integer  function getKXC(self)
  implicit none
  type(InputParams), intent(in) :: self

  getKxc = self%KXC
end function

!-------------------------------------------------------------------------------
!> Returns KFORCE.
!> @param[in] self InputParams object
integer  function getKFORCE(self)
  implicit none
  type(InputParams), intent(in) :: self

  getKforce = self%KFORCE
end function

!-------------------------------------------------------------------------------
!> Returns JIJ.
!> @param[in] self InputParams object
logical  function getJIJ(self)
  implicit none
  type(InputParams), intent(in) :: self

  getJij = self%JIJ
end function

!-------------------------------------------------------------------------------
!> Returns LDAU.
!> @param[in] self InputParams object
logical  function getLDAU(self)
  implicit none
  type(InputParams), intent(in) :: self

  getLdau = self%LDAU
end function

!-------------------------------------------------------------------------------
!> Returns QMRBOUND.
!> @param[in] self InputParams object
double precision  function getQMRBOUND(self)
  implicit none
  type(InputParams), intent(in) :: self

  getQmrbound = self%QMRBOUND
end function

!-------------------------------------------------------------------------------
!> Returns ICST.
!> @param[in] self InputParams object
integer  function getICST(self)
  implicit none
  type(InputParams), intent(in) :: self

  getIcst = self%ICST
end function

!-------------------------------------------------------------------------------
!> Returns NSRA.
!> @param[in] self InputParams object
integer  function getNSRA(self)
  implicit none
  type(InputParams), intent(in) :: self

  getNsra = self%NSRA
end function

!-------------------------------------------------------------------------------
!> Returns NCLS.
!> @param[in] self InputParams object
integer  function getNCLS(self)
  implicit none
  type(InputParams), intent(in) :: self

  getNcls = self%NCLS
end function

!-------------------------------------------------------------------------------
!> Returns NREF.
!> @param[in] self InputParams object
integer  function getNREF(self)
  implicit none
  type(InputParams), intent(in) :: self

  getNref = self%NREF
end function

!-------------------------------------------------------------------------------
!> Returns IELAST.
!> @param[in] self InputParams object
integer  function getIELAST(self)
  implicit none
  type(InputParams), intent(in) :: self

  getIelast = self%IELAST
end function

!-------------------------------------------------------------------------------
!> Returns NSYMAT.
!> @param[in] self InputParams object
integer  function getNSYMAT(self)
  implicit none
  type(InputParams), intent(in) :: self

  getNsymat = self%NSYMAT
end function

!-------------------------------------------------------------------------------
!> Returns MAXMESH.
!> @param[in] self InputParams object
integer  function getMAXMESH(self)
  implicit none
  type(InputParams), intent(in) :: self

  getMaxmesh = self%MAXMESH
end function

!-------------------------------------------------------------------------------
!> Returns NR.
!> @param[in] self InputParams object
integer  function getNR(self)
  implicit none
  type(InputParams), intent(in) :: self

  getNr = self%NR
end function

end module
