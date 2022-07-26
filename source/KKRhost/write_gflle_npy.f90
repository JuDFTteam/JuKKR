module mod_write_gflle 

  implicit none

contains

  !-------------------------------------------------------------------------------
  !> Summary: Write gflle file out in npy format
  !> Author: Philipp Rüßmann
  !> Category: writeout
  !> Deprecated: False 
  !> Creates one file per atom and energy, otherwise file can be very large which
  !> might be problematic in post-processing
  !-------------------------------------------------------------------------------
  subroutine write_gflle_to_npy(lmmaxd, ielast, nqdos, i1, gflle)

    use mod_datatypes, only: dp
    use m_npy, only: save_npy
    implicit none
    integer, intent(in) :: lmmaxd, ielast, nqdos, i1
    complex (kind=dp) :: gflle(lmmaxd,lmmaxd,ielast,nqdos)
    character (len=100) :: filename
    integer :: ie
    do ie = 1, ielast
      write(filename, "(A,1I0.3,A,1I0.3,A)") "gllke.", I1, ".", IE, ".npy"
      call save_npy(trim(filename), gflle(:,:,ie, :))
    end do

  end subroutine write_gflle_to_npy

end module mod_write_gflle 
