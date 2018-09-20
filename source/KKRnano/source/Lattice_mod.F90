module Lattice_mod
!-------------------------------------------------------------------------------
!> Summary: Computes the inverse of the Bravais matrix
!> Author: Paul F Baumeister, Elias Rabel
!> Category: KKRnano, geometry, initialization
!-------------------------------------------------------------------------------
  implicit none
  private
  
  public :: lattix99
  
  contains

  !*==lattix99.f    processed by spag 6.05rc at 17:56 on 17 may 2004
  subroutine lattix99(alat, bravais, recbv, volume, output)
    use VectorMath_mod, only: ddet33, spatpr, cross
    use Constants_mod, only: pi
    ! **********************************************************************
    ! *                                                                    *
    ! * lattix99 generates the real space and reciprocal lattices.         *
    ! * bravais(i,j) are basis vectors, with i=x,y,z and j=a,b,c           *
    ! *              input - normalised vectors                            *
    ! * reciprocal space vectors are in units of 2*pi/alatc - output       *
    ! * rr are the direct space vectors - output                           *
    ! * nr+1 is the number of direct space vectors created - output        *
    ! * (structure dependent output).                                      *
    ! *                                                                    *
    ! **********************************************************************
    logical, intent(in) :: output
    double precision, intent(in) :: alat
    double precision, intent(out) :: volume
    double precision, intent(in) :: bravais(3,3) !> real space bravais vectors. read in normalised to alat
    double precision, intent(out) :: recbv(3,3) !> reciprocal lattice vectors in 2*pi/alat

    integer :: i!, j
    double precision :: voluc, det, tpia
    character(len=*), parameter :: F99001="(9x,32(1h-),6x,32(1h-))", F99002="(5x,a2,i1,':',3f10.6,8x,3f10.6)"

    tpia = 2.d0*pi/alat

    ! ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo output
    if (output) then
      write(6,'(79(1h=))')
      write(6,'(23x,a)') '  LATTIX99: bulk geometry mode'
      write(6,'(79(1H=))')
      write(6,*)
      write(6,'(5x,a,f12.8,4x,a,f12.8,/)') 'Lattice constants :  ALAT =',alat,' 2*PI/ALAT =',tpia
      
      write(6,'(5x,a,/)') 'Direct lattice cell vectors :'
      write(6,'(9x,a,21x,a)') 'normalised (ALAT)','a.u.'
      write(6, fmt=F99001)
      do i = 1, 3
        write(6, fmt=F99002) 'a_',i, bravais(1:3,i), bravais(1:3,i)*alat
      enddo ! i
      write(6, fmt=F99001)
      write(6,*)
    endif ! output

    ! ----------------------------------------------------------------------
    ! now generate the reciprocal lattice unit-vectors,
    ! and calculate the unit-cell volume in units au**3.
    ! ----------------------------------------------------------------------

    det = ddet33(bravais)
    if ( abs(det) < 1d-8 ) stop ' ERROR: 3d Bravais vectors are linearly dependent'

    recbv(1:3,1) = cross(bravais(1:3,2), bravais(1:3,3))
    recbv(1:3,2) = cross(bravais(1:3,3), bravais(1:3,1))
    recbv(1:3,3) = cross(bravais(1:3,1), bravais(1:3,2))

    voluc = abs(spatpr(bravais(1:3,2), bravais(1:3,3), bravais(1:3,1)))
    
    recbv(:,:) = recbv(:,:)/voluc
    ! ----------------------------------------------------------------------

    ! --> test on volume unit cell:

    if (voluc < 1.0d-5) stop ' ERROR: Unit-cell volume suspiciously small ( < 1D-5)'

    if (output) then
      write(6,fmt='(5x,a,f8.4,a,f14.8,a,/)') 'Unit cell volume :  V =',voluc,' (ALAT**3) = ', voluc*(alat**3),' (a.u.**3)'
    endif ! output

    ! --> check volume of unit cell vs. average WS-radius

    volume = voluc * alat**3

    ! ----------------------------------------------------------------------
    !  Reciprocal lattice unit-vectors and unit-cell volume calculated
    ! ----------------------------------------------------------------------

    if (output) then
      write(6,'(5x,a,/)') 'Reciprocal lattice cell vectors :'
      write(6,'(9x,a,16x,a)') 'normalised (2*PI/ALAT)','1/a.u.'
      write(6, fmt=F99001)
      do i = 1, 3
        write(6, fmt=F99002) 'b_', i, recbv(1:3,i), recbv(1:3,i)*tpia
      enddo ! i
      write(6, fmt=F99001)
      write(6,*)
    endif ! output

  endsubroutine ! lattix99

endmodule Lattice_mod      
