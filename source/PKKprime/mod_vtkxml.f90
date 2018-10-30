!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module mod_vtkxml

  implicit none

  character(len=*), parameter :: vtkfmt_ivtkfile = '("<VTKFile type=""PolyData"" version=""0.1"" byte_order=""LittleEndian"">")',&
                               & vtkfmt_fvtkfile = '("</VTKFile>")',&
                               & vtkfmt_ipolydata = '(2X,"<PolyData>")',&
                               & vtkfmt_fpolydata = '(2X,"</PolyData>")',&
                               & vtkfmt_ipiece = '(4X,"<Piece NumberOfPoints=""",I0,""" NumberOfVerts=""0""&
                                                 & NumberOfLines=""0"" NumberOfStrips=""0"" NumberOfPolys=""",I0,""">")',&
                               & vtkfmt_fpiece = '(4X,"</Piece>")',&
                               & vtkfmt_ipoints    = '(6X,"<Points>")',&
                               & vtkfmt_fpoints    = '(6X,"</Points>")',&
                               & vtkfmt_ipolys     = '(6X,"<Polys>")',&
                               & vtkfmt_fpolys     = '(6X,"</Polys>")',&
                               & vtkfmt_ipointdata = '(6X,"<PointData ",(A)," ",(A),">")',&
                               & vtkfmt_fpointdata = '(6X,"</PointData>")',&
                               & vtkfmt_icelldata = '(6X,"<CellData ",(A)," ",(A),">")',&
                               & vtkfmt_fcelldata = '(6X,"</CellData>")',&
                               & vtkfmt_idata_points       = '(8X,"<DataArray NumberOfComponents=""3"" type=""Float32"" format=""ascii"">")',&
                               & vtkfmt_idata_connectivity = '(8X,"<DataArray type=""Int32"" Name=""connectivity"" format=""ascii"">")',&
                               & vtkfmt_idata_offsets      = '(8X,"<DataArray type=""Int32"" Name=""offsets"" format=""ascii"">")',&
                               & vtkfmt_idata_general      = '(8X,"<DataArray type=""",(A),""" Name=""",(A),"""&
                                                             & NumberOfComponents=""",I0,""" format=""",(A),""">")',&
                               & vtkfmt_fdata              = '(8X,"</DataArray>")'

  integer, parameter :: vtkfmxXdata=10

contains

  subroutine write_IBZ_rot(filename,npoints,points,nfaces,nfaceverts,ifaceverts,nsym,rotmat,isym)
    ! This subroutine writes a vtk-file, containing the wedge of the IBZ and
    !  its equvalents by symmetry. As color-code, the number of the ith wedge
    !  (related to the original wedge by symmetry operation i) is given.

    implicit none

    character(len=*), intent(in) :: filename
    integer, intent(in) :: npoints, nfaces, nsym
    integer, intent(in) :: nfaceverts(nfaces), ifaceverts(:,:), isym(nsym)
    double precision, intent(in) :: points(3,npoints), rotmat(64,3,3)

    integer :: ipoint, iface, isy, partsum(nfaces)
    character(len=50)  :: fmtstr
!   character(len=80)  :: filename

    integer, parameter :: ifile= 36

    do iface=1,nfaces
      partsum(iface) = sum(nfaceverts(1:iface))
    end do

    open(unit=ifile, file=trim(filename), form='formatted', action='write')

    !*** START: write start of header
    write(ifile,vtkfmt_ivtkfile)
      write(ifile,vtkfmt_ipolydata)
        write(ifile,vtkfmt_ipiece) nsym*npoints, nsym*nfaces

          !write the points
          write(ifile,vtkfmt_ipoints)
            write(ifile,vtkfmt_idata_points)
              write(fmtstr,'("(",I0,"X,3ES14.5)")') vtkfmxXdata
              do isy=1,nsym
                do ipoint=1,npoints
                  write(ifile,fmtstr)  rotmat(isym(isy),1,:)*points(1,ipoint) &
                                   & + rotmat(isym(isy),2,:)*points(2,ipoint) &
                                   & + rotmat(isym(isy),3,:)*points(3,ipoint)
                end do!ipoint
              end do!isy
            write(ifile,vtkfmt_fdata)
          write(ifile,vtkfmt_fpoints)

          !*** START: write the geometry
          write(ifile,vtkfmt_ipolys)

            !write the connectivity
            write(ifile,vtkfmt_idata_connectivity)
              do isy=1,nsym
                do iface=1,nfaces
                  write(fmtstr,'("(",I0,"X,",I0,"(I0,X))")') vtkfmxXdata, nfaceverts(iface)
                  write(ifile,fmtstr) ifaceverts(1:nfaceverts(iface),iface)-1 + (isy-1)*npoints
                end do!iface
              end do!isy
            write(ifile,vtkfmt_fdata)

            !write the offsets
            write(ifile,vtkfmt_idata_offsets)
              write(fmtstr,'("(",I0,"X,",I0,"(I0,X))")') vtkfmxXdata, 30
              do isy=1,nsym
                write(ifile,fmtstr) partsum + (isy-1)*sum(nfaceverts)
              end do!isy
            write(ifile,vtkfmt_fdata)

          write(ifile,vtkfmt_fpolys)
          !*** END: write the geometry

          !*** START: write the color code
          write(ifile,vtkfmt_icelldata) 'Scalars="isym"', ''
            write(ifile,vtkfmt_idata_general) 'Int32', 'isym', 1, 'ascii'
              write(fmtstr,'("(",I0,"X,",I0,"I)")') vtkfmxXdata, nfaces
              do isy=1,nsym
                write(ifile,fmtstr) ( isy ,iface=1,nfaces )
              end do!isy
            write(ifile,vtkfmt_fdata)
          write(ifile,vtkfmt_fcelldata)
          !*** END: write the color code


        write(ifile,vtkfmt_fpiece)
      write(ifile,vtkfmt_fpolydata)
    write(ifile,vtkfmt_fvtkfile)
    !*** END: write end of header

    close(ifile)

  end subroutine write_IBZ_rot


  subroutine write_pointdata_rot( filename,npoints,points,         &
                                & nscal,scalardata,scalarstring,   &
                                & nvect,vectordata,vectorstring,   &
                                & nsym,rotmat,isym,nall_in,kpt2irr )

    implicit none

    character(len=*), intent(in) :: filename
    integer, intent(in) :: npoints, nscal, nvect, nsym
    integer, intent(in) :: isym(nsym)
    double precision, intent(in) :: points(3,npoints), rotmat(64,3,3)
    double precision, allocatable, intent(in) :: scalardata(:,:), vectordata(:,:,:)
    character(len=*), intent(in) :: scalarstring(:), vectorstring(:)
    integer, intent(in), optional :: nall_in, kpt2irr(:)

    character(len=256) :: str1, str2, fmtstr
    integer :: ii, isy, ipoint, nall

    open(unit=123894, file=trim(filename), action='write', form='formatted')

    nall = npoints
    if(present(nall_in)) nall = nall_in

    !write header
    write(123894, FMT=vtkfmt_ivtkfile)
      write(123894, FMT=vtkfmt_ipolydata)
        write(123894, FMT=vtkfmt_ipiece) nsym*npoints, nsym*nall/3

          !write points
          write(123894, FMT=vtkfmt_ipoints)
            write(123894, FMT=vtkfmt_idata_points)

              write(fmtstr,'("(",I0,"X,9ES14.6)")') vtkfmxXdata
              do isy=1,nsym
                do ipoint=1,npoints
                  write(123894,fmtstr)   rotmat(isym(isy),1,:)*points(1,ipoint) &
                                     & + rotmat(isym(isy),2,:)*points(2,ipoint) &
                                     & + rotmat(isym(isy),3,:)*points(3,ipoint)
                end do!ipoint
              end do!isy

            write(123894, FMT=vtkfmt_fdata)
          write(123894, FMT=vtkfmt_fpoints)


          !write polys
          write(123894, FMT=vtkfmt_ipolys)
            write(123894, FMT=vtkfmt_idata_connectivity)
              do isy=1,nsym
                if(present(kpt2irr))then
                    write(123894, FMT='(10X,4I8)' ) kpt2irr-1 + (isy-1)*npoints
                else!present
                  do ii=1,npoints
                    write(123894, FMT='(10X,4I8)' ) ii-1 + (isy-1)*npoints
                  end do!icub
                end if!present
              end do!isy
            write(123894, FMT=vtkfmt_fdata)
            write(123894, FMT=vtkfmt_idata_offsets)
              write(fmtstr, '("(",I0,"X,6(I0,X))")') vtkfmxXdata
              do isy=1,nsym
                do ii=1,nall/3
                  write(123894, FMT='(10X,4I8)' ) ii*3 + (isy-1)*nall
                end do!icub
              end do!isy
            write(123894, FMT=vtkfmt_fdata)
          write(123894, FMT=vtkfmt_fpolys)


          !=======================!
          !=== write pointdata ===!
          !=======================!
          !gather header information
          if(nscal>0) then
            write(str1,'("Scalars=""",(A),"""")') trim(scalarstring(1))
          else
            str1=''
          end if
          if(nvect>0) then
            write(str2,'("Vectors=""",(A),"""")') trim(vectorstring(1))
          else
            str2=''
          end if
          write(123894, FMT=vtkfmt_ipointdata) trim(str1), trim(str2)

          !write scalar data
          do ii=1,nscal
            write(123894, FMT=vtkfmt_idata_general) 'Float32', trim(scalarstring(ii)), 1, 'ascii'
              write(fmtstr,'("(",I0,"X,10ES14.6)")') vtkfmxXdata
              do isy=1,nsym
                write(123894,fmtstr) scalardata(:,ii)
              end do!isy
            write(123894, FMT=vtkfmt_fdata)
          end do!ii

          !write vector data
          do ii=1,nvect
            write(123894, FMT=vtkfmt_idata_general) 'Float32', trim(vectorstring(ii)), 3, 'ascii'
              write(fmtstr,'("(",I0,"X,3ES14.6)")') vtkfmxXdata
              do isy=1,nsym
                do ipoint=1,npoints
                  write(123894,fmtstr) rotmat(isym(isy),1,:)*vectordata(1,ipoint,ii) &
                                   & + rotmat(isym(isy),2,:)*vectordata(2,ipoint,ii) &
                                   & + rotmat(isym(isy),3,:)*vectordata(3,ipoint,ii)
                end do!ipoint
              end do!isy
            write(123894, FMT=vtkfmt_fdata)
          end do!ii

          write(123894, FMT=vtkfmt_fpointdata)


        !write tail
        write(123894, FMT=vtkfmt_fpiece)
      write(123894, FMT=vtkfmt_fpolydata)
    write(123894, FMT=vtkfmt_fvtkfile)

    close(123894)



  end subroutine write_pointdata_rot

end module mod_vtkxml
