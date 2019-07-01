module mod_read_angle

contains

subroutine read_angle(natom,my_rank,density)
use type_density
use mod_config, only: config_runflag
use mod_version_info
use mod_types, only: t_inc
implicit none
!interface
integer  :: natom
integer  :: my_rank
type(density_type)  :: density(natom)
!local
integer :: iatom,ios,ierror
character(len=200)  ::string1
integer,save :: first=1

if (first==1) then
  density(1)%nangleconfigur=0
  if ( config_runflag('force_angles') ) then
    call read_numbofangles(density(1)%nangleconfigur,natom)
    if (my_rank==0) then 
      write(*,*) '  ###############################################'
      write(*,*) 'Using magnetic configuration mode                 '
      write(*,*) 'different magnetic directions are read in in'
      write(*,*) 'in each iteration'
      write(*,*) '  Number of configutations is',density(1)%nangleconfigur
      write(*,*) '  ###############################################'
    end if
  end if

  open(unit=33952084, file='kkrflex_angle', status='old', iostat=ierror)
  if (ierror/=0) then
    if (my_rank==0) then 
      write(*,*) '[read_potential] angle file does not exist'
      write(*,*) '                 setting all starting angles to zero'
      write(*,*) '  ###############################################'
      write(*,*) 'Initial angles for non-collinear calculation'
      write(*,*) '  ###############################################'
      write(*,*) '        iatom      theta               phi                  moment fixed'
      if (t_inc%i_write>0) write(1337,*) 'Initial angles for non-collinear calculation'
      if (t_inc%i_write>0) write(1337,*) 'iatom  theta                    phi  moment fixed'
    end if
    do iatom=1,natom
      density(iatom)%theta=0.0D0
      density(iatom)%phi  =0.0D0
      density(iatom)%magmomentfixed=1
      if (t_inc%i_write>0) write(1337,*) iatom,density(iatom)%theta,density(iatom)%phi,density(iatom)%magmomentfixed
    end do
    return
  else
    call version_check_header(33952084)
  end if

end if



if (my_rank==0) then 
  write(*,*) '  ###############################################'
  write(*,*) 'Initial angles for non-collinear calculation'
  write(*,*) '  ###############################################'
  write(*,*) '        iatom      theta               phi                  moment fixed'
end if
if (t_inc%i_write>0) write(1337,*) 'Initial angles for non-collinear calculation'
if (t_inc%i_write>0) write(1337,*) 'iatom  theta                    phi  moment fixed'
do iatom=1,natom
   string1=this_readline(33952084,ios)
   if (ios/=0) stop '[read_atominfo] Error reading atom info2'
   if (ios==-1) stop '[read_atominfo] EOF'
   read(string1,*) density(iatom)%theta,density(iatom)%phi,density(iatom)%magmomentfixed
   if (t_inc%i_write>0) write(1337,*) iatom,density(iatom)%theta,density(iatom)%phi,density(iatom)%magmomentfixed
   if (my_rank==0) write(*,*) iatom,density(iatom)%theta,density(iatom)%phi,density(iatom)%magmomentfixed
   density(iatom)%theta = density(iatom)%theta/360.0D0*8.0D0*datan(1.0D0)
   density(iatom)%phi   = density(iatom)%phi  /360.0D0*8.0D0*datan(1.0D0)
end do
! close(33952084)
first=0
end subroutine read_angle

subroutine read_numbofangles(nangleconfigur,natom)
use mod_version_info
implicit none
integer :: nangleconfigur, natom
integer :: ios,linecount, ierror
character(len=200)  ::string1

open(unit=345345318, file='kkrflex_angle', status='old', iostat=ierror)
call version_check_header(345345318)
ios=0
linecount = -1
do while (ios/=-1)
   string1=this_readline(345345318,ios)
   linecount=linecount+1
end do

nangleconfigur=linecount/natom
if (nangleconfigur*natom/=linecount) then
  print *,'[read_angle] number of angles given in file kkrflex_angles'
  print *,'             is not correct'
  print *,'nangleconfigur',nangleconfigur
  print *,'linecount',linecount
  print *,'natom',natom
  stop
end if
close(345345318)
end subroutine read_numbofangles

subroutine check_angle(density)
use type_density
use mod_config, only: config_runflag
implicit none
!interface
type(density_type)  :: density(:)
!local
integer  :: natom
integer :: iatom
! double precision,parameter :: pi = 4.d0*datan(1.d0)
double precision :: theta,phi
natom=ubound(density,1)
do iatom=1,natom
  theta = density(iatom)%theta
  phi   = density(iatom)%phi
  if (theta<0.0D0 .or. theta>pi+1.0D-12) then
    write(*,*) 'angle theta of atom ',iatom, 'has value',theta
    write(*,*) 'but is supposed to be between 0 and 180 deg'
    stop
  end if
end do !iatom

end subroutine check_angle







function this_readline(ifile,ios)
!--------------------------------------------------------
!--  reads the next line in file unit ifile            --
!--------------------------------------------------------
!--  files starting with a dash (#) are treated as     --
!--  a comment !!!                                     --
!--  OUTPUT: next line which is not commented out      --
!--          error variable IOS (should be zero)       --
!--------------------------------------------------------
! input variables
  implicit none
integer,intent(in)               :: ifile
integer,intent(out)              :: ios
! local variables
character(len=200)  ::this_readline
do
  read(unit=ifile,fmt='(A)', iostat=ios) this_readline
  if (ios/=0 .or. this_readline(1:1)/='#') exit
end do
end function this_readline

end module mod_read_angle
