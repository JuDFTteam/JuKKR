  write(*,'(a)') "0 -9", "0 9" !! y-axis
  write(*,*)
  write(*,'(a)') "-9 0", "9 0" !! x-axis
  write(*,*)
  do i = -600, 600
    x = 0.01d0*i
    write(*,*) x, (.5d0*x - 1)**2 * (.25d0*x + 1)
    if (i == 0 .or. i == 200) then
      write(*,*)
      write(*,*) x, (.5d0*x - 1)**2 * (.25d0*x + 1)
    endif
  enddo
  write(*,*)
  write(*,'(a)') "2 0", "9 0" !! zero-continuation
  write(*,*)
  write(*,'(a)') "0 1", "2 0" !! critical dots
  
  end
