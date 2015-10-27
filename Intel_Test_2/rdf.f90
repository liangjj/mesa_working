program rdf
  implicit none
  integer, parameter :: wp = selected_real_kind(12)
  integer, save :: nuh=7
  integer       :: n, i, jstart, j, ios
  real(wp), allocatable :: eval(:), x(:)

  open (unit=nuh, file='dc_vecs.001', form='unformatted',         &
       status='old', access='sequential', action='read')
  read (nuh) n
  write (*,*) ' n = ', n
  allocate (eval(n), x(n))
  read (nuh) eval
  do i = 1, n
     write (*,*) i, ' eval = ', eval(i)
  end do
  jst: do j = 1, n
     read (nuh,*, iostat=ios) jstart
     if (ios /= 0) then
        write (*,*) '***, j, ios = ', j, ios
        write (*,*) 'jstart = ', jstart
        stop
     end if
     write (*,*) 'jstart = ', jstart
     read (nuh, *, iostat=ios) x
     if (ios /= 0) then
        write (*,*) '***** j, jstart = ', j, jstart, ' ios = ', ios
        write (*,*) x
        stop
     end if
     if (j == n) write (*,*) 'x = ', x
     if (j > 10) cycle
     write (*,*) 'x = ', x
  end do jst
  write (*,*) 'eigenvectors done'
end program rdf
