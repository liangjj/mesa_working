program rdf
  implicit none
  integer, parameter :: wp = selected_real_kind(12)
  integer, save :: nuh=7
  integer       :: n, i, jstart, j
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
     read (nuh) jstart
     write (*,*) 'jstart = ', jstart
     read (nuh) x
     if (j == n) write (*,*) 'x = ', x
     if (j > 10) cycle
     write (*,*) 'x = ', x
  end do jst
  write (*,*) 'eigenvectors done'
end program rdf
