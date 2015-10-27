program whittaker_w_test
! Test run of exponentially decaying negative-energy Coulomb function
! routine whittaker_w
  use precisn, only: wp
  use io_units, only: fi, fo
  use whittaker_w, only: coulomb_whittaker
  implicit none
  real(wp)            :: eta, rho, fx, fxp
  real(wp), pointer   :: eta_vals(:), rho_vals(:)
  integer, pointer    :: l_vals(:)
  integer             :: j, i, k, m, ikm, l, status, sf
  integer             :: num_l, num_eta, num_rho

  write (fo,'(10x,a,/)') 'Whittaker_w test run output.'
  write (fo,'(a,i10)') 'Machine arithmetic: precision = ', precision(fx)
  write (fo,'(19x,a,i10,/)') ' range     = ', range(fx) 
  cases: do j = 1, 1000
     read (fi,'(3i10)') num_l, num_eta, num_rho
     if (num_l <= 0) exit
     write (fo,'(4(a,i6),/)') 'set ', j, ' : num_l = ',               &
          num_l, ' num_eta = ', num_eta, ' num_rho = ', num_rho
     write (fo,'(2(3x,a1),2(7x,a3),17x,a1,16x,a2,4x,a2,/)') 'i', 'l', &
          'eta', 'rho', 'w', 'wp', 'sf'
     allocate (l_vals(num_l), eta_vals(num_eta), rho_vals(num_rho),   &
          stat=status)
     if (status /= 0) then
        write (fo,'(a,i3)') 'allocation error, status = ', status
        stop
     end if
     read (fi,'(5i10)') l_vals
     read (fi,'(5e14.6)') eta_vals
     read (fi,'(5e14.6)') rho_vals
     ikm = 0
     l_values: do i = 1, num_l
        l = l_vals(i)
        eta_values: do k = 1, num_eta
           eta = eta_vals(k)
           rho_values: do m = 1, num_rho
              rho = rho_vals(m)
              ikm = ikm + 1
              call coulomb_whittaker (eta, l, rho, fx, fxp, sf)
              write (fo,'(2i4,2f10.4,2e18.10,i6)') ikm, l, eta,      &
                   rho, fx, fxp, sf
           end do rho_values
           write (fo,'(1x)')
        end do eta_values
     end do l_values
     deallocate (l_vals, eta_vals, rho_vals, stat=status)
     if (status /= 0) then
        write (fo,'(a,i3)') 'deallocation error, status = ', status
        stop
     end if
  end do cases
  stop
end program whittaker_w_test
