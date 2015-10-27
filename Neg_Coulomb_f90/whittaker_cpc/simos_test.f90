program simos_test
! Test run of negative-energy Coulomb function integration
! routine simos 
  use precisn, only: wp
  use io_units, only: fi, fo
  use whittaker_w, only: coulomb_whittaker, scale_real
  use simos, only: simos_cint, simos_cder2
  implicit none

  real(wp)               :: ws1, ws2, w, wd
  real(wp)               :: ws(2)
  real(wp)               :: rho, rhop, rhoi, eta, hi
  real(wp), pointer      :: rho_vals(:)
  integer                :: ns, l, sf, sf1, sf2, sft, num_rho
  integer                :: status, ikm, m, j

  write (fo,'(10x,a,/)') 'simos test run output'
  cases: do j = 1, 1000
     read (fi,'(i10,e14.6,i10)') l, eta, num_rho
     if (l < 0) exit
     write (fo,'(2(a,i6),a,e14.8,a,i6,/)') 'set ', j, ' : l = ', l,   &
          ' eta = ', eta, ' num_rho = ', num_rho
     write (fo,'(2(3x,a1),2(7x,a3),17x,a1,16x,a2,4x,a2,/)') 'i', 'l', &
          'eta', 'rho', 'w', 'wp', 'sf'
     allocate (rho_vals(num_rho), stat=status)
     if (status /= 0) then
        write (fo,'(a,i3)') 'allocation error, status = ', status
        stop
     end if
     read (fi,'(5e14.6)') rho_vals
     ikm = 0
     rho_values: do m = 2, num_rho
        rhoi = rho_vals(m-1)
        rho = rho_vals(m)
        ikm = ikm + 1
        ns = 512 * MAX(INT(ABS(rho - rhoi)), 1)
        hi = (rho - rhoi) / REAL(ns,wp)
        call coulomb_whittaker (eta, l, rhoi, ws1, wd, sf1)
        rhop = rhoi - hi
        call coulomb_whittaker (eta, l, rhop, ws2, wd, sf2)
        if (sf1 /= sf2) ws2 = ws2 * 10.0_wp**(sf2-sf1)
        ws = (/ws1, ws2/)
        call simos_cint (rhoi, rho, eta, l, hi, ws)
        ws1 = ws(1)
        call simos_cder2 (rho, hi, ws, wd)
        call scale_real (ws1, w, sft)
        wd = wd * 10.0_wp**(-sft)
        sf = sft + sf1
        write (fo,'(2i4,2f10.4,2e18.10,i6)') ikm, l, eta,      &
             rho, w, wd, sf
     end do rho_values
     write (fo,'(1x)')
  end do cases
  stop
end program simos_test
