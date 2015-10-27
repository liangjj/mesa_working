program pgd
! Solution of Generalized Eigenproblem using BSR data provided by
! Oleg Zatsarinny
! Time-stamp: "2008-10-29 06:37:36 cjn"
  use precisn, only: wp
  use io_units, only: fo
  use blacs, only: nprocs, iam, ictxt, io_processor, p, q, ctxt, myrow,&
       mycol, p0, q0, set_cur_ctxt, p_error, nblock
  use pg_diag, only: pgdiag
  use nmlst, only: rd_nmlst
  implicit none
  integer                     :: i, j, ij, status
  real(wp)                    :: t0, t1
  integer                     :: imycol, imyrow, ip, iq
  integer                     :: c0, c1, cr
  integer                     :: parms(3)

  call cpu_time (t0)
  call system_clock (count=c0)

  call BLACS_PINFO (iam, nprocs) ! find process #, total # processors
  call BLACS_GET (-1, 0, ictxt)  ! find default context, ictxt

! define blacs grid to include all processes to broadcast information
  call BLACS_GRIDINIT (ictxt, 'Row-major', 1, nprocs)
  call BLACS_GRIDINFO (ictxt, ip, iq, imyrow, imycol)
  io_processor = (imycol == 0)
  call set_cur_ctxt (1)

  if (io_processor) then ! This is the i/o processor
     write (fo,'(//,15x,a)') '============'
     write (fo,'(15x,a)') 'Program PGD'
     write (fo,'(15x,a,//)') '============'
     write (fo,'(a)') 'BSR Two-Electron Matrix: Parallel General Eigen&
          &solution'
     write (fo,'(a,i6)') 'Number of processors = ', nprocs

     call rd_nmlst    ! read namelist input
     write (fo,'(a,2i6)') 'Blacs process dimensions, p, q = ', p, q
     write (fo,'(a,2i6)') 'Blacs grid size, nblock        = ', nblock

! send grid parameters to processors:
     parms = (/nblock, p, q/)
     call igebs2d (ictxt, 'all', ' ', 3, 1, parms,3)
  else   ! receiving processors
     call igebr2d (ictxt, 'all', ' ', 3, 1, parms, 3, 0, 0)
     nblock = parms(1)
     p = parms(2);     q = parms(3)
  end if

  if (io_processor) io_processor = .false.

  call BLACS_GRIDEXIT (ictxt)   ! kill current mesh
  call BLACS_GET (-1, 0, ctxt)  ! find default context, ctxt

  call BLACS_GRIDINIT (ctxt, 'Row-major', p, q)
  call BLACS_GRIDINFO (ctxt, ip, iq, myrow, mycol)
! now ip,ip are context ctxt dimensions
!     myrow, mycol are processor coords in context

! processors not included in the grid skip to the end
  if (myrow >= 0 .and. myrow < p .and. mycol >= 0 .and. mycol < q) then
     p0 = 0;     q0 = 0
     io_processor = (myrow == p0 .and. mycol == q0)
     call set_cur_ctxt (2)

     call BLACS_BARRIER (ctxt, 'all')   ! temporary precaution
     call pgdiag (p0, q0)   ! perform eigensolution

     if (io_processor) then
        write (fo,'(/,a,/)') 'end of PGD'
        call cpu_time (t1)
        write (fo,'(a,f16.4,a)') 'CPU time     = ', t1 - t0, ' secs'
        call system_clock (count=c1, count_rate=cr)
        write (fo,'(a,f16.4,a)') 'Elapsed time = ', REAL(c1-c0,wp) / &
             REAL(cr,wp), ' secs'
     end if
     call BLACS_GRIDEXIT (ctxt)
  end if
  call BLACS_BARRIER (ictxt, 'All')
  call BLACS_EXIT()
  stop
end program pgd
