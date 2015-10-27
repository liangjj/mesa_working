module blacs
! store blacs context and processor data
! Time-stamp: "2008-10-29 06:32:14 cjn"
  implicit none

! define parameters used to descibe BLACS distributed arrays:
  integer, parameter :: dlen_ = 9
  integer, parameter :: block_cyclic_2d = 1
  integer, parameter :: ctxt_ = 2
  integer, parameter :: m_ = 3
  integer, parameter :: n_ = 4
  integer, parameter :: mb_ = 5
  integer, parameter :: nb_ = 6
  integer, parameter :: rsrc_ = 7
  integer, parameter :: csrc_ = 8
  integer, parameter :: lld_ = 9

! initial linear context:
  integer, save      :: ictxt  ! initial blacs context handle
  integer, save      :: iam    ! processor id in initial context
  integer, save      :: nprocs ! total # processors

! context for main calculation: p * q <= nprocs
  integer, save      :: p=1, q=1   ! # processors in grid
  integer, save      :: ctxt   ! blacs context handle
  integer, save      :: mycol, myrow ! processor coordinates
  logical, save      :: io_processor   ! i/o processor flag
  integer, save      :: p0, q0 ! i/o processor coordinates
  integer, save      :: nblock = 64 ! blacs blocking factor
  integer, save      :: mynumrows   ! Hamiltonian local # rows
  integer, save      :: mynumcols   ! Hamiltonian local # cols

  integer, save      :: cur_ctxt=-5 ! current contxt for p_error

  private
  public p_error, set_cur_ctxt
  public dlen_, block_cyclic_2d, ctxt_, m_, n_, mb_, nb_
  public rsrc_, csrc_, lld_
  public ictxt, nprocs, iam
  public ctxt, io_processor, p, q, mycol, myrow, p0, q0
  public nblock, mynumrows, mynumcols

contains

  subroutine set_cur_ctxt (n)
! Define the Blacs context handle to be used by p_error
    integer, intent(in)    :: n
    select case (n)
    case (1)
       cur_ctxt = ictxt
    case (2)
       cur_ctxt = ctxt
    case default
       cur_ctxt = ictxt
    end select
  end subroutine set_cur_ctxt

  subroutine p_error (err, msg)
!    use io_units, only: fo
    integer, parameter :: fo = 6
    integer, intent(in)          :: err   ! flag
    character(len=*), intent(in) :: msg   ! error message
    integer                      :: no

    if (cur_ctxt /= -5) then
       no = 2
       call igsum2d (cur_ctxt, 'all', ' ', 1, 1, err, 1, -1, -1)
       if (err /= 0) then
          if (iam == 0) write (fo, fmt='(a)') msg
          call BLACS_ABORT (cur_ctxt, no)
       end if
    else   ! no context set
       no = 1
       if (err /= 0) then
          write (fo, fmt='(a)') msg
          call BLACS_GET (-1, 0, cur_ctxt) ! find default context
          call BLACS_ABORT (cur_ctxt, no)
       end if
    end if
  end subroutine p_error

end module blacs
program pgd
! Solution of Generalized Eigenproblem using BSR data provided by
! Oleg Zatsarinny
! Time-stamp: "2008-10-30 13:21:58 cjn"
!!  use precisn, only: wp
!!  use io_units, only: fo
  use blacs, only: nprocs, iam, ictxt, io_processor, p, q, ctxt, myrow,&
       mycol, p0, q0, set_cur_ctxt, p_error, nblock
!!  use pg_diag, only: pgdiag
!!  use nmlst, only: rd_nmlst
  implicit none
  integer, parameter          :: wp = selected_real_kind(12)
  integer, parameter          :: fo = 8
  integer                     :: i, j, ij, status
  real(wp)                    :: t0, t1
  integer                     :: imycol, imyrow, ip, iq
  integer                     :: c0, c1, cr
  integer                     :: parms(3)
  real*8                      :: test, ddot

  call cpu_time (t0)
  call system_clock (count=c0)
  call pdpotrf
  call pdsyngst
  call pdsyevd
  test=ddot
  call BLACS_PINFO (iam, nprocs) ! find process #, total # processors
  call BLACS_GET (-1, 0, ictxt)  ! find default context, ictxt

! define blacs grid to include all processes to broadcast information
  call BLACS_GRIDINIT (ictxt, 'Row-major', 1, nprocs)
  call BLACS_GRIDINFO (ictxt, ip, iq, imyrow, imycol)
  io_processor = (imycol == 0)
!!  call set_cur_ctxt (1)

  if (io_processor) then ! This is the i/o processor
     write (fo,'(//,15x,a)') '============'
     write (fo,'(15x,a)') 'Program PGD'
     write (fo,'(15x,a,//)') '============'
     write (fo,'(a)') 'BSR Two-Electron Matrix: Parallel General Eigen&
          &solution'
     write (fo,'(a,i6)') 'Number of processors = ', nprocs

!!     call rd_nmlst    ! read namelist input
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
!!     call set_cur_ctxt (2)

     call BLACS_BARRIER (ctxt, 'all')   ! temporary precaution
!!     call pgdiag (p0, q0)   ! perform eigensolution

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


