module blacs
! store blacs context and processor data
! Time-stamp: "2008-11-02 03:59:26 Clifford Noble"
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
! Blacs shutdown on error condition
    use io_units, only: fo
    integer, intent(in)          :: err   ! flag
    character(len=*), intent(in) :: msg   ! error message

    if (err /= 0) then
       write (fo, fmt='(a)') msg
       call BLACS_GET (-1, 0, cur_ctxt) ! find default context
       write (fo, fmt='(a,i6,a)') 'node ', iam, ': calling BLACS_ABORT'
       call BLACS_ABORT (cur_ctxt, err)
    end if
  end subroutine p_error

end module blacs
