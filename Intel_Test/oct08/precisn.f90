module precisn
! define working precision kind value
! Time-stamp: "2008-10-23 14:40:32 cjn"
  implicit none

! Kind number for at least IEEE 754 double precision
  integer, parameter :: dp = selected_real_kind(15, 307)

! Kind number for quad precision, fall back on double if quad is not available
  integer, parameter, private :: qp_preferred = selected_real_kind(30,1000)
!  integer, parameter          :: qp = (1+SIGN(1,qp_preferred))/2 * &
!       qp_preferred + (1-SIGN(1,qp_preferred))/2 * dp
  integer, parameter :: qp = qp_preferred

! Kind numbers for real type variables and parameters
  integer, parameter :: sp = selected_real_kind(7)
  integer, parameter :: wp = dp
  integer, parameter :: ep = qp

! kind value for 32-bit integer:
  integer, parameter :: long    = selected_int_kind(9)   ! Long  integer
! kind value for 64-bit integer
  integer, parameter :: llong_t = selected_int_kind(18)  ! LLong integer
! kind value for 64-bit integer, if one available:
!  integer, parameter :: llong   = (((1+SIGN(1,LLong_t))/2)*&
!       LLong_t) + (((1-SIGN(1,LLong_t))/2)*Long)
  integer, parameter :: llong = llong_t

  private
  public long, llong
  public sp, wp, ep
  public wrt_precisn

contains

  subroutine wrt_precisn
    use io_units, only: fo
    integer    :: i

    write (fo,'(/,a)') 'Precision definitions:'
    write (fo,'(a)') 'IEEE 754 double precision -> working &
         &precision (wp).'
    write (fo,'(a, i6)') 'wp = selected_real_kind(15, 307) = ', wp
    if (wp == 0) then
       write (fo,'(a)') 'Required precision not available'
       stop
    end if
    write (fo,'(a,i6)') 'Extended precision, ep = ', ep
    write (fo,'(a,i6)') 'Default integer, long = ', long
    write (fo,'(a,i6)') 'Long integer, llong = ', llong
    if (kind(i) /= long) write (fo,'(a,i6)') 'integer kind = ', kind(i)
    write (fo,'(/)')
  end subroutine wrt_precisn

end module precisn
