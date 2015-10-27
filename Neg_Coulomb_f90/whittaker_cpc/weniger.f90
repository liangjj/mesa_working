module weniger
! Time-stamp: "2001-03-27 14:58:36 cjn"
  use precisn, only: wp
  implicit none
  private
  public levin_weniger

  integer, save               :: minn = 30 ! 60 30
  integer, save               :: maxn = 50 ! 80 50

contains
  subroutine levin_weniger (a, f, n, err)
! wrapper for levin_wen
    real(wp), intent(in)     :: a(0:maxn+1) ! elements of series
    real(wp), intent(out)    :: f        ! final value returned
    integer, intent(out)     :: n        ! # of terms giving 'best' result
    real(wp), intent(out)    :: err      ! estimated uncertainty of result
    integer                  :: l, m
    l = minn
    m = maxn
    call levin_wen (a, l, m, f, n, err)
  end subroutine levin_weniger

  subroutine levin_wen (a, l, m, f, n, err)
! computes the accelerated sum of a series or limit of a sequence.
    integer, intent(in)      :: m        ! most # of terms to be used
    real(wp), intent(in)     :: a(0:m+1) ! elements of series
    integer, intent(in)      :: L        ! least # of terms to be used
    real(wp), intent(out)    :: f        ! final value returned
    integer, intent(out)     :: n        ! # of terms giving 'best' result
    real(wp), intent(out)    :: err      ! estimated uncertainty of result
    logical                  :: improved, start_testing
    logical                  :: previous_iteration
    integer                  :: i
    real(wp)                 :: result(0:m+1), trunc(m+1)
    real(wp)                 :: w1(0:m+1), w2(0:m+1)
    real(wp)                 :: s, test
    real(wp), parameter      :: small = 0.01_wp
    real(wp), parameter      :: huge = 1.0e75_wp

    call wen0 (a, L, w1, w2, result(L), s) ! no test of initial terms
    trunc(L) = 0.0_wp
    improved = .false.
    start_testing = .false.
    previous_iteration = .false.
    test = huge
    n = L
    i_loop: do i = L+1, M
       call wenn (a, i, w1, w2, result(i), s) ! accelerate i-th term
       trunc(i) = ABS(result(i)-result(i-1))  ! truncation error
       improved = (trunc(i) < trunc(i-1)) .or. trunc(i) < small * &
            ABS(result(i))
       start_testing = start_testing .or. (improved .and. previous_iteration)
       previous_iteration = improved
       if (.NOT.start_testing) cycle
       if (trunc(i) < test) then ! find term giving best result
          n = i
          test = trunc(i)
       else
          if (n /= i-1) test = (test + trunc(i)) / 2.0_wp
       end if
    end do i_loop
    f = result(n)
    err = trunc(n)
  end subroutine levin_wen

  subroutine wen0 (a, n, qnum, qden, value, s)
! Weniger algorithm for accelerating a series
    integer, intent(in)      :: n         ! # of elements in a
    real(wp), intent(in)     :: a(0:n+1)  ! terms of series
    real(wp), intent(out)    :: qnum(0:n) ! backward diag of num array
    real(wp), intent(out)    :: qden(0:n) ! backward diag of den array
    real(wp), intent(out)    :: value     ! accelerated value of sum
    real(wp), intent(out)    :: s         ! partial sum of series
    integer                  :: next, j, nmj
    real(wp)                 :: term, termp, c, num

    term = a(0)
    termp = 1.0_wp / a(1)
    s = a(0)
    qden(0) = termp
    qnum(0) = s * termp

    if (n > 0.0_wp) then
       term = a(1)
       termp = 1.0_wp / a(2)
       s = a(1) + s
       qden(1) = termp
       qnum(1) = s * termp
       qden(0) = qden(1) - qden(0)
       qnum(0) = qnum(1) -  qnum(0)
    end if

    do next = 2, n
       term = a(next)
       termp = 1.0_wp / a(next+1)
       s = a(next) + s
       qden(next) = termp
       qnum(next) = s * termp

       num = REAL(next*(next-1),wp)
       j_loop: do j = 1, next
          nmj = next - j
          c = num / REAL((next + j - 1) * (next + j - 2),wp)
          qden(nmj) = qden(nmj+1) - c * qden(nmj)
          qnum(nmj) = qnum(nmj+1) - c *  qnum(nmj)
       end do j_loop
    end do

    value = qnum(0) / qden(0)
  end subroutine wen0

  subroutine wenn (a, n, qnum, qden, value, s)
! Weniger algorithm for accelerating a series
! get n-terms result from result of n-1 terms
    integer, intent(in)      :: n         ! # of elements in a
    real(wp), intent(in)     :: a(0:n+1)  ! elements of series
    real(wp), intent(out)    :: qnum(0:n) ! backward diag of num array
    real(wp), intent(out)    :: qden(0:n) ! backward diag of den array
    real(wp), intent(out)    :: value     ! accelerated value of sum
    real(wp), intent(inout)  :: s         ! partial sum of series
    integer                  :: next, j, nmj
    real(wp)                 :: term, termp, c, num

    next = n
    term = a(next)
    termp = 1.0_wp / a(next+1)
    s = a(next) + s
    qden(next) = termp
    qnum(next) = s * termp

    if (next > 1) then
       num = REAL(next * (next - 1), wp)
       j_loop: do j = 1, next
          nmj = next - j
          c = num / REAL((next + j - 1) * (next + j - 2), wp)
          qden(nmj) = qden(nmj+1) - c * qden(nmj)
          qnum(nmj) = qnum(nmj+1) - c * qnum(nmj)
       end do j_loop
    else if (next == 1) then    ! unlikely case
       qden(0) = qden(1) - qden(0)
       qnum(0) = qnum(1) - qnum(0)
    end if

    value = qnum(0) / qden(0)
  end subroutine wenn
end module weniger

