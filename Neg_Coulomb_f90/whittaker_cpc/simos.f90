module simos
! Simos-Numerov integration of Negative Energy Coulomb equation
! Time-stamp: "2001-05-10 15:23:31 cjn"
  use precisn, only: wp
  use io_units, only: fo
  implicit none

  real(wp), save     :: eta    ! Coulomb parameter
  integer, save      :: l      ! orbital a.m.
  integer, save      :: range_type
  logical, save      :: osc    ! region type (oscillatory/exponential)

  private
  public simos_int, simos_cint, simos_cder2

contains
  subroutine simos_cint (ri, rf, eta_t, l_t, hi, y)
! integrate Negative Energy Coulomb equation from rho = ri to rho = rf
    real(wp), intent(in)     :: ri     ! initial rho value
    real(wp), intent(in)     :: rf     ! final rho value
    real(wp), intent(in)     :: eta_t  ! Coulomb parameter (eta < 0)
    integer, intent(in)      :: l_t    ! orbital ang. momentum
    real(wp), intent(in)     :: hi     ! step-size
    real(wp), intent(inout)  :: y(2)   ! Coulomb fn values
    real(wp)                 :: phi, rt, hip, r1, r2, rb, rb1, rb2
    integer                  :: nsp
    real(wp), parameter      :: eps = 1.0e-12_wp
 
    eta = eta_t              ! set global variables ...
    l = l_t

    call ranges (ri, rf, r1, r2)
    if (SIGN(1.0_wp,hi) /= SIGN(1.0_wp,rf-ri)) then
       write (fo,'(a,2e16.8)') 'hi, sgn = ', hi, sign(1.0_wp,hi)
       write (fo,'(a,2e16.8)') 'ri-rf, sgn = ', ri-rf, sign(1.0_wp,ri-rf)
       write (fo,'(a,e16.8)') 'sn_int: hi sign inconsistent with ri, rf'
       stop
    end if

    select case (range_type)
    case (1)            ! one exponential region in range
       osc = .false.
       phi = 0.5_wp * (SQRT(ABS(yd(rf))) + SQRT(ABS(yd(ri))))
       call simos_int (y, hi, ri, rf, phi)
    case (4)          ! one oscillatory region in range
       osc = .true.
       rt = REAL(l*(l+1),wp) / ABS(eta)
       if (ABS(ABS(rt-ri) + ABS(rt-rf) - ABS(rf-ri)) <= eps) then
          phi = (SQRT(ABS(yd(rf))) + SQRT(ABS(yd(ri))) + &
               SQRT(ABS(yd(rt)))) / 3.0_wp
       else
          phi = 0.5_wp * (SQRT(ABS(yd(rf))) + SQRT(ABS(yd(ri))))
       end if
       call simos_int (y, hi, ri, rf, phi)
    case (2,5)          ! two regions bounded by rb in range
       if (range_type == 2) then
          osc = ri > rf
          rb = r1
       else
          osc = rf > ri
          rb = r2
       end if
       rt = REAL(l*(l+1),wp) / ABS(eta)
       if (ABS(ABS(rt-ri) + ABS(rt-rb) - ABS(rb-ri)) <= eps) then
          phi = 0.5_wp * SQRT(ABS(yd(rt)))
       else
          phi = 0.5_wp * SQRT(ABS(yd(ri)))
       end if
       nsp = MAX(NINT((rb - ri) / hi), 1)
       hip = (rb - ri) / REAL(nsp,wp)
       call simos_int (y, hip, ri, rb, phi)
       osc = .NOT.osc
       if (ABS(ABS(rt-rb) + ABS(rt-rf) - ABS(rb-rf)) <= eps) then
          phi = 0.5_wp * SQRT(ABS(yd(rt)))
       else
          phi = 0.5_wp * SQRT(ABS(yd(rf)))
       end if
       nsp = MAX(NINT((rf - rb) / hi), 1)
       hip = (rf - rb) / REAL(nsp,wp)
       call simos_int (y, hip, rb, rf, phi)
    case (3)           ! three regions in range (boundaries rb1, rb2)
       osc = .false.
       if (rf > ri) then
          rb1 = r1
          rb2 = r2
       else
          rb1 = r2
          rb2 = r1
       end if
       phi = 0.5_wp * SQRT(ABS(yd(ri)))
       nsp = MAX(NINT((rb1 - ri) / hi), 1)
       hip = (rb1 - ri) / REAL(nsp,wp)
       call simos_int (y, hip, ri, rb1, phi)
       osc = .true.
       rt = REAL(l*(l+1),wp) / ABS(eta)
       phi = 0.5_wp * SQRT(ABS(yd(rt)))
       nsp = MAX(NINT((rb2 - rb1) / hi), 1)
       hip = (rb2 - rb1) / REAL(nsp,wp)
       call simos_int (y, hip, rb1, rb2, phi)
       osc = .false.
       phi = 0.5_wp * SQRT(ABS(yd(rf)))
       nsp = MAX(NINT((rf - rb2) / hi), 1)
       hip = (rf - rb2) / REAL(nsp,wp)
       call simos_int (y, hip, rb2, rf, phi)
    end select
  end subroutine simos_cint

  subroutine ranges (ri, rf, r1, r2)
! determine the behaviour (oscillatory, exponential) over the
! integration range. Set global variable range_type accordingly
    real(wp), intent(in)  :: ri ! initial integration point
    real(wp), intent(in)  :: rf ! final integration point
    real(wp), intent(out) :: r1 ! boundary oscillatory/exponential
    real(wp), intent(out) :: r2 ! boundary oscillatory/exponential
 
    call coulomb_region (r1, r2)
    if (r2 == 0.0_wp) then
       range_type = 1        ! single exponential region
    else
       if (ri >= r2) then ! outer exponential region
          if (rf >= r2) then
             range_type = 1
          else if (rf < r2 .and. rf >= r1) then
             range_type = 5
          else if (rf < r1) then
             range_type = 3
          end if
       else if (ri >= r1) then ! oscillatory region
          if (rf > r2) then
             range_type = 5
          else if (rf < r1) then
             range_type = 2
          else
             range_type = 4
          end if
       else              ! inner exponential region
          if (rf <= r1) then
             range_type = 1
          else if (rf > r1 .and. rf <= r2) then
             range_type = 2
          else if (rf > r2) then
             range_type = 3
          end if
       end if
    end if
  end subroutine ranges

  subroutine coulomb_region (r1, r2)
! Define character of negative-energy Coulomb solution
    real(wp), intent(out)  :: r1  ! boundary of inner exponential region
    real(wp), intent(out)  :: r2  ! boundary of outer exponential region
    real(wp)               :: sq
    integer                :: l_crit
 
    r1 = 0.0_wp
    r2 = 0.0_wp
    l_crit = CEILING(SQRT(eta * eta + 0.25_wp) - 0.5_wp)
    if (eta >= 0.0_wp .or. l >= l_crit) return
    sq = SQRT(eta * eta - REAL(l*(l+1),wp))
    r1 = ABS(eta) - sq
    r2 = ABS(eta) + sq
  end subroutine coulomb_region

  function yd (r)
    real(wp)                :: yd
    real(wp), intent(in)    :: r

    yd = (1.0_wp + (2.0_wp * eta + REAL(l*(l+1),wp) / r) / r)
  end function yd

  subroutine simos_int (y, hi, ri, rf, phi)
! Simos integration of negative energy Coulomb equation
! Integration range entirely in the oscillatory region
    real(wp), intent(in)    :: ri   ! initial radius
    real(wp), intent(in)    :: rf   ! final radius
    real(wp), intent(in)    :: hi   ! stepsize
    real(wp), intent(inout) :: y(2) ! fn at ri -> fn at rf
    real(wp), intent(in)    :: phi  ! phase
    real(wp)                :: h, h2, yn, ym, yp, r, yb, yb1
    real(wp)                :: ah2, b0h2, b1h2, dn, dm, db1, db
    real(wp)                :: a, a1, b0, b1
    integer                 :: ns, i

    ns = MAX(NINT((rf - ri) / hi), 1)
    if (ABS(ri + ns * hi - rf) >= 1.0e-6_wp * ABS(hi)) then
       write (fo,'(a)') 'simos: inconsistent input data'
       write (fo,'(a,e16.8)') 'rf - ri = ', rf - ri
       write (fo,'(a,i6)') 'nsteps = ', ns
       write (fo,'(a,e16.8)') 'nsteps * hi = ', ns * hi
       stop
    end if
    h = (rf - ri) / REAL(ns,wp)

    if (osc) then   ! fitted sin, cos behaviour
       call trig_consts (a, a1, b0, b1, h, phi)
    else            ! fitted exponential behaviour
       call exp_consts (a, a1, b0, b1, h, phi)
    end if

    r = ri
    h2 = h * h
    ah2 = a * h2
    b0h2 = b0 * h2
    b1h2 = b1 * h2
    yn = y(1)
    ym = y(2)

    dn = yd(ri) * yn
    dm = yd(ri-hi) * ym
    steps: do i = 1, ns
       r = r + h
       yb1 = 2.0_wp * yn - ym + h2 * dn
       db1 = yd(r) * yb1
       yb = yn - ah2 * (db1 - 2.0_wp * dn + dm)
       db = yd(r-h) * yb
       yp = b0h2 * (db1 + dm) + b1h2 * db - a1 * yn - ym
       ym = yn
       dm = dn
       yn = yp
       dn = yd(r) * yp
    end do steps
    y = (/yn, ym/)
  end subroutine simos_int

  subroutine simos_cder2 (rf, h, y, yp)
! two-point formula for first derivative of Coulomb function
    real(wp), intent(in)     :: rf    ! radius
    real(wp), intent(in)     :: h     ! step size
    real(wp), intent(in)     :: y(2)  ! y(rf), y(rf-h)
    real(wp), intent(out)    :: yp    ! y'(rf)
    real(wp)                 :: on, om, a, eps, r2, v

    v = SQRT(ABS(yd(rf - 0.5_wp * h)))
    if (osc) then
       r2 = rf - h
       eps = ATAN2(y(2) * COS(v * rf) - y(1) * COS(v * r2), &
            y(2) * SIN(v * rf) - y(1) * SIN(v * r2))
       a = y(1) / COS(v * rf + eps)
       yp = - a * v * SIN(v * rf + eps)
    else
       on = EXP(v * rf)
       om = EXP(v * (rf - h))
       a = (y(1) * on - y(2) * om) / ((on - om) * (on + om))
       eps = on * (y(1) - a * on)
       yp = v * (a * on - eps / on)
    end if
  end subroutine simos_cder2

  subroutine trig_consts (a, a1, b0, b1, h, phi)
! integration coefficients
    real(wp), intent(in)    :: h
    real(wp), intent(in)    :: phi
    real(wp), intent(out)   :: a
    real(wp), intent(out)   :: a1
    real(wp), intent(out)   :: b0
    real(wp), intent(out)   :: b1
    real(wp)                :: w, c, s, w2, w4, w6, w8, w10, w12
    real(wp), parameter     :: wsmall = 0.2_wp

    w = h * phi
    if (w > wsmall) then
       c = COS(w)
       s = SIN(w)
       w2 = w * w
       w4 = w2 * w2
       a1 = (3.0_wp * (3.0_wp * w2 - 16.0_wp) * C + w *              &
            (w2 - 33.0_wp) * s) / 24.0_wp
       b0 = ((5.0_wp - w2) * S - 5.0_wp * w * C) / (8.0_wp * w * w2)
       b1 = (w * (10.0_wp - 7.0_wp * w2) * C - (w4 - 17.0_wp * w2 +  &
            10.0_wp) * S) / (8.0_wp * w * W2)
       a = (3.0_wp * w * C + (w2 - 3.0_wp) * S) /                    &
            (((w4 - 17.0_wp * w2 + 10.0_wp) * S + w * (7.0_wp * w2 - &
            10.0_wp) * C) * (-3.0_wp * w2))
    else
       w2 = w * w
       w4 = w2 * w2
       w6 = w2 * w4
       w8 = w2 * w6
       w10 = w2 * w8
       w12 = w2 * w10
       a1 = -2.0_wp + w8 / 20160.0_wp - w10 / 453600.0_wp +          &
            w12 / 23950080.0_wp
       b0 = 1.0_wp / 12.0_wp - w4 / 3360.0_wp + w6 / 90720.0_wp -    &
            w8 / 5322240.0_wp + w10 / 518918400.0_wp -               &
            w12 / 74724249600.0_wp
       b1 = 5.0_wp / 6.0_wp + w4 / 1680.0_wp - w6 / 4536.0_wp +      &
            23.0_wp * w8 / 2661120.0_wp - w10 / 6486480.0_wp +       &
            61.0_wp * w12 / 37362124800.0_wp
       a = - 1.0_wp / 300.0_wp + w2 / 4200.0_wp - w4 / 236250.0_wp - &
            277.0_wp * w6 / 291060000.0_wp + 18847.0_wp * w8 /       &
            189189000000.0_wp - 31463.0_wp * w10 /                   &
            8939180250000.0_wp - 13888499.0_wp * w12 /               &
            60786425700000000.0_wp
    end if
  end subroutine trig_consts

  subroutine exp_consts (a, a1, b0, b1, h, phi)
! integration coefficients - exponential phases
    real(wp), intent(in)    :: h
    real(wp), intent(in)    :: phi
    real(wp), intent(out)   :: a
    real(wp), intent(out)   :: a1
    real(wp), intent(out)   :: b0
    real(wp), intent(out)   :: b1
    real(wp)                :: w, c, s, w2, w4, w6, w8, w10, w12
    real(wp), parameter     :: wsmall = 0.2_wp

    w = h * phi
    if (w > wsmall) then
       c = COSH(w)
       s = SINH(w)
       w2 = w * w
       w4 = w2 * w2
       a1 = (-3.0_wp * (3.0_wp * w2 + 16.0_wp) * C + w * &
            (w2 + 33.0_wp) * s) / 24.0_wp
       b0 = (-(5.0_wp + w2) * S + 5.0_wp * w * C) / (8.0_wp * w * w2)
       b1 = (-w * (10.0_wp + 7.0_wp * w2) * C + (w4 + 17.0_wp * w2 + &
            10.0_wp) * S) / (8.0_wp * w * W2)
       a = (3.0_wp * w * C - (w2 + 3.0_wp) * S) / &
            (((w4 + 17.0_wp * w2 + 10.0_wp) * S - w * (7.0_wp * w2 + &
            10.0_wp) * C) * (3.0_wp * w2))
    else
       w2 = w * w
       w4 = w2 * w2
       w6 = w2 * w4
       w8 = w2 * w6
       w10 = w2 * w8
       w12 = w2 * w10
       a1 = -2.0_wp + w8 / 20160.0_wp + w10 / 453600.0_wp + &
            w12 / 23950080.0_wp
       b0 = 1.0_wp / 12.0_wp - w4 / 3360.0_wp - w6 / 90720.0_wp -    &
            w8 / 5322240.0_wp - w10 / 518918400.0_wp - &
            w12 / 74724249600.0_wp
       b1 = 5.0_wp / 6.0_wp + w4 / 1680.0_wp + w6 / 4536.0_wp +      &
            23.0_wp * w8 / 2661120.0_wp + w10 / 6486480.0_wp +       &
            61.0_wp * w12 / 37362124800.0_wp
       a = - 1.0_wp / 300.0_wp - w2 / 4200.0_wp - w4 / 236250.0_wp + &
            277.0_wp * w6 / 291060000.0_wp + 18847.0_wp * w8 /       &
            189189000000.0_wp + 31463.0_wp * w10 /                   &
            8939180250000.0_wp - 13888499.0_wp * w12 /               &
            60786425700000000.0_wp
    end if
  end subroutine exp_consts
end module simos
