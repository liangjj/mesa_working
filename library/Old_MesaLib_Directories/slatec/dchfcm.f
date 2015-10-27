*deck dchfcm
      integer function dchfcm (d1, d2, delta)
c***begin prologue  dchfcm
c***subsidiary
c***purpose  check a single cubic for monotonicity.
c***library   slatec (pchip)
c***type      double precision (chfcm-s, dchfcm-d)
c***author  fritsch, f. n., (llnl)
c***description
c
c *usage:
c
c        double precision  d1, d2, delta
c        integer  ismon, dchfcm
c
c        ismon = dchfcm (d1, d2, delta)
c
c *arguments:
c
c     d1,d2:in  are the derivative values at the ends of an interval.
c
c     delta:in  is the data slope over that interval.
c
c *function return values:
c     ismon : indicates the monotonicity of the cubic segment:
c             ismon = -3  if function is probably decreasing;
c             ismon = -1  if function is strictly decreasing;
c             ismon =  0  if function is constant;
c             ismon =  1  if function is strictly increasing;
c             ismon =  2  if function is non-monotonic;
c             ismon =  3  if function is probably increasing.
c           if abs(ismon)=3, the derivative values are too close to the
c           boundary of the monotonicity region to declare monotonicity
c           in the presence of roundoff error.
c
c *description:
c
c          dchfcm:  cubic hermite function -- check monotonicity.
c
c    called by  dpchcm  to determine the monotonicity properties of the
c    cubic with boundary derivative values d1,d2 and chord slope delta.
c
c *cautions:
c     this is essentially the same as old dchfmc, except that a
c     new output value, -3, was added february 1989.  (formerly, -3
c     and +3 were lumped together in the single value 3.)  codes that
c     flag nonmonotonicity by "if (ismon.eq.2)" need not be changed.
c     codes that check via "if (ismon.ge.3)" should change the test to
c     "if (iabs(ismon).ge.3)".  codes that declare monotonicity via
c     "if (ismon.le.1)" should change to "if (iabs(ismon).le.1)".
c
c   refer to  dpchcm
c
c***routines called  d1mach
c***revision history  (yymmdd)
c   820518  date written
c   820805  converted to slatec library version.
c   831201  changed from  isign  to sign  to correct bug that
c           produced wrong sign when -1 .lt. delta .lt. 0 .
c   890206  added save statements.
c   890209  added sign to returned value ismon=3 and corrected
c           argument description accordingly.
c   890306  added caution about changed output.
c   890407  changed name from dchfmc to dchfcm, as requested at the
c           march 1989 slatec cml meeting, and made a few other
c           minor modifications necessitated by this change.
c   890407  converted to new slatec format.
c   890407  modified description to ldoc format.
c   891214  moved save statements.  (wrb)
c***end prologue  dchfcm
c
c  fortran intrinsics used:  dsign.
c  other routines used:  d1mach.
c
c ----------------------------------------------------------------------
c
c  programming notes:
c
c     ten is actually a tuning parameter, which determines the width of
c     the fuzz around the elliptical boundary.
c
c     to produce a single precision version, simply:
c        a. change dchfcm to chfcm wherever it occurs,
c        b. change the double precision declarations to real, and
c        c. change the constants zero, one, ... to single precision.
c
c  declare arguments.
c
      double precision  d1, d2, delta, d1mach
c
c  declare local variables.
c
      integer ismon, itrue
      double precision  a, b, eps, four, one, phi, ten, three, two,
     * zero
      save zero, one, two, three, four
      save ten
c
c  initialize.
c
      data zero /0.d0/, one/1.d0/, two/2.d0/, three/3.d0/, four/4.d0/,
     1      ten /10.d0/
c
c        machine-dependent parameter -- should be about 10*uround.
c***first executable statement  dchfcm
      eps = ten*d1mach(4)
c
c  make the check.
c
      if (delta .eq. zero)  then
c        case of constant data.
         if ((d1.eq.zero) .and. (d2.eq.zero))  then
            ismon = 0
         else
            ismon = 2
         endif
      else
c        data is not constant -- pick up sign.
         itrue = dsign (one, delta)
         a = d1/delta
         b = d2/delta
         if ((a.lt.zero) .or. (b.lt.zero))  then
            ismon = 2
         else if ((a.le.three-eps) .and. (b.le.three-eps))  then
c           inside square (0,3)x(0,3)  implies   ok.
            ismon = itrue
         else if ((a.gt.four+eps) .and. (b.gt.four+eps))  then
c           outside square (0,4)x(0,4)  implies   nonmonotonic.
            ismon = 2
         else
c           must check against boundary of ellipse.
            a = a - two
            b = b - two
            phi = ((a*a + b*b) + a*b) - three
            if (phi .lt. -eps)  then
               ismon = itrue
            else if (phi .gt. eps)  then
               ismon = 2
            else
c              to close to boundary to tell,
c                  in the presence of round-off errors.
               ismon = 3*itrue
            endif
         endif
      endif
c
c  return value.
c
      dchfcm = ismon
      return
c------------- last line of dchfcm follows -----------------------------
      end
