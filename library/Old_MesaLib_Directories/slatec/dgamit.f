*deck dgamit
      double precision function dgamit (a, x)
c***begin prologue  dgamit
c***purpose  calculate tricomi's form of the incomplete gamma function.
c***library   slatec (fnlib)
c***category  c7e
c***type      double precision (gamit-s, dgamit-d)
c***keywords  complementary incomplete gamma function, fnlib,
c             special functions, tricomi
c***author  fullerton, w., (lanl)
c***description
c
c   evaluate tricomi's incomplete gamma function defined by
c
c   dgamit = x**(-a)/gamma(a) * integral from 0 to x of exp(-t) *
c              t**(a-1.)
c
c   for a .gt. 0.0 and by analytic continuation for a .le. 0.0.
c   gamma(x) is the complete gamma function of x.
c
c   dgamit is evaluated for arbitrary real values of a and for non-
c   negative values of x (even though dgamit is defined for x .lt.
c   0.0), except that for x = 0 and a .le. 0.0, dgamit is infinite,
c   which is a fatal error.
c
c   the function and both arguments are double precision.
c
c   a slight deterioration of 2 or 3 digits accuracy will occur when
c   dgamit is very large or very small in absolute value, because log-
c   arithmic variables are used.  also, if the parameter  a  is very
c   close to a negative integer (but not a negative integer), there is
c   a loss of accuracy, which is reported if the result is less than
c   half machine precision.
c
c***references  w. gautschi, a computational procedure for incomplete
c                 gamma functions, acm transactions on mathematical
c                 software 5, 4 (december 1979), pp. 466-481.
c               w. gautschi, incomplete gamma functions, algorithm 542,
c                 acm transactions on mathematical software 5, 4
c                 (december 1979), pp. 482-489.
c***routines called  d1mach, d9gmit, d9lgic, d9lgit, dgamr, dlgams,
c                    dlngam, xerclr, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920528  description and references sections revised.  (wrb)
c***end prologue  dgamit
      double precision a, x, aeps, ainta, algap1, alneps, alng, alx,
     1  bot, h, sga, sgngam, sqeps, t, d1mach, dgamr, d9gmit, d9lgit,
     2  dlngam, d9lgic
      logical first
      save alneps, sqeps, bot, first
      data first /.true./
c***first executable statement  dgamit
      if (first) then
         alneps = -log (d1mach(3))
         sqeps = sqrt(d1mach(4))
         bot = log (d1mach(1))
      endif
      first = .false.
c
      if (x .lt. 0.d0) call xermsg ('slatec', 'dgamit', 'x is negative'
     +   , 2, 2)
c
      if (x.ne.0.d0) alx = log (x)
      sga = 1.0d0
      if (a.ne.0.d0) sga = sign (1.0d0, a)
      ainta = aint (a + 0.5d0*sga)
      aeps = a - ainta
c
      if (x.gt.0.d0) go to 20
      dgamit = 0.0d0
      if (ainta.gt.0.d0 .or. aeps.ne.0.d0) dgamit = dgamr(a+1.0d0)
      return
c
 20   if (x.gt.1.d0) go to 30
      if (a.ge.(-0.5d0) .or. aeps.ne.0.d0) call dlgams (a+1.0d0, algap1,
     1  sgngam)
      dgamit = d9gmit (a, x, algap1, sgngam, alx)
      return
c
 30   if (a.lt.x) go to 40
      t = d9lgit (a, x, dlngam(a+1.0d0))
      if (t.lt.bot) call xerclr
      dgamit = exp (t)
      return
c
 40   alng = d9lgic (a, x, alx)
c
c evaluate dgamit in terms of log (dgamic (a, x))
c
      h = 1.0d0
      if (aeps.eq.0.d0 .and. ainta.le.0.d0) go to 50
c
      call dlgams (a+1.0d0, algap1, sgngam)
      t = log (abs(a)) + alng - algap1
      if (t.gt.alneps) go to 60
c
      if (t.gt.(-alneps)) h = 1.0d0 - sga * sgngam * exp(t)
      if (abs(h).gt.sqeps) go to 50
c
      call xerclr
      call xermsg ('slatec', 'dgamit', 'result lt half precision', 1,
     +   1)
c
 50   t = -a*alx + log(abs(h))
      if (t.lt.bot) call xerclr
      dgamit = sign (exp(t), h)
      return
c
 60   t = t - a*alx
      if (t.lt.bot) call xerclr
      dgamit = -sga * sgngam * exp(t)
      return
c
      end
