*deck gamit
      real function gamit (a, x)
c***begin prologue  gamit
c***purpose  calculate tricomi's form of the incomplete gamma function.
c***library   slatec (fnlib)
c***category  c7e
c***type      single precision (gamit-s, dgamit-d)
c***keywords  complementary incomplete gamma function, fnlib,
c             special functions, tricomi
c***author  fullerton, w., (lanl)
c***description
c
c   evaluate tricomi's incomplete gamma function defined by
c
c   gamit = x**(-a)/gamma(a) * integral from 0 to x of exp(-t) *
c             t**(a-1.)
c
c   for a .gt. 0.0 and by analytic continuation for a .le. 0.0.
c   gamma(x) is the complete gamma function of x.
c
c   gamit is evaluated for arbitrary real values of a and for non-
c   negative values of x (even though gamit is defined for x .lt.
c   0.0), except that for x = 0 and a .le. 0.0, gamit is infinite,
c   which is a fatal error.
c
c   the function and both arguments are real.
c
c   a slight deterioration of 2 or 3 digits accuracy will occur when
c   gamit is very large or very small in absolute value, because log-
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
c***routines called  algams, alngam, gamr, r1mach, r9gmit, r9lgic,
c                    r9lgit, xerclr, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920528  description and references sections revised.  (wrb)
c***end prologue  gamit
      logical first
      save alneps, sqeps, bot, first
      data first /.true./
c***first executable statement  gamit
      if (first) then
         alneps = -log(r1mach(3))
         sqeps = sqrt(r1mach(4))
         bot = log(r1mach(1))
      endif
      first = .false.
c
      if (x .lt. 0.0) call xermsg ('slatec', 'gamit', 'x is negative',
     +   2, 2)
c
      if (x.ne.0.0) alx = log(x)
      sga = 1.0
      if (a.ne.0.0) sga = sign (1.0, a)
      ainta = aint (a+0.5*sga)
      aeps = a - ainta
c
      if (x.gt.0.0) go to 20
      gamit = 0.0
      if (ainta.gt.0.0 .or. aeps.ne.0.0) gamit = gamr(a+1.0)
      return
c
 20   if (x.gt.1.0) go to 40
      if (a.ge.(-0.5) .or. aeps.ne.0.0) call algams (a+1.0, algap1,
     1  sgngam)
      gamit = r9gmit (a, x, algap1, sgngam, alx)
      return
c
 40   if (a.lt.x) go to 50
      t = r9lgit (a, x, alngam(a+1.0))
      if (t.lt.bot) call xerclr
      gamit = exp(t)
      return
c
 50   alng = r9lgic (a, x, alx)
c
c evaluate gamit in terms of log(gamic(a,x))
c
      h = 1.0
      if (aeps.eq.0.0 .and. ainta.le.0.0) go to 60
      call algams (a+1.0, algap1, sgngam)
      t = log(abs(a)) + alng - algap1
      if (t.gt.alneps) go to 70
      if (t.gt.(-alneps)) h = 1.0 - sga*sgngam*exp(t)
      if (abs(h).gt.sqeps) go to 60
      call xerclr
      call xermsg ('slatec', 'gamit', 'result lt half precision', 1, 1)
c
 60   t = -a*alx + log(abs(h))
      if (t.lt.bot) call xerclr
      gamit = sign (exp(t), h)
      return
c
 70   t = t - a*alx
      if (t.lt.bot) call xerclr
      gamit = -sga*sgngam*exp(t)
      return
c
      end
