*deck dgamic
      double precision function dgamic (a, x)
c***begin prologue  dgamic
c***purpose  calculate the complementary incomplete gamma function.
c***library   slatec (fnlib)
c***category  c7e
c***type      double precision (gamic-s, dgamic-d)
c***keywords  complementary incomplete gamma function, fnlib,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c   evaluate the complementary incomplete gamma function
c
c   dgamic = integral from x to infinity of exp(-t) * t**(a-1.)  .
c
c   dgamic is evaluated for arbitrary real values of a and for non-
c   negative values of x (even though dgamic is defined for x .lt.
c   0.0), except that for x = 0 and a .le. 0.0, dgamic is undefined.
c
c   dgamic, a, and x are double precision.
c
c   a slight deterioration of 2 or 3 digits accuracy will occur when
c   dgamic is very large or very small in absolute value, because log-
c   arithmic variables are used.  also, if the parameter a is very close
c   to a negative integer (but not a negative integer), there is a loss
c   of accuracy, which is reported if the result is less than half
c   machine precision.
c
c***references  w. gautschi, a computational procedure for incomplete
c                 gamma functions, acm transactions on mathematical
c                 software 5, 4 (december 1979), pp. 466-481.
c               w. gautschi, incomplete gamma functions, algorithm 542,
c                 acm transactions on mathematical software 5, 4
c                 (december 1979), pp. 482-489.
c***routines called  d1mach, d9gmic, d9gmit, d9lgic, d9lgit, dlgams,
c                    dlngam, xerclr, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920528  description and references sections revised.  (wrb)
c***end prologue  dgamic
      double precision a, x, aeps, ainta, algap1, alneps, alngs, alx,
     1  bot, e, eps, gstar, h, sga, sgng, sgngam, sgngs, sqeps, t,
     2  d1mach, dlngam, d9gmic, d9gmit, d9lgic, d9lgit
      logical first
      save eps, sqeps, alneps, bot, first
      data first /.true./
c***first executable statement  dgamic
      if (first) then
         eps = 0.5d0*d1mach(3)
         sqeps = sqrt(d1mach(4))
         alneps = -log (d1mach(3))
         bot = log (d1mach(1))
      endif
      first = .false.
c
      if (x .lt. 0.d0) call xermsg ('slatec', 'dgamic', 'x is negative'
     +   , 2, 2)
c
      if (x.gt.0.d0) go to 20
      if (a .le. 0.d0) call xermsg ('slatec', 'dgamic',
     +   'x = 0 and a le 0 so dgamic is undefined', 3, 2)
c
      dgamic = exp (dlngam(a+1.d0) - log(a))
      return
c
 20   alx = log (x)
      sga = 1.0d0
      if (a.ne.0.d0) sga = sign (1.0d0, a)
      ainta = aint (a + 0.5d0*sga)
      aeps = a - ainta
c
      izero = 0
      if (x.ge.1.0d0) go to 40
c
      if (a.gt.0.5d0 .or. abs(aeps).gt.0.001d0) go to 30
      e = 2.0d0
      if (-ainta.gt.1.d0) e = 2.d0*(-ainta+2.d0)/(ainta*ainta-1.0d0)
      e = e - alx * x**(-0.001d0)
      if (e*abs(aeps).gt.eps) go to 30
c
      dgamic = d9gmic (a, x, alx)
      return
c
 30   call dlgams (a+1.0d0, algap1, sgngam)
      gstar = d9gmit (a, x, algap1, sgngam, alx)
      if (gstar.eq.0.d0) izero = 1
      if (gstar.ne.0.d0) alngs = log (abs(gstar))
      if (gstar.ne.0.d0) sgngs = sign (1.0d0, gstar)
      go to 50
c
 40   if (a.lt.x) dgamic = exp (d9lgic(a, x, alx))
      if (a.lt.x) return
c
      sgngam = 1.0d0
      algap1 = dlngam (a+1.0d0)
      sgngs = 1.0d0
      alngs = d9lgit (a, x, algap1)
c
c evaluation of dgamic(a,x) in terms of tricomi-s incomplete gamma fn.
c
 50   h = 1.d0
      if (izero.eq.1) go to 60
c
      t = a*alx + alngs
      if (t.gt.alneps) go to 70
      if (t.gt.(-alneps)) h = 1.0d0 - sgngs*exp(t)
c
      if (abs(h).lt.sqeps) call xerclr
      if (abs(h) .lt. sqeps) call xermsg ('slatec', 'dgamic',
     +   'result lt half precision', 1, 1)
c
 60   sgng = sign (1.0d0, h) * sga * sgngam
      t = log(abs(h)) + algap1 - log(abs(a))
      if (t.lt.bot) call xerclr
      dgamic = sgng * exp(t)
      return
c
 70   sgng = -sgngs * sga * sgngam
      t = t + algap1 - log(abs(a))
      if (t.lt.bot) call xerclr
      dgamic = sgng * exp(t)
      return
c
      end
