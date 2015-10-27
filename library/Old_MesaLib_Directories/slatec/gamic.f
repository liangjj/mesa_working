*deck gamic
      real function gamic (a, x)
c***begin prologue  gamic
c***purpose  calculate the complementary incomplete gamma function.
c***library   slatec (fnlib)
c***category  c7e
c***type      single precision (gamic-s, dgamic-d)
c***keywords  complementary incomplete gamma function, fnlib,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c   evaluate the complementary incomplete gamma function
c
c   gamic = integral from x to infinity of exp(-t) * t**(a-1.)  .
c
c   gamic is evaluated for arbitrary real values of a and for non-
c   negative values of x (even though gamic is defined for x .lt.
c   0.0), except that for x = 0 and a .le. 0.0, gamic is undefined.
c
c   gamic, a, and x are real.
c
c   a slight deterioration of 2 or 3 digits accuracy will occur when
c   gamic is very large or very small in absolute value, because log-
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
c***routines called  algams, alngam, r1mach, r9gmic, r9gmit, r9lgic,
c                    r9lgit, xerclr, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920528  description and references sections revised.  (wrb)
c***end prologue  gamic
      logical first
      save eps, sqeps, alneps, bot, first
      data first /.true./
c***first executable statement  gamic
      if (first) then
         eps = 0.5*r1mach(3)
         sqeps = sqrt(r1mach(4))
         alneps = -log(r1mach(3))
         bot = log(r1mach(1))
      endif
      first = .false.
c
      if (x .lt. 0.0) call xermsg ('slatec', 'gamic', 'x is negative',
     +   2, 2)
c
      if (x.gt.0.0) go to 20
      if (a .le. 0.0) call xermsg ('slatec', 'gamic',
     +   'x = 0 and a le 0 so gamic is undefined', 3, 2)
c
      gamic = exp (alngam(a+1.0) - log(a))
      return
c
 20   alx = log(x)
      sga = 1.0
      if (a.ne.0.0) sga = sign (1.0, a)
      ma = a + 0.5*sga
      aeps = a - ma
c
      izero = 0
      if (x.ge.1.0) go to 60
c
      if (a.gt.0.5 .or. abs(aeps).gt.0.001) go to 50
      fm = -ma
      e = 2.0
      if (fm.gt.1.0) e = 2.0*(fm+2.0)/(fm*fm-1.0)
      e = e - alx*x**(-0.001)
      if (e*abs(aeps).gt.eps) go to 50
c
      gamic = r9gmic (a, x, alx)
      return
c
 50   call algams (a+1.0, algap1, sgngam)
      gstar = r9gmit (a, x, algap1, sgngam, alx)
      if (gstar.eq.0.0) izero = 1
      if (gstar.ne.0.0) alngs = log (abs(gstar))
      if (gstar.ne.0.0) sgngs = sign (1.0, gstar)
      go to 70
c
 60   if (a.lt.x) gamic = exp (r9lgic(a, x, alx))
      if (a.lt.x) return
c
      sgngam = 1.0
      algap1 = alngam (a+1.0)
      sgngs = 1.0
      alngs = r9lgit (a, x, algap1)
c
c evaluation of gamic(a,x) in terms of tricomi-s incomplete gamma fn.
c
 70   h = 1.0
      if (izero.eq.1) go to 80
c
      t = a*alx + alngs
      if (t.gt.alneps) go to 90
      if (t.gt.(-alneps)) h = 1.0 - sgngs*exp(t)
c
      if (abs(h).lt.sqeps) call xerclr
      if (abs(h) .lt. sqeps) call xermsg ('slatec', 'gamic',
     +   'result lt half precision', 1, 1)
c
 80   sgng = sign (1.0, h) * sga * sgngam
      t = log(abs(h)) + algap1 - log(abs(a))
      if (t.lt.bot) call xerclr
      gamic = sgng * exp(t)
      return
c
 90   sgng = -sgngs * sga * sgngam
      t = t + algap1 - log(abs(a))
      if (t.lt.bot) call xerclr
      gamic = sgng * exp(t)
      return
c
      end
