*deck r9lgit
      function r9lgit (a, x, algap1)
c***begin prologue  r9lgit
c***subsidiary
c***purpose  compute the logarithm of tricomi's incomplete gamma
c            function with perron's continued fraction for large x and
c            a .ge. x.
c***library   slatec (fnlib)
c***category  c7e
c***type      single precision (r9lgit-s, d9lgit-d)
c***keywords  fnlib, incomplete gamma function, logarithm,
c             perron's continued fraction, special functions, tricomi
c***author  fullerton, w., (lanl)
c***description
c
c compute the log of tricomi's incomplete gamma function with perron's
c continued fraction for large x and for a .ge. x.
c
c***references  (none)
c***routines called  r1mach, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c***end prologue  r9lgit
      save eps, sqeps
      data eps, sqeps / 2*0.0 /
c***first executable statement  r9lgit
      if (eps.eq.0.0) eps = 0.5*r1mach(3)
      if (sqeps.eq.0.0) sqeps = sqrt(r1mach(4))
c
      if (x .le. 0.0 .or. a .lt. x) call xermsg ('slatec', 'r9lgit',
     +   'x should be gt 0.0 and le a', 2, 2)
c
      ax = a + x
      a1x = ax + 1.0
      r = 0.0
      p = 1.0
      s = p
      do 20 k=1,200
        fk = k
        t = (a+fk)*x*(1.0+r)
        r = t/((ax+fk)*(a1x+fk)-t)
        p = r*p
        s = s + p
        if (abs(p).lt.eps*s) go to 30
 20   continue
      call xermsg ('slatec', 'r9lgit',
     +   'no convergence in 200 terms of continued fraction', 3, 2)
c
 30   hstar = 1.0 - x*s/a1x
      if (hstar .lt. sqeps) call xermsg ('slatec', 'r9lgit',
     +   'result less than half precision', 1, 1)
c
      r9lgit = -x - algap1 - log(hstar)
c
      return
      end
