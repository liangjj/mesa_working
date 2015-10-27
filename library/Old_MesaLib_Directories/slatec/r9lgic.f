*deck r9lgic
      function r9lgic (a, x, alx)
c***begin prologue  r9lgic
c***subsidiary
c***purpose  compute the log complementary incomplete gamma function
c            for large x and for a .le. x.
c***library   slatec (fnlib)
c***category  c7e
c***type      single precision (r9lgic-s, d9lgic-d)
c***keywords  complementary incomplete gamma function, fnlib, large x,
c             logarithm, special functions
c***author  fullerton, w., (lanl)
c***description
c
c compute the log complementary incomplete gamma function for large x
c and for a .le. x.
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
c***end prologue  r9lgic
      save eps
      data eps / 0.0 /
c***first executable statement  r9lgic
      if (eps.eq.0.0) eps = 0.5*r1mach(3)
c
      xpa = x + 1.0 - a
      xma = x - 1.0 - a
c
      r = 0.0
      p = 1.0
      s = p
      do 10 k=1,200
        fk = k
        t = fk*(a-fk)*(1.0+r)
        r = -t/((xma+2.0*fk)*(xpa+2.0*fk)+t)
        p = r*p
        s = s + p
        if (abs(p).lt.eps*s) go to 20
 10   continue
      call xermsg ('slatec', 'r9lgic',
     +   'no convergence in 200 terms of continued fraction', 1, 2)
c
 20   r9lgic = a*alx - x + log(s/xpa)
c
      return
      end
