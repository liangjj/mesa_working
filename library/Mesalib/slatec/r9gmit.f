*deck r9gmit
      function r9gmit (a, x, algap1, sgngam, alx)
c***begin prologue  r9gmit
c***subsidiary
c***purpose  compute tricomi's incomplete gamma function for small
c            arguments.
c***library   slatec (fnlib)
c***category  c7e
c***type      single precision (r9gmit-s, d9gmit-d)
c***keywords  complementary incomplete gamma function, fnlib, small x,
c             special functions, tricomi
c***author  fullerton, w., (lanl)
c***description
c
c compute tricomi's incomplete gamma function for small x.
c
c***references  (none)
c***routines called  alngam, r1mach, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c***end prologue  r9gmit
      save eps, bot
      data eps, bot / 2*0.0 /
c***first executable statement  r9gmit
      if (eps.eq.0.0) eps = 0.5*r1mach(3)
      if (bot.eq.0.0) bot = log(r1mach(1))
c
      if (x .le. 0.0) call xermsg ('slatec', 'r9gmit',
     +   'x should be gt 0', 1, 2)
c
      ma = a + 0.5
      if (a.lt.0.0) ma = a - 0.5
      aeps = a - ma
c
      ae = a
      if (a.lt.(-0.5)) ae = aeps
c
      t = 1.0
      te = ae
      s = t
      do 20 k=1,200
        fk = k
        te = -x*te/fk
        t = te/(ae+fk)
        s = s + t
        if (abs(t).lt.eps*abs(s)) go to 30
 20   continue
      call xermsg ('slatec', 'r9gmit',
     +   'no convergence in 200 terms of taylor-s series', 2, 2)
c
 30   if (a.ge.(-0.5)) algs = -algap1 + log(s)
      if (a.ge.(-0.5)) go to 60
c
      algs = -alngam(1.0+aeps) + log(s)
      s = 1.0
      m = -ma - 1
      if (m.eq.0) go to 50
      t = 1.0
      do 40 k=1,m
        t = x*t/(aeps-m-1+k)
        s = s + t
        if (abs(t).lt.eps*abs(s)) go to 50
 40   continue
c
 50   r9gmit = 0.0
      algs = -ma*log(x) + algs
      if (s.eq.0.0 .or. aeps.eq.0.0) go to 60
c
      sgng2 = sgngam*sign(1.0,s)
      alg2 = -x - algap1 + log(abs(s))
c
      if (alg2.gt.bot) r9gmit = sgng2*exp(alg2)
      if (algs.gt.bot) r9gmit = r9gmit + exp(algs)
      return
c
 60   r9gmit = exp(algs)
      return
c
      end
