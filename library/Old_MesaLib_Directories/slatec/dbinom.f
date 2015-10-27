*deck dbinom
      double precision function dbinom (n, m)
c***begin prologue  dbinom
c***purpose  compute the binomial coefficients.
c***library   slatec (fnlib)
c***category  c1
c***type      double precision (binom-s, dbinom-d)
c***keywords  binomial coefficients, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c dbinom(n,m) calculates the double precision binomial coefficient
c for integer arguments n and m.  the result is (n!)/((m!)(n-m)!).
c
c***references  (none)
c***routines called  d1mach, d9lgmc, dlnrel, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dbinom
      double precision corr, fintmx, sq2pil, xk, xn, xnk, d9lgmc,
     1  dlnrel, d1mach, bilnmx
      logical first
      save sq2pil, bilnmx, fintmx, first
      data sq2pil / 0.9189385332 0467274178 0329736405 62 d0 /
      data first /.true./
c***first executable statement  dbinom
      if (first) then
         bilnmx = log(d1mach(2)) - 0.0001d0
         fintmx = 0.9d0/d1mach(3)
      endif
      first = .false.
c
      if (n .lt. 0 .or. m .lt. 0) call xermsg ('slatec', 'dbinom',
     +   'n or m lt zero', 1, 2)
      if (n .lt. m) call xermsg ('slatec', 'dbinom', 'n lt m', 2, 2)
c
      k = min (m, n-m)
      if (k.gt.20) go to 30
      if (k*log(amax0(n,1)).gt.bilnmx) go to 30
c
      dbinom = 1.0d0
      if (k.eq.0) return
      do 20 i=1,k
        xn = n - i + 1
        xk = i
        dbinom = dbinom * (xn/xk)
 20   continue
c
      if (dbinom.lt.fintmx) dbinom = aint (dbinom+0.5d0)
      return
c
c if k.lt.9, approx is not valid and answer is close to the overflow lim
 30   if (k .lt. 9) call xermsg ('slatec', 'dbinom',
     +   'result overflows because n and/or m too big', 3, 2)
c
      xn = n + 1
      xk = k + 1
      xnk = n - k + 1
c
      corr = d9lgmc(xn) - d9lgmc(xk) - d9lgmc(xnk)
      dbinom = xk*log(xnk/xk) - xn*dlnrel(-(xk-1.0d0)/xn)
     1  -0.5d0*log(xn*xnk/xk) + 1.0d0 - sq2pil + corr
c
      if (dbinom .gt. bilnmx) call xermsg ('slatec', 'dbinom',
     +   'result overflows because n and/or m too big', 3, 2)
c
      dbinom = exp (dbinom)
      if (dbinom.lt.fintmx) dbinom = aint (dbinom+0.5d0)
c
      return
      end
