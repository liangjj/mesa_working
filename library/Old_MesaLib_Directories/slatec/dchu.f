*deck dchu
      double precision function dchu (a, b, x)
c***begin prologue  dchu
c***purpose  compute the logarithmic confluent hypergeometric function.
c***library   slatec (fnlib)
c***category  c11
c***type      double precision (chu-s, dchu-d)
c***keywords  fnlib, logarithmic confluent hypergeometric function,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c dchu(a,b,x) calculates the double precision logarithmic confluent
c hypergeometric function u(a,b,x) for double precision arguments
c a, b, and x.
c
c this routine is not valid when 1+a-b is close to zero if x is small.
c
c***references  (none)
c***routines called  d1mach, d9chu, dexprl, dgamma, dgamr, dpoch,
c                    dpoch1, xermsg
c***revision history  (yymmdd)
c   770801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900727  added external statement.  (wrb)
c***end prologue  dchu
      double precision a, b, x, aintb, alnx, a0, beps, b0, c0, eps,
     1  factor, gamri1, gamrni, pch1ai, pch1i, pi, pochai, sum, t,
     2  xeps1, xi, xi1, xn, xtoeps,  d1mach, dpoch, dgamma, dgamr,
     3  dpoch1, dexprl, d9chu
      external dgamma
      save pi, eps
      data pi / 3.1415926535 8979323846 2643383279 503 d0 /
      data eps / 0.0d0 /
c***first executable statement  dchu
      if (eps.eq.0.0d0) eps = d1mach(3)
c
      if (x .eq. 0.0d0) call xermsg ('slatec', 'dchu',
     +   'x is zero so dchu is infinite', 1, 2)
      if (x .lt. 0.0d0) call xermsg ('slatec', 'dchu',
     +   'x is negative, use cchu', 2, 2)
c
      if (max(abs(a),1.0d0)*max(abs(1.0d0+a-b),1.0d0).lt.
     1  0.99d0*abs(x)) go to 120
c
c the ascending series will be used, because the descending rational
c approximation (which is based on the asymptotic series) is unstable.
c
      if (abs(1.0d0+a-b) .lt. sqrt(eps)) call xermsg ('slatec', 'dchu',
     +   'algorithmis bad when 1+a-b is near zero for small x', 10, 2)
c
      if (b.ge.0.0d0) aintb = aint(b+0.5d0)
      if (b.lt.0.0d0) aintb = aint(b-0.5d0)
      beps = b - aintb
      n = aintb
c
      alnx = log(x)
      xtoeps = exp (-beps*alnx)
c
c evaluate the finite sum.     -----------------------------------------
c
      if (n.ge.1) go to 40
c
c consider the case b .lt. 1.0 first.
c
      sum = 1.0d0
      if (n.eq.0) go to 30
c
      t = 1.0d0
      m = -n
      do 20 i=1,m
        xi1 = i - 1
        t = t*(a+xi1)*x/((b+xi1)*(xi1+1.0d0))
        sum = sum + t
 20   continue
c
 30   sum = dpoch(1.0d0+a-b, -a)*sum
      go to 70
c
c now consider the case b .ge. 1.0.
c
 40   sum = 0.0d0
      m = n - 2
      if (m.lt.0) go to 70
      t = 1.0d0
      sum = 1.0d0
      if (m.eq.0) go to 60
c
      do 50 i=1,m
        xi = i
        t = t * (a-b+xi)*x/((1.0d0-b+xi)*xi)
        sum = sum + t
 50   continue
c
 60   sum = dgamma(b-1.0d0) * dgamr(a) * x**(1-n) * xtoeps * sum
c
c next evaluate the infinite sum.     ----------------------------------
c
 70   istrt = 0
      if (n.lt.1) istrt = 1 - n
      xi = istrt
c
      factor = (-1.0d0)**n * dgamr(1.0d0+a-b) * x**istrt
      if (beps.ne.0.0d0) factor = factor * beps*pi/sin(beps*pi)
c
      pochai = dpoch (a, xi)
      gamri1 = dgamr (xi+1.0d0)
      gamrni = dgamr (aintb+xi)
      b0 = factor * dpoch(a,xi-beps) * gamrni * dgamr(xi+1.0d0-beps)
c
      if (abs(xtoeps-1.0d0).gt.0.5d0) go to 90
c
c x**(-beps) is close to 1.0d0, so we must be careful in evaluating the
c differences.
c
      pch1ai = dpoch1 (a+xi, -beps)
      pch1i = dpoch1 (xi+1.0d0-beps, beps)
      c0 = factor * pochai * gamrni * gamri1 * (
     1  -dpoch1(b+xi,-beps) + pch1ai - pch1i + beps*pch1ai*pch1i)
c
c xeps1 = (1.0 - x**(-beps))/beps = (x**(-beps) - 1.0)/(-beps)
      xeps1 = alnx*dexprl(-beps*alnx)
c
      dchu = sum + c0 + xeps1*b0
      xn = n
      do 80 i=1,1000
        xi = istrt + i
        xi1 = istrt + i - 1
        b0 = (a+xi1-beps)*b0*x/((xn+xi1)*(xi-beps))
        c0 = (a+xi1)*c0*x/((b+xi1)*xi)
     1    - ((a-1.0d0)*(xn+2.d0*xi-1.0d0) + xi*(xi-beps)) * b0
     2    / (xi*(b+xi1)*(a+xi1-beps))
        t = c0 + xeps1*b0
        dchu = dchu + t
        if (abs(t).lt.eps*abs(dchu)) go to 130
 80   continue
      call xermsg ('slatec', 'dchu',
     +   'no convergence in 1000 terms of the ascending series', 3, 2)
c
c x**(-beps) is very different from 1.0, so the straightforward
c formulation is stable.
c
 90   a0 = factor * pochai * dgamr(b+xi) * gamri1 / beps
      b0 = xtoeps * b0 / beps
c
      dchu = sum + a0 - b0
      do 100 i=1,1000
        xi = istrt + i
        xi1 = istrt + i - 1
        a0 = (a+xi1)*a0*x/((b+xi1)*xi)
        b0 = (a+xi1-beps)*b0*x/((aintb+xi1)*(xi-beps))
        t = a0 - b0
        dchu = dchu + t
        if (abs(t).lt.eps*abs(dchu)) go to 130
 100  continue
      call xermsg ('slatec', 'dchu',
     +   'no convergence in 1000 terms of the ascending series', 3, 2)
c
c use luke-s rational approximation in the asymptotic region.
c
 120  dchu = x**(-a) * d9chu(a,b,x)
c
 130  return
      end
