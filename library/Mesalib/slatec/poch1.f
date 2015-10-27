*deck poch1
      function poch1 (a, x)
c***begin prologue  poch1
c***purpose  calculate a generalization of pochhammer's symbol starting
c            from first order.
c***library   slatec (fnlib)
c***category  c1, c7a
c***type      single precision (poch1-s, dpoch1-d)
c***keywords  first order, fnlib, pochhammer, special functions
c***author  fullerton, w., (lanl)
c***description
c
c evaluate a generalization of pochhammer's symbol for special
c situations that require especially accurate values when x is small in
c        poch1(a,x) = (poch(a,x)-1)/x
c                   = (gamma(a+x)/gamma(a) - 1.0)/x .
c this specification is particularly suited for stably computing
c expressions such as
c        (gamma(a+x)/gamma(a) - gamma(b+x)/gamma(b))/x
c             = poch1(a,x) - poch1(b,x)
c note that poch1(a,0.0) = psi(a)
c
c when abs(x) is so small that substantial cancellation will occur if
c the straightforward formula is used, we  use an expansion due
c to fields and discussed by y. l. luke, the special functions and their
c approximations, vol. 1, academic press, 1969, page 34.
c
c the ratio poch(a,x) = gamma(a+x)/gamma(a) is written by luke as
c        (a+(x-1)/2)**x * polynomial in (a+(x-1)/2)**(-2) .
c in order to maintain significance in poch1, we write for positive a
c        (a+(x-1)/2)**x = exp(x*log(a+(x-1)/2)) = exp(q)
c                       = 1.0 + q*exprel(q) .
c likewise the polynomial is written
c        poly = 1.0 + x*poly1(a,x) .
c thus,
c        poch1(a,x) = (poch(a,x) - 1) / x
c                   = exprel(q)*(q/x + q*poly1(a,x)) + poly1(a,x)
c
c***references  (none)
c***routines called  cot, exprel, poch, psi, r1mach, xermsg
c***revision history  (yymmdd)
c   770801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900727  added external statement.  (wrb)
c***end prologue  poch1
      dimension bern(9), gbern(10)
      logical first
      external cot
      save bern, pi, sqtbig, alneps, first
      data bern( 1) /   .8333333333 3333333e-01 /
      data bern( 2) /  -.1388888888 8888889e-02 /
      data bern( 3) /   .3306878306 8783069e-04 /
      data bern( 4) /  -.8267195767 1957672e-06 /
      data bern( 5) /   .2087675698 7868099e-07 /
      data bern( 6) /  -.5284190138 6874932e-09 /
      data bern( 7) /   .1338253653 0684679e-10 /
      data bern( 8) /  -.3389680296 3225829e-12 /
      data bern( 9) /   .8586062056 2778446e-14 /
      data pi / 3.1415926535 8979324 e0 /
      data first /.true./
c***first executable statement  poch1
      if (first) then
         sqtbig = 1.0/sqrt(24.0*r1mach(1))
         alneps = log(r1mach(3))
      endif
      first = .false.
c
      if (x.eq.0.0) poch1 = psi(a)
      if (x.eq.0.0) return
c
      absx = abs(x)
      absa = abs(a)
      if (absx.gt.0.1*absa) go to 70
      if (absx*log(max(absa,2.0)).gt.0.1) go to 70
c
      bp = a
      if (a.lt.(-0.5)) bp = 1.0 - a - x
      incr = 0
      if (bp.lt.10.0) incr = 11.0 - bp
      b = bp + incr
c
      var = b + 0.5*(x-1.0)
      alnvar = log(var)
      q = x*alnvar
c
      poly1 = 0.0
      if (var.ge.sqtbig) go to 40
      var2 = (1.0/var)**2
c
      rho = 0.5*(x+1.0)
      gbern(1) = 1.0
      gbern(2) = -rho/12.0
      term = var2
      poly1 = gbern(2)*term
c
      nterms = -0.5*alneps/alnvar + 1.0
      if (nterms .gt. 9) call xermsg ('slatec', 'poch1',
     +   'nterms is too big, maybe r1mach(3) is bad', 1, 2)
      if (nterms.lt.2) go to 40
c
      do 30 k=2,nterms
        gbk = 0.0
        do 20 j=1,k
          ndx = k - j + 1
          gbk = gbk + bern(ndx)*gbern(j)
 20     continue
        gbern(k+1) = -rho*gbk/k
c
        term = term * (2*k-2.-x)*(2*k-1.-x)*var2
        poly1 = poly1 + gbern(k+1)*term
 30   continue
c
 40   poly1 = (x-1.0)*poly1
      poch1 = exprel(q)*(alnvar + q*poly1) + poly1
c
      if (incr.eq.0) go to 60
c
c we have poch1(b,x).  but bp is small, so we use backwards recursion
c to obtain poch1(bp,x).
c
      do 50 ii=1,incr
        i = incr - ii
        binv = 1.0/(bp+i)
        poch1 = (poch1-binv)/(1.0+x*binv)
 50   continue
c
 60   if (bp.eq.a) return
c
c we have poch1(bp,x), but a is lt -0.5.  we therefore use a reflection
c formula to obtain poch1(a,x).
c
      sinpxx = sin(pi*x)/x
      sinpx2 = sin(0.5*pi*x)
      trig = sinpxx*cot(pi*b) - 2.0*sinpx2*(sinpx2/x)
c
      poch1 = trig + (1.0 + x*trig) * poch1
      return
c
 70   poch1 = (poch(a,x) - 1.0) / x
      return
c
      end
