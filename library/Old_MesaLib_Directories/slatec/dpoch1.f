*deck dpoch1
      double precision function dpoch1 (a, x)
c***begin prologue  dpoch1
c***purpose  calculate a generalization of pochhammer's symbol starting
c            from first order.
c***library   slatec (fnlib)
c***category  c1, c7a
c***type      double precision (poch1-s, dpoch1-d)
c***keywords  first order, fnlib, pochhammer, special functions
c***author  fullerton, w., (lanl)
c***description
c
c evaluate a double precision generalization of pochhammer's symbol
c for double precision a and x for special situations that require
c especially accurate values when x is small in
c        poch1(a,x) = (poch(a,x)-1)/x
c                   = (gamma(a+x)/gamma(a) - 1.0)/x .
c this specification is particularly suited for stably computing
c expressions such as
c        (gamma(a+x)/gamma(a) - gamma(b+x)/gamma(b))/x
c             = poch1(a,x) - poch1(b,x)
c note that poch1(a,0.0) = psi(a)
c
c when abs(x) is so small that substantial cancellation will occur if
c the straightforward formula is used, we use an expansion due
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
c***routines called  d1mach, dcot, dexprl, dpoch, dpsi, xermsg
c***revision history  (yymmdd)
c   770801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900727  added external statement.  (wrb)
c***end prologue  dpoch1
      double precision a, x, absa, absx, alneps, alnvar, b, bern(20),
     1  binv, bp, gbern(21), gbk, pi, poly1, q, rho, sinpxx, sinpx2,
     2  sqtbig, term, trig, var, var2, d1mach, dpsi, dexprl, dcot, dpoch
      logical first
      external dcot
      save bern, pi, sqtbig, alneps, first
      data bern  (  1) / +.8333333333 3333333333 3333333333 333 d-1    /
      data bern  (  2) / -.1388888888 8888888888 8888888888 888 d-2    /
      data bern  (  3) / +.3306878306 8783068783 0687830687 830 d-4    /
      data bern  (  4) / -.8267195767 1957671957 6719576719 576 d-6    /
      data bern  (  5) / +.2087675698 7868098979 2100903212 014 d-7    /
      data bern  (  6) / -.5284190138 6874931848 4768220217 955 d-9    /
      data bern  (  7) / +.1338253653 0684678832 8269809751 291 d-10   /
      data bern  (  8) / -.3389680296 3225828668 3019539124 944 d-12   /
      data bern  (  9) / +.8586062056 2778445641 3590545042 562 d-14   /
      data bern  ( 10) / -.2174868698 5580618730 4151642386 591 d-15   /
      data bern  ( 11) / +.5509002828 3602295152 0265260890 225 d-17   /
      data bern  ( 12) / -.1395446468 5812523340 7076862640 635 d-18   /
      data bern  ( 13) / +.3534707039 6294674716 9322997780 379 d-20   /
      data bern  ( 14) / -.8953517427 0375468504 0261131811 274 d-22   /
      data bern  ( 15) / +.2267952452 3376830603 1095073886 816 d-23   /
      data bern  ( 16) / -.5744724395 2026452383 4847971943 400 d-24   /
      data bern  ( 17) / +.1455172475 6148649018 6626486727 132 d-26   /
      data bern  ( 18) / -.3685994940 6653101781 8178247990 866 d-28   /
      data bern  ( 19) / +.9336734257 0950446720 3255515278 562 d-30   /
      data bern  ( 20) / -.2365022415 7006299345 5963519636 983 d-31   /
      data pi / 3.1415926535 8979323846 2643383279 503 d0 /
      data first /.true./
c***first executable statement  dpoch1
      if (first) then
         sqtbig = 1.0d0/sqrt(24.0d0*d1mach(1))
         alneps = log(d1mach(3))
      endif
      first = .false.
c
      if (x.eq.0.0d0) dpoch1 = dpsi(a)
      if (x.eq.0.0d0) return
c
      absx = abs(x)
      absa = abs(a)
      if (absx.gt.0.1d0*absa) go to 70
      if (absx*log(max(absa,2.0d0)).gt.0.1d0) go to 70
c
      bp = a
      if (a.lt.(-0.5d0)) bp = 1.0d0 - a - x
      incr = 0
      if (bp.lt.10.0d0) incr = 11.0d0 - bp
      b = bp + incr
c
      var = b + 0.5d0*(x-1.0d0)
      alnvar = log(var)
      q = x*alnvar
c
      poly1 = 0.0d0
      if (var.ge.sqtbig) go to 40
      var2 = (1.0d0/var)**2
c
      rho = 0.5d0*(x+1.0d0)
      gbern(1) = 1.0d0
      gbern(2) = -rho/12.0d0
      term = var2
      poly1 = gbern(2)*term
c
      nterms = -0.5d0*alneps/alnvar + 1.0d0
      if (nterms .gt. 20) call xermsg ('slatec', 'dpoch1',
     +   'nterms is too big, maybe d1mach(3) is bad', 1, 2)
      if (nterms.lt.2) go to 40
c
      do 30 k=2,nterms
        gbk = 0.0d0
        do 20 j=1,k
          ndx = k - j + 1
          gbk = gbk + bern(ndx)*gbern(j)
 20     continue
        gbern(k+1) = -rho*gbk/k
c
        term = term * (2*k-2-x)*(2*k-1-x)*var2
        poly1 = poly1 + gbern(k+1)*term
 30   continue
c
 40   poly1 = (x-1.0d0)*poly1
      dpoch1 = dexprl(q)*(alnvar+q*poly1) + poly1
c
      if (incr.eq.0) go to 60
c
c we have dpoch1(b,x), but bp is small, so we use backwards recursion
c to obtain dpoch1(bp,x).
c
      do 50 ii=1,incr
        i = incr - ii
        binv = 1.0d0/(bp+i)
        dpoch1 = (dpoch1 - binv) / (1.0d0 + x*binv)
 50   continue
c
 60   if (bp.eq.a) return
c
c we have dpoch1(bp,x), but a is lt -0.5.  we therefore use a reflection
c formula to obtain dpoch1(a,x).
c
      sinpxx = sin(pi*x)/x
      sinpx2 = sin(0.5d0*pi*x)
      trig = sinpxx*dcot(pi*b) - 2.0d0*sinpx2*(sinpx2/x)
c
      dpoch1 = trig + (1.0d0 + x*trig)*dpoch1
      return
c
 70   dpoch1 = (dpoch(a,x) - 1.0d0) / x
      return
c
      end
