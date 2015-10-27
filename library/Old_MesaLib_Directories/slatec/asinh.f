*deck asinh
      function asinh (x)
c***begin prologue  asinh
c***purpose  compute the arc hyperbolic sine.
c***library   slatec (fnlib)
c***category  c4c
c***type      single precision (asinh-s, dasinh-d, casinh-c)
c***keywords  arc hyperbolic sine, asinh, elementary functions, fnlib,
c             inverse hyperbolic sine
c***author  fullerton, w., (lanl)
c***description
c
c asinh(x) computes the arc hyperbolic sine of x.
c
c series for asnh       on the interval  0.          to  1.00000d+00
c                                        with weighted error   2.19e-17
c                                         log weighted error  16.66
c                               significant figures required  15.60
c                                    decimal places required  17.31
c
c***references  (none)
c***routines called  csevl, inits, r1mach
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  asinh
      dimension asnhcs(20)
      logical first
      save aln2, asnhcs, nterms, xmax, sqeps, first
      data aln2 /0.6931471805 5994530942e0/
      data asnhcs( 1) /   -.1282003991 1738186e0 /
      data asnhcs( 2) /   -.0588117611 89951768e0 /
      data asnhcs( 3) /    .0047274654 32212481e0 /
      data asnhcs( 4) /   -.0004938363 16265361e0 /
      data asnhcs( 5) /    .0000585062 07058557e0 /
      data asnhcs( 6) /   -.0000074669 98328931e0 /
      data asnhcs( 7) /    .0000010011 69358355e0 /
      data asnhcs( 8) /   -.0000001390 35438587e0 /
      data asnhcs( 9) /    .0000000198 23169483e0 /
      data asnhcs(10) /   -.0000000028 84746841e0 /
      data asnhcs(11) /    .0000000004 26729654e0 /
      data asnhcs(12) /   -.0000000000 63976084e0 /
      data asnhcs(13) /    .0000000000 09699168e0 /
      data asnhcs(14) /   -.0000000000 01484427e0 /
      data asnhcs(15) /    .0000000000 00229037e0 /
      data asnhcs(16) /   -.0000000000 00035588e0 /
      data asnhcs(17) /    .0000000000 00005563e0 /
      data asnhcs(18) /   -.0000000000 00000874e0 /
      data asnhcs(19) /    .0000000000 00000138e0 /
      data asnhcs(20) /   -.0000000000 00000021e0 /
      data first /.true./
c***first executable statement  asinh
      if (first) then
         nterms = inits (asnhcs, 20, 0.1*r1mach(3))
         sqeps = sqrt (r1mach(3))
         xmax = 1.0/sqeps
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.1.0) go to 20
c
      asinh = x
      if (y.gt.sqeps) asinh = x*(1.0 + csevl (2.*x*x-1., asnhcs,nterms))
      return
c
 20   if (y.lt.xmax) asinh = log (y + sqrt(y**2+1.))
      if (y.ge.xmax) asinh = aln2 + log(y)
      asinh = sign (asinh, x)
c
      return
      end
