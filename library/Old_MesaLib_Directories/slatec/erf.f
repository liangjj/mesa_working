*deck erf
      function erf (x)
c***begin prologue  erf
c***purpose  compute the error function.
c***library   slatec (fnlib)
c***category  c8a, l5a1e
c***type      single precision (erf-s, derf-d)
c***keywords  erf, error function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c erf(x) calculates the single precision error function for
c single precision argument x.
c
c series for erf        on the interval  0.          to  1.00000d+00
c                                        with weighted error   7.10e-18
c                                         log weighted error  17.15
c                               significant figures required  16.31
c                                    decimal places required  17.71
c
c***references  (none)
c***routines called  csevl, erfc, inits, r1mach
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900727  added external statement.  (wrb)
c   920618  removed space from variable name.  (rwc, wrb)
c***end prologue  erf
      dimension erfcs(13)
      logical first
      external erfc
      save erfcs, sqrtpi, nterf, xbig, sqeps, first
      data erfcs( 1) /   -.0490461212 34691808e0 /
      data erfcs( 2) /   -.1422612051 0371364e0 /
      data erfcs( 3) /    .0100355821 87599796e0 /
      data erfcs( 4) /   -.0005768764 69976748e0 /
      data erfcs( 5) /    .0000274199 31252196e0 /
      data erfcs( 6) /   -.0000011043 17550734e0 /
      data erfcs( 7) /    .0000000384 88755420e0 /
      data erfcs( 8) /   -.0000000011 80858253e0 /
      data erfcs( 9) /    .0000000000 32334215e0 /
      data erfcs(10) /   -.0000000000 00799101e0 /
      data erfcs(11) /    .0000000000 00017990e0 /
      data erfcs(12) /   -.0000000000 00000371e0 /
      data erfcs(13) /    .0000000000 00000007e0 /
      data sqrtpi /1.772453850 9055160e0/
      data first /.true./
c***first executable statement  erf
      if (first) then
         nterf = inits (erfcs, 13, 0.1*r1mach(3))
         xbig = sqrt(-log(sqrtpi*r1mach(3)))
         sqeps = sqrt(2.0*r1mach(3))
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.1.) go to 20
c
c erf(x) = 1. - erfc(x) for -1. .le. x .le. 1.
c
      if (y.le.sqeps) erf = 2.0*x/sqrtpi
      if (y.gt.sqeps) erf = x*(1.0 + csevl(2.*x**2-1., erfcs, nterf))
      return
c
c erf(x) = 1. - erfc(x) for  abs(x) .gt. 1.
c
 20   if (y.le.xbig) erf = sign (1.0-erfc(y), x)
      if (y.gt.xbig) erf = sign (1.0, x)
c
      return
      end
