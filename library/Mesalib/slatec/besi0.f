*deck besi0
      function besi0 (x)
c***begin prologue  besi0
c***purpose  compute the hyperbolic bessel function of the first kind
c            of order zero.
c***library   slatec (fnlib)
c***category  c10b1
c***type      single precision (besi0-s, dbesi0-d)
c***keywords  first kind, fnlib, hyperbolic bessel function,
c             modified bessel function, order zero, special functions
c***author  fullerton, w., (lanl)
c***description
c
c besi0(x) computes the modified (hyperbolic) bessel function
c of the first kind of order zero and real argument x.
c
c series for bi0        on the interval  0.          to  9.00000d+00
c                                        with weighted error   2.46e-18
c                                         log weighted error  17.61
c                               significant figures required  17.90
c                                    decimal places required  18.15
c
c***references  (none)
c***routines called  besi0e, csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  besi0
      dimension bi0cs(12)
      logical first
      save bi0cs, nti0, xsml, xmax, first
      data bi0cs( 1) /   -.0766054725 2839144951e0 /
      data bi0cs( 2) /   1.9273379539 93808270e0 /
      data bi0cs( 3) /    .2282644586 920301339e0 /
      data bi0cs( 4) /    .0130489146 6707290428e0 /
      data bi0cs( 5) /    .0004344270 9008164874e0 /
      data bi0cs( 6) /    .0000094226 5768600193e0 /
      data bi0cs( 7) /    .0000001434 0062895106e0 /
      data bi0cs( 8) /    .0000000016 1384906966e0 /
      data bi0cs( 9) /    .0000000000 1396650044e0 /
      data bi0cs(10) /    .0000000000 0009579451e0 /
      data bi0cs(11) /    .0000000000 0000053339e0 /
      data bi0cs(12) /    .0000000000 0000000245e0 /
      data first /.true./
c***first executable statement  besi0
      if (first) then
         nti0 = inits (bi0cs, 12, 0.1*r1mach(3))
         xsml = sqrt (4.5*r1mach(3))
         xmax = log (r1mach(2))
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.3.0) go to 20
c
      besi0 = 1.0
      if (y.gt.xsml) besi0 = 2.75 + csevl (y*y/4.5-1.0, bi0cs, nti0)
      return
c
 20   if (y .gt. xmax) call xermsg ('slatec', 'besi0',
     +   'abs(x) so big i0 overflows', 1, 2)
c
      besi0 = exp(y) * besi0e(x)
c
      return
      end
