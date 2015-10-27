*deck besj0
      function besj0 (x)
c***begin prologue  besj0
c***purpose  compute the bessel function of the first kind of order
c            zero.
c***library   slatec (fnlib)
c***category  c10a1
c***type      single precision (besj0-s, dbesj0-d)
c***keywords  bessel function, first kind, fnlib, order zero,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c besj0(x) calculates the bessel function of the first kind of
c order zero for real argument x.
c
c series for bj0        on the interval  0.          to  1.60000d+01
c                                        with weighted error   7.47e-18
c                                         log weighted error  17.13
c                               significant figures required  16.98
c                                    decimal places required  17.68
c
c series for bm0        on the interval  0.          to  6.25000d-02
c                                        with weighted error   4.98e-17
c                                         log weighted error  16.30
c                               significant figures required  14.97
c                                    decimal places required  16.96
c
c series for bth0       on the interval  0.          to  6.25000d-02
c                                        with weighted error   3.67e-17
c                                         log weighted error  16.44
c                               significant figures required  15.53
c                                    decimal places required  17.13
c
c***references  (none)
c***routines called  csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890210  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  besj0
      dimension bj0cs(13), bm0cs(21), bth0cs(24)
      logical first
      save bj0cs, bm0cs, bth0cs, pi4, ntj0, ntm0, ntth0, xsml, xmax,
     1   first
      data bj0cs( 1) /    .1002541619 68939137e0 /
      data bj0cs( 2) /   -.6652230077 64405132e0 /
      data bj0cs( 3) /    .2489837034 98281314e0 /
      data bj0cs( 4) /   -.0332527231 700357697e0 /
      data bj0cs( 5) /    .0023114179 304694015e0 /
      data bj0cs( 6) /   -.0000991127 741995080e0 /
      data bj0cs( 7) /    .0000028916 708643998e0 /
      data bj0cs( 8) /   -.0000000612 108586630e0 /
      data bj0cs( 9) /    .0000000009 838650793e0 /
      data bj0cs(10) /   -.0000000000 124235515e0 /
      data bj0cs(11) /    .0000000000 001265433e0 /
      data bj0cs(12) /   -.0000000000 000010619e0 /
      data bj0cs(13) /    .0000000000 000000074e0 /
      data bm0cs( 1) /    .0928496163 7381644e0 /
      data bm0cs( 2) /   -.0014298770 7403484e0 /
      data bm0cs( 3) /    .0000283057 9271257e0 /
      data bm0cs( 4) /   -.0000014330 0611424e0 /
      data bm0cs( 5) /    .0000001202 8628046e0 /
      data bm0cs( 6) /   -.0000000139 7113013e0 /
      data bm0cs( 7) /    .0000000020 4076188e0 /
      data bm0cs( 8) /   -.0000000003 5399669e0 /
      data bm0cs( 9) /    .0000000000 7024759e0 /
      data bm0cs(10) /   -.0000000000 1554107e0 /
      data bm0cs(11) /    .0000000000 0376226e0 /
      data bm0cs(12) /   -.0000000000 0098282e0 /
      data bm0cs(13) /    .0000000000 0027408e0 /
      data bm0cs(14) /   -.0000000000 0008091e0 /
      data bm0cs(15) /    .0000000000 0002511e0 /
      data bm0cs(16) /   -.0000000000 0000814e0 /
      data bm0cs(17) /    .0000000000 0000275e0 /
      data bm0cs(18) /   -.0000000000 0000096e0 /
      data bm0cs(19) /    .0000000000 0000034e0 /
      data bm0cs(20) /   -.0000000000 0000012e0 /
      data bm0cs(21) /    .0000000000 0000004e0 /
      data bth0cs( 1) /   -.2463916377 4300119e0 /
      data bth0cs( 2) /    .0017370983 07508963e0 /
      data bth0cs( 3) /   -.0000621836 33402968e0 /
      data bth0cs( 4) /    .0000043680 50165742e0 /
      data bth0cs( 5) /   -.0000004560 93019869e0 /
      data bth0cs( 6) /    .0000000621 97400101e0 /
      data bth0cs( 7) /   -.0000000103 00442889e0 /
      data bth0cs( 8) /    .0000000019 79526776e0 /
      data bth0cs( 9) /   -.0000000004 28198396e0 /
      data bth0cs(10) /    .0000000001 02035840e0 /
      data bth0cs(11) /   -.0000000000 26363898e0 /
      data bth0cs(12) /    .0000000000 07297935e0 /
      data bth0cs(13) /   -.0000000000 02144188e0 /
      data bth0cs(14) /    .0000000000 00663693e0 /
      data bth0cs(15) /   -.0000000000 00215126e0 /
      data bth0cs(16) /    .0000000000 00072659e0 /
      data bth0cs(17) /   -.0000000000 00025465e0 /
      data bth0cs(18) /    .0000000000 00009229e0 /
      data bth0cs(19) /   -.0000000000 00003448e0 /
      data bth0cs(20) /    .0000000000 00001325e0 /
      data bth0cs(21) /   -.0000000000 00000522e0 /
      data bth0cs(22) /    .0000000000 00000210e0 /
      data bth0cs(23) /   -.0000000000 00000087e0 /
      data bth0cs(24) /    .0000000000 00000036e0 /
      data pi4 / 0.7853981633 9744831e0 /
      data first /.true./
c***first executable statement  besj0
      if (first) then
         ntj0 = inits (bj0cs, 13, 0.1*r1mach(3))
         ntm0 = inits (bm0cs, 21, 0.1*r1mach(3))
         ntth0 = inits (bth0cs, 24, 0.1*r1mach(3))
c
         xsml = sqrt (8.0*r1mach(3))
         xmax = 1.0/r1mach(4)
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.4.0) go to 20
c
      besj0 = 1.0
      if (y.gt.xsml) besj0 = csevl (.125*y*y-1., bj0cs, ntj0)
      return
c
 20   if (y .gt. xmax) call xermsg ('slatec', 'besj0',
     +   'no precision because abs(x) is too big', 1, 2)
c
      z = 32.0/y**2 - 1.0
      ampl = (0.75 + csevl (z, bm0cs, ntm0)) / sqrt(y)
      theta = y - pi4 + csevl (z, bth0cs, ntth0) / y
      besj0 = ampl * cos (theta)
c
      return
      end
