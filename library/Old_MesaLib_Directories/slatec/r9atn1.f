*deck r9atn1
      function r9atn1 (x)
c***begin prologue  r9atn1
c***subsidiary
c***purpose  evaluate atan(x) from first order relative accuracy so that
c            atan(x) = x + x**3*r9atn1(x).
c***library   slatec (fnlib)
c***category  c4a
c***type      single precision (r9atn1-s, d9atn1-d)
c***keywords  arc tangent, elementary functions, first order, fnlib,
c             trigonometric
c***author  fullerton, w., (lanl)
c***description
c
c evaluate  atan(x)  from first order, that is, evaluate
c (atan(x)-x)/x**3  with relative error accuracy so that
c        atan(x) = x + x**3*r9atn1(x).
c
c series for atn1       on the interval  0.          to  1.00000d+00
c                                        with weighted error   2.21e-17
c                                         log weighted error  16.66
c                               significant figures required  15.44
c                                    decimal places required  17.32
c
c***references  (none)
c***routines called  csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   780401  date written
c   890206  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c***end prologue  r9atn1
      dimension atn1cs(21)
      logical first
      save atn1cs, ntatn1, xsml, xbig, xmax, first
      data atn1cs( 1) /   -.0328399753 5355202e0 /
      data atn1cs( 2) /    .0583343234 3172412e0 /
      data atn1cs( 3) /   -.0074003696 9671964e0 /
      data atn1cs( 4) /    .0010097841 9933728e0 /
      data atn1cs( 5) /   -.0001439787 1635652e0 /
      data atn1cs( 6) /    .0000211451 2648992e0 /
      data atn1cs( 7) /   -.0000031723 2107425e0 /
      data atn1cs( 8) /    .0000004836 6203654e0 /
      data atn1cs( 9) /   -.0000000746 7746546e0 /
      data atn1cs(10) /    .0000000116 4800896e0 /
      data atn1cs(11) /   -.0000000018 3208837e0 /
      data atn1cs(12) /    .0000000002 9019082e0 /
      data atn1cs(13) /   -.0000000000 4623885e0 /
      data atn1cs(14) /    .0000000000 0740552e0 /
      data atn1cs(15) /   -.0000000000 0119135e0 /
      data atn1cs(16) /    .0000000000 0019240e0 /
      data atn1cs(17) /   -.0000000000 0003118e0 /
      data atn1cs(18) /    .0000000000 0000506e0 /
      data atn1cs(19) /   -.0000000000 0000082e0 /
      data atn1cs(20) /    .0000000000 0000013e0 /
      data atn1cs(21) /   -.0000000000 0000002e0 /
      data first /.true./
c***first executable statement  r9atn1
      if (first) then
         eps = r1mach(3)
         ntatn1 = inits (atn1cs, 21, 0.1*eps)
c
         xsml = sqrt (0.1*eps)
         xbig = 1.571/sqrt(eps)
         xmax = 1.571/eps
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.1.0) go to 20
c
      if (y.le.xsml) r9atn1 = -1.0/3.0
      if (y.le.xsml) return
c
      r9atn1 = -0.25 + csevl (2.0*y*y-1., atn1cs, ntatn1)
      return
c
 20   if (y .gt. xmax) call xermsg ('slatec', 'r9atn1',
     +   'no precision in answer because x is too big', 2, 2)
      if (y .gt. xbig) call xermsg ('slatec', 'r9atn1',
     +   'answer lt half precision because x is too big', 1, 1)
c
      r9atn1 = (atan(x) - x) / x**3
      return
c
      end
