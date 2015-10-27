*deck r9lgmc
      function r9lgmc (x)
c***begin prologue  r9lgmc
c***subsidiary
c***purpose  compute the log gamma correction factor so that
c            log(gamma(x)) = log(sqrt(2*pi)) + (x-.5)*log(x) - x
c            + r9lgmc(x).
c***library   slatec (fnlib)
c***category  c7e
c***type      single precision (r9lgmc-s, d9lgmc-d, c9lgmc-c)
c***keywords  complete gamma function, correction term, fnlib,
c             log gamma, logarithm, special functions
c***author  fullerton, w., (lanl)
c***description
c
c compute the log gamma correction factor for x .ge. 10.0 so that
c  log (gamma(x)) = log(sqrt(2*pi)) + (x-.5)*log(x) - x + r9lgmc(x)
c
c series for algm       on the interval  0.          to  1.00000d-02
c                                        with weighted error   3.40e-16
c                                         log weighted error  15.47
c                               significant figures required  14.39
c                                    decimal places required  15.86
c
c***references  (none)
c***routines called  csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c***end prologue  r9lgmc
      dimension algmcs(6)
      logical first
      save algmcs, nalgm, xbig, xmax, first
      data algmcs( 1) /    .1666389480 45186e0 /
      data algmcs( 2) /   -.0000138494 817606e0 /
      data algmcs( 3) /    .0000000098 108256e0 /
      data algmcs( 4) /   -.0000000000 180912e0 /
      data algmcs( 5) /    .0000000000 000622e0 /
      data algmcs( 6) /   -.0000000000 000003e0 /
      data first /.true./
c***first executable statement  r9lgmc
      if (first) then
         nalgm = inits (algmcs, 6, r1mach(3))
         xbig = 1.0/sqrt(r1mach(3))
         xmax = exp (min(log(r1mach(2)/12.0), -log(12.0*r1mach(1))) )
      endif
      first = .false.
c
      if (x .lt. 10.0) call xermsg ('slatec', 'r9lgmc',
     +   'x must be ge 10', 1, 2)
      if (x.ge.xmax) go to 20
c
      r9lgmc = 1.0/(12.0*x)
      if (x.lt.xbig) r9lgmc = csevl (2.0*(10./x)**2-1., algmcs, nalgm)/x
      return
c
 20   r9lgmc = 0.0
      call xermsg ('slatec', 'r9lgmc', 'x so big r9lgmc underflows', 2,
     +   1)
      return
c
      end
