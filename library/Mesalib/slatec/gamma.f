*deck gamma
      function gamma (x)
c***begin prologue  gamma
c***purpose  compute the complete gamma function.
c***library   slatec (fnlib)
c***category  c7a
c***type      single precision (gamma-s, dgamma-d, cgamma-c)
c***keywords  complete gamma function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c gamma computes the gamma function at x, where x is not 0, -1, -2, ....
c gamma and x are single precision.
c
c***references  (none)
c***routines called  csevl, gamlim, inits, r1mach, r9lgmc, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  gamma
      dimension gcs(23)
      logical first
      save gcs, pi, sq2pil, ngcs, xmin, xmax, dxrel, first
      data gcs   ( 1) / .0085711955 90989331e0/
      data gcs   ( 2) / .0044153813 24841007e0/
      data gcs   ( 3) / .0568504368 1599363e0/
      data gcs   ( 4) /-.0042198353 96418561e0/
      data gcs   ( 5) / .0013268081 81212460e0/
      data gcs   ( 6) /-.0001893024 529798880e0/
      data gcs   ( 7) / .0000360692 532744124e0/
      data gcs   ( 8) /-.0000060567 619044608e0/
      data gcs   ( 9) / .0000010558 295463022e0/
      data gcs   (10) /-.0000001811 967365542e0/
      data gcs   (11) / .0000000311 772496471e0/
      data gcs   (12) /-.0000000053 542196390e0/
      data gcs   (13) / .0000000009 193275519e0/
      data gcs   (14) /-.0000000001 577941280e0/
      data gcs   (15) / .0000000000 270798062e0/
      data gcs   (16) /-.0000000000 046468186e0/
      data gcs   (17) / .0000000000 007973350e0/
      data gcs   (18) /-.0000000000 001368078e0/
      data gcs   (19) / .0000000000 000234731e0/
      data gcs   (20) /-.0000000000 000040274e0/
      data gcs   (21) / .0000000000 000006910e0/
      data gcs   (22) /-.0000000000 000001185e0/
      data gcs   (23) / .0000000000 000000203e0/
      data pi /3.14159 26535 89793 24e0/
c sq2pil is log (sqrt (2.*pi) )
      data sq2pil /0.91893 85332 04672 74e0/
      data first /.true./
c
c lanl dependent code removed 81.02.04
c
c***first executable statement  gamma
      if (first) then
c
c ---------------------------------------------------------------------
c initialize.  find legal bounds for x, and determine the number of
c terms in the series required to attain an accuracy ten times better
c than machine precision.
c
         ngcs = inits (gcs, 23, 0.1*r1mach(3))
c
         call gamlim (xmin, xmax)
         dxrel = sqrt (r1mach(4))
c
c ---------------------------------------------------------------------
c finish initialization.  start evaluating gamma(x).
c
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.10.0) go to 50
c
c compute gamma(x) for abs(x) .le. 10.0.  reduce interval and
c find gamma(1+y) for 0. .le. y .lt. 1. first of all.
c
      n = x
      if (x.lt.0.) n = n - 1
      y = x - n
      n = n - 1
      gamma = 0.9375 + csevl(2.*y-1., gcs, ngcs)
      if (n.eq.0) return
c
      if (n.gt.0) go to 30
c
c compute gamma(x) for x .lt. 1.
c
      n = -n
      if (x .eq. 0.) call xermsg ('slatec', 'gamma', 'x is 0', 4, 2)
      if (x .lt. 0. .and. x+n-2 .eq. 0.) call xermsg ('slatec', 'gamma'
     1, 'x is a negative integer', 4, 2)
      if (x .lt. (-0.5) .and. abs((x-aint(x-0.5))/x) .lt. dxrel) call
     1xermsg ( 'slatec', 'gamma',
     2'answer lt half precision because x too near negative integer'
     3, 1, 1)
c
      do 20 i=1,n
        gamma = gamma / (x+i-1)
 20   continue
      return
c
c gamma(x) for x .ge. 2.
c
 30   do 40 i=1,n
        gamma = (y+i)*gamma
 40   continue
      return
c
c compute gamma(x) for abs(x) .gt. 10.0.  recall y = abs(x).
c
 50   if (x .gt. xmax) call xermsg ('slatec', 'gamma',
     +   'x so big gamma overflows', 3, 2)
c
      gamma = 0.
      if (x .lt. xmin) call xermsg ('slatec', 'gamma',
     +   'x so small gamma underflows', 2, 1)
      if (x.lt.xmin) return
c
      gamma = exp((y-0.5)*log(y) - y + sq2pil + r9lgmc(y) )
      if (x.gt.0.) return
c
      if (abs((x-aint(x-0.5))/x) .lt. dxrel) call xermsg ('slatec',
     +   'gamma',
     +   'answer lt half precision, x too near negative integer', 1, 1)
c
      sinpiy = sin (pi*y)
      if (sinpiy .eq. 0.) call xermsg ('slatec', 'gamma',
     +   'x is a negative integer', 4, 2)
c
      gamma = -pi / (y*sinpiy*gamma)
c
      return
      end
