*deck sdzro
      subroutine sdzro (ae, f, h, n, nq, iroot, re, t, yh, uround, b, c,
     8   fb, fc, y)
c***begin prologue  sdzro
c***subsidiary
c***purpose  sdzro searches for a zero of a function f(n, t, y, iroot)
c            between the given values b and c until the width of the
c            interval (b, c) has collapsed to within a tolerance
c            specified by the stopping criterion,
c              abs(b - c) .le. 2.*(rw*abs(b) + ae).
c***library   slatec (sdrive)
c***type      single precision (sdzro-s, ddzro-d, cdzro-c)
c***author  kahaner, d. k., (nist)
c             national institute of standards and technology
c             gaithersburg, md  20899
c           sutherland, c. d., (lanl)
c             mail stop d466
c             los alamos national laboratory
c             los alamos, nm  87545
c***description
c
c     this is a special purpose version of zeroin, modified for use with
c     the sdriv package.
c
c     sandia mathematical program library
c     mathematical computing services division 5422
c     sandia laboratories
c     p. o. box 5800
c     albuquerque, new mexico  87115
c     control data 6600 version 4.5, 1 november 1971
c
c     parameters
c        f     - name of the external function, which returns a
c                real result.  this name must be in an
c                external statement in the calling program.
c        b     - one end of the interval (b, c).  the value returned for
c                b usually is the better approximation to a zero of f.
c        c     - the other end of the interval (b, c).
c        re    - relative error used for rw in the stopping criterion.
c                if the requested re is less than machine precision,
c                then rw is set to approximately machine precision.
c        ae    - absolute error used in the stopping criterion.  if the
c                given interval (b, c) contains the origin, then a
c                nonzero value should be chosen for ae.
c
c***references  l. f. shampine and h. a. watts, zeroin, a root-solving
c                 routine, sc-tm-70-631, sept 1970.
c               t. j. dekker, finding a zero by means of successive
c                 linear interpolation, constructive aspects of the
c                 fundamental theorem of algebra, edited by b. dejon
c                 and p. henrici, 1969.
c***routines called  sdntp
c***revision history  (yymmdd)
c   790601  date written
c   900329  initial submission to slatec.
c***end prologue  sdzro
      integer ic, iroot, kount, n, nq
      real a, acbs, acmb, ae, b, c, cmb, er, f, fa, fb, fc,
     8     h, p, q, re, rw, t, tol, uround, y(*), yh(n,*)
c***first executable statement  sdzro
      er = 4.e0*uround
      rw = max(re, er)
      ic = 0
      acbs = abs(b - c)
      a = c
      fa = fc
      kount = 0
c                                                    perform interchange
 10   if (abs(fc) .lt. abs(fb)) then
        a = b
        fa = fb
        b = c
        fb = fc
        c = a
        fc = fa
      end if
      cmb = 0.5e0*(c - b)
      acmb = abs(cmb)
      tol = rw*abs(b) + ae
c                                                test stopping criterion
      if (acmb .le. tol) return
      if (kount .gt. 50) return
c                                    calculate new iterate implicitly as
c                                    b + p/q, where we arrange p .ge. 0.
c                         the implicit form is used to prevent overflow.
      p = (b - a)*fb
      q = fa - fb
      if (p .lt. 0.e0) then
        p = -p
        q = -q
      end if
c                          update a and check for satisfactory reduction
c                          in the size of our bounding interval.
      a = b
      fa = fb
      ic = ic + 1
      if (ic .ge. 4) then
        if (8.e0*acmb .ge. acbs) then
c                                                                 bisect
          b = 0.5e0*(c + b)
          go to 20
        end if
        ic = 0
      end if
      acbs = acmb
c                                            test for too small a change
      if (p .le. abs(q)*tol) then
c                                                 increment by tolerance
        b = b + sign(tol, cmb)
c                                               root ought to be between
c                                               b and (c + b)/2.
      else if (p .lt. cmb*q) then
c                                                            interpolate
        b = b + p/q
      else
c                                                                 bisect
        b = 0.5e0*(c + b)
      end if
c                                             have completed computation
c                                             for new iterate b.
 20   call sdntp (h, 0, n, nq, t, b, yh,  y)
      fb = f(n, b, y, iroot)
      if (n .eq. 0) return
      if (fb .eq. 0.e0) return
      kount = kount + 1
c
c             decide whether next step is interpolation or extrapolation
c
      if (sign(1.0e0, fb) .eq. sign(1.0e0, fc)) then
        c = a
        fc = fa
      end if
      go to 10
      end
