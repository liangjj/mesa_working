*deck dpchsw
      subroutine dpchsw (dfmax, iextrm, d1, d2, h, slope, ierr)
c***begin prologue  dpchsw
c***subsidiary
c***purpose  limits excursion from data for dpchcs
c***library   slatec (pchip)
c***type      double precision (pchsw-s, dpchsw-d)
c***author  fritsch, f. n., (llnl)
c***description
c
c         dpchsw:  dpchcs switch excursion limiter.
c
c     called by  dpchcs  to adjust d1 and d2 if necessary to insure that
c     the extremum on this interval is not further than dfmax from the
c     extreme data value.
c
c ----------------------------------------------------------------------
c
c  calling sequence:
c
c        integer  iextrm, ierr
c        double precision  dfmax, d1, d2, h, slope
c
c        call  dpchsw (dfmax, iextrm, d1, d2, h, slope, ierr)
c
c   parameters:
c
c     dfmax -- (input) maximum allowed difference between f(iextrm) and
c           the cubic determined by derivative values d1,d2.  (assumes
c           dfmax.gt.0.)
c
c     iextrm -- (input) index of the extreme data value.  (assumes
c           iextrm = 1 or 2 .  any value .ne.1 is treated as 2.)
c
c     d1,d2 -- (input) derivative values at the ends of the interval.
c           (assumes d1*d2 .le. 0.)
c          (output) may be modified if necessary to meet the restriction
c           imposed by dfmax.
c
c     h -- (input) interval length.  (assumes  h.gt.0.)
c
c     slope -- (input) data slope on the interval.
c
c     ierr -- (output) error flag.  should be zero.
c           if ierr=-1, assumption on d1 and d2 is not satisfied.
c           if ierr=-2, quadratic equation locating extremum has
c                       negative discriminant (should never occur).
c
c    -------
c    warning:  this routine does no validity-checking of arguments.
c    -------
c
c  fortran intrinsics used:  abs, sign, sqrt.
c
c***see also  dpchcs
c***routines called  d1mach, xermsg
c***revision history  (yymmdd)
c   820218  date written
c   820805  converted to slatec library version.
c   870707  corrected xerror calls for d.p. name(s).
c   870707  replaced data statement for small with a use of d1mach.
c   870813  minor cosmetic changes.
c   890206  corrected xerror calls.
c   890411  1. added save statements (vers. 3.2).
c           2. added double precision declaration for d1mach.
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   910408  updated author and date written sections in prologue.  (wrb)
c   920526  eliminated possible divide by zero problem.  (fnf)
c   930503  improved purpose.  (fnf)
c***end prologue  dpchsw
c
c**end
c
c  declare arguments.
c
      integer  iextrm, ierr
      double precision  dfmax, d1, d2, h, slope
c
c  declare local variables.
c
      double precision  cp, fact, hphi, lambda, nu, one, phi, radcal,
     *                  rho, sigma, small, that, third, three, two, zero
      save zero, one, two, three, fact
      save third
      double precision  d1mach
c
      data  zero /0.d0/,  one /1.d0/,  two /2.d0/, three /3.d0/,
     *      fact /100.d0/
c        third should be slightly less than 1/3.
      data  third /0.33333d0/
c
c  notation and general remarks.
c
c     rho is the ratio of the data slope to the derivative being tested.
c     lambda is the ratio of d2 to d1.
c     that = t-hat(rho) is the normalized location of the extremum.
c     phi is the normalized value of p(x)-f1 at x = xhat = x-hat(rho),
c           where  that = (xhat - x1)/h .
c        that is, p(xhat)-f1 = d*h*phi,  where d=d1 or d2.
c     similarly,  p(xhat)-f2 = d*h*(phi-rho) .
c
c      small should be a few orders of magnitude greater than macheps.
c***first executable statement  dpchsw
      small = fact*d1mach(4)
c
c  do main calculation.
c
      if (d1 .eq. zero)  then
c
c        special case -- d1.eq.zero .
c
c          if d2 is also zero, this routine should not have been called.
         if (d2 .eq. zero)  go to 5001
c
         rho = slope/d2
c          extremum is outside interval when rho .ge. 1/3 .
         if (rho .ge. third)  go to 5000
         that = (two*(three*rho-one)) / (three*(two*rho-one))
         phi = that**2 * ((three*rho-one)/three)
c
c          convert to distance from f2 if iextrm.ne.1 .
         if (iextrm .ne. 1)  phi = phi - rho
c
c          test for exceeding limit, and adjust accordingly.
         hphi = h * abs(phi)
         if (hphi*abs(d2) .gt. dfmax)  then
c           at this point, hphi.gt.0, so divide is ok.
            d2 = sign (dfmax/hphi, d2)
         endif
      else
c
         rho = slope/d1
         lambda = -d2/d1
         if (d2 .eq. zero)  then
c
c           special case -- d2.eq.zero .
c
c             extremum is outside interval when rho .ge. 1/3 .
            if (rho .ge. third)  go to 5000
            cp = two - three*rho
            nu = one - two*rho
            that = one / (three*nu)
         else
            if (lambda .le. zero)  go to 5001
c
c           normal case -- d1 and d2 both nonzero, opposite signs.
c
            nu = one - lambda - two*rho
            sigma = one - rho
            cp = nu + sigma
            if (abs(nu) .gt. small)  then
               radcal = (nu - (two*rho+one))*nu + sigma**2
               if (radcal .lt. zero)  go to 5002
               that = (cp - sqrt(radcal)) / (three*nu)
            else
               that = one/(two*sigma)
            endif
         endif
         phi = that*((nu*that - cp)*that + one)
c
c          convert to distance from f2 if iextrm.ne.1 .
         if (iextrm .ne. 1)  phi = phi - rho
c
c          test for exceeding limit, and adjust accordingly.
         hphi = h * abs(phi)
         if (hphi*abs(d1) .gt. dfmax)  then
c           at this point, hphi.gt.0, so divide is ok.
            d1 = sign (dfmax/hphi, d1)
            d2 = -lambda*d1
         endif
      endif
c
c  normal return.
c
 5000 continue
      ierr = 0
      return
c
c  error returns.
c
 5001 continue
c     d1 and d2 both zero, or both nonzero and same sign.
      ierr = -1
      call xermsg ('slatec', 'dpchsw', 'd1 and/or d2 invalid', ierr, 1)
      return
c
 5002 continue
c     negative value of radical (should never occur).
      ierr = -2
      call xermsg ('slatec', 'dpchsw', 'negative radical', ierr, 1)
      return
c------------- last line of dpchsw follows -----------------------------
      end
