*deck betai
      real function betai (x, pin, qin)
c***begin prologue  betai
c***purpose  calculate the incomplete beta function.
c***library   slatec (fnlib)
c***category  c7f
c***type      single precision (betai-s, dbetai-d)
c***keywords  fnlib, incomplete beta function, special functions
c***author  fullerton, w., (lanl)
c***description
c
c   betai calculates the real incomplete beta function.
c
c   the incomplete beta function ratio is the probability that a
c   random variable from a beta distribution having parameters pin and
c   qin will be less than or equal to x.
c
c     -- input arguments -- all arguments are real.
c   x      upper limit of integration.  x must be in (0,1) inclusive.
c   pin    first beta distribution parameter.  pin must be .gt. 0.0.
c   qin    second beta distribution parameter.  qin must be .gt. 0.0.
c
c***references  nancy e. bosten and e. l. battiste, remark on algorithm
c                 179, communications of the acm 17, 3 (march 1974),
c                 pp. 156.
c***routines called  albeta, r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920528  description and references sections revised.  (wrb)
c***end prologue  betai
      logical first
      save eps, alneps, sml, alnsml, first
      data first /.true./
c***first executable statement  betai
      if (first) then
         eps = r1mach(3)
         alneps = log(eps)
         sml = r1mach(1)
         alnsml = log(sml)
      endif
      first = .false.
c
      if (x .lt. 0. .or. x .gt. 1.0) call xermsg ('slatec', 'betai',
     +   'x is not in the range (0,1)', 1, 2)
      if (pin .le. 0. .or. qin .le. 0.) call xermsg ('slatec', 'betai',
     +   'p and/or q is le zero', 2, 2)
c
      y = x
      p = pin
      q = qin
      if (q.le.p .and. x.lt.0.8) go to 20
      if (x.lt.0.2) go to 20
      y = 1.0 - y
      p = qin
      q = pin
c
 20   if ((p+q)*y/(p+1.).lt.eps) go to 80
c
c evaluate the infinite sum first.
c term will equal y**p/beta(ps,p) * (1.-ps)i * y**i / fac(i)
c
      ps = q - aint(q)
      if (ps.eq.0.) ps = 1.0
      xb = p*log(y) -  albeta(ps, p) - log(p)
      betai = 0.0
      if (xb.lt.alnsml) go to 40
c
      betai = exp (xb)
      term = betai*p
      if (ps.eq.1.0) go to 40
c
      n = max (alneps/log(y), 4.0e0)
      do 30 i=1,n
        term = term*(i-ps)*y/i
        betai = betai + term/(p+i)
 30   continue
c
c now evaluate the finite sum, maybe.
c
 40   if (q.le.1.0) go to 70
c
      xb = p*log(y) + q*log(1.0-y) - albeta(p,q) - log(q)
      ib = max (xb/alnsml, 0.0e0)
      term = exp (xb - ib*alnsml)
      c = 1.0/(1.0-y)
      p1 = q*c/(p+q-1.)
c
      finsum = 0.0
      n = q
      if (q.eq.real(n)) n = n - 1
      do 50 i=1,n
        if (p1.le.1.0 .and. term/eps.le.finsum) go to 60
        term = (q-i+1)*c*term/(p+q-i)
c
        if (term.gt.1.0) ib = ib - 1
        if (term.gt.1.0) term = term*sml
c
        if (ib.eq.0) finsum = finsum + term
 50   continue
c
 60   betai = betai + finsum
 70   if (y.ne.x .or. p.ne.pin) betai = 1.0 - betai
      betai = max (min (betai, 1.0), 0.0)
      return
c
 80   betai = 0.0
      xb = p*log(max(y,sml)) - log(p) - albeta(p,q)
      if (xb.gt.alnsml .and. y.ne.0.) betai = exp (xb)
      if (y.ne.x .or. p.ne.pin) betai = 1.0 - betai
      return
c
      end
