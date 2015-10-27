*deck dbetai
      double precision function dbetai (x, pin, qin)
c***begin prologue  dbetai
c***purpose  calculate the incomplete beta function.
c***library   slatec (fnlib)
c***category  c7f
c***type      double precision (betai-s, dbetai-d)
c***keywords  fnlib, incomplete beta function, special functions
c***author  fullerton, w., (lanl)
c***description
c
c   dbetai calculates the double precision incomplete beta function.
c
c   the incomplete beta function ratio is the probability that a
c   random variable from a beta distribution having parameters pin and
c   qin will be less than or equal to x.
c
c     -- input arguments -- all arguments are double precision.
c   x      upper limit of integration.  x must be in (0,1) inclusive.
c   pin    first beta distribution parameter.  pin must be .gt. 0.0.
c   qin    second beta distribution parameter.  qin must be .gt. 0.0.
c
c***references  nancy e. bosten and e. l. battiste, remark on algorithm
c                 179, communications of the acm 17, 3 (march 1974),
c                 pp. 156.
c***routines called  d1mach, dlbeta, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920528  description and references sections revised.  (wrb)
c***end prologue  dbetai
      double precision x, pin, qin, alneps, alnsml, c, eps, finsum, p,
     1  ps, q, sml, term, xb, xi, y, d1mach, dlbeta, p1
      logical first
      save eps, alneps, sml, alnsml, first
      data first /.true./
c***first executable statement  dbetai
      if (first) then
         eps = d1mach(3)
         alneps = log (eps)
         sml = d1mach(1)
         alnsml = log (sml)
      endif
      first = .false.
c
      if (x .lt. 0.d0 .or. x .gt. 1.d0) call xermsg ('slatec', 'dbetai',
     +   'x is not in the range (0,1)', 1, 2)
      if (pin .le. 0.d0 .or. qin .le. 0.d0) call xermsg ('slatec',
     +   'dbetai', 'p and/or q is le zero', 2, 2)
c
      y = x
      p = pin
      q = qin
      if (q.le.p .and. x.lt.0.8d0) go to 20
      if (x.lt.0.2d0) go to 20
      y = 1.0d0 - y
      p = qin
      q = pin
c
 20   if ((p+q)*y/(p+1.d0).lt.eps) go to 80
c
c evaluate the infinite sum first.  term will equal
c y**p/beta(ps,p) * (1.-ps)-sub-i * y**i / fac(i) .
c
      ps = q - aint(q)
      if (ps.eq.0.d0) ps = 1.0d0
      xb = p*log(y) - dlbeta(ps,p) - log(p)
      dbetai = 0.0d0
      if (xb.lt.alnsml) go to 40
c
      dbetai = exp (xb)
      term = dbetai*p
      if (ps.eq.1.0d0) go to 40
      n = max (alneps/log(y), 4.0d0)
      do 30 i=1,n
        xi = i
        term = term * (xi-ps)*y/xi
        dbetai = dbetai + term/(p+xi)
 30   continue
c
c now evaluate the finite sum, maybe.
c
 40   if (q.le.1.0d0) go to 70
c
      xb = p*log(y) + q*log(1.0d0-y) - dlbeta(p,q) - log(q)
      ib = max (xb/alnsml, 0.0d0)
      term = exp(xb - ib*alnsml)
      c = 1.0d0/(1.d0-y)
      p1 = q*c/(p+q-1.d0)
c
      finsum = 0.0d0
      n = q
      if (q.eq.dble(n)) n = n - 1
      do 50 i=1,n
        if (p1.le.1.0d0 .and. term/eps.le.finsum) go to 60
        xi = i
        term = (q-xi+1.0d0)*c*term/(p+q-xi)
c
        if (term.gt.1.0d0) ib = ib - 1
        if (term.gt.1.0d0) term = term*sml
c
        if (ib.eq.0) finsum = finsum + term
 50   continue
c
 60   dbetai = dbetai + finsum
 70   if (y.ne.x .or. p.ne.pin) dbetai = 1.0d0 - dbetai
      dbetai = max (min (dbetai, 1.0d0), 0.0d0)
      return
c
 80   dbetai = 0.0d0
      xb = p*log(max(y,sml)) - log(p) - dlbeta(p,q)
      if (xb.gt.alnsml .and. y.ne.0.0d0) dbetai = exp(xb)
      if (y.ne.x .or. p.ne.pin) dbetai = 1.0d0 - dbetai
c
      return
      end
