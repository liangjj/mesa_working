*deck polfit
      subroutine polfit (n, x, y, w, maxdeg, ndeg, eps, r, ierr, a)
c***begin prologue  dpolft
c***purpose  fit discrete data in a least squares sense by polynomials
c            in one variable.
c***library   slatec
c***category  k1a1a2
c***type      double precision (polfit-s, dpolft-d)
c***keywords  curve fitting, data fitting, least squares, polynomial fit
c***author  shampine, l. f., (snla)
c           davenport, s. m., (snla)
c           huddleston, r. e., (snll)
c***description
c
c     abstract
c
c     given a collection of points x(i) and a set of values y(i) which
c     correspond to some function or measurement at each of the x(i),
c     subroutine  dpolft  computes the weighted least-squares polynomial
c     fits of all degrees up to some degree either specified by the user
c     or determined by the routine.  the fits thus obtained are in
c     orthogonal polynomial form.  subroutine  dp1vlu  may then be
c     called to evaluate the fitted polynomials and any of their
c     derivatives at any point.  the subroutine  dpcoef  may be used to
c     express the polynomial fits as powers of (x-c) for any specified
c     point c.
c
c     the parameters for  dpolft  are
c
c     input -- all type real variables are double precision
c         n -      the number of data points.  the arrays x, y and w
c                  must be dimensioned at least  n  (n .ge. 1).
c         x -      array of values of the independent variable.  these
c                  values may appear in any order and need not all be
c                  distinct.
c         y -      array of corresponding function values.
c         w -      array of positive values to be used as weights.  if
c                  w(1) is negative,  dpolft  will set all the weights
c                  to 1.0, which means unweighted least squares error
c                  will be minimized.  to minimize relative error, the
c                  user should set the weights to:  w(i) = 1.0/y(i)**2,
c                  i = 1,...,n .
c         maxdeg - maximum degree to be allowed for polynomial fit.
c                  maxdeg  may be any non-negative integer less than  n.
c                  note -- maxdeg  cannot be equal to  n-1  when a
c                  statistical test is to be used for degree selection,
c                  i.e., when input value of  eps  is negative.
c         eps -    specifies the criterion to be used in determining
c                  the degree of fit to be computed.
c                  (1)  if  eps  is input negative,  dpolft  chooses the
c                       degree based on a statistical f test of
c                       significance.  one of three possible
c                       significance levels will be used:  .01, .05 or
c                       .10.  if  eps=-1.0 , the routine will
c                       automatically select one of these levels based
c                       on the number of data points and the maximum
c                       degree to be considered.  if  eps  is input as
c                       -.01, -.05, or -.10, a significance level of
c                       .01, .05, or .10, respectively, will be used.
c                  (2)  if  eps  is set to 0.,  dpolft  computes the
c                       polynomials of degrees 0 through  maxdeg .
c                  (3)  if  eps  is input positive,  eps  is the rms
c                       error tolerance which must be satisfied by the
c                       fitted polynomial.  dpolft  will increase the
c                       degree of fit until this criterion is met or
c                       until the maximum degree is reached.
c
c     output -- all type real variables are double precision
c         ndeg -   degree of the highest degree fit computed.
c         eps -    rms error of the polynomial of degree  ndeg .
c         r -      vector of dimension at least ndeg containing values
c                  of the fit of degree  ndeg  at each of the  x(i) .
c                  except when the statistical test is used, these
c                  values are more accurate than results from subroutine
c                  dp1vlu  normally are.
c         ierr -   error flag with the following possible values.
c             1 -- indicates normal execution, i.e., either
c                  (1)  the input value of  eps  was negative, and the
c                       computed polynomial fit of degree  ndeg
c                       satisfies the specified f test, or
c                  (2)  the input value of  eps  was 0., and the fits of
c                       all degrees up to  maxdeg  are complete, or
c                  (3)  the input value of  eps  was positive, and the
c                       polynomial of degree  ndeg  satisfies the rms
c                       error requirement.
c             2 -- invalid input parameter.  at least one of the input
c                  parameters has an illegal value and must be corrected
c                  before  dpolft  can proceed.  valid input results
c                  when the following restrictions are observed
c                       n .ge. 1
c                       0 .le. maxdeg .le. n-1  for  eps .ge. 0.
c                       0 .le. maxdeg .le. n-2  for  eps .lt. 0.
c                       w(1)=-1.0  or  w(i) .gt. 0., i=1,...,n .
c             3 -- cannot satisfy the rms error requirement with a
c                  polynomial of degree no greater than  maxdeg .  best
c                  fit found is of degree  maxdeg .
c             4 -- cannot satisfy the test for significance using
c                  current value of  maxdeg .  statistically, the
c                  best fit found is of order  nord .  (in this case,
c                  ndeg will have one of the values:  maxdeg-2,
c                  maxdeg-1, or maxdeg).  using a higher value of
c                  maxdeg  may result in passing the test.
c         a -      work and output array having at least 3n+3maxdeg+3
c                  locations
c
c     note - dpolft  calculates all fits of degrees up to and including
c            ndeg .  any or all of these fits can be evaluated or
c            expressed as powers of (x-c) using  dp1vlu  and  dpcoef
c            after just one call to  dpolft .
c
c***references  l. f. shampine, s. m. davenport and r. e. huddleston,
c                 curve fitting by polynomials in one variable, report
c                 sla-74-0270, sandia laboratories, june 1974.
c***routines called  dp1vlu, xermsg
c***revision history  (yymmdd)
c   740601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900911  added variable yp to double precision declaration.  (wrb)
c   920501  reformatted the references section.  (wrb)
c   920527  corrected erroneous statements in description.  (wrb)
c***end prologue  dpolft
      integer i,idegf,ierr,j,jp1,jpas,k1,k1pj,k2,k2pj,k3,k3pi,k4,
     * k4pi,k5,k5pi,ksig,m,maxdeg,mop1,ndeg,nder,nfail
      real*8 temd1,temd2
      real*8 a(*),degf,den,eps,etst,f,fcrit,r(*),sig,sigj,
     * sigjm1,sigpas,temp,x(*),xm,y(*),yp,w(*),w1,w11
      real*8 co(4,3)
      save co
      data  co(1,1), co(2,1), co(3,1), co(4,1), co(1,2), co(2,2),
     1      co(3,2), co(4,2), co(1,3), co(2,3), co(3,3),
     2  co(4,3)/-13.086850d0,-2.4648165d0,-3.3846535d0,-1.2973162d0,
     3          -3.3381146d0,-1.7812271d0,-3.2578406d0,-1.6589279d0,
     4          -1.6282703d0,-1.3152745d0,-3.2640179d0,-1.9829776d0/
c***first executable statement  dpolft
      m = abs(n)
      if (m .eq. 0) go to 30
      if (maxdeg .lt. 0) go to 30
      a(1) = maxdeg
      mop1 = maxdeg + 1
      if (m .lt. mop1) go to 30
      if (eps .lt. 0.0d0 .and.  m .eq. mop1) go to 30
      xm = m
      etst = eps*eps*xm
      if (w(1) .lt. 0.0d0) go to 2
      do 1 i = 1,m
        if (w(i) .le. 0.0d0) go to 30
 1      continue
      go to 4
 2    do 3 i = 1,m
 3      w(i) = 1.0d0
 4    if (eps .ge. 0.0d0) go to 8
c
c determine significance level index to be used in statistical test for
c choosing degree of polynomial fit
c
      if (eps .gt. (-.55d0)) go to 5
      idegf = m - maxdeg - 1
      ksig = 1
      if (idegf .lt. 10) ksig = 2
      if (idegf .lt. 5) ksig = 3
      go to 8
 5    ksig = 1
      if (eps .lt. (-.03d0)) ksig = 2
      if (eps .lt. (-.07d0)) ksig = 3
c
c initialize indexes and coefficients for fitting
c
 8    k1 = maxdeg + 1
      k2 = k1 + maxdeg
      k3 = k2 + maxdeg + 2
      k4 = k3 + m
      k5 = k4 + m
      do 9 i = 2,k4
 9      a(i) = 0.0d0
      w11 = 0.0d0
      if (n .lt. 0) go to 11
c
c unconstrained case
c
      do 10 i = 1,m
        k4pi = k4 + i
        a(k4pi) = 1.0d0
 10     w11 = w11 + w(i)
      go to 13
c
c constrained case
c
 11   do 12 i = 1,m
        k4pi = k4 + i
 12     w11 = w11 + w(i)*a(k4pi)**2
c
c compute fit of degree zero
c
 13   temd1 = 0.0d0
      do 14 i = 1,m
        k4pi = k4 + i
        temd1 = temd1 + w(i)*y(i)*a(k4pi)
 14     continue
      temd1 = temd1/w11
      a(k2+1) = temd1
      sigj = 0.0d0
      do 15 i = 1,m
        k4pi = k4 + i
        k5pi = k5 + i
        temd2 = temd1*a(k4pi)
        r(i) = temd2
        a(k5pi) = temd2 - r(i)
 15     sigj = sigj + w(i)*((y(i)-r(i)) - a(k5pi))**2
      j = 0
c
c see if polynomial of degree 0 satisfies the degree selection criterion
c
      if (eps) 24,26,27
c
c increment degree
c
 16   j = j + 1
      jp1 = j + 1
      k1pj = k1 + j
      k2pj = k2 + j
      sigjm1 = sigj
c
c compute new b coefficient except when j = 1
c
      if (j .gt. 1) a(k1pj) = w11/w1
c
c compute new a coefficient
c
      temd1 = 0.0d0
      do 18 i = 1,m
        k4pi = k4 + i
        temd2 = a(k4pi)
        temd1 = temd1 + x(i)*w(i)*temd2*temd2
 18     continue
      a(jp1) = temd1/w11
c
c evaluate orthogonal polynomial at data points
c
      w1 = w11
      w11 = 0.0d0
      do 19 i = 1,m
        k3pi = k3 + i
        k4pi = k4 + i
        temp = a(k3pi)
        a(k3pi) = a(k4pi)
        a(k4pi) = (x(i)-a(jp1))*a(k3pi) - a(k1pj)*temp
 19     w11 = w11 + w(i)*a(k4pi)**2
c
c get new orthogonal polynomial coefficient using partial double
c precision
c
      temd1 = 0.0d0
      do 20 i = 1,m
        k4pi = k4 + i
        k5pi = k5 + i
        temd2 = w(i)*((y(i)-r(i))-a(k5pi))*a(k4pi)
 20     temd1 = temd1 + temd2
      temd1 = temd1/w11
      a(k2pj+1) = temd1
c
c update polynomial evaluations at each of the data points, and
c accumulate sum of squares of errors.  the polynomial evaluations are
c computed and stored in extended precision.  for the i-th data point,
c the most significant bits are stored in  r(i) , and the least
c significant bits are in  a(k5pi) .
c
      sigj = 0.0d0
      do 21 i = 1,m
        k4pi = k4 + i
        k5pi = k5 + i
        temd2 = r(i) + a(k5pi) + temd1*a(k4pi)
        r(i) = temd2
        a(k5pi) = temd2 - r(i)
 21     sigj = sigj + w(i)*((y(i)-r(i)) - a(k5pi))**2
c
c see if degree selection criterion has been satisfied or if degree
c maxdeg  has been reached
c
      if (eps) 23,26,27
c
c compute f statistics  (input eps .lt. 0.)
c
 23   if (sigj .eq. 0.0d0) go to 29
      degf = m - j - 1
      den = (co(4,ksig)*degf + 1.0d0)*degf
      fcrit = (((co(3,ksig)*degf) + co(2,ksig))*degf + co(1,ksig))/den
      fcrit = fcrit*fcrit
      f = (sigjm1 - sigj)*degf/sigj
      if (f .lt. fcrit) go to 25
c
c polynomial of degree j satisfies f test
c
 24   sigpas = sigj
      jpas = j
      nfail = 0
      if (maxdeg .eq. j) go to 32
      go to 16
c
c polynomial of degree j fails f test.  if there have been three
c successive failures, a statistically best degree has been found.
c
 25   nfail = nfail + 1
      if (nfail .ge. 3) go to 29
      if (maxdeg .eq. j) go to 32
      go to 16
c
c raise the degree if degree  maxdeg  has not yet been reached  (input
c eps = 0.)
c
 26   if (maxdeg .eq. j) go to 28
      go to 16
c
c see if rms error criterion is satisfied  (input eps .gt. 0.)
c
 27   if (sigj .le. etst) go to 28
      if (maxdeg .eq. j) go to 31
      go to 16
c
c returns
c
 28   ierr = 1
      ndeg = j
      sig = sigj
      go to 33
 29   ierr = 1
      ndeg = jpas
      sig = sigpas
      go to 33
 30   ierr = 2
      call lnkerr ('invalid input parameter in polfit')
      go to 37
 31   ierr = 3
      ndeg = maxdeg
      sig = sigj
      go to 33
 32   ierr = 4
      ndeg = jpas
      sig = sigpas
c
 33   a(k3) = ndeg
c
c when statistical test has been used, evaluate the best polynomial at
c all the data points if  r  does not already contain these values
c
      if(eps .ge. 0.0  .or.  ndeg .eq. maxdeg) go to 36
      nder = 0
      do 35 i = 1,m
        call polvlu (ndeg,nder,x(i),r(i),yp,a)
 35     continue
 36   eps = sqrt(sig/xm)
 37   return
      end
