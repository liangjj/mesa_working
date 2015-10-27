*deck qawce
      subroutine qawce (f, a, b, c, epsabs, epsrel, limit, result,
     +   abserr, neval, ier, alist, blist, rlist, elist, iord, last)
c***begin prologue  qawce
c***purpose  the routine calculates an approximation result to a
c            cauchy principal value i = integral of f*w over (a,b)
c            (w(x) = 1/(x-c), (c.ne.a, c.ne.b), hopefully satisfying
c            following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i))
c***library   slatec (quadpack)
c***category  h2a2a1, j4
c***type      single precision (qawce-s, dqawce-d)
c***keywords  automatic integrator, cauchy principal value,
c             clenshaw-curtis method, quadpack, quadrature,
c             special-purpose
c***author  piessens, robert
c             applied mathematics and programming division
c             k. u. leuven
c           de doncker, elise
c             applied mathematics and programming division
c             k. u. leuven
c***description
c
c        computation of a cauchy principal value
c        standard fortran subroutine
c        real version
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            b      - real
c                     upper limit of integration
c
c            c      - real
c                     parameter in the weight function, c.ne.a, c.ne.b
c                     if c = a or c = b, the routine will end with
c                     ier = 6.
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     in the partition of (a,b), limit.ge.1
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine
c                             the estimates for integral and error are
c                             less reliable. it is assumed that the
c                             requested accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more sub-
c                             divisions by increasing the value of
c                             limit. however, if this yields no
c                             improvement it is advised to analyze the
c                             the integrand, in order to determine the
c                             the integration difficulties. if the
c                             position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling
c                             appropriate integrators on the subranges.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour
c                             occurs at some interior points of
c                             the integration interval.
c                         = 6 the input is invalid, because
c                             c = a or c = b or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.1.
c                             result, abserr, neval, rlist(1), elist(1),
c                             iord(1) and last are set to zero. alist(1)
c                             and blist(1) are set to a and b
c                             respectively.
c
c            alist   - real
c                      vector of dimension at least limit, the first
c                       last  elements of which are the left
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            blist   - real
c                      vector of dimension at least limit, the first
c                       last  elements of which are the right
c                      end points of the subintervals in the partition
c                      of the given integration range (a,b)
c
c            rlist   - real
c                      vector of dimension at least limit, the first
c                       last  elements of which are the integral
c                      approximations on the subintervals
c
c            elist   - real
c                      vector of dimension limit, the first  last
c                      elements of which are the moduli of the absolute
c                      error estimates on the subintervals
c
c            iord    - integer
c                      vector of dimension at least limit, the first k
c                      elements of which are pointers to the error
c                      estimates over the subintervals, so that
c                      elist(iord(1)), ..., elist(iord(k)) with k = last
c                      if last.le.(limit/2+2), and k = limit+1-last
c                      otherwise, form a decreasing sequence
c
c            last    - integer
c                      number of subintervals actually produced in
c                      the subdivision process
c
c***references  (none)
c***routines called  qc25c, qpsrt, r1mach
c***revision history  (yymmdd)
c   800101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  qawce
c
      real a,aa,abserr,alist,area,area1,area12,area2,a1,a2,b,bb,blist,
     1  b1,b2,c,r1mach,elist,epmach,epsabs,epsrel,errbnd,errmax,error1,
     2  error2,errsum,f,result,rlist,uflow
      integer ier,iord,iroff1,iroff2,k,krule,last,limit,maxerr,nev,
     1  neval,nrmax
c
      dimension alist(*),blist(*),rlist(*),elist(*),
     1  iord(*)
c
      external f
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest
c                       error estimate
c           errmax    - elist(maxerr)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  qawce
      epmach = r1mach(4)
      uflow = r1mach(1)
c
c
c           test on validity of parameters
c           ------------------------------
c
      ier = 6
      neval = 0
      last = 0
      alist(1) = a
      blist(1) = b
      rlist(1) = 0.0e+00
      elist(1) = 0.0e+00
      iord(1) = 0
      result = 0.0e+00
      abserr = 0.0e+00
      if (c.eq.a.or.c.eq.b.or.(epsabs.le.0.0e+00.and.
     1    epsrel.lt.max(0.5e+02*epmach,0.5e-14))) go to 999
c
c           first approximation to the integral
c           -----------------------------------
c
      aa=a
      bb=b
      if (a.le.b) go to 10
      aa=b
      bb=a
10    ier=0
      krule = 1
      call qc25c(f,aa,bb,c,result,abserr,krule,neval)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      alist(1) = a
      blist(1) = b
c
c           test on accuracy
c
      errbnd = max(epsabs,epsrel*abs(result))
      if(limit.eq.1) ier = 1
      if(abserr.lt.min(0.1e-01*abs(result),errbnd)
     1  .or.ier.eq.1) go to 70
c
c           initialization
c           --------------
c
      alist(1) = aa
      blist(1) = bb
      rlist(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      nrmax = 1
      iroff1 = 0
      iroff2 = 0
c
c           main do-loop
c           ------------
c
      do 40 last = 2,limit
c
c           bisect the subinterval with nrmax-th largest
c           error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5e+00*(alist(maxerr)+blist(maxerr))
        b2 = blist(maxerr)
        if(c.le.b1.and.c.gt.a1) b1 = 0.5e+00*(c+b2)
        if(c.gt.b1.and.c.lt.b2) b1 = 0.5e+00*(a1+c)
        a2 = b1
        krule = 2
        call qc25c(f,a1,b1,c,area1,error1,krule,nev)
        neval = neval+nev
        call qc25c(f,a2,b2,c,area2,error2,krule,nev)
        neval = neval+nev
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(abs(rlist(maxerr)-area12).lt.0.1e-04*abs(area12)
     1    .and.erro12.ge.0.99e+00*errmax.and.krule.eq.0)
     2    iroff1 = iroff1+1
        if(last.gt.10.and.erro12.gt.errmax.and.krule.eq.0)
     1    iroff2 = iroff2+1
        rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = max(epsabs,epsrel*abs(area))
        if(errsum.le.errbnd) go to 15
c
c           test for roundoff error and eventually
c           set error flag.
c
        if(iroff1.ge.6.and.iroff2.gt.20) ier = 2
c
c           set error flag in the case that number of interval
c           bisections exceeds limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at a point of the integration range.
c
        if(max(abs(a1),abs(b2)).le.(0.1e+01+0.1e+03*epmach)
     1    *(abs(a2)+0.1e+04*uflow)) ier = 3
c
c           append the newly-created intervals to the list.
c
   15   if(error2.gt.error1) go to 20
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 30
   20   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine qpsrt to maintain the descending ordering
c           in the list of error estimates and select the
c           subinterval with nrmax-th largest error estimate (to be
c           bisected next).
c
   30    call qpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
c ***jump out of do-loop
        if(ier.ne.0.or.errsum.le.errbnd) go to 50
   40 continue
c
c           compute final result.
c           ---------------------
c
   50 result = 0.0e+00
      do 60 k=1,last
        result = result+rlist(k)
   60 continue
      abserr = errsum
   70 if (aa.eq.b) result=-result
  999 return
      end
