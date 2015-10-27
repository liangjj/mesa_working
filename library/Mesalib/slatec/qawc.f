*deck qawc
      subroutine qawc (f, a, b, c, epsabs, epsrel, result, abserr,
     +   neval, ier, limit, lenw, last, iwork, work)
c***begin prologue  qawc
c***purpose  the routine calculates an approximation result to a
c            cauchy principal value i = integral of f*w over (a,b)
c            (w(x) = 1/((x-c), c.ne.a, c.ne.b), hopefully satisfying
c            following claim for accuracy
c            abs(i-result).le.max(epsabe,epsrel*abs(i)).
c***library   slatec (quadpack)
c***category  h2a2a1, j4
c***type      single precision (qawc-s, dqawc-d)
c***keywords  automatic integrator, cauchy principal value,
c             clenshaw-curtis method, globally adaptive, quadpack,
c             quadrature, special-purpose
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
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            a      - real
c                     under limit of integration
c
c            b      - real
c                     upper limit of integration
c
c            c      - parameter in the weight function, c.ne.a, c.ne.b.
c                     if c = a or c = b, the routine will end with
c                     ier = 6 .
c
c            epsabs - real
c                     absolute accuracy requested
c            epsrel - real
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c         on return
c            result - real
c                     approximation to the integral
c
c            abserr - real
c                     estimate or the modulus of the absolute error,
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
c                             divisions by increasing the value of limit
c                             (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties.
c                             if the position of a local difficulty
c                             can be determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling
c                             appropriate integrators on the subranges.
c                         = 2 the occurrence of roundoff error is detec-
c                             ted, which prevents the requested
c                             tolerance from being achieved.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 6 the input is invalid, because
c                             c = a or c = b or
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                             or limit.lt.1 or lenw.lt.limit*4.
c                             result, abserr, neval, last are set to
c                             zero.  except when lenw or limit is
c                             invalid, iwork(1), work(limit*2+1) and
c                             work(limit*3+1) are set to zero, work(1)
c                             is set to a and work(limit+1) to b.
c
c         dimensioning parameters
c            limit - integer
c                    dimensioning parameter for iwork
c                    limit determines the maximum number of subintervals
c                    in the partition of the given integration interval
c                    (a,b), limit.ge.1.
c                    if limit.lt.1, the routine will end with ier = 6.
c
c           lenw   - integer
c                    dimensioning parameter for work
c                    lenw must be at least limit*4.
c                    if lenw.lt.limit*4, the routine will end with
c                    ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdivision process, which
c                    determines the number of significant elements
c                    actually in the work arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least limit, the first k
c                    elements of which contain pointers
c                    to the error estimates over the subintervals,
c                    such that work(limit*3+iwork(1)), ... ,
c                    work(limit*3+iwork(k)) form a decreasing
c                    sequence, with k = last if last.le.(limit/2+2),
c                    and k = limit+1-last otherwise
c
c            work  - real
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left
c                     end points of the subintervals in the
c                     partition of (a,b),
c                    work(limit+1), ..., work(limit+last) contain
c                     the right end points,
c                    work(limit*2+1), ..., work(limit*2+last) contain
c                     the integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3+last)
c                     contain the error estimates.
c
c***references  (none)
c***routines called  qawce, xermsg
c***revision history  (yymmdd)
c   800101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  qawc
c
      real a,abserr,b,c,epsabs,epsrel,f,result,work
      integer ier,iwork,lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(*),work(*)
c
      external f
c
c         check validity of limit and lenw.
c
c***first executable statement  qawc
      ier = 6
      neval = 0
      last = 0
      result = 0.0e+00
      abserr = 0.0e+00
      if(limit.lt.1.or.lenw.lt.limit*4) go to 10
c
c         prepare call for qawce.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
      call qawce(f,a,b,c,epsabs,epsrel,limit,result,abserr,neval,ier,
     1  work(1),work(l1),work(l2),work(l3),iwork,last)
c
c         call error handler if necessary.
c
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if (ier .ne. 0) call xermsg ('slatec', 'qawc',
     +   'abnormal return', ier, lvl)
      return
      end
