*deck qawfe
      subroutine qawfe (f, a, omega, integr, epsabs, limlst, limit,
     +   maxp1, result, abserr, neval, ier, rslst, erlst, ierlst, lst,
     +   alist, blist, rlist, elist, iord, nnlog, chebmo)
c***begin prologue  qawfe
c***purpose  the routine calculates an approximation result to a
c            given fourier integral
c            i = integral of f(x)*w(x) over (a,infinity)
c             where w(x) = cos(omega*x) or w(x) = sin(omega*x),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.epsabs.
c***library   slatec (quadpack)
c***category  h2a3a1
c***type      single precision (qawfe-s, dqawfe-d)
c***keywords  automatic integrator, convergence acceleration,
c             fourier integrals, integration between zeros, quadpack,
c             quadrature, special-purpose integral
c***author  piessens, robert
c             applied mathematics and programming division
c             k. u. leuven
c           de doncker, elise
c             applied mathematics and programming division
c             k. u. leuven
c***description
c
c        computation of fourier integrals
c        standard fortran subroutine
c        real version
c
c        parameters
c         on entry
c            f      - real
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to
c                     be declared e x t e r n a l in the driver program.
c
c            a      - real
c                     lower limit of integration
c
c            omega  - real
c                     parameter in the weight function
c
c            integr - integer
c                     indicates which weight function is used
c                     integr = 1      w(x) = cos(omega*x)
c                     integr = 2      w(x) = sin(omega*x)
c                     if integr.ne.1.and.integr.ne.2, the routine will
c                     end with ier = 6.
c
c            epsabs - real
c                     absolute accuracy requested, epsabs.gt.0
c                     if epsabs.le.0, the routine will end with ier = 6.
c
c            limlst - integer
c                     limlst gives an upper bound on the number of
c                     cycles, limlst.ge.1.
c                     if limlst.lt.3, the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     allowed in the partition of each cycle, limit.ge.1
c                     each cycle, limit.ge.1.
c
c            maxp1  - integer
c                     gives an upper bound on the number of
c                     chebyshev moments which can be stored, i.e.
c                     for the intervals of lengths abs(b-a)*2**(-l),
c                     l=0,1, ..., maxp1-2, maxp1.ge.1
c
c         on return
c            result - real
c                     approximation to the integral x
c
c            abserr - real
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - ier = 0 normal and reliable termination of
c                             the routine. it is assumed that the
c                             requested accuracy has been achieved.
c                     ier.gt.0 abnormal termination of the routine. the
c                             estimates for integral and error are less
c                             reliable. it is assumed that the requested
c                             accuracy has not been achieved.
c            error messages
c                    if omega.ne.0
c                     ier = 1 maximum number of  cycles  allowed
c                             has been achieved., i.e. of subintervals
c                             (a+(k-1)c,a+kc) where
c                             c = (2*int(abs(omega))+1)*pi/abs(omega),
c                             for k = 1, 2, ..., lst.
c                             one can allow more cycles by increasing
c                             the value of limlst (and taking the
c                             according dimension adjustments into
c                             account).
c                             examine the array iwork which contains
c                             the error flags on the cycles, in order to
c                             look for eventual local integration
c                             difficulties. if the position of a local
c                             difficulty can be determined (e.g.
c                             singularity, discontinuity within the
c                             interval) one will probably gain from
c                             splitting up the interval at this point
c                             and calling appropriate integrators on
c                             the subranges.
c                         = 4 the extrapolation table constructed for
c                             convergence acceleration of the series
c                             formed by the integral contributions over
c                             the cycles, does not converge to within
c                             the requested accuracy. as in the case of
c                             ier = 1, it is advised to examine the
c                             array iwork which contains the error
c                             flags on the cycles.
c                         = 6 the input is invalid because
c                             (integr.ne.1 and integr.ne.2) or
c                              epsabs.le.0 or limlst.lt.3.
c                              result, abserr, neval, lst are set
c                              to zero.
c                         = 7 bad integrand behaviour occurs within one
c                             or more of the cycles. location and type
c                             of the difficulty involved can be
c                             determined from the vector ierlst. here
c                             lst is the number of cycles actually
c                             needed (see below).
c                             ierlst(k) = 1 the maximum number of
c                                           subdivisions (= limit) has
c                                           been achieved on the k th
c                                           cycle.
c                                       = 2 occurrence of roundoff error
c                                           is detected and prevents the
c                                           tolerance imposed on the
c                                           k th cycle, from being
c                                           achieved.
c                                       = 3 extremely bad integrand
c                                           behaviour occurs at some
c                                           points of the k th cycle.
c                                       = 4 the integration procedure
c                                           over the k th cycle does
c                                           not converge (to within the
c                                           required accuracy) due to
c                                           roundoff in the
c                                           extrapolation procedure
c                                           invoked on this cycle. it
c                                           is assumed that the result
c                                           on this interval is the
c                                           best which can be obtained.
c                                       = 5 the integral over the k th
c                                           cycle is probably divergent
c                                           or slowly convergent. it
c                                           must be noted that
c                                           divergence can occur with
c                                           any other value of
c                                           ierlst(k).
c                    if omega = 0 and integr = 1,
c                    the integral is calculated by means of dqagie
c                    and ier = ierlst(1) (with meaning as described
c                    for ierlst(k), k = 1).
c
c            rslst  - real
c                     vector of dimension at least limlst
c                     rslst(k) contains the integral contribution
c                     over the interval (a+(k-1)c,a+kc) where
c                     c = (2*int(abs(omega))+1)*pi/abs(omega),
c                     k = 1, 2, ..., lst.
c                     note that, if omega = 0, rslst(1) contains
c                     the value of the integral over (a,infinity).
c
c            erlst  - real
c                     vector of dimension at least limlst
c                     erlst(k) contains the error estimate corresponding
c                     with rslst(k).
c
c            ierlst - integer
c                     vector of dimension at least limlst
c                     ierlst(k) contains the error flag corresponding
c                     with rslst(k). for the meaning of the local error
c                     flags see description of output parameter ier.
c
c            lst    - integer
c                     number of subintervals needed for the integration
c                     if omega = 0 then lst is set to 1.
c
c            alist, blist, rlist, elist - real
c                     vector of dimension at least limit,
c
c            iord, nnlog - integer
c                     vector of dimension at least limit, providing
c                     space for the quantities needed in the subdivision
c                     process of each cycle
c
c            chebmo - real
c                     array of dimension at least (maxp1,25), providing
c                     space for the chebyshev moments needed within the
c                     cycles
c
c***references  (none)
c***routines called  qagie, qawoe, qelg, r1mach
c***revision history  (yymmdd)
c   800101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891009  removed unreferenced variable.  (wrb)
c   891009  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  qawfe
c
      real a,abseps,abserr,alist,blist,chebmo,correc,cycle,
     1  c1,c2,dl,drl,elist,ep,eps,epsa,epsabs,erlst,
     2  errsum,fact,omega,p,pi,p1,psum,reseps,result,res3la,rlist,rslst
     3  ,r1mach,uflow
      integer ier,ierlst,integr,iord,ktmin,l,lst,limit,ll,maxp1,
     1    nev,neval,nnlog,nres,numrl2
c
      dimension alist(*),blist(*),chebmo(maxp1,25),elist(*),
     1  erlst(*),ierlst(*),iord(*),nnlog(*),psum(52),
     2  res3la(3),rlist(*),rslst(*)
c
      external f
c
c
c            the dimension of  psum  is determined by the value of
c            limexp in subroutine qelg (psum must be
c            of dimension (limexp+2) at least).
c
c           list of major variables
c           -----------------------
c
c           c1, c2    - end points of subinterval (of length
c                       cycle)
c           cycle     - (2*int(abs(omega))+1)*pi/abs(omega)
c           psum      - vector of dimension at least (limexp+2)
c                       (see routine qelg)
c                       psum contains the part of the epsilon
c                       table which is still needed for further
c                       computations.
c                       each element of psum is a partial sum of
c                       the series which should sum to the value of
c                       the integral.
c           errsum    - sum of error estimates over the
c                       subintervals, calculated cumulatively
c           epsa      - absolute tolerance requested over current
c                       subinterval
c           chebmo    - array containing the modified chebyshev
c                       moments (see also routine qc25f)
c
      save p, pi
      data p/0.9e+00/,pi/0.31415926535897932e+01/
c
c           test on validity of parameters
c           ------------------------------
c
c***first executable statement  qawfe
      result = 0.0e+00
      abserr = 0.0e+00
      neval = 0
      lst = 0
      ier = 0
      if((integr.ne.1.and.integr.ne.2).or.epsabs.le.0.0e+00.or.
     1  limlst.lt.3) ier = 6
      if(ier.eq.6) go to 999
      if(omega.ne.0.0e+00) go to 10
c
c           integration by qagie if omega is zero
c           --------------------------------------
c
      if(integr.eq.1) call qagie(f,a,1,epsabs,0.0e+00,limit,
     1  result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      rslst(1) = result
      erlst(1) = abserr
      ierlst(1) = ier
      lst = 1
      go to 999
c
c           initializations
c           ---------------
c
   10 l = abs(omega)
      dl = 2*l+1
      cycle = dl*pi/abs(omega)
      ier = 0
      ktmin = 0
      neval = 0
      numrl2 = 0
      nres = 0
      c1 = a
      c2 = cycle+a
      p1 = 0.1e+01-p
      eps = epsabs
      uflow = r1mach(1)
      if(epsabs.gt.uflow/p1) eps = epsabs*p1
      ep = eps
      fact = 0.1e+01
      correc = 0.0e+00
      abserr = 0.0e+00
      errsum = 0.0e+00
c
c           main do-loop
c           ------------
c
      do 50 lst = 1,limlst
c
c           integrate over current subinterval.
c
        epsa = eps*fact
        call qawoe(f,c1,c2,omega,integr,epsa,0.0e+00,limit,lst,maxp1,
     1  rslst(lst),erlst(lst),nev,ierlst(lst),last,alist,blist,rlist,
     2  elist,iord,nnlog,momcom,chebmo)
        neval = neval+nev
        fact = fact*p
        errsum = errsum+erlst(lst)
        drl = 0.5e+02*abs(rslst(lst))
c
c           test on accuracy with partial sum
c
        if(errsum+drl.le.epsabs.and.lst.ge.6) go to 80
        correc = max(correc,erlst(lst))
        if(ierlst(lst).ne.0) eps = max(ep,correc*p1)
        if(ierlst(lst).ne.0) ier = 7
        if(ier.eq.7.and.(errsum+drl).le.correc*0.1e+02.and.
     1  lst.gt.5) go to 80
        numrl2 = numrl2+1
        if(lst.gt.1) go to 20
        psum(1) = rslst(1)
        go to 40
   20   psum(numrl2) = psum(ll)+rslst(lst)
        if(lst.eq.2) go to 40
c
c           test on maximum number of subintervals
c
        if(lst.eq.limlst) ier = 1
c
c           perform new extrapolation
c
        call qelg(numrl2,psum,reseps,abseps,res3la,nres)
c
c           test whether extrapolated result is influenced by
c           roundoff
c
        ktmin = ktmin+1
        if(ktmin.ge.15.and.abserr.le.0.1e-02*(errsum+drl)) ier = 4
        if(abseps.gt.abserr.and.lst.ne.3) go to 30
        abserr = abseps
        result = reseps
        ktmin = 0
c
c           if ier is not 0, check whether direct result (partial
c           sum) or extrapolated result yields the best integral
c           approximation
c
        if((abserr+0.1e+02*correc).le.epsabs.or.
     1  (abserr.le.epsabs.and.0.1e+02*correc.ge.epsabs)) go to 60
   30   if(ier.ne.0.and.ier.ne.7) go to 60
   40   ll = numrl2
        c1 = c2
        c2 = c2+cycle
   50 continue
c
c         set final result and error estimate
c         -----------------------------------
c
   60 abserr = abserr+0.1e+02*correc
      if(ier.eq.0) go to 999
      if(result.ne.0.0e+00.and.psum(numrl2).ne.0.0e+00) go to 70
      if(abserr.gt.errsum) go to 80
      if(psum(numrl2).eq.0.0e+00) go to 999
   70 if(abserr/abs(result).gt.(errsum+drl)/abs(psum(numrl2)))
     1  go to 80
      if(ier.ge.1.and.ier.ne.7) abserr = abserr+drl
      go to 999
   80 result = psum(numrl2)
      abserr = errsum+drl
  999 return
      end
