*deck qpdoc
      subroutine qpdoc
c***begin prologue  qpdoc
c***purpose  documentation for quadpack, a package of subprograms for
c            automatic evaluation of one-dimensional definite integrals.
c***library   slatec (quadpack)
c***category  h2, z
c***type      all (qpdoc-a)
c***keywords  documentation, guidelines for selection, quadpack,
c             quadrature, survey of integrators
c***author  piessens, robert
c             applied mathematics and programming division
c             k. u. leuven
c           de doncker, elise
c             applied mathematics and programming division
c             k. u. leuven
c           kahaner, d. k., (nbs)
c***description
c
c 1. introduction
c    ------------
c    quadpack is a fortran subroutine package for the numerical
c    computation of definite one-dimensional integrals. it originated
c    from a joint project of r. piessens and e. de doncker (appl.
c    math. and progr. div.- k.u.leuven, belgium), c. ueberhuber (inst.
c    fuer math.- techn. u. wien, austria), and d. kahaner (national
c    bureau of standards- washington d.c., u.s.a.).
c
c    documentation routine qpdoc describes the package in the form it
c    was released from a.m.p.d.- leuven, for adherence to the slatec
c    library in may 1981. apart from a survey of the integrators, some
c    guidelines will be given in order to help the quadpack user with
c    selecting an appropriate routine or a combination of several
c    routines for handling his problem.
c
c    in the long description of qpdoc it is demonstrated how to call
c    the integrators, by means of small example calling programs.
c
c    for precise guidelines involving the use of each routine in
c    particular, we refer to the extensive introductory comments
c    within each routine.
c
c 2. survey
c    ------
c    the following list gives an overview of the quadpack integrators.
c    the routine names for the double precision versions are preceded
c    by the letter d.
c
c    - qng  : is a simple non-adaptive automatic integrator, based on
c             a sequence of rules with increasing degree of algebraic
c             precision (patterson, 1968).
c
c    - qag  : is a simple globally adaptive integrator using the
c             strategy of aind (piessens, 1973). it is possible to
c             choose between 6 pairs of gauss-kronrod quadrature
c             formulae for the rule evaluation component. the pairs
c             of high degree of precision are suitable for handling
c             integration difficulties due to a strongly oscillating
c             integrand.
c
c    - qags : is an integrator based on globally adaptive interval
c             subdivision in connection with extrapolation (de doncker,
c             1978) by the epsilon algorithm (wynn, 1956).
c
c    - qagp : serves the same purposes as qags, but also allows
c             for eventual user-supplied information, i.e. the
c             abscissae of internal singularities, discontinuities
c             and other difficulties of the integrand function.
c             the algorithm is a modification of that in qags.
c
c    - qagi : handles integration over infinite intervals. the
c             infinite range is mapped onto a finite interval and
c             then the same strategy as in qags is applied.
c
c    - qawo : is a routine for the integration of cos(omega*x)*f(x)
c             or sin(omega*x)*f(x) over a finite interval (a,b).
c             omega is is specified by the user
c             the rule evaluation component is based on the
c             modified clenshaw-curtis technique.
c             an adaptive subdivision scheme is used connected with
c             an extrapolation procedure, which is a modification
c             of that in qags and provides the possibility to deal
c             even with singularities in f.
c
c    - qawf : calculates the fourier cosine or fourier sine
c             transform of f(x), for user-supplied interval (a,
c             infinity), omega, and f. the procedure of qawo is
c             used on successive finite intervals, and convergence
c             acceleration by means of the epsilon algorithm (wynn,
c             1956) is applied to the series of the integral
c             contributions.
c
c    - qaws : integrates w(x)*f(x) over (a,b) with a.lt.b finite,
c             and   w(x) = ((x-a)**alfa)*((b-x)**beta)*v(x)
c             where v(x) = 1 or log(x-a) or log(b-x)
c                            or log(x-a)*log(b-x)
c             and   alfa.gt.(-1), beta.gt.(-1).
c             the user specifies a, b, alfa, beta and the type of
c             the function v.
c             a globally adaptive subdivision strategy is applied,
c             with modified clenshaw-curtis integration on the
c             subintervals which contain a or b.
c
c    - qawc : computes the cauchy principal value of f(x)/(x-c)
c             over a finite interval (a,b) and for
c             user-determined c.
c             the strategy is globally adaptive, and modified
c             clenshaw-curtis integration is used on the subranges
c             which contain the point x = c.
c
c  each of the routines above also has a "more detailed" version
c    with a name ending in e, as qage.  these provide more
c    information and control than the easier versions.
c
c
c   the preceding routines are all automatic.  that is, the user
c      inputs his problem and an error tolerance.  the routine
c      attempts to perform the integration to within the requested
c      absolute or relative error.
c   there are, in addition, a number of non-automatic integrators.
c      these are most useful when the problem is such that the
c      user knows that a fixed rule will provide the accuracy
c      required.  typically they return an error estimate but make
c      no attempt to satisfy any particular input error request.
c
c      qk15
c      qk21
c      qk31
c      qk41
c      qk51
c      qk61
c           estimate the integral on [a,b] using 15, 21,..., 61
c           point rule and return an error estimate.
c      qk15i 15 point rule for (semi)infinite interval.
c      qk15w 15 point rule for special singular weight functions.
c      qc25c 25 point rule for cauchy principal values
c      qc25f 25 point rule for sin/cos integrand.
c      qmomo integrates k-th degree chebyshev polynomial times
c            function with various explicit singularities.
c
c 3. guidelines for the use of quadpack
c    ----------------------------------
c    here it is not our purpose to investigate the question when
c    automatic quadrature should be used. we shall rather attempt
c    to help the user who already made the decision to use quadpack,
c    with selecting an appropriate routine or a combination of
c    several routines for handling his problem.
c
c    for both quadrature over finite and over infinite intervals,
c    one of the first questions to be answered by the user is
c    related to the amount of computer time he wants to spend,
c    versus his -own- time which would be needed, for example, for
c    manual subdivision of the interval or other analytic
c    manipulations.
c
c    (1) the user may not care about computer time, or not be
c        willing to do any analysis of the problem. especially when
c        only one or a few integrals must be calculated, this attitude
c        can be perfectly reasonable. in this case it is clear that
c        either the most sophisticated of the routines for finite
c        intervals, qags, must be used, or its analogue for infinite
c        intervals, gagi. these routines are able to cope with
c        rather difficult, even with improper integrals.
c        this way of proceeding may be expensive. but the integrator
c        is supposed to give you an answer in return, with additional
c        information in the case of a failure, through its error
c        estimate and flag. yet it must be stressed that the programs
c        cannot be totally reliable.
c        ------
c
c    (2) the user may want to examine the integrand function.
c        if bad local difficulties occur, such as a discontinuity, a
c        singularity, derivative singularity or high peak at one or
c        more points within the interval, the first advice is to
c        split up the interval at these points. the integrand must
c        then be examined over each of the subintervals separately,
c        so that a suitable integrator can be selected for each of
c        them. if this yields problems involving relative accuracies
c        to be imposed on -finite- subintervals, one can make use of
c        qagp, which must be provided with the positions of the local
c        difficulties. however, if strong singularities are present
c        and a high accuracy is requested, application of qags on the
c        subintervals may yield a better result.
c
c        for quadrature over finite intervals we thus dispose of qags
c        and
c        - qng for well-behaved integrands,
c        - qag for functions with an oscillating behaviour of a non
c          specific type,
c        - qawo for functions, eventually singular, containing a
c          factor cos(omega*x) or sin(omega*x) where omega is known,
c        - qaws for integrands with algebraico-logarithmic end point
c          singularities of known type,
c        - qawc for cauchy principal values.
c
c        remark
c        ------
c        on return, the work arrays in the argument lists of the
c        adaptive integrators contain information about the interval
c        subdivision process and hence about the integrand behaviour:
c        the end points of the subintervals, the local integral
c        contributions and error estimates, and eventually other
c        characteristics. for this reason, and because of its simple
c        globally adaptive nature, the routine qag in particular is
c        well-suited for integrand examination. difficult spots can
c        be located by investigating the error estimates on the
c        subintervals.
c
c        for infinite intervals we provide only one general-purpose
c        routine, qagi. it is based on the qags algorithm applied
c        after a transformation of the original interval into (0,1).
c        yet it may eventuate that another type of transformation is
c        more appropriate, or one might prefer to break up the
c        original interval and use qagi only on the infinite part
c        and so on. these kinds of actions suggest a combined use of
c        different quadpack integrators. note that, when the only
c        difficulty is an integrand singularity at the finite
c        integration limit, it will in general not be necessary to
c        break up the interval, as qagi deals with several types of
c        singularity at the boundary point of the integration range.
c        it also handles slowly convergent improper integrals, on
c        the condition that the integrand does not oscillate over
c        the entire infinite interval. if it does we would advise
c        to sum succeeding positive and negative contributions to
c        the integral -e.g. integrate between the zeros- with one
c        or more of the finite-range integrators, and apply
c        convergence acceleration eventually by means of quadpack
c        subroutine qelg which implements the epsilon algorithm.
c        such quadrature problems include the fourier transform as
c        a special case. yet for the latter we have an automatic
c        integrator available, qawf.
c
c *long description:
c
c 4. example programs
c    ----------------
c 4.1. calling program for qng
c      -----------------------
c
c            real a,abserr,b,f,epsabs,epsrel,result
c            integer ier,neval
c            external f
c            a = 0.0e0
c            b = 1.0e0
c            epsabs = 0.0e0
c            epsrel = 1.0e-3
c            call qng(f,a,b,epsabs,epsrel,result,abserr,neval,ier)
c      c  include write statements
c            stop
c            end
c      c
c            real function f(x)
c            real x
c            f = exp(x)/(x*x+0.1e+01)
c            return
c            end
c
c 4.2. calling program for qag
c      -----------------------
c
c            real a,abserr,b,epsabs,epsrel,f,result,work
c            integer ier,iwork,key,last,lenw,limit,neval
c            dimension iwork(100),work(400)
c            external f
c            a = 0.0e0
c            b = 1.0e0
c            epsabs = 0.0e0
c            epsrel = 1.0e-3
c            key = 6
c            limit = 100
c            lenw = limit*4
c            call qag(f,a,b,epsabs,epsrel,key,result,abserr,neval,
c           *  ier,limit,lenw,last,iwork,work)
c      c  include write statements
c            stop
c            end
c      c
c            real function f(x)
c            real x
c            f = 2.0e0/(2.0e0+sin(31.41592653589793e0*x))
c            return
c            end
c
c 4.3. calling program for qags
c      ------------------------
c
c            real a,abserr,b,epsabs,epsrel,f,result,work
c            integer ier,iwork,last,lenw,limit,neval
c            dimension iwork(100),work(400)
c            external f
c            a = 0.0e0
c            b = 1.0e0
c            epsabs = 0.0e0
c            epsrel = 1.0e-3
c            limit = 100
c            lenw = limit*4
c            call qags(f,a,b,epsabs,epsrel,result,abserr,neval,ier,
c           *  limit,lenw,last,iwork,work)
c      c  include write statements
c            stop
c            end
c      c
c            real function f(x)
c            real x
c            f = 0.0e0
c            if(x.gt.0.0e0) f = 1.0e0/sqrt(x)
c            return
c            end
c
c 4.4. calling program for qagp
c      ------------------------
c
c            real a,abserr,b,epsabs,epsrel,f,points,result,work
c            integer ier,iwork,last,leniw,lenw,limit,neval,npts2
c            dimension iwork(204),points(4),work(404)
c            external f
c            a = 0.0e0
c            b = 1.0e0
c            npts2 = 4
c            points(1) = 1.0e0/7.0e0
c            points(2) = 2.0e0/3.0e0
c            limit = 100
c            leniw = limit*2+npts2
c            lenw = limit*4+npts2
c            call qagp(f,a,b,npts2,points,epsabs,epsrel,result,abserr,
c           *  neval,ier,leniw,lenw,last,iwork,work)
c      c  include write statements
c            stop
c            end
c      c
c            real function f(x)
c            real x
c            f = 0.0e+00
c            if(x.ne.1.0e0/7.0e0.and.x.ne.2.0e0/3.0e0) f =
c           *  abs(x-1.0e0/7.0e0)**(-0.25e0)*
c           *  abs(x-2.0e0/3.0e0)**(-0.55e0)
c            return
c            end
c
c 4.5. calling program for qagi
c      ------------------------
c
c            real abserr,boun,epsabs,epsrel,f,result,work
c            integer ier,inf,iwork,last,lenw,limit,neval
c            dimension iwork(100),work(400)
c            external f
c            boun = 0.0e0
c            inf = 1
c            epsabs = 0.0e0
c            epsrel = 1.0e-3
c            limit = 100
c            lenw = limit*4
c            call qagi(f,boun,inf,epsabs,epsrel,result,abserr,neval,
c           *  ier,limit,lenw,last,iwork,work)
c      c  include write statements
c            stop
c            end
c      c
c            real function f(x)
c            real x
c            f = 0.0e0
c            if(x.gt.0.0e0) f = sqrt(x)*log(x)/
c           *                   ((x+1.0e0)*(x+2.0e0))
c            return
c            end
c
c 4.6. calling program for qawo
c      ------------------------
c
c            real a,abserr,b,epsabs,epsrel,f,result,omega,work
c            integer ier,integr,iwork,last,leniw,lenw,limit,maxp1,neval
c            dimension iwork(200),work(925)
c            external f
c            a = 0.0e0
c            b = 1.0e0
c            omega = 10.0e0
c            integr = 1
c            epsabs = 0.0e0
c            epsrel = 1.0e-3
c            limit = 100
c            leniw = limit*2
c            maxp1 = 21
c            lenw = limit*4+maxp1*25
c            call qawo(f,a,b,omega,integr,epsabs,epsrel,result,abserr,
c           *  neval,ier,leniw,maxp1,lenw,last,iwork,work)
c      c  include write statements
c            stop
c            end
c      c
c            real function f(x)
c            real x
c            f = 0.0e0
c            if(x.gt.0.0e0) f = exp(-x)*log(x)
c            return
c            end
c
c 4.7. calling program for qawf
c      ------------------------
c
c            real a,abserr,epsabs,f,result,omega,work
c            integer ier,integr,iwork,last,leniw,lenw,limit,limlst,
c           *  lst,maxp1,neval
c            dimension iwork(250),work(1025)
c            external f
c            a = 0.0e0
c            omega = 8.0e0
c            integr = 2
c            epsabs = 1.0e-3
c            limlst = 50
c            limit = 100
c            leniw = limit*2+limlst
c            maxp1 = 21
c            lenw = leniw*2+maxp1*25
c            call qawf(f,a,omega,integr,epsabs,result,abserr,neval,
c           *  ier,limlst,lst,leniw,maxp1,lenw,iwork,work)
c      c  include write statements
c            stop
c            end
c      c
c            real function f(x)
c            real x
c            if(x.gt.0.0e0) f = sin(50.0e0*x)/(x*sqrt(x))
c            return
c            end
c
c 4.8. calling program for qaws
c      ------------------------
c
c            real a,abserr,alfa,b,beta,epsabs,epsrel,f,result,work
c            integer ier,integr,iwork,last,lenw,limit,neval
c            dimension iwork(100),work(400)
c            external f
c            a = 0.0e0
c            b = 1.0e0
c            alfa = -0.5e0
c            beta = -0.5e0
c            integr = 1
c            epsabs = 0.0e0
c            epsrel = 1.0e-3
c            limit = 100
c            lenw = limit*4
c            call qaws(f,a,b,alfa,beta,integr,epsabs,epsrel,result,
c           *  abserr,neval,ier,limit,lenw,last,iwork,work)
c      c  include write statements
c            stop
c            end
c      c
c            real function f(x)
c            real x
c            f = sin(10.0e0*x)
c            return
c            end
c
c 4.9. calling program for qawc
c      ------------------------
c
c            real a,abserr,b,c,epsabs,epsrel,f,result,work
c            integer ier,iwork,last,lenw,limit,neval
c            dimension iwork(100),work(400)
c            external f
c            a = -1.0e0
c            b = 1.0e0
c            c = 0.5e0
c            epsabs = 0.0e0
c            epsrel = 1.0e-3
c            limit = 100
c            lenw = limit*4
c            call qawc(f,a,b,c,epsabs,epsrel,result,abserr,neval,
c           *  ier,limit,lenw,last,iwork,work)
c      c  include write statements
c            stop
c            end
c      c
c            real function f(x)
c            real x
c            f = 1.0e0/(x*x+1.0e-4)
c            return
c            end
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   810401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900723  purpose section revised.  (wrb)
c***end prologue  qpdoc
c***first executable statement  qpdoc
      return
      end
