*deck xnrmp
      subroutine xnrmp (nu, mu1, mu2, sarg, mode, spn, ipn, isig,
     1   ierror)
c***begin prologue  xnrmp
c***purpose  compute normalized legendre polynomials.
c***library   slatec
c***category  c3a2, c9
c***type      single precision (xnrmp-s, dxnrmp-d)
c***keywords  legendre functions
c***author  lozier, daniel w., (national bureau of standards)
c           smith, john m., (nbs and george mason university)
c***description
c
c        subroutine to calculate normalized legendre polynomials
c        (dxnrmp is double-precision version)
c        xnrmp calculates normalized legendre polynomials of varying
c        order and fixed argument and degree. the order mu and degree
c        nu are non-negative integers and the argument is real. because
c        the algorithm requires the use of numbers outside the normal
c        machine range, this subroutine employs a special arithmetic
c        called extended-range arithmetic. see j.m. smith, f.w.j. olver,
c        and d.w. lozier, extended-range arithmetic and normalized
c        legendre polynomials, acm transactions on mathematical soft-
c        ware, 93-105, march 1981, for a complete description of the
c        algorithm and special arithmetic. also see program comments
c        in xset.
c
c        the normalized legendre polynomials are multiples of the
c        associated legendre polynomials of the first kind where the
c        normalizing coefficients are chosen so as to make the integral
c        from -1 to 1 of the square of each function equal to 1. see
c        e. jahnke, f. emde and f. losch, tables of higher functions,
c        mcgraw-hill, new york, 1960, p. 121.
c
c        the input values to xnrmp are nu, mu1, mu2, sarg, and mode.
c        these must satisfy
c          1. nu .ge. 0 specifies the degree of the normalized legendre
c             polynomial that is wanted.
c          2. mu1 .ge. 0 specifies the lowest-order normalized legendre
c             polynomial that is wanted.
c          3. mu2 .ge. mu1 specifies the highest-order normalized leg-
c             endre polynomial that is wanted.
c         4a. mode = 1 and -1.0 .le. sarg .le. 1.0 specifies that
c             normalized legendre(nu, mu, sarg) is wanted for mu = mu1,
c             mu1 + 1, ..., mu2.
c         4b. mode = 2 and -3.14159... .lt. sarg .lt. 3.14159... spec-
c             ifies that normalized legendre(nu, mu, cos(sarg)) is want-
c             ed for mu = mu1, mu1 + 1, ..., mu2.
c
c        the output of xnrmp consists of the two vectors spn and ipn
c        and the error estimate isig. the computed values are stored as
c        extended-range numbers such that
c             (spn(1),ipn(1))=normalized legendre(nu,mu1,x)
c             (spn(2),ipn(2))=normalized legendre(nu,mu1+1,x)
c                .
c                .
c             (spn(k),ipn(k))=normalized legendre(nu,mu2,x)
c        where k = mu2 - mu1 + 1 and x = sarg or cos(sarg) according
c        to whether mode = 1 or 2. finally, isig is an estimate of the
c        number of decimal digits lost through rounding errors in the
c        computation. for example if sarg is accurate to 12 significant
c        decimals, then the computed function values are accurate to
c        12 - isig significant decimals (except in neighborhoods of
c        zeros).
c
c        the interpretation of (spn(i),ipn(i)) is spn(i)*(ir**ipn(i))
c        where ir is the internal radix of the computer arithmetic. when
c        ipn(i) = 0 the value of the normalized legendre polynomial is
c        contained entirely in spn(i) and subsequent single-precision
c        computations can be performed without further consideration of
c        extended-range arithmetic. however, if ipn(i) .ne. 0 the corre-
c        sponding value of the normalized legendre polynomial cannot be
c        represented in single-precision because of overflow or under-
c        flow. the user must test ipn(i) in his/her program. in the case
c        that ipn(i) is nonzero, the user should try using double pre-
c        cision if it has a wider exponent range. if double precision
c        fails, the user could rewrite his/her program to use extended-
c        range arithmetic.
c
c        the interpretation of (spn(i),ipn(i)) can be changed to
c        spn(i)*(10**ipn(i)) by calling the extended-range subroutine
c        xcon. this should be done before printing the computed values.
c        as an example of usage, the fortran coding
c              j = k
c              do 20 i = 1, k
c              call xcon(spn(i), ipn(i),ierror)
c              if (ierror.ne.0) return
c              print 10, spn(i), ipn(i)
c           10 format(1x, e30.8 , i15)
c              if ((ipn(i) .eq. 0) .or. (j .lt. k)) go to 20
c              j = i - 1
c           20 continue
c        will print all computed values and determine the largest j
c        such that ipn(1) = ipn(2) = ... = ipn(j) = 0. because of the
c        change of representation caused by calling xcon, (spn(i),
c        ipn(i)) for i = j+1, j+2, ... cannot be used in subsequent
c        extended-range computations.
c
c        ierror is an error indicator. if no errors are detected,
c        ierror=0 when control returns to the calling routine. if
c        an error is detected, ierror is returned as nonzero. the
c        calling routine must check the value of ierror.
c
c        if ierror=112 or 113, invalid input was provided to xnrmp.
c        if ierror=101,102,103, or 104, invalid input was provided
c        to xset.
c        if ierror=105 or 106, an internal consistency error occurred
c        in xset (probably due to a software malfunction in the
c        library routine i1mach).
c        if ierror=107, an overflow or underflow of an extended-range
c        number was detected in xadj.
c        if ierror=108, an overflow or underflow of an extended-range
c        number was detected in xc210.
c
c***see also  xset
c***references  smith, olver and lozier, extended-range arithmetic and
c                 normalized legendre polynomials, acm trans on math
c                 softw, v 7, n 1, march 1981, pp 93--105.
c***routines called  xadd, xadj, xermsg, xred, xset
c***revision history  (yymmdd)
c   820712  date written
c   881020  revised to meet slatec cml recommendations.  (dwl and jms)
c   901019  revisions to prologue.  (dwl and wrb)
c   901106  changed all specific intrinsics to generic.  (wrb)
c           corrected order of sections in prologue and added type
c           section.  (wrb)
c           calls to xerror changed to calls to xermsg.  (wrb)
c   920127  revised purpose section of prologue.  (dwl)
c***end prologue  xnrmp
      integer nu, mu1, mu2, mode, ipn, isig
      real sarg, spn
      dimension spn(*), ipn(*)
      real c1,c2,p,p1,p2,p3,s,sx,t,tx,x,rk
c call xset to initialize extended-range arithmetic (see xset
c listing for details)
c***first executable statement  xnrmp
      ierror=0
      call xset (0, 0, 0.0, 0,ierror)
      if (ierror.ne.0) return
c
c        test for proper input values.
c
      if (nu.lt.0) go to 110
      if (mu1.lt.0) go to 110
      if (mu1.gt.mu2) go to 110
      if (nu.eq.0) go to 90
      if (mode.lt.1 .or. mode.gt.2) go to 110
      go to (10, 20), mode
   10 if (abs(sarg).gt.1.0) go to 120
      if (abs(sarg).eq.1.0) go to 90
      x = sarg
      sx = sqrt((1.0+abs(x))*((0.5-abs(x))+0.5))
      tx = x/sx
      isig = log10(2.0*nu*(5.0+tx**2))
      go to 30
   20 if (abs(sarg).gt.4.0*atan(1.0)) go to 120
      if (sarg.eq.0.0) go to 90
      x = cos(sarg)
      sx = abs(sin(sarg))
      tx = x/sx
      isig = log10(2.0*nu*(5.0+abs(sarg*tx)))
c
c        begin calculation
c
   30 mu = mu2
      i = mu2 - mu1 + 1
c
c        if mu.gt.nu, normalized legendre(nu,mu,x)=0.
c
   40 if (mu.le.nu) go to 50
      spn(i) = 0.0
      ipn(i) = 0
      i = i - 1
      mu = mu - 1
      if (i .gt. 0) go to 40
      isig = 0
      go to 160
   50 mu = nu
c
c        p1 = 0. = normalized legendre(nu,nu+1,x)
c
      p1 = 0.0
      ip1 = 0
c
c        calculate p2 = normalized legendre(nu,nu,x)
c
      p2 = 1.0
      ip2 = 0
      p3 = 0.5
      rk = 2.0
      do 60 j=1,nu
        p3 = ((rk+1.0)/rk)*p3
        p2 = p2*sx
        call xadj(p2, ip2,ierror)
        if (ierror.ne.0) return
        rk = rk + 2.0
   60 continue
      p2 = p2*sqrt(p3)
      call xadj(p2, ip2,ierror)
      if (ierror.ne.0) return
      s = 2.0*tx
      t = 1.0/nu
      if (mu2.lt.nu) go to 70
      spn(i) = p2
      ipn(i) = ip2
      i = i - 1
      if (i .eq. 0) go to 140
c
c        recurrence process
c
   70 p = mu*t
      c1 = 1.0/sqrt((1.0-p+t)*(1.0+p))
      c2 = s*p*c1*p2
      c1 = -sqrt((1.0+p+t)*(1.0-p))*c1*p1
      call xadd(c2, ip2, c1, ip1, p, ip,ierror)
      if (ierror.ne.0) return
      mu = mu - 1
      if (mu.gt.mu2) go to 80
c
c        store in array spn for return to calling routine.
c
      spn(i) = p
      ipn(i) = ip
      i = i - 1
      if (i .eq. 0) go to 140
   80 p1 = p2
      ip1 = ip2
      p2 = p
      ip2 = ip
      if (mu.le.mu1) go to 140
      go to 70
c
c        special case when x=-1 or +1, or nu=0.
c
   90 k = mu2 - mu1 + 1
      do 100 i=1,k
        spn(i) = 0.0
        ipn(i) = 0
  100 continue
      isig = 0
      if (mu1.gt.0) go to 160
      isig = 1
      spn(1) = sqrt(nu+0.5)
      ipn(1) = 0
      if (mod(nu,2).eq.0) go to 160
      if (mode.eq.1 .and. sarg.eq.1.0) go to 160
      if (mode.eq.2) go to 160
      spn(1) = -spn(1)
      go to 160
c
c          error printouts and termination.
c
  110 call xermsg ('slatec', 'xnrmp', 'nu, mu1, mu2 or mode not valid',
     +             112, 1)
      ierror=112
      return
  120 call xermsg ('slatec', 'xnrmp', 'sarg out of range', 113, 1)
      ierror=113
      return
c
c        return to calling program
c
  140 k = mu2 - mu1 + 1
      do 150 i=1,k
        call xred(spn(i),ipn(i),ierror)
        if (ierror.ne.0) return
  150 continue
  160 return
      end
