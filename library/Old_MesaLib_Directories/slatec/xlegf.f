*deck xlegf
      subroutine xlegf (dnu1, nudiff, mu1, mu2, theta, id, pqa, ipqa,
     1   ierror)
c***begin prologue  xlegf
c***purpose  compute normalized legendre polynomials and associated
c            legendre functions.
c***library   slatec
c***category  c3a2, c9
c***type      single precision (xlegf-s, dxlegf-d)
c***keywords  legendre functions
c***author  smith, john m., (nbs and george mason university)
c***description
c
c   xlegf: extended-range single-precision legendre functions
c
c   a feature of the xlegf subroutine for legendre functions is
c the use of extended-range arithmetic, a software extension of
c ordinary floating-point arithmetic that greatly increases the
c exponent range of the representable numbers. this avoids the
c need for scaling the solutions to lie within the exponent range
c of the most restrictive manufacturer's hardware. the increased
c exponent range is achieved by allocating an integer storage
c location together with each floating-point storage location.
c
c   the interpretation of the pair (x,i) where x is floating-point
c and i is integer is x*(ir**i) where ir is the internal radix of
c the computer arithmetic.
c
c   this subroutine computes one of the following vectors:
c
c 1. legendre function of the first kind of negative order, either
c    a. p(-mu1,nu,x), p(-mu1-1,nu,x), ..., p(-mu2,nu,x) or
c    b. p(-mu,nu1,x), p(-mu,nu1+1,x), ..., p(-mu,nu2,x)
c 2. legendre function of the second kind, either
c    a. q(mu1,nu,x), q(mu1+1,nu,x), ..., q(mu2,nu,x) or
c    b. q(mu,nu1,x), q(mu,nu1+1,x), ..., q(mu,nu2,x)
c 3. legendre function of the first kind of positive order, either
c    a. p(mu1,nu,x), p(mu1+1,nu,x), ..., p(mu2,nu,x) or
c    b. p(mu,nu1,x), p(mu,nu1+1,x), ..., p(mu,nu2,x)
c 4. normalized legendre polynomials, either
c    a. pn(mu1,nu,x), pn(mu1+1,nu,x), ..., pn(mu2,nu,x) or
c    b. pn(mu,nu1,x), pn(mu,nu1+1,x), ..., pn(mu,nu2,x)
c
c where x = cos(theta).
c
c   the input values to xlegf are dnu1, nudiff, mu1, mu2, theta,
c and id. these must satisfy
c
c    dnu1 is real and greater than or equal to -0.5;
c    nudiff is integer and non-negative;
c    mu1 is integer and non-negative;
c    mu2 is integer and greater than or equal to mu1;
c    theta is real and in the half-open interval (0,pi/2];
c    id is integer and equal to 1, 2, 3 or 4;
c
c and  additionally either nudiff = 0 or mu2 = mu1.
c
c   if id=1 and nudiff=0, a vector of type 1a above is computed
c with nu=dnu1.
c
c   if id=1 and mu1=mu2, a vector of type 1b above is computed
c with nu1=dnu1, nu2=dnu1+nudiff and mu=mu1.
c
c   if id=2 and nudiff=0, a vector of type 2a above is computed
c with nu=dnu1.
c
c   if id=2 and mu1=mu2, a vector of type 2b above is computed
c with nu1=dnu1, nu2=dnu1+nudiff and mu=mu1.
c
c   if id=3 and nudiff=0, a vector of type 3a above is computed
c with nu=dnu1.
c
c   if id=3 and mu1=mu2, a vector of type 3b above is computed
c with nu1=dnu1, nu2=dnu1+nudiff and mu=mu1.
c
c   if id=4 and nudiff=0, a vector of type 4a above is computed
c with nu=dnu1.
c
c   if id=4 and mu1=mu2, a vector of type 4b above is computed
c with nu1=dnu1, nu2=dnu1+nudiff and mu=mu1.
c
c   in each case the vector of computed legendre function values
c is returned in the extended-range vector (pqa(i),ipqa(i)). the
c length of this vector is either mu2-mu1+1 or nudiff+1.
c
c   where possible, xlegf returns ipqa(i) as zero. in this case the
c value of the legendre function is contained entirely in pqa(i),
c so it can be used in subsequent computations without further
c consideration of extended-range arithmetic. if ipqa(i) is nonzero,
c then the value of the legendre function is not representable in
c floating-point because of underflow or overflow. the program that
c calls xlegf must test ipqa(i) to ensure correct usage.
c
c   ierror is an error indicator. if no errors are detected, ierror=0
c when control returns to the calling routine. if an error is detected,
c ierror is returned as nonzero. the calling routine must check the
c value of ierror.
c
c   if ierror=110 or 111, invalid input was provided to xlegf.
c   if ierror=101,102,103, or 104, invalid input was provided to xset.
c   if ierror=105 or 106, an internal consistency error occurred in
c xset (probably due to a software malfunction in the library routine
c i1mach).
c   if ierror=107, an overflow or underflow of an extended-range number
c was detected in xadj.
c   if ierror=108, an overflow or underflow of an extended-range number
c was detected in xc210.
c
c***see also  xset
c***references  olver and smith, associated legendre functions on the
c                 cut, j comp phys, v 51, n 3, sept 1983, pp 502--518.
c               smith, olver and lozier, extended-range arithmetic and
c                 normalized legendre polynomials, acm trans on math
c                 softw, v 7, n 1, march 1981, pp 93--105.
c***routines called  xermsg, xpmu, xpmup, xpnrm, xpqnu, xqmu, xqnu,
c                    xred, xset
c***revision history  (yymmdd)
c   820728  date written
c   890126  revised to meet slatec cml recommendations.  (dwl and jms)
c   901019  revisions to prologue.  (dwl and wrb)
c   901106  changed all specific intrinsics to generic.  (wrb)
c           corrected order of sections in prologue and added type
c           section.  (wrb)
c           calls to xerror changed to calls to xermsg.  (wrb)
c   920127  revised purpose section of prologue.  (dwl)
c***end prologue  xlegf
      real pqa,dnu1,dnu2,sx,theta,x,pi2
      dimension pqa(*),ipqa(*)
c
c***first executable statement  xlegf
      ierror=0
      call xset (0, 0, 0.0, 0,ierror)
      if (ierror.ne.0) return
      pi2=2.*atan(1.)
c
c        zero output arrays
c
      l=(mu2-mu1)+nudiff+1
      do 290 i=1,l
      pqa(i)=0.
  290 ipqa(i)=0
c
c        check for valid input values
c
      if(nudiff.lt.0) go to 400
      if(dnu1.lt.-.5) go to 400
      if(mu2.lt.mu1) go to 400
      if(mu1.lt.0) go to 400
      if(theta.le.0..or.theta.gt.pi2) go to 420
      if(id.lt.1.or.id.gt.4) go to 400
      if((mu1.ne.mu2).and.(nudiff.gt.0)) go to 400
c
c        if dnu1 is not an integer, normalized p(mu,dnu,x)
c        cannot be calculated.  if dnu1 is an integer and
c        mu1.gt.dnu2 then all values of p(+mu,dnu,x) and
c        normalized p(mu,nu,x) will be zero.
c
      dnu2=dnu1+nudiff
      if((id.eq.3).and.(mod(dnu1,1.).ne.0.)) go to 295
      if((id.eq.4).and.(mod(dnu1,1.).ne.0.)) go to 400
      if((id.eq.3.or.id.eq.4).and.mu1.gt.dnu2) return
  295 continue
c
      x=cos(theta)
      sx=1./sin(theta)
      if(id.eq.2) go to 300
      if(mu2-mu1.le.0) go to 360
c
c        fixed nu, variable mu
c        call xpmu to calculate p(-mu1,nu,x),....,p(-mu2,nu,x)
c
      call xpmu(dnu1,dnu2,mu1,mu2,theta,x,sx,id,pqa,ipqa,ierror)
      if (ierror.ne.0) return
      go to 380
c
  300 if(mu2.eq.mu1) go to 320
c
c        fixed nu, variable mu
c        call xqmu to calculate q(mu1,nu,x),....,q(mu2,nu,x)
c
      call xqmu(dnu1,dnu2,mu1,mu2,theta,x,sx,id,pqa,ipqa,ierror)
      if (ierror.ne.0) return
      go to 390
c
c        fixed mu, variable nu
c        call xqnu to calculate q(mu,dnu1,x),....,q(mu,dnu2,x)
c
  320 call xqnu(dnu1,dnu2,mu1,theta,x,sx,id,pqa,ipqa,ierror)
      if (ierror.ne.0) return
      go to 390
c
c        fixed mu, variable nu
c        call xpqnu to calculate p(-mu,dnu1,x),....,p(-mu,dnu2,x)
c
  360 call xpqnu(dnu1,dnu2,mu1,theta,id,pqa,ipqa,ierror)
      if (ierror.ne.0) return
c
c        if id = 3, transform p(-mu,nu,x) vector into
c        p(mu,nu,x) vector.
c
  380 if(id.eq.3) call xpmup(dnu1,dnu2,mu1,mu2,pqa,ipqa,ierror)
      if (ierror.ne.0) return
c
c        if id = 4, transform p(-mu,nu,x) vector into
c        normalized p(mu,nu,x) vector.
c
      if(id.eq.4) call xpnrm(dnu1,dnu2,mu1,mu2,pqa,ipqa,ierror)
      if (ierror.ne.0) return
c
c        place results in reduced form if possible
c        and return to main program.
c
  390 do 395 i=1,l
      call xred(pqa(i),ipqa(i),ierror)
      if (ierror.ne.0) return
  395 continue
      return
c
c        *****     error termination     *****
c
  400 call xermsg ('slatec', 'xlegf',
     +             'dnu1, nudiff, mu1, mu2, or id not valid', 110, 1)
      ierror=110
      return
  420 call xermsg ('slatec', 'xlegf', 'theta out of range', 111, 1)
      ierror=111
      return
      end
