*deck exprel
      function exprel (x)
c***begin prologue  exprel
c***purpose  calculate the relative error exponential (exp(x)-1)/x.
c***library   slatec (fnlib)
c***category  c4b
c***type      single precision (exprel-s, dexprl-d, cexprl-c)
c***keywords  elementary functions, exponential, first order, fnlib
c***author  fullerton, w., (lanl)
c***description
c
c evaluate  exprel(x) = (exp(x) - 1.0) / x.   for small abs(x) the
c taylor series is used.  if x is negative, the reflection formula
c         exprel(x) = exp(x) * exprel(abs(x))
c may be used.  this reflection formula will be of use when the
c evaluation for small abs(x) is done by chebyshev series rather than
c taylor series.  exprel and x are single precision.
c
c***references  (none)
c***routines called  r1mach
c***revision history  (yymmdd)
c   770801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  exprel
      logical first
      save nterms, xbnd, first
      data first /.true./
c***first executable statement  exprel
      if (first) then
         alneps = log(r1mach(3))
         xn = 3.72 - 0.3*alneps
         xln = log((xn+1.0)/1.36)
         nterms = xn - (xn*xln+alneps)/(xln+1.36) + 1.5
         xbnd = r1mach(3)
      endif
      first = .false.
c
      absx = abs(x)
      if (absx.gt.0.5) exprel = (exp(x) - 1.0) / x
      if (absx.gt.0.5) return
c
      exprel = 1.0
      if (absx.lt.xbnd) return
c
      exprel = 0.0
      do 20 i=1,nterms
        exprel = 1.0 + exprel*x/(nterms+2-i)
 20   continue
c
      return
      end
