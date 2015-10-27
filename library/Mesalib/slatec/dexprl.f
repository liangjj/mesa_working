*deck dexprl
      double precision function dexprl (x)
c***begin prologue  dexprl
c***purpose  calculate the relative error exponential (exp(x)-1)/x.
c***library   slatec (fnlib)
c***category  c4b
c***type      double precision (exprel-s, dexprl-d, cexprl-c)
c***keywords  elementary functions, exponential, first order, fnlib
c***author  fullerton, w., (lanl)
c***description
c
c evaluate  exprel(x) = (exp(x) - 1.0) / x.   for small abs(x) the
c taylor series is used.  if x is negative the reflection formula
c         exprel(x) = exp(x) * exprel(abs(x))
c may be used.  this reflection formula will be of use when the
c evaluation for small abs(x) is done by chebyshev series rather than
c taylor series.
c
c***references  (none)
c***routines called  d1mach
c***revision history  (yymmdd)
c   770801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  dexprl
      double precision x, absx, alneps, xbnd, xln, xn,  d1mach
      logical first
      save nterms, xbnd, first
      data first /.true./
c***first executable statement  dexprl
      if (first) then
         alneps = log(d1mach(3))
         xn = 3.72d0 - 0.3d0*alneps
         xln = log((xn+1.0d0)/1.36d0)
         nterms = xn - (xn*xln+alneps)/(xln+1.36d0) + 1.5d0
         xbnd = d1mach(3)
      endif
      first = .false.
c
      absx = abs(x)
      if (absx.gt.0.5d0) dexprl = (exp(x)-1.0d0)/x
      if (absx.gt.0.5d0) return
c
      dexprl = 1.0d0
      if (absx.lt.xbnd) return
c
      dexprl = 0.0d0
      do 20 i=1,nterms
        dexprl = 1.0d0 + dexprl*x/(nterms+2-i)
 20   continue
c
      return
      end
