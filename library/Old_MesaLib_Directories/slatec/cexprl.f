*deck cexprl
      complex function cexprl (z)
c***begin prologue  cexprl
c***purpose  calculate the relative error exponential (exp(x)-1)/x.
c***library   slatec (fnlib)
c***category  c4b
c***type      complex (exprel-s, dexprl-d, cexprl-c)
c***keywords  elementary functions, exponential, first order, fnlib
c***author  fullerton, w., (lanl)
c***description
c
c evaluate  (exp(z)-1)/z .  for small abs(z), we use the taylor
c series.  we could instead use the expression
c        cexprl(z) = (exp(x)*exp(i*y)-1)/z
c                  = (x*exprel(x) * (1 - 2*sin(y/2)**2) - 2*sin(y/2)**2
c                                    + i*sin(y)*(1+x*exprel(x))) / z
c
c***references  (none)
c***routines called  r1mach
c***revision history  (yymmdd)
c   770801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  cexprl
      complex z
      logical first
      save nterms, rbnd, first
      data first / .true. /
c***first executable statement  cexprl
      if (first) then
         alneps = log(r1mach(3))
         xn = 3.72 - 0.3*alneps
         xln = log((xn+1.0)/1.36)
         nterms = xn - (xn*xln+alneps)/(xln+1.36) + 1.5
         rbnd = r1mach(3)
      endif
      first = .false.
c
      r = abs(z)
      if (r.gt.0.5) cexprl = (exp(z) - 1.0) / z
      if (r.gt.0.5) return
c
      cexprl = (1.0, 0.0)
      if (r.lt.rbnd) return
c
      cexprl = (0.0, 0.0)
      do 20 i=1,nterms
        cexprl = 1.0 + cexprl*z/(nterms+2-i)
 20   continue
c
      return
      end
