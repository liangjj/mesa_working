*deck cbeta
      complex function cbeta (a, b)
c***begin prologue  cbeta
c***purpose  compute the complete beta function.
c***library   slatec (fnlib)
c***category  c7b
c***type      complex (beta-s, dbeta-d, cbeta-c)
c***keywords  complete beta function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c cbeta computes the complete beta function of complex parameters a
c and b.
c input parameters:
c       a   complex and the real part of a positive
c       b   complex and the real part of b positive
c
c***references  (none)
c***routines called  cgamma, clbeta, gamlim, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890206  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   900727  added external statement.  (wrb)
c***end prologue  cbeta
      complex a, b, cgamma, clbeta
      external cgamma
      save xmax
      data xmax / 0.0 /
c***first executable statement  cbeta
      if (xmax.eq.0.0) then
         call gamlim (xmin, xmaxt)
         xmax = xmaxt
      endif
c
      if (real(a) .le. 0.0 .or. real(b) .le. 0.0) call xermsg ('slatec',
     +   'cbeta', 'real part of both arguments must be gt 0', 1, 2)
c
      if (real(a)+real(b).lt.xmax) cbeta = cgamma(a) * (cgamma(b)/
     1  cgamma(a+b) )
      if (real(a)+real(b).lt.xmax) return
c
      cbeta = exp (clbeta(a, b))
c
      return
      end
