*deck dbeta
      double precision function dbeta (a, b)
c***begin prologue  dbeta
c***purpose  compute the complete beta function.
c***library   slatec (fnlib)
c***category  c7b
c***type      double precision (beta-s, dbeta-d, cbeta-c)
c***keywords  complete beta function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c dbeta(a,b) calculates the double precision complete beta function
c for double precision arguments a and b.
c
c***references  (none)
c***routines called  d1mach, dgamlm, dgamma, dlbeta, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900727  added external statement.  (wrb)
c***end prologue  dbeta
      double precision a, b, alnsml, xmax, xmin, dlbeta, dgamma, d1mach
      logical first
      external dgamma
      save xmax, alnsml, first
      data first /.true./
c***first executable statement  dbeta
      if (first) then
         call dgamlm (xmin, xmax)
         alnsml = log (d1mach(1))
      endif
      first = .false.
c
      if (a .le. 0.d0 .or. b .le. 0.d0) call xermsg ('slatec', 'dbeta',
     +   'both arguments must be gt 0', 2, 2)
c
      if (a+b.lt.xmax) dbeta = dgamma(a)*dgamma(b)/dgamma(a+b)
      if (a+b.lt.xmax) return
c
      dbeta = dlbeta (a, b)
      if (dbeta.lt.alnsml) go to 20
      dbeta = exp (dbeta)
      return
c
 20   dbeta = 0.d0
      call xermsg ('slatec', 'dbeta',
     +   'a and/or b so big beta underflows', 1, 1)
      return
c
      end
