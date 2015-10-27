*deck inits
      function inits (os, nos, eta)
c***begin prologue  inits
c***purpose  determine the number of terms needed in an orthogonal
c            polynomial series so that it meets a specified accuracy.
c***library   slatec (fnlib)
c***category  c3a2
c***type      single precision (inits-s, initds-d)
c***keywords  chebyshev, fnlib, initialize, orthogonal polynomial,
c             orthogonal series, special functions
c***author  fullerton, w., (lanl)
c***description
c
c  initialize the orthogonal series, represented by the array os, so
c  that inits is the number of terms needed to insure the error is no
c  larger than eta.  ordinarily, eta will be chosen to be one-tenth
c  machine precision.
c
c             input arguments --
c   os     single precision array of nos coefficients in an orthogonal
c          series.
c   nos    number of coefficients in os.
c   eta    single precision scalar containing requested accuracy of
c          series.
c
c***references  (none)
c***routines called  xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890831  modified array declarations.  (wrb)
c   891115  modified error message.  (wrb)
c   891115  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  inits
      real os(*)
c***first executable statement  inits
      if (nos .lt. 1) call xermsg ('slatec', 'inits',
     +   'number of coefficients is less than 1', 2, 1)
c
      err = 0.
      do 10 ii = 1,nos
        i = nos + 1 - ii
        err = err + abs(os(i))
        if (err.gt.eta) go to 20
   10 continue
c
   20 if (i .eq. nos) call xermsg ('slatec', 'inits',
     +   'chebyshev series too short for specified accuracy', 1, 1)
      inits = i
c
      return
      end
