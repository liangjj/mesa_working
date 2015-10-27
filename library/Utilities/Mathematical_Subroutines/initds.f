*deck initds
      function initds (os, nos, eta)
c***begin prologue  initds
c***purpose  determine the number of terms needed in an orthogonal
c            polynomial series so that it meets a specified accuracy.
c***library   slatec (fnlib)
c***category  c3a2
c***type      double precision (inits-s, initds-d)
c***keywords  chebyshev, fnlib, initialize, orthogonal polynomial,
c             orthogonal series, special functions
c***author  fullerton, w., (lanl)
c***description
c
c  initialize the orthogonal series, represented by the array os, so
c  that initds is the number of terms needed to insure the error is no
c  larger than eta.  ordinarily, eta will be chosen to be one-tenth
c  machine precision.
c
c             input arguments --
c   os     double precision array of nos coefficients in an orthogonal
c          series.
c   nos    number of coefficients in os.
c   eta    single precision scalar containing requested accuracy of
c          series.
c
c***references  (none)
c***routines called  xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891115  modified error message.  (wrb)
c   891115  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  initds
      real*8 os(*)
c***first executable statement  initds
      if (nos .lt. 1) call lnkerr ('initds:no. coef. less than 1')
c
      err = 0.
      do 10 ii = 1,nos
        i = nos + 1 - ii
        err = err + abs(real(os(i)))
        if (err.gt.eta) go to 20
   10 continue
c
   20 if (i .eq. nos) call lnkerr ('initds:cheby series too short')
      initds = i
c
      return
      end
