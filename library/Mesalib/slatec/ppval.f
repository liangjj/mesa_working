*deck ppval
      function ppval (ldc, c, xi, lxi, k, ideriv, x, inppv)
c***begin prologue  ppval
c***purpose  calculate the value of the ideriv-th derivative of the
c            b-spline from the pp-representation.
c***library   slatec
c***category  e3, k6
c***type      single precision (ppval-s, dppval-d)
c***keywords  b-spline, data fitting, interpolation, splines
c***author  amos, d. e., (snla)
c***description
c
c     written by carl de boor and modified by d. e. amos
c
c     abstract
c         ppval is the ppvalu function of the reference.
c
c         ppval calculates (at x) the value of the ideriv-th
c         derivative of the b-spline from the pp-representation
c         (c,xi,lxi,k).  the taylor expansion about xi(j) for x in
c         the interval xi(j) .le. x .lt. xi(j+1) is evaluated, j=1,lxi.
c         right limiting values at x=xi(j) are obtained.  ppval will
c         extrapolate beyond xi(1) and xi(lxi+1).
c
c         to obtain left limiting values (left derivatives) at xi(j),
c         replace lxi by j-1 and set x=xi(j),j=2,lxi+1.
c
c     description of arguments
c         input
c          ldc     - leading dimension of c matrix, ldc .ge. k
c          c       - matrix of dimension at least (k,lxi) containing
c                    right derivatives at break points xi(*).
c          xi      - break point vector of length lxi+1
c          lxi     - number of polynomial pieces
c          k       - order of b-spline, k .ge. 1
c          ideriv  - order of the derivative, 0 .le. ideriv .le. k-1
c                    ideriv=0 gives the b-spline value
c          x       - argument, xi(1) .le. x .le. xi(lxi+1)
c          inppv   - an initialization parameter which must be set
c                    to 1 the first time ppval is called.
c
c         output
c          inppv   - inppv contains information for efficient process-
c                    ing after the initial call and inppv must not
c                    be changed by the user.  distinct splines require
c                    distinct inppv parameters.
c          ppval   - value of the ideriv-th derivative at x
c
c     error conditions
c         improper input is a fatal error
c
c***references  carl de boor, package for calculating with b-splines,
c                 siam journal on numerical analysis 14, 3 (june 1977),
c                 pp. 441-472.
c***routines called  intrv, xermsg
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920501  reformatted the references section.  (wrb)
c***end prologue  ppval
c
      integer i, ideriv, inppv, j, k, ldc, lxi, ndummy
      real c, dx, fltk, x, xi
      dimension xi(*), c(ldc,*)
c***first executable statement  ppval
      ppval = 0.0e0
      if(k.lt.1) go to 90
      if(ldc.lt.k) go to 80
      if(lxi.lt.1) go to 85
      if(ideriv.lt.0 .or. ideriv.ge.k) go to 95
      i = k - ideriv
      fltk = i
      call intrv(xi, lxi, x, inppv, i, ndummy)
      dx = x - xi(i)
      j = k
   10 ppval = (ppval/fltk)*dx + c(j,i)
      j = j - 1
      fltk = fltk - 1.0e0
      if (fltk.gt.0.0e0) go to 10
      return
c
c
   80 continue
      call xermsg ('slatec', 'ppval', 'ldc does not satisfy ldc.ge.k',
     +   2, 1)
      return
   85 continue
      call xermsg ('slatec', 'ppval', 'lxi does not satisfy lxi.ge.1',
     +   2, 1)
      return
   90 continue
      call xermsg ('slatec', 'ppval', 'k does not satisfy k.ge.1', 2,
     +   1)
      return
   95 continue
      call xermsg ('slatec', 'ppval',
     +   'ideriv does not satisfy 0.le.ideriv.lt.k', 2, 1)
      return
      end
