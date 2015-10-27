*deck dbsppp
      subroutine dbsppp (t, a, n, k, ldc, c, xi, lxi, work)
c***begin prologue  dbsppp
c***purpose  convert the b-representation of a b-spline to the piecewise
c            polynomial (pp) form.
c***library   slatec
c***category  e3, k6
c***type      double precision (bsppp-s, dbsppp-d)
c***keywords  b-spline, piecewise polynomial
c***author  amos, d. e., (snla)
c***description
c
c     written by carl de boor and modified by d. e. amos
c
c     abstract    **** a double precision routine ****
c         dbsppp is the bsplpp routine of the reference.
c
c         dbsppp converts the b-representation (t,a,n,k) to the
c         piecewise polynomial (pp) form (c,xi,lxi,k) for use with
c         dppval.  here xi(*), the break point array of length lxi, is
c         the knot array t(*) with multiplicities removed.  the columns
c         of the matrix c(i,j) contain the right taylor derivatives
c         for the polynomial expansion about xi(j) for the intervals
c         xi(j) .le. x .le. xi(j+1), i=1,k, j=1,lxi.  function dppval
c         makes this evaluation at a specified point x in
c         xi(1) .le. x .le. xi(lxi+1)
c
c     description of arguments
c
c         input      t,a are double precision
c          t       - knot vector of length n+k
c          a       - b-spline coefficient vector of length n
c          n       - number of b-spline coefficients
c                    n = sum of knot multiplicities-k
c          k       - order of the b-spline, k .ge. 1
c          ldc     - leading dimension of c, ldc .ge. k
c
c         output     c,xi,work are double precision
c          c       - matrix of dimension at least (k,lxi) containing
c                    right derivatives at break points
c          xi      - xi break point vector of length lxi+1
c          lxi     - number of break points, lxi .le. n-k+1
c          work    - work vector of length k*(n+3)
c
c     error conditions
c         improper input is a fatal error
c
c***references  carl de boor, package for calculating with b-splines,
c                 siam journal on numerical analysis 14, 3 (june 1977),
c                 pp. 441-472.
c***routines called  dbspdr, dbspev, xermsg
c***revision history  (yymmdd)
c   800901  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dbsppp
c
      integer ileft, inev, k, ldc, lxi, n, nk
      double precision a, c, t, work, xi
c     dimension t(n+k),xi(lxi+1),c(ldc,*)
c     here, * = the final value of the output parameter lxi.
      dimension t(*), a(*), work(*), xi(*), c(ldc,*)
c***first executable statement  dbsppp
      if(k.lt.1) go to 100
      if(n.lt.k) go to 105
      if(ldc.lt.k) go to 110
      call dbspdr(t, a, n, k, k, work)
      lxi = 0
      xi(1) = t(k)
      inev = 1
      nk = n*k + 1
      do 10 ileft=k,n
        if (t(ileft+1).eq.t(ileft)) go to 10
        lxi = lxi + 1
        xi(lxi+1) = t(ileft+1)
        call dbspev(t,work(1),n,k, k,xi(lxi),inev,c(1,lxi),work(nk))
   10 continue
      return
  100 continue
      call xermsg ('slatec', 'dbsppp', 'k does not satisfy k.ge.1', 2,
     +   1)
      return
  105 continue
      call xermsg ('slatec', 'dbsppp', 'n does not satisfy n.ge.k', 2,
     +   1)
      return
  110 continue
      call xermsg ('slatec', 'dbsppp', 'ldc does not satisfy ldc.ge.k',
     +   2, 1)
      return
      end
