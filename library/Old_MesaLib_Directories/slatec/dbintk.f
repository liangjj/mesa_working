*deck dbintk
      subroutine dbintk (x, y, t, n, k, bcoef, q, work)
c***begin prologue  dbintk
c***purpose  compute the b-representation of a spline which interpolates
c            given data.
c***library   slatec
c***category  e1a
c***type      double precision (bintk-s, dbintk-d)
c***keywords  b-spline, data fitting, interpolation
c***author  amos, d. e., (snla)
c***description
c
c     written by carl de boor and modified by d. e. amos
c
c     abstract    **** a double precision routine ****
c
c         dbintk is the splint routine of the reference.
c
c         dbintk produces the b-spline coefficients, bcoef, of the
c         b-spline of order k with knots t(i), i=1,...,n+k, which
c         takes on the value y(i) at x(i), i=1,...,n.  the spline or
c         any of its derivatives can be evaluated by calls to dbvalu.
c
c         the i-th equation of the linear system a*bcoef = b for the
c         coefficients of the interpolant enforces interpolation at
c         x(i), i=1,...,n.  hence, b(i) = y(i), for all i, and a is
c         a band matrix with 2k-1 bands if a is invertible.  the matrix
c         a is generated row by row and stored, diagonal by diagonal,
c         in the rows of q, with the main diagonal going into row k.
c         the banded system is then solved by a call to dbnfac (which
c         constructs the triangular factorization for a and stores it
c         again in q), followed by a call to dbnslv (which then
c         obtains the solution bcoef by substitution).  dbnfac does no
c         pivoting, since the total positivity of the matrix a makes
c         this unnecessary.  the linear system to be solved is
c         (theoretically) invertible if and only if
c                 t(i) .lt. x(i) .lt. t(i+k),        for all i.
c         equality is permitted on the left for i=1 and on the right
c         for i=n when k knots are used at x(1) or x(n).  otherwise,
c         violation of this condition is certain to lead to an error.
c
c     description of arguments
c
c         input       x,y,t are double precision
c           x       - vector of length n containing data point abscissa
c                     in strictly increasing order.
c           y       - corresponding vector of length n containing data
c                     point ordinates.
c           t       - knot vector of length n+k
c                     since t(1),..,t(k) .le. x(1) and t(n+1),..,t(n+k)
c                     .ge. x(n), this leaves only n-k knots (not nec-
c                     essarily x(i) values) interior to (x(1),x(n))
c           n       - number of data points, n .ge. k
c           k       - order of the spline, k .ge. 1
c
c         output      bcoef,q,work are double precision
c           bcoef   - a vector of length n containing the b-spline
c                     coefficients
c           q       - a work vector of length (2*k-1)*n, containing
c                     the triangular factorization of the coefficient
c                     matrix of the linear system being solved.  the
c                     coefficients for the interpolant of an
c                     additional data set (x(i),yy(i)), i=1,...,n
c                     with the same abscissa can be obtained by loading
c                     yy into bcoef and then executing
c                         call dbnslv (q,2k-1,n,k-1,k-1,bcoef)
c           work    - work vector of length 2*k
c
c     error conditions
c         improper input is a fatal error
c         singular system of equations is a fatal error
c
c***references  d. e. amos, computation with splines and b-splines,
c                 report sand78-1968, sandia laboratories, march 1979.
c               carl de boor, package for calculating with b-splines,
c                 siam journal on numerical analysis 14, 3 (june 1977),
c                 pp. 441-472.
c               carl de boor, a practical guide to splines, applied
c                 mathematics series 27, springer-verlag, new york,
c                 1978.
c***routines called  dbnfac, dbnslv, dbspvn, xermsg
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dbintk
c
      integer iflag, iwork, k, n, i, ilp1mx, j, jj, km1, kpkm2, left,
     1 lenq, np1
      double precision bcoef(*), y(*), q(*), t(*), x(*), xi, work(*)
c     dimension q(2*k-1,n), t(n+k)
c***first executable statement  dbintk
      if(k.lt.1) go to 100
      if(n.lt.k) go to 105
      jj = n - 1
      if(jj.eq.0) go to 6
      do 5 i=1,jj
      if(x(i).ge.x(i+1)) go to 110
    5 continue
    6 continue
      np1 = n + 1
      km1 = k - 1
      kpkm2 = 2*km1
      left = k
c                zero out all entries of q
      lenq = n*(k+km1)
      do 10 i=1,lenq
        q(i) = 0.0d0
   10 continue
c
c  ***   loop over i to construct the  n  interpolation equations
      do 50 i=1,n
        xi = x(i)
        ilp1mx = min(i+k,np1)
c        *** find  left  in the closed interval (i,i+k-1) such that
c                t(left) .le. x(i) .lt. t(left+1)
c        matrix is singular if this is not possible
        left = max(left,i)
        if (xi.lt.t(left)) go to 80
   20   if (xi.lt.t(left+1)) go to 30
        left = left + 1
        if (left.lt.ilp1mx) go to 20
        left = left - 1
        if (xi.gt.t(left+1)) go to 80
c        *** the i-th equation enforces interpolation at xi, hence
c        a(i,j) = b(j,k,t)(xi), all j. only the  k  entries with  j =
c        left-k+1,...,left actually might be nonzero. these  k  numbers
c        are returned, in  bcoef (used for temp. storage here), by the
c        following
   30   call dbspvn(t, k, k, 1, xi, left, bcoef, work, iwork)
c        we therefore want  bcoef(j) = b(left-k+j)(xi) to go into
c        a(i,left-k+j), i.e., into  q(i-(left+j)+2*k,(left+j)-k) since
c        a(i+j,j)  is to go into  q(i+k,j), all i,j,  if we consider  q
c        as a two-dim. array , with  2*k-1  rows (see comments in
c        dbnfac). in the present program, we treat  q  as an equivalent
c        one-dimensional array (because of fortran restrictions on
c        dimension statements) . we therefore want  bcoef(j) to go into
c        entry
c            i -(left+j) + 2*k + ((left+j) - k-1)*(2*k-1)
c                   =  i-left+1 + (left -k)*(2*k-1) + (2*k-2)*j
c        of  q .
        jj = i - left + 1 + (left-k)*(k+km1)
        do 40 j=1,k
          jj = jj + kpkm2
          q(jj) = bcoef(j)
   40   continue
   50 continue
c
c     ***obtain factorization of  a  , stored again in  q.
      call dbnfac(q, k+km1, n, km1, km1, iflag)
      go to (60, 90), iflag
c     *** solve  a*bcoef = y  by backsubstitution
   60 do 70 i=1,n
        bcoef(i) = y(i)
   70 continue
      call dbnslv(q, k+km1, n, km1, km1, bcoef)
      return
c
c
   80 continue
      call xermsg ('slatec', 'dbintk',
     +   'some abscissa was not in the support of the corresponding ' //
     +   'basis function and the system is singular.', 2, 1)
      return
   90 continue
      call xermsg ('slatec', 'dbintk',
     +   'the system of solver detects a singular system although ' //
     +   'the theoretical conditions for a solution were satisfied.',
     +   8, 1)
      return
  100 continue
      call xermsg ('slatec', 'dbintk', 'k does not satisfy k.ge.1', 2,
     +   1)
      return
  105 continue
      call xermsg ('slatec', 'dbintk', 'n does not satisfy n.ge.k', 2,
     +   1)
      return
  110 continue
      call xermsg ('slatec', 'dbintk',
     +   'x(i) does not satisfy x(i).lt.x(i+1) for some i', 2, 1)
      return
      end
