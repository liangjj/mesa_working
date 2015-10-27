*deck dbspdr
      subroutine dbspdr (t, a, n, k, nderiv, ad)
c***begin prologue  dbspdr
c***purpose  use the b-representation to construct a divided difference
c            table preparatory to a (right) derivative calculation.
c***library   slatec
c***category  e3, k6
c***type      double precision (bspdr-s, dbspdr-d)
c***keywords  b-spline, data fitting, differentiation of splines,
c             interpolation
c***author  amos, d. e., (snla)
c***description
c
c     written by carl de boor and modified by d. e. amos
c
c     abstract     **** a double precision routine ****
c         dbspdr is the bspldr routine of the reference.
c
c         dbspdr uses the b-representation (t,a,n,k) to construct a
c         divided difference table adif preparatory to a (right)
c         derivative calculation in dbspev.  the lower triangular matrix
c         adif is stored in vector ad by columns.  the arrays are
c         related by
c
c           adif(i,j) = ad(i-j+1 + (2*n-j+2)*(j-1)/2)
c
c         i = j,n   ,   j=1,nderiv.
c
c     description of arguments
c
c         input      t,a are double precision
c          t       - knot vector of length n+k
c          a       - b-spline coefficient vector of length n
c          n       - number of b-spline coefficients
c                    n = sum of knot multiplicities-k
c          k       - order of the spline, k .ge. 1
c          nderiv  - number of derivatives, 1 .le. nderiv .le. k.
c                    nderiv=1 gives the zero-th derivative =
c                    function value
c
c         output     ad is double precision
c          ad      - table of differences in a vector of length
c                    (2*n-nderiv+1)*nderiv/2 for input to dbspev
c
c     error conditions
c         improper input is a fatal error
c
c***references  carl de boor, package for calculating with b-splines,
c                 siam journal on numerical analysis 14, 3 (june 1977),
c                 pp. 441-472.
c***routines called  xermsg
c***revision history  (yymmdd)
c   800901  date written
c   890831  modified array declarations.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dbspdr
c
c
      integer i, id, ii, ipkmid, jj, jm, k, kmid, n, nderiv
      double precision a, ad, diff, fkmid, t
c     dimension t(n+k), ad((2*n-nderiv+1)*nderiv/2)
      dimension t(*), a(*), ad(*)
c***first executable statement  dbspdr
      if(k.lt.1) go to 100
      if(n.lt.k) go to 105
      if(nderiv.lt.1 .or. nderiv.gt.k) go to 110
      do 10 i=1,n
        ad(i) = a(i)
   10 continue
      if (nderiv.eq.1) return
      kmid = k
      jj = n
      jm = 0
      do 30 id=2,nderiv
        kmid = kmid - 1
        fkmid = kmid
        ii = 1
        do 20 i=id,n
          ipkmid = i + kmid
          diff = t(ipkmid) - t(i)
          if (diff.ne.0.0d0) ad(ii+jj) = (ad(ii+jm+1)-ad(ii+jm))/
     1     diff*fkmid
          ii = ii + 1
   20   continue
        jm = jj
        jj = jj + n - id + 1
   30 continue
      return
c
c
  100 continue
      call xermsg ('slatec', 'dbspdr', 'k does not satisfy k.ge.1', 2,
     +   1)
      return
  105 continue
      call xermsg ('slatec', 'dbspdr', 'n does not satisfy n.ge.k', 2,
     +   1)
      return
  110 continue
      call xermsg ('slatec', 'dbspdr',
     +   'nderiv does not satisfy 1.le.nderiv.le.k', 2, 1)
      return
      end
