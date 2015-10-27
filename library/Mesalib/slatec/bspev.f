*deck bspev
      subroutine bspev (t, ad, n, k, nderiv, x, inev, svalue, work)
c***begin prologue  bspev
c***purpose  calculate the value of the spline and its derivatives from
c            the b-representation.
c***library   slatec
c***category  e3, k6
c***type      single precision (bspev-s, dbspev-d)
c***keywords  b-spline, data fitting, interpolation, splines
c***author  amos, d. e., (snla)
c***description
c
c     written by carl de boor and modified by d. e. amos
c
c     abstract
c         bspev is the bsplev routine of the reference.
c
c         bspev calculates the value of the spline and its derivatives
c         at x from the b-representation (t,a,n,k) and returns them
c         in svalue(i),i=1,nderiv, t(k) .le. x .le. t(n+1).  ad(i) can
c         be the b-spline coefficients a(i), i=1,n if nderiv=1.  other-
c         wise ad must be computed before hand by a call to bspdr (t,a,
c         n,k,nderiv,ad).  if x=t(i),i=k,n, right limiting values are
c         obtained.
c
c         to compute left derivatives or left limiting values at a
c         knot t(i), replace n by i-1 and set x=t(i), i=k+1,n+1.
c
c         bspev calls intrv, bspvn
c
c     description of arguments
c         input
c          t       - knot vector of length n+k
c          ad      - vector of length (2*n-nderiv+1)*nderiv/2 containing
c                    the difference table from bspdr.
c          n       - number of b-spline coefficients
c                    n = sum of knot multiplicities-k
c          k       - order of the b-spline, k .ge. 1
c          nderiv  - number of derivatives, 1 .le. nderiv .le. k.
c                    nderiv=1 gives the zero-th derivative = function
c                    value
c          x       - argument, t(k) .le. x .le. t(n+1)
c          inev    - an initialization parameter which must be set
c                    to 1 the first time bspev is called.
c
c         output
c          inev    - inev contains information for efficient process-
c                    ing after the initial call and inev must not
c                    be changed by the user.  distinct splines require
c                    distinct inev parameters.
c          svalue  - vector of length nderiv containing the spline
c                    value in svalue(1) and the nderiv-1 derivatives
c                    in the remaining components.
c          work    - work vector of length 3*k
c
c     error conditions
c         improper input is a fatal error.
c
c***references  carl de boor, package for calculating with b-splines,
c                 siam journal on numerical analysis 14, 3 (june 1977),
c                 pp. 441-472.
c***routines called  bspvn, intrv, xermsg
c***revision history  (yymmdd)
c   800901  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  bspev
c
      integer i,id,inev,iwork,jj,k,kp1,kp1mn,l,left,ll,mflag,
     1 n, nderiv
      real ad, svalue, sum, t, work, x
c     dimension t(n+k)
      dimension t(*), ad(*), svalue(*), work(*)
c***first executable statement  bspev
      if(k.lt.1) go to 100
      if(n.lt.k) go to 105
      if(nderiv.lt.1 .or. nderiv.gt.k) go to 115
      id = nderiv
      call intrv(t, n+1, x, inev, i, mflag)
      if (x.lt.t(k)) go to 110
      if (mflag.eq.0) go to 30
      if (x.gt.t(i)) go to 110
   20 if (i.eq.k) go to 120
      i = i - 1
      if (x.eq.t(i)) go to 20
c
c *i* has been found in (k,n) so that t(i) .le. x .lt. t(i+1)
c     (or .le. t(i+1), if t(i) .lt. t(i+1) = t(n+1) ).
   30 kp1mn = k + 1 - id
      kp1 = k + 1
      call bspvn(t, kp1mn, k, 1, x, i, work(1),work(kp1),iwork)
      jj = (n+n-id+2)*(id-1)/2
c     adif(leftpl,id) = ad(leftpl-id+1 + (2*n-id+2)*(id-1)/2)
c     leftpl = left + l
   40 left = i - kp1mn
      sum = 0.0e0
      ll = left + jj + 2 - id
      do 50 l=1,kp1mn
        sum = sum + work(l)*ad(ll)
        ll = ll + 1
   50 continue
      svalue(id) = sum
      id = id - 1
      if (id.eq.0) go to 60
      jj = jj-(n-id+1)
      kp1mn = kp1mn + 1
      call bspvn(t, kp1mn, k, 2, x, i, work(1), work(kp1),iwork)
      go to 40
c
   60 return
c
c
  100 continue
      call xermsg ('slatec', 'bspev', 'k does not satisfy k.ge.1', 2,
     +   1)
      return
  105 continue
      call xermsg ('slatec', 'bspev', 'n does not satisfy n.ge.k', 2,
     +   1)
      return
  110 continue
      call xermsg ('slatec', 'bspev', 'x is not in t(k).le.x.le.t(n+1)'
     +   , 2, 1)
      return
  115 continue
      call xermsg ('slatec', 'bspev',
     +   'nderiv does not satisfy 1.le.nderiv.le.k', 2, 1)
      return
  120 continue
      call xermsg ('slatec', 'bspev',
     +   'a left limiting value cannot be obtained at t(k)', 2, 1)
      return
      end
